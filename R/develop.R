# Automation functions (scripts) for developing models over a given period using
# the standard stuff (AirNow/AirSis, PurpleAir, MAIAC AOD, BlueSky, NARR)

# This function requires MAIAC and NARR data to be downloaded already
develop_model <- function(dt1, dt2, states, pa_cutoff = 100000, seed = 1977,
                          bluesky_special = NULL, pa_data = NULL) {

  # Get and prep AirNow and AirSIS data
  print("AirNow and AirSIS data...")
  an_ws <- get_airnow_daterange(dt1, dt2, states) %>%
    recast_monitors()
  as_ws <- get_airsis_daterange(dt1, dt2, states) %>%
    recast_monitors()
  mon <- rbind(an_ws, as_ws)

  an_vg <- create_airnow_variograms(mon)
  # trying a larger fraction of test data to make a more complete RF model
  mon_split <- split_airnow_data(mon, test_fraction = 0.3)
  an_ok <- krige_airnow(mon_split, an_vg)

  ank_complete <- krige_airnow_all(mon, mon, an_vg) %>%
    select(monitorID, Day, PM25_log_ANK, PM25_log_var) %>%
    distinct()

  # Prep Bluesky For 2020, the grid changed midway through october - will need
  # to process them separately. I added this hack to address it
  print("BlueSky data...")
  if (!is.null(bluesky_special)) {
    if (bluesky_special == "2020") {
      bluesky_stack1 <- stack_bluesky_archive(dt1, as.Date("2020-10-09"))
      bluesky_stack2 <- stack_bluesky_archive(as.Date("2020-10-10"), dt2)
      mon1 <- mon[mon@data$Day <= as.POSIXct("2020-10-09"),]
      mon2 <- mon[mon@data$Day > as.POSIXct("2020-10-09"),]
      bluesky1 <- preprocessed_bluesky_at_airnow(bluesky_stack1, mon1)
      bluesky2 <- preprocessed_bluesky_at_airnow(bluesky_stack2, mon2)
      bluesky <- bind_rows(bluesky1, bluesky2)
    } else if (bluesky_special == "HAQAST") {
      #Use the HAQAST S1 CMAQ run instead of the standard Bluesky archive
      bluesky_stack <- stack_haqast_archive(dt1, dt2)
      bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
    } else {
      stop(paste("bluesky_special:", period, "not supported"))
    }
  } else {
    bluesky_stack <- stack_bluesky_archive(dt1, dt2)
    bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
  }

  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(mon)

  # Prep NARR
  print("Preparing NARR...")
  narr <- narr_at_airnow(mon)

  # Get and prep PurpleAir data For dates prior to April 5, 2019, the data must
  # be acquired by create_purpleair_archive before running this script
  print("PurpleAir data...")
  if (is.null(pa_data)) {
    pa_data <- get_purpleair_daterange(dt1, dt2, states)
  } else {
    pa_data <- readRDS(pa_data) %>%
      filter(Date >= dt1,
             Date <=dt2)
  }
  pas <- purpleair_spatial(pa_data)
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)

  # Krige purple air data for all locations in the monitor data set
  pa_ok <- krige_purpleair(pa_clean, mon, pa_vg)

  # Now need to put altogether in prep for RF modeling
  print("Combining inputs for modeling...")
  train <- mon_split$train@data %>%
    mutate(Split = "Train")
  test <- mon_split$test@data %>%
    mutate(Split = "Test")

  # These have the measured PM25 - the target variable
  model_in <- bind_rows(train, test)

  # All the inputs
  # AirNow Kriged
  ank_in <- an_ok@data %>%
    select(monitorID, Day, PM25_log_ANK = var1.pred)
  # PurpleAir Kriged
  pak_in <- pa_ok@data %>%
    select(monitorID, Day, PM25_log_PAK = purp.pred)
  # NARR
  narr_in <- narr %>%
    select(-PM25, -Hours, -PM25_log)
  # Bluesky
  bluesky_in <- bluesky %>%
    select(monitorID, Day, PM25_bluesky)
  # MAIAC AOD
  maiac_in <- maiac %>%
    select(monitorID, Day, MAIAC_AOD)

  model_in <- model_in %>%
    left_join(ank_in, by = c("monitorID", "Day")) %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(narr_in, by = c("monitorID", "Day")) %>%
    left_join(bluesky_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day"))

  # Check for missing values
  model_in <- model_in %>%
    filter(!is.na(PM25_bluesky))

  # Replace missing and non-finite values with overall median
  model_in <- model_in %>%
    mutate(across(PM25_log_PAK:MAIAC_AOD,
                  ~if_else(is.finite(.x), .x, median(.x, na.rm = TRUE))))

  # Possibly save the expensive portions
  # saveRDS(model_in, "model_in_Aug2021.RDS")
  # saveRDS(an_vg, "an_vg_Aug2021.RDS")

  # Train the model
  set.seed(seed)
  train_control <- caret::trainControl(method = "cv", number = 10)
  tune_grid <- data.frame(mtry = c(2,3,4,5))

  model_in_test <- filter(model_in, Split == "Test")

  model <- caret::train(PM25_log ~ PM25_log_ANK + PM25_log_PAK + PM25_bluesky + MAIAC_AOD +
                          air.2m + uwnd.10m + vwnd.10m + rhum.2m + apcp + hpbl,
                          data = model_in_test, tuneGrid = tune_grid, do.trace = 100,
                          ntree = 500, trControl = train_control, method = "rf",
                          importance = TRUE)

  print(model)

  # Now, krige all the monitor sites and make a prediction at each. Then return
  # both the model and the data frame including the inputs and prediction.
  output <- bind_rows(train, test) %>%
    left_join(ank_complete, by = c("monitorID", "Day")) %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(narr_in, by = c("monitorID", "Day")) %>%
    left_join(bluesky_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day")) %>%
    filter(!is.na(PM25_bluesky)) %>%
    filter(if_all(where(is.numeric), is.finite))

  output$PM25_log_RF <- predict(model$finalModel, output)
  list(model = model, output=output)

}

