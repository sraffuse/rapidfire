# Automation functions (scripts) for developing models over a given period using
# the standard stuff (AirNow/AirSis, PurpleAir, MAIAC AOD, BlueSky, NARR)

# This function requires MAIAC and NARR data to be downloaded already
develop_model <- function(dt1, dt2, states, pa_cutoff = 100000, seed = 1977,
                          period = NULL) {

  # Get and prep AirNow and AirSIS data
  an_ws <- get_airnow_daterange(dt1, dt2, states) %>%
    recast_monitors()
  as_ws <- get_airsis_daterange(dt1, dt2, states) %>%
    recast_monitors()
  mon <- rbind(an_ws, as_ws)

  an_vg <- create_airnow_variograms(mon)
  mon_split <- split_airnow_data(mon)
  an_ok <- krige_airnow(mon_split, an_vg)


  # Prep Bluesky For 2020, the grid changed midway through october - will need
  # to process them separately. I added this hack to address it
  if (period == "2020") {
    bluesky_stack1 <- stack_bluesky_archive(as.Date("2020-08-01"), as.Date("2020-10-09"))
    bluesky_stack2 <- stack_bluesky_archive(as.Date("2020-10-10"), as.Date("2020-10-31"))
    mon1 <- mon[mon@data$Day <= as.POSIXct("2020-10-09"),]
    mon2 <- mon[mon@data$Day > as.POSIXct("2020-10-09"),]
    ## YOU ARE HERE - The extract in here is not working properly
    bluesky1 <- preprocessed_bluesky_at_airnow(bluesky_stack1, mon1)
    bluesky2 <- preprocessed_bluesky_at_airnow(bluesky_stack2, mon2)
    bluesky <- bind_rows(bluesky1, bluesky2)
  } else {
    bluesky_stack <- stack_bluesky_archive(dt1, dt2)
    bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
  }

  # Prep MAIAC AOD
  maiac <- maiac_at_airnow(mon)

  # Prep NARR
  narr <- narr_at_airnow(mon)

  # Get and prep PurpleAir data
  # For dates prior to 2019, the data must first be acquired by create_purpleair_archive
  pa <- get_purpleair_daterange(dt1, dt2, states)
  pas <- purpleair_spatial(pa)
  pa_clean <- purpleair_clean_spatial_outliers(pas)
  pa_vg <- create_purpleair_variograms(pa_clean, cutoff = pa_cutoff)

  # Krige purple air data for all locations in the monitor data set
  pa_ok <- krige_purpleair(pa_clean, mon, pa_vg)

  # Now need to put altogether in prep for RF modeling
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

  # possibly save these model results?
  # saveRDS(model, "./data/models/xxxxx")
  return(model)

}

