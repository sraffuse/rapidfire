# Use an existing developed model to predict values at given new locations and dates

predict_locs <- function(dt1, dt2, states = "CA", model, locations,
                         pa_cutoff = 100000, bluesky_special = NULL,
                         pa_data = NULL) {

  # Get and prep AirNow and AirSIS data
  print("AirNow and AirSIS data...")
  an_ws <- get_airnow_daterange(dt1, dt2, states) %>%
    recast_monitors()
  as_ws <- get_airsis_daterange(dt1, dt2, states) %>%
    recast_monitors()
  mon <- rbind(an_ws, as_ws)

  an_vg <- create_airnow_variograms(mon)
  ank <- krige_airnow_sitedates(mon, locations, an_vg)

  # Prep Bluesky For 2020, the grid changed midway through october - will need
  # to process them separately. I added this hack to address it
  print("BlueSky data...")
  if (!is.null(bluesky_special)) {
    if (bluesky_special == "2020") {
      bluesky <- bluesky_archive_at_locs(locations)
    } else if (bluesky_special == "HAQAST") {
      #Use the HAQAST S1 CMAQ run instead of the standard Bluesky archive
      stop("Need to implement bluesky_haqast_at_locs")
      bluesky_stack <- stack_haqast_archive(dt1, dt2)
      bluesky <- preprocessed_bluesky_at_airnow(bluesky_stack, mon)
    } else {
      stop(paste("bluesky_special:", period, "not supported"))
    }
  } else {
    bluesky <- bluesky_archive_at_locs(locations)
  }

  # Prep MAIAC AOD
  print("Preparing MAIAC...")
  maiac <- maiac_at_airnow(locations)

  # Prep NARR
  print("Preparing NARR...")
  narr <- narr_at_airnow(locations)

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
  pa_ok <- krige_purpleair_sitedates(pa_clean, locations, pa_vg)

  # All the inputs
  # PurpleAir Kriged
  pak_in <- pa_ok %>%
    select(monitorID, Day, PM25_log_PAK)
  # NARR
  narr_in <- narr %>%
    select(-one_of("PM25"), -one_of("Source"))
  # Bluesky
  bluesky_in <- bluesky %>%
    select(monitorID, Day, PM25_bluesky)
  # MAIAC AOD
  maiac_in <- maiac %>%
    select(monitorID, Day, MAIAC_AOD)

  results <- ank %>%
    left_join(pak_in, by = c("monitorID", "Day")) %>%
    left_join(narr_in, by = c("monitorID", "Day")) %>%
    left_join(bluesky_in, by = c("monitorID", "Day")) %>%
    left_join(maiac_in, by = c("monitorID", "Day")) %>%
  # Replace missing and non-finite values with overall median
    mutate(across(where(is.numeric),
                  ~ifelse(is.finite(.x), .x, median(.x, na.rm = TRUE))))

  results$PM25_log_RF <- predict(model, results)
  return(results)
}
