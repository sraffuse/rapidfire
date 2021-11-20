
# PurpleAir will be one of the inputs to the RF model, after a small bit of
# scrubbing and interpolation by ordinary kriging

# pa <- get_purpleair_daterange(dt, dt, c("CA", "NV"))
get_purpleair_daterange <- function(start, end, states) {
  AirSensor::setArchiveBaseUrl("https://airfire-data-exports.s3-us-west-2.amazonaws.com/PurpleAir/v1")

  get_one_day <- function(dt, states) {

    # only include sensors that reported recently (hr 20 on this day)
    valid_datetime <- as.POSIXct(dt) + (20 * 60 * 60)
    AirSensor::pas_load(strftime(dt, format = "%Y%m%d"),
                        timezone = "America/Los_Angeles") %>%
      filter(stateCode %in% states,
             DEVICE_LOCATIONTYPE == "outside",
             lastSeenDate > valid_datetime,
             statsLastModifiedDate > valid_datetime) %>%
      mutate(Date = dt)
  }

  dates <- seq.Date(from = start, to = end, by = "1 day")
  purrr::map_dfr(dates, get_one_day, states)

}

# Get log of daily mean, select only one sensor per location, and convert to
# spatialpointsdataframe
# pas <- purpleair_spatial(pa)
purpleair_spatial <- function(df) {
  df <- df %>%
    filter(pm25_1day > 0) %>%
    mutate(Purple_log = log(pm25_1day)) %>%
    select(Date, deviceID, latitude, longitude, pm25_1day, Purple_log) %>%
    group_by(Date, latitude, longitude) %>%
    summarise(deviceID = first(deviceID),
              pm25_1day = first(pm25_1day),
              Purple_log = first(Purple_log),
              .groups = "drop")

  sp::coordinates(df) <- ~longitude+latitude
  sp::proj4string(df) <- sp::CRS("+init=epsg:4326")
  sp::spTransform(df, sp::CRS("+init=epsg:3395"))

}

# Scrub spatial outliers - remove any sensor that is > 2 sd away from the median
# within a 10 km radius
# pa_clean <- purplueair_clean_spatial_outliers(pas)
purpleair_clean_spatial_outliers <- function(spdf) {


  region_stats <- function(ind, var) {
    med <- median(var[ind])
    std <- sd(var[ind])
    n <- length(ind)
    data.frame(med, std, n)
  }

  # Need to convert to sf to use the operator we need
  sfs <- sf::st_as_sf(spdf)

  # Could not get grouping to work, so do this one date at a time manually
  datelist <- unique(sfs$Date)

  calc_daily_stats <- function(dt, sfs, spdf) {
    sfs <- filter(sfs, Date == dt)
    within <- sf::st_is_within_distance(sfs, dist = 10000)
    stats <- purrr::map_dfr(within, region_stats, sfs$pm25_1day)
    sfs <- bind_cols(sfs, stats) %>%
      tidyr::replace_na(list(std = 0)) %>%
      mutate(AbsDiff = abs(pm25_1day - med),
             Outlier = if_else(AbsDiff > 2 * std, 1, 0))

    good_devices <- as.data.frame(sfs) %>%
      filter(Outlier == 0) %>%
      .$deviceID
    clean <- spdf[spdf$deviceID %in% good_devices,]
    clean <- clean[clean$Date == dt,]

  }

  clean_spdf <- purrr::map(datelist, calc_daily_stats, sfs, spdf)
  merged <- do.call(rbind, clean_spdf)
}


# Create binned variograms using cleaned dataset, one day at a time
# pa_vgms <- create_purpleair_variograms(pa_clean, cutoff = 100000)
create_purpleair_variograms <- function(df, cutoff = NULL) {

  # If fit does not converge, just pass null
  try_fit <- purrr::possibly(gstat::fit.variogram, otherwise = NULL)

  # intialize a Gaussian variagram model to fit with a nugget of 0.02
  vgm_mod <- gstat::vgm(model = "Gau", nugget = 0.02)

  daily_variogram <- function(day, cutoff) {
    df <- df[df$Date == day,]
    vgm_data <- gstat::variogram(Purple_log ~  1, data = df, cutoff = cutoff, width = 15000) %>%
      mutate(Date = day)
    mod_fit <- try_fit(vgm_data, vgm_mod)
    return(list(d=vgm_data, m=mod_fit))
  }

  dates <- unique(df$Date)
  purrr:::map(dates, daily_variogram, cutoff)

}


# Now run OK day-by-day outputting to the airnow test and training sites
# Run Ordinary Kriging on the training data at the test locations using the
# models created in create_airnow_variograms
krige_purpleair <- function(pa_data, out_locs, vgms) {

  dates <- unique(pa_data$Date)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Date)

  process_one_ok <- function(date, out_locs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Date == date) %>%
      distinct() %>%
      .$row
    test <- out_locs[out_locs$Day == date,]
    mod <- vgms[[row]]$m
    train <- pa_data[pa_data$Date == date,]
    print(date)
    ok <- gstat::krige(Purple_log ~ 1, locations = train, newdata = test,
                       model = mod)

    # Attach to measured values
    test$purp.pred <- ok$var1.pred
    test$purp.var <- ok$var1.var
    test

  }

  all <- purrr::map(dates, process_one_ok, out_locs, vgms, rows)
  # Combine together - this looks like it will be slow
  merged <- do.call(rbind, all)

}

krige_purpleair_all <- function(pa_data, out_locs, vgms) {
  dates <- unique(pa_data$Date)
  rows <- purrr::map_dfr(vgms, ~.x$d, .id = "row") %>%
    mutate(row = as.numeric(row)) %>%
    select(row, Date)

  process_one_ok <- function(date, out_locs, vgms, rows) {
    # extract the correct model data
    row <- filter(rows, Date == date) %>%
      distinct() %>%
      .$row
    mod <- vgms[[row]]$m
    train <- pa_data[pa_data$Date == date,]
    print(date)
    ok <- gstat::krige(Purple_log ~ 1, locations = train, newdata = out_locs,
                       model = mod)

    # Attach to measured values
    output <- out_locs@data
    output$Day <- date
    output$PM25_log_PAK <- ok$var1.pred
    output$PM25_log_PAvar <- ok$var1.var
    output

  }

  all <- purrr::map_dfr(dates, process_one_ok, out_locs, vgms, rows)

}

# Find avaiable purple air sensors for creating an archive (pre-canned archive
# only goes back to April 5, 2019)
create_purpleair_archive <- function(dt1, dt2, states = "CA", pas = NULL,
                                     path = "./data/purpleair") {
  # Get the oldest possible snapshot
  AirSensor::setArchiveBaseUrl("http://data.mazamascience.com/PurpleAir/v1")
  if (is.null(pas)) {
    pas <- AirSensor::pas_load("20190406", archival = TRUE) %>%
      filter(stateCode %in% states)
  }

  ids <- AirSensor::pas_getDeviceDeploymentIDs(pas, isOutside = TRUE)

  pb <- progress::progress_bar$new(total = length(ids))
  pas <- purrr::map_dfr(ids, pat_to_paslike, dt1, dt2, pas, pb)

}

# For a given device and time range, acquire available data raw data and return
# 24-hr average results in a pas-like format
pat_to_paslike <- function(deviceDeploymentID, dt1, dt2, pas, pb) {

  pb$tick()
  try_pat <- purrr::possibly(AirSensor::pat_createNew, otherwise = NULL)
  pat <- try_pat(id = deviceDeploymentID, pas = pas, startdate = dt1,
                 enddate = dt2, verbose = TRUE)
  if (is.null(pat)) {
    return(NULL)
  }

  m <- pat$meta

  # Keep only data when the two channels are tracking
  # Keep only values that either have abs(SRD) < 0.5 or Avg < 2
  df <- pat$data %>%
    mutate(datetime = lubridate::with_tz(datetime, tzone = m$timezone[1]),
           Date = lubridate::floor_date(datetime, unit = "days"),
           SAD = (pm25_A - pm25_B) / sqrt(2),
           Avg = (pm25_A + pm25_B) / 2,
           SRD = SAD/Avg) %>%
    filter(abs(SRD) <= 0.5 | Avg < 2) %>%
    group_by(Day) %>%
    summarise(across(pm25_A:pm10_atm_B, ~ mean(.x, na.rm = TRUE)),
              Count = n()) %>%
    mutate(Fraction = Count / 1440,
           pm25_1day = (pm25_A + pm25_B) / 2,
           ID = m$ID)

  df <- right_join(m, df, by = "ID")

}

