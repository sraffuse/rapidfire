# Build PurpleAir Archive for CA using new API
#

read_key <- "insert read key here"

nwlng <- -124.41
nwlat <- 42.01
selng <- -114.13
selat <- 32.53

sensors_in_box <- pa_find_sensors(nwlng, nwlat, selng, selat, key = read_key)

# Limit to within state
ca <- sf::st_read("U:/SharedData/GIS/states_hi_res.shp") %>%
  filter(STATE_NAME == "California")

sens_sf <- sf::st_as_sf(sensors_in_box, coords = c("longitude", "latitude"),
                        crs = 4326)
ints <- sf::st_intersects(sens_sf, ca, sparse = FALSE)
pa_ca <- sensors_in_box[ints[,1],]
saveRDS(pa_ca, "../rapidfire_project/data/purpleair/sensors/sensor_locs.RDS")

runs <- pa_ca %>%
  select(sensor_index, start_date=date_created, end_date=last_seen) %>%
  mutate(key = read_key)

run_test <- slice(runs, 1:4)

archive <- purrr::pmap_dfr(run_test, pa_sensor_history, .progress = "TRUE")

# Actaully, let's save each monitor to RDS in case this fails at any time, so we
# can pick up where we left off
save_pa <- function(sensor_index, start_date, end_date, key) {
  df <- pa_sensor_history(sensor_index, start_date, end_date, key)
  saveRDS(df, glue::glue("../rapidfire_project/data/purpleair/sensors/{sensor_index}.RDS"))
}

# Limit to 2017-2021 (should still have most of the sites)
runs2 <- runs %>%
  mutate(start_date = if_else(start_date < as.Date("2017-01-01"),
                              as.POSIXct("2017-01-01", origin = "1970-01-01"),
                              start_date),
         end_date = if_else(end_date > as.Date("2022-01-01"),
                            as.POSIXct("2022-01-01", origin = "1970-01-01"),
                            end_date)) %>%
  filter(start_date < end_date)

# Okay, here goes nothing!
purrr::pwalk(runs2, save_pa, .progress = TRUE)

