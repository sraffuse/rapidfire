# Test the prediction developed in test_processing - just do one month because R
# crashed when I tried to krige purpleair over two months on a 4 km grid

library(rapidfire)

# Now, we need to create an output grid and calculate all of the variables on
# that grid, using the full set of AN data, then we can output results.
# Cannot do finer grids yet because cannot allocate vector. Can we get around this?

dt1 <- as.Date("2020-09-01")
dt2 <- as.Date("2020-09-30")

# Get and prep AirNow data
an_ws <- get_airnow_daterange(dt1, dt2, "CA")
an <- recast_monitors(an_ws)
an_vg <- create_airnow_variograms(an)

# Get and prep PurpleAir data
pa <- get_purpleair_daterange(dt1, dt2, "CA")
pas <- purpleair_spatial(pa)
pa_clean <- purpleair_clean_spatial_outliers(pas)
pa_vg <- create_purpleair_variograms(pa_clean, cutoff = 100000)


grid_4k <- make_grid(pa_clean, cellsize = 4000)

an_4k <- krige_airnow_all(an, grid_4k, an_vg)

# These are pretty fast
narr_4k <- narr_at_grid(dt1, dt2, grid_4k)
bluesky_4k <- bluesky_at_grid(dt1, dt2, grid_4k)

# This takes many minutes per day (6-8 hours for a month)
pa_4k <- krige_purpleair_all(pa_clean, grid_4k, pa_vg)
saveRDS(pa_4k, "data/pa_4k.RDS")


# Prepare for RF prediction using premade model
results_4k <- an_4k %>%
  inner_join(pa_4k, by = c("Id", "Day")) %>%
  inner_join(narr_4k, by = c("Id", "Day")) %>%
  inner_join(bluesky_4k, by = c("Id", "Day"))

# Get the results of the regression trees
results_4k$PM25_log_RF <- predict(mod1.cv$finalModel, results_4k)

# Make some sanity check plots



# Would be good to use the Airsis data as a validation tool, as they are 100%
# independent at this point
