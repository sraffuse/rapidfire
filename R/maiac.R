
# Data pool location for MAIAC MCD19A2
# https://e4ftl01.cr.usgs.gov/MOTA/MCD19A2.006/


# For a given bounding box, return the MODIS tiles needed for coverage using the
# MODIS siusoidal grid. Definition found at
# https://modis-land.gsfc.nasa.gov/pdf/sn_bound_10deg.txt
modis_tile_pick <- function(llat, ulat, llon, rlon) {

  grid <- read.table("./data/MAIAC/modis_grid/sn_bound_10deg.txt", header = TRUE,
                   row.names = NULL, sep = "", skip = 6, nrow = 648)
  # Get tile for each corner
  ur <- grid %>%
    filter(lon_min < rlon, lon_max > rlon,
           lat_min < ulat, lat_max > ulat)
  lr <- grid %>%
    filter(lon_min < rlon, lon_max > rlon,
           lat_min < llat, lat_max > llat)
  ul <- grid %>%
    filter(lon_min < llon, lon_max > llon,
           lat_min < ulat, lat_max > ulat)
  ll <- grid %>%
    filter(lon_min < llon, lon_max > llon,
           lat_min < llat, lat_max > llat)

  all <- bind_rows(ur, lr, ul, ll) %>%
    distinct()

}

maiac_aod <- function(fname) {
  # sds <- MODIS::getSds(fname)
  # tile <- raster::raster(rgdal::readGDAL(sds$SDS4gdal[6], as.is = TRUE))
  sds <- terra::sds(fname)
  stack <- sds$Optical_Depth_047
  # There are multiple overpasses in the file, so we combine them with avg
  r <- terra::mean(stack, na.rm = TRUE)

}

maiac_mosaic <- function(tiles) {

  purrr::reduce(tiles, terra::mosaic, fun = "mean")

}

# Use raster::focal to fill in some missing values by interpolating neighbors,
# Then replace the remaining missing values with the median value

maiac_fill_gaps <- function(maiac, window = 7) {

  md <- terra::global(maiac, fun = median, na.rm = TRUE) %>%
    .$global
  sr <- terra::focal(maiac, w = window, fun = "mean", na.rm = TRUE, na.only = TRUE)

  blank_space <- sr == 0
  fill <- blank_space * md
  filled <- sr + fill

}


extract_maiac <- function(maiac, locs) {

  sin_proj <- sp::CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  loc_sin <- sp::spTransform(locs, sin_proj)
  vals <- terra::extract(maiac, loc_sin@coords)

}

maiac_one_day <- function(dt, input_path = "./data/MAIAC/",
                          tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  date_str <- strftime(dt, format = "%Y%j")
  file_prefix <- paste0("MCD19A2.A", date_str)

  all_files <- fs::path_file(fs::dir_ls(input_path))

  pattern <- paste0(file_prefix, ".(?:", paste(tiles_needed, collapse = "|"), ").*hdf")

  files <- stringr::str_subset(all_files, pattern = pattern)
  paths <- paste0(input_path, files)

  tiles <- purrr::map(paths, maiac_aod)
  maiac <- maiac_mosaic(tiles) %>%
    maiac_fill_gaps()

}

maiac_at_airnow <- function(an, maiac_path = "./data/MAIAC/",
                            tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  daily_extract <- function(dt) {
    print(dt)
    maiac_aod <- maiac_one_day(dt, input_path = maiac_path, tiles_needed = tiles_needed)
    i <- an$Day == dt
    locs = an[i, ]
    e <- extract_maiac(maiac_aod, locs)
    df <- locs@data %>%
      mutate(MAIAC_AOD = e$mean)
  }

  dates <- unique(an$Day)
  purrr::map_dfr(dates, daily_extract)

}

maiac_at_grid <- function(start, end, grid, maiac_path = "./data/MAIAC/",
                          tiles_needed = c("h08v04", "h08v05", "h09v04")) {

  daily_extract <- function(dt) {
    print(dt)
    maiac_aod <- maiac_one_day(dt, input_path = maiac_path, tiles_needed = tiles_needed)
    r <- raster::raster(maiac_aod)

    on_grid <- extract_maiac(maiac_aod, grid)

    df <- grid@data %>%
      mutate(MAIAC_AOD = on_grid$mean,
             Day = dt)
  }

  dates <- seq.Date(from = start, to = end, by = "1 day")
  purrr::map_dfr(dates, daily_extract)


}
