
# Tools for downloading and extracting BlueSky output

grab_bluesky <- function(start, end, output_path = "./data/bluesky/NAM/",
                         model = "NAM-3km") {

  # Each of the models seems stored a bit differently - adding one at a time as
  # needed
  if (model != "NAM-3km") {
    stop("this model instance isn't supported yet")
  }

  path <- "http://haze.airfire.org/bluesky-daily/output/standard/NAM84-0.15deg/"

  dates <- seq.Date(from = start, to = end, by = "1 day")

  get_one <- function(dt) {
    full_path <- paste0(path, strftime(dt, format = "%Y%m%d"),
                        "00/combined/data/smoke_dispersion.nc")
    output <- fs::path_join(c(output_path, paste0(strftime(dt, format = "%Y%m%d"),
                                                  "00zNAM3km.nc")))
    download.file(full_path, output, mode = "wb", quiet = TRUE)
  }

  purrr::walk(dates, get_one)

}

bluesky_at_airnow <- function(an, bluesky_path = "./data/bluesky/NAM",
                              debug = FALSE) {

  dates <- unique(an$Day)

  one_day <- function(dt, an) {

    filename <- paste0(strftime(dt, format = "%Y%m%d"), "00zNAM3km.nc")
    full_path <- fs::path_join(c(bluesky_path, filename))

    nc <- ncdf4::nc_open(full_path)
    g <- ncdf4::ncatt_get(nc, 0)
    pm <- ncdf4::ncvar_get(nc, "PM25")

    xmn <- g$XORIG
    ymn <- g$YORIG
    xmx <- xmn + g$XCELL * g$NCOLS
    ymx <- ymn + g$YCELL * g$NROWS

    pmb <- raster::brick(pm, xmn = xmn, ymn = ymn, xmx = xmx, ymx = ymx,
                         transpose = TRUE,  crs = "+proj=longlat +datum=WGS84")
    pmb <- raster::flip(pmb, "y")

    # Calculate the 24-hr avg
    pm24 <- raster::mean(pmb[[8:31]], na.rm = TRUE)
    if (debug) {

      rdf <- raster::as.data.frame(pm24, xy = TRUE)
      usa <- rgdal::readOGR("U:/SharedData/GIS", "states_lo_res")
      ca <- usa[usa$STATE_ABBR == "CA",]
      ca <- sp::spTransform(ca, sp::CRS("+proj=longlat +datum=WGS84"))
      g <- ggplot(rdf) +
        geom_raster(aes(x = x, y = y, fill = layer)) +
        geom_polygon(data = ca, aes(x = long, y = lat, group = group),
                     fill = NA, colour = "black") +
        scale_fill_viridis_c()
      g
    }

    i <- an$Day == dt
    locs = an[i, ]
    e <- raster::extract(pm24, locs)

    # The CMAQ output is not really lognormal, so probably should use the linear version
    df <- locs@data %>%
      mutate(PM25_bluesky = e,
             bluesky_log = log(PM25_bluesky))

  }

  purrr::map_dfr(dates, one_day, an)


}


# Take a raster stack of 24-hr avg PM2.5 and extract values at a set of
# locations on multiple days
preprocessed_bluesky_at_airnow <- function(stack, an) {

  dates <- unique(an$Day)
  layers <- seq(1, raster::nlayers(stack), by = 1)

  if (length(layers) != length(dates)) {
    stop("Stack and dates do not match!")
  }


  one_day <- function(dt, layer, an) {


    i <- an$Day == dt
    locs = an[i, ]
    e <- raster::extract(stack[[layer]], locs)

    # The CMAQ output is not really lognormal, so probably should use the linear version
    df <- locs@data %>%
      mutate(PM25_bluesky = e,
             bluesky_log = log(PM25_bluesky))

  }

  purrr::map2_dfr(dates, layers, one_day, an)

}

bluesky_at_grid <- function(start, end, grid, bluesky_path = "./data/bluesky/NAM") {

  dates <- seq.Date(from = start, to = end, by = "1 day")
  one_day <- function(dt, grid) {
    print(dt)
    filename <- paste0(strftime(dt, format = "%Y%m%d"), "00zNAM3km.nc")
    full_path <- fs::path_join(c(bluesky_path, filename))

    nc <- ncdf4::nc_open(full_path)
    g <- ncdf4::ncatt_get(nc, 0)
    pm <- ncdf4::ncvar_get(nc, "PM25")

    xmn <- g$XORIG
    ymn <- g$YORIG
    xmx <- xmn + g$XCELL * g$NCOLS
    ymx <- ymn + g$YCELL * g$NROWS

    pmb <- raster::brick(pm, xmn = xmn, ymn = ymn, xmx = xmx, ymx = ymx,
                         transpose = TRUE,  crs = "+proj=longlat +datum=WGS84")
    pmb <- raster::flip(pmb, "y")

    # Calculate the 24-hr avg
    pm24 <- raster::mean(pmb[[8:31]], na.rm = TRUE)

    e <- raster::extract(pm24, grid)

    # The CMAQ output is not really lognormal, so probably should use the linear version
    df <- grid@data %>%
      mutate(PM25_bluesky = e,
             Day = dt)

  }

  purrr::map_dfr(dates, one_day, grid)


}

# Take a raster stack of 24-hr avg PM2.5 and extract values at all grid
# locations on multiple days. Note that the raster stack must start on the same
# date as start
preprocessed_bluesky_at_grid <- function(start, end, stack, grid) {

  dates <- seq.Date(from = start, to = end, by = "1 day")

  one_day <- function(dt, layer, grid) {

    e <- raster::extract(stack[[layer]], grid)

    # The CMAQ output is not really lognormal, so probably should use the linear version
    df <- grid@data %>%
      mutate(PM25_bluesky = e,
             Day = dt)

  }

  layer <- 1:length(dates)
  purrr::map2_dfr(dates, layer, one_day, grid)

}
