# Tools to create an output grid for predection and pretty pictures. We can also
# just do the predictions at participant locations, but that would require
# running kriging and extraction every time a location changed.

# Make a grid for output predictions
make_grid <- function(spdf, cellsize) {
  box <- spdf@bbox
  xrange <- box[3] - box[1]
  yrange <- box[4] - box[2]
  xcells <- ceiling(xrange / cellsize)
  ycells <- ceiling(yrange / cellsize)

  x <- seq(1, xcells, by = 1) * cellsize + box[1]
  y <- seq(1, ycells, by = 1) * cellsize + box[2]

  df <- data.frame(Id = seq(1, length(x) * length(y), by = 1))

  grid <- sp::SpatialPointsDataFrame(cbind(rep(x, length(y)), rep(y, each = length(x))),
                        proj4string = sp::CRS("+init=epsg:3395"), data = df)
  sp::gridded(grid) <- TRUE
  return(grid)
}


