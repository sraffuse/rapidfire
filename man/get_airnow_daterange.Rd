% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monitors.R
\name{get_airnow_daterange}
\alias{get_airnow_daterange}
\title{get_airnow_daterange}
\usage{
get_airnow_daterange(start, end, states)
}
\arguments{
\item{start}{Date The earliest date to search}

\item{end}{Date The latest date to search. Must be in the same calendar year
as \code{start}.}

\item{states}{character The states to include as a vector of two-character
state abbreviations}
}
\value{
AirNow data for the specified dates and locations as SpatialPointsDataFrame
}
\description{
For the specified date range and US states, downloads AirNow data.
}
\examples{
dt1 <- as.Date("2018-11-01")
  dt2 <- as.Date("2018-11-30")
  an_ws <- get_airnow_daterange(dt1, dt2, c("CA", "NV"))
}
