# func/wind.R
# Comparisons of wind forcing against plume size


# Libraries ---------------------------------------------------------------

library(tidyverse)


# Meta-data ---------------------------------------------------------------

river_mouths <- data.frame()


# Functions ---------------------------------------------------------------

spatial_wind_calc <- function(river_mouth){
  
  # Start with lon/lat coords of river mouth
  
  # Load and subset wind data to + 0.5 N and E from mouth
  
  # Calculate spatial average of wind vectors by day
  
  # Calculate wind speed and direction
  
  # Determine simple upwelling/downwelling index based on coastal direction (based on river mouth name)
  
  # Load panache time series based on river mouth name
  
  # Compare panache size against wind speed, against direction
  
  # Two separate comparisons based on upwelling or downwelling times
  
  # Probably in a separate analysis the tides (neap <-> spring) will be folded in
  
}
