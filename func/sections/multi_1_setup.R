# func/multi.R
# Loads all drivers of plume size, performs stats, plots results


# Analysis ideas ----------------------------------------------------------

# Treat each pixel like its own time series and see what is happening with the forces when the pixel is triggered as a panache
## also how high SPM is while all this is happening
## show primary wind direction when pixel is triggered
## also relationship with SPM and tide range or category
## number of times pixel is flagged related to the size of the total panache when it is flagged
### Would need to relate wind with time lag to this as well
## could also tally the shape of the panache whenever pixel is flagged

# GIFS
## Animate mean offshore distance of centroid
## Wuld be interesting to have the centroid visualised as a 21 dot
### it could leave a trail of 1 dots behind it per day
### or fill colour could be left to show tide, wind, etc. on that day

# More ideas
## nmds of mean characteristics during plume events
## extreme event analysis
### need to establish a reasonable baseline, and go from there
### could be interesting to use the time varying X11 seasonality
### but then is it relevant to calculate event stats from the median signal?
  ### one could modify the analysis to always take from the base, being zero
### how would one establish the 90th percentile? Or rather, just always take everything over the seasonal signal
## EMD - see Vincent email
### percent contribution of each component to time series
### changes over time as well
## like a CTD cast, figure out a way to measure from when a plume peaks and then goes down below a certain threshold as a way of determining if it is an individual event
### and then get statistics from that
### also add up number of days with onshore wind, neap tide, etc.
### the spatial threshold could be based on the percentile of the total plume area over the full time series
### start searching by creating a contour plot of every 10th percentile
## also account for the lon/lat of the centroid and how that relates to the other drivers
# be able to say if the plume is increasing or decreasing in size so that the drivers on that day can be categorised under what the plume is doing

## Ideas for driver decomposition
# A SOM analysis of the spatial footprint during given drivers may work.
## Average the primary drivers over the timespan of the given large plume
## Then see if the SOM organises them by surface signature. I.e. size, centroid, shape, direction, duration
# Or PCA/ordination of days of panache given certain values for variables
# Use rle() to determine contiguous events temporarily

## Temperature analysis
# Panache should be colder than coastal water
# It should be possible to a priori define the colder surface temperatures based on the panache pixels
# Use the ODATIS-MR SST for this


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ncdf4)
library(seasonal) # For X11 analysis (currently not used)
library(patchwork)
library(sandwich) # For HAC covariance tests (driver_plume_trend)
library(lmtest) # For more detailed linear model tests (driver_plume_trend)

# Common function
source("func/util.R")

# Zones, north to south (see func/util.R::ZONE_ORDER/order_zones())
zones <- ZONE_ORDER


# Zone / gauge metadata -----------------------------------------------------

# Canonical river mouth -> zone -> tide gauge lookup. Replaces the identical
# if/else block that was previously copy-pasted in flow_comp(), flow_trend(),
# flow_plume_trend_plus() (all three in the old flow.R), tide_calc() (tide.R),
# spatial_wind_calc() (wind.R), and surface_plot() (surface.R).
zone_meta <- river_mouths |>
  dplyr::mutate(
    zone = dplyr::case_when(
      mouth_name == "Seine"       ~ "BAY_OF_SEINE",
      mouth_name == "Gironde"     ~ "BAY_OF_BISCAY",
      mouth_name == "Loire"       ~ "SOUTHERN_BRITTANY",
      mouth_name == "Grand Rhone" ~ "GULF_OF_LION",
      TRUE ~ NA_character_
    ),
    gauge = dplyr::case_when(
      zone == "BAY_OF_SEINE"      ~ "LE_HAVRE",
      zone == "BAY_OF_BISCAY"     ~ "PORT-BLOC",
      zone == "SOUTHERN_BRITTANY" ~ "SAINT-NAZAIRE",
      zone == "GULF_OF_LION"      ~ "MARSEILLE",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::arrange(match(zone, ZONE_ORDER))

# Look up zone metadata by either the river mouth name (as used in
# river_mouths/zone_meta) or the zone code (as used in zones_bbox / output/
# paths). Exactly one of mouth_name / zone_name should be supplied.
# get_zone_meta(mouth_name = "Seine")
# get_zone_meta(zone_name = "GULF_OF_LION")
get_zone_meta <- function(mouth_name = NULL, zone_name = NULL){
  if(!is.null(mouth_name)){
    out <- dplyr::filter(zone_meta, mouth_name == !!mouth_name)
  } else if(!is.null(zone_name)){
    out <- dplyr::filter(zone_meta, zone == !!zone_name)
  } else {
    stop("Supply either mouth_name or zone_name to get_zone_meta().")
  }
  if(nrow(out) != 1) stop("Zone/mouth not recognised in zone_meta.")
  return(out)
}

# Human-facing zone display name (e.g. "BAY_OF_SEINE" -> "Bay of Seine")
zone_display_names <- c(
  BAY_OF_SEINE      = "Bay of Seine",
  BAY_OF_BISCAY     = "Gironde shelf",
  SOUTHERN_BRITTANY = "Southern Brittany",
  GULF_OF_LION      = "Rhône shelf"
)
zone_title <- function(zone_name){
  zone_name <- as.character(zone_name)
  out <- unname(zone_display_names[zone_name])
  if(anyNA(out) && !anyNA(zone_name)) stop("zone_title(): unrecognised zone code(s): ",
                                            paste(setdiff(zone_name, names(zone_display_names)), collapse = ", "))
  out
}


