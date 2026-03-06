# func/dl.R

# R code for running downloads

# TODO: Export as much of this to 0_download_data.py as possible


# Setup -------------------------------------------------------------------

# Load tidyverse
library(tidyverse)

# Set multi-core
library(doParallel); registerDoParallel(cores = detectCores()-2)

# Get internal functions, meta-data, etc.
source("func/util.R")

# Repo available here:
# https://github.com/RiOMar-projet/sat_access
source("~/sat_access/sat_access_script.R")

# AVISO+ credentials
aviso_plus_cred <- read.csv("~/pCloudDrive/Documents/info/aviso_plus_pswd.csv")

# Time steps to walk through
vec_time_steps <- c("daily")#, "weekly")#, "monthly")
# NB: Not enough HDD space for the weekly data
# NB: Monthly is too coarse for validation and TS analyses

# The variables we want in the two different ODATIS-MR correction types
vars_nirswir <- c("SPM", "CHL", "TUR", "CDOM", "RRS", "SST")
vars_polymer <- c("SPM", "CHL", "TUR", "CDOM", "RRS")


# Functions ---------------------------------------------------------------

# Convenience wrapper to help plyr::m_ply
download_nc_ply <- function(username, password, 
                            dl_product, dl_sensor, dl_correction, 
                            date_start, date_end, time_step, 
                            zone, lon_min, lon_max, lat_min, lat_max, 
                            dl_var){
  download_nc(
    dl_var = dl_var,
    dl_dates = c(date_start, date_end), dl_time_step = time_step,
    dl_product = dl_product, dl_sensor = dl_sensor, dl_correction = dl_correction, 
    dl_bbox = c(lon_min, lon_max, lat_min, lat_max),
    username = username, password = password,
    output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/", dl_sensor, zone, time_step),
    overwrite = FALSE)
}


# Download all SEXTANT ----------------------------------------------------

# Unnecessary as this is all run via the main python script
# Though I would like to remove the file structure complexity in the project generally...


# Download ODATIS-MR ------------------------------------------------------

## MODIS ------------------------------------------------------------------

# Create a stack of arguments for download_nc_ply that can be fed via plyr::m_ply()
modis_ply_df <- data.frame(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                           dl_product = "ODATIS-MR", dl_sensor = "MODIS", dl_correction = "nirswir", 
                           date_start = "2002-07-04", date_end = "2024-12-31", time_step = "daily", # NB: Eventually we may want the 8-day data as well
                           zone = rep(zones_bbox$zone, each = 6),
                           lon_min = rep(zones_bbox$lon_min, each = 6),
                           lon_max = rep(zones_bbox$lon_max, each = 6),
                           lat_min = rep(zones_bbox$lat_min, each = 6),
                           lat_max = rep(zones_bbox$lat_max, each = 6),
                           dl_var = vars_nirswir)

# Run them all
plyr::m_ply(.data = modis_ply_df[1:2,], .fun = download_nc_ply, .parallel = TRUE)


# NB: Realistically this can't all be run in one go
# It takes multiple days and will require caching multiple hundreds of Gigs of RAM
# But one can stop the process whenever necessary and relaunch without issue
for(i in 1:nrow(zones_bbox)){
  # Get the target zone
  zone <- zones_bbox[3:4,][i,]
  for(j in 1:length(vec_time_steps)){
    # Pick one time step
    time_step <- vec_time_steps[j]
    for(k in 1:length(vars_nirswir[1])){
      # Download the data
      download_nc(
        dl_var = vars_nirswir[k],
        dl_dates = c("2002-07-04", "2024-12-31"), dl_time_step = time_step,
        dl_product = "ODATIS-MR", dl_sensor = "MODIS", dl_correction = "nirswir", 
        dl_bbox = c(zone$lon_min, zone$lon_max, zone$lat_min, zone$lat_max),
        username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
        output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS", zone$zone, time_step),
        overwrite = FALSE
      ); gc()
    }
  }
}


## MERIS -------------------------------------------------------------------

for(i in 1:nrow(zones_bbox)){
  # Get the target zone
  zone <- zones_bbox[i,]
  for(j in 1:length(vec_time_steps)){
    # Pick one time step
    time_step <- vec_time_steps[j]
    for(k in 1:length(vars_polymer[1])){
      # Download the data
      download_nc(
        dl_var = vars_polymer[k],
        dl_dates = c("2002-06-19", "2012-04-08"), dl_time_step = time_step,
        dl_product = "ODATIS-MR", dl_sensor = "MERIS", dl_correction = "polymer", 
        dl_bbox = c(zone$lon_min, zone$lon_max, zone$lat_min, zone$lat_max),
        username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
        output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS", zone$zone, time_step),
        overwrite = FALSE
      ); gc()
    }
  }
}


## OLCI-A ------------------------------------------------------------------

for(i in 1:nrow(zones_bbox)){
  # Get the target zone
  zone <- zones_bbox[i,]
  for(j in 1:length(vec_time_steps)){
    # Pick one time step
    time_step <- vec_time_steps[j]
    for(k in 1:length(vars_polymer[1])){
      # Download the data
      download_nc(
        dl_var = vars_polymer[k],
        dl_dates = c("2016-04-26", "2024-12-31"), dl_time_step = time_step,
        dl_product = "ODATIS-MR", dl_sensor = "OLCI-A", dl_correction = "polymer", 
        dl_bbox = c(zone$lon_min, zone$lon_max, zone$lat_min, zone$lat_max),
        username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
        output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-A", zone$zone, time_step),
        overwrite = FALSE
      ); gc()
    }
  }
}


## OLCI-B ------------------------------------------------------------------

for(i in 1:nrow(zones_bbox)){
  # Get the target zone
  zone <- zones_bbox[i,]
  for(j in 1:length(vec_time_steps)){
    # Pick one time step
    time_step <- vec_time_steps[j]
    for(k in 1:length(vars_polymer[1])){
      # Download the data
      download_nc(
        dl_var = vars_polymer[k],
        dl_dates = c("2018-05-15", "2024-12-31"), dl_time_step = time_step,
        dl_product = "ODATIS-MR", dl_sensor = "OLCI-B", dl_correction = "polymer", 
        dl_bbox = c(zone$lon_min, zone$lon_max, zone$lat_min, zone$lat_max),
        username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
        output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-B", zone$zone, time_step),
        overwrite = FALSE
      ); gc()
    }
  }
}

