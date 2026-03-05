# func/dl.R

# R code for running downloads

# TODO: Export as much of this to 0_download_data.py as possible


# Setup -------------------------------------------------------------------

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
vars_nirswir <- c("CDOM", "CHL", "SPM", "TUR", "RRS", "SST")
vars_polymer <- c("CDOM", "CHL", "SPM", "TUR", "RRS")


# Download all SEXTANT ----------------------------------------------------

# Unnecessary as this is all run via the main python script
# Though I would like to remove the file structure complexity in the project generally...


# Download ODATIS-MR ------------------------------------------------------

## MODIS ------------------------------------------------------------------

# NB: Realistically this can't all be run in one go
# It takes multiple days and will require caching multiple hndreds of Gigs of RAM
# But one can stop the process whenever necessary and relaunch without issue
for(i in 1:nrow(zones_bbox)){
  # Get the target zone
  zone <- zones_bbox[i,]
  for(j in 1:length(vec_time_steps)){
    # Pick one time step
    time_step <- vec_time_steps[j]
    for(k in 1:length(vars_nirswir)){
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
    for(k in 1:length(vars_polymer)){
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
    for(k in 1:length(vars_polymer)){
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
    for(k in 1:length(vars_polymer)){
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

