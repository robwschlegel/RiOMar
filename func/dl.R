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
                            time_step, zone, lon_min, lon_max, lat_min, lat_max,
                            dl_var, date_start, date_end){
  download_nc(
    dl_var = dl_var,
    dl_dates = c(date_start, date_end), dl_time_step = time_step,
    dl_product = dl_product, dl_sensor = dl_sensor, dl_correction = dl_correction, 
    dl_bbox = c(lon_min, lon_max, lat_min, lat_max),
    username = username, password = password,
    output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/", dl_sensor, zone, time_step),
    overwrite = FALSE)
}

# Another wrapper that calls this one... probably a better way to do this
# NB: It is assumed that all arguments are one value except zone and dl_var
# username = aviso_plus_cred$usrname; password = aviso_plus_cred$psswrd 
# dl_product = "ODATIS-MR"; dl_sensor = "MODIS"; dl_correction = "nirswir" 
# date_start = "2002-07-04"; date_end = "2024-12-31"; time_step = "daily"
# zone_info = zones_bbox[4,]; dl_var = vars_nirswir
download_study_area <- function(username, password, 
                                dl_product, dl_sensor, dl_correction, 
                                date_start, date_end, time_step, 
                                zone_info, dl_var){
  
  # Get a dataframe of download dates by year
  year_range <- year(date_start):year(date_end)
  date_df <- data.frame(date_start = paste0(year_range,"-01-01"),
                        date_end = paste0(year_range,"-12-31"))
  
  # Correct start and end dates accordingly
  date_df$date_start[1] <- ifelse(date_df$date_start[1] < date_start, date_start, date_df$date_start[1])
  date_df$date_end[nrow(date_df)] <- ifelse(date_df$date_end[nrow(date_df)] > date_end, date_end, date_df$date_end[nrow(date_df)])
  
  # Create the full ply dataframe
  ply_df <- data.frame(username = username, password = password,
                       dl_product = dl_product, dl_sensor = dl_sensor, dl_correction = dl_correction,
                       # date_start = date_start, date_end = date_end, 
                       time_step = time_step,
                       zone = rep(zone_info$zone, each = length(dl_var)),
                       lon_min = rep(zone_info$lon_min, each = length(dl_var)),
                       lon_max = rep(zone_info$lon_max, each = length(dl_var)),
                       lat_min = rep(zone_info$lat_min, each = length(dl_var)),
                       lat_max = rep(zone_info$lat_max, each = length(dl_var)),
                       dl_var = dl_var)
  
  # Bigger grid
  ply_date_df <- expand_grid(ply_df, date_df)
  
  # Ply it
  plyr::m_ply(.data = ply_date_df, .fun = download_nc_ply, .parallel = TRUE); gc()
}


# Download all SEXTANT ----------------------------------------------------

# Unnecessary as this is all run via the main python script
# Though I would like to remove the file structure complexity in the project generally...


# Download ODATIS-MR ------------------------------------------------------

## MODIS ------------------------------------------------------------------

# Run a loop across all sites
# NB: The multi-core is targeted at running one core for each year of data in the subset by variable
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "MODIS", dl_correction = "nirswir",
                      date_start = "2002-07-04", date_end = "2024-12-31", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_nirswir[1]) # NB: Just getting SPM for the moment
}


## MERIS -------------------------------------------------------------------

# Run a loop across all sites
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "MERIS", dl_correction = "polymer",
                      date_start = "2002-06-19", date_end = "2012-04-08", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_polymer[1]) # NB: Just getting SPM for the moment
}


## OLCI-A ------------------------------------------------------------------

# Run a loop across all sites
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "OLCI-A", dl_correction = "polymer",
                      date_start = "2016-04-26", date_end = "2024-12-31", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_polymer[1]) # NB: Just getting SPM for the moment
}


## OLCI-B ------------------------------------------------------------------

# Run a loop across all sites
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "OLCI-B", dl_correction = "polymer",
                      date_start = "2018-05-15", date_end = "2024-12-31", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_polymer[1]) # NB: Just getting SPM for the moment
}

