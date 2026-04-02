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

# Create FRANCE bounding box with same structure as zones_bbox
france_bbox <- data.frame(zone = "FRANCE",
                         lon_min = c(-7.8),
                         lon_max = c(10.3),
                         lat_min  = c(41.2),
                         lat_max = c(51.5)) 

# Credentials
aviso_plus_cred <- read.csv("~/pCloudDrive/Documents/info/aviso_plus_pswd.csv")
odatis_mr_expert_cred <- read.csv("~/pCloudDrive/Documents/info/odatis_mr_expert_pswrd.csv")

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
                            dl_product, dl_sensor, dl_correction, dl_processing,
                            time_step, zone, lon_min, lon_max, lat_min, lat_max,
                            dl_var, date_start, date_end){
  download_nc(
    dl_var = dl_var,
    dl_dates = c(date_start, date_end), dl_time_step = time_step,
    dl_product = dl_product, dl_sensor = dl_sensor, 
    dl_correction = dl_correction, dl_processing = dl_processing,
    dl_bbox = c(lon_min, lon_max, lat_min, lat_max),
    username = username, password = password,
    # output_dir = file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/", dl_sensor, zone, time_step),
    output_dir = file.path("~/data/ODATIS-MR/", dl_sensor, zone, time_step),
    overwrite = FALSE)
}

# Another wrapper that calls this one... probably a better way to do this
# NB: It is assumed that all arguments are one value except zone and dl_var
# username = aviso_plus_cred$usrname; password = aviso_plus_cred$psswrd 
# dl_product = "ODATIS-MR"; dl_sensor = "MODIS"; dl_correction = "nirswir" 
# date_start = "2002-07-04"; date_end = "2024-12-31"; time_step = "daily"
# zone_info = zones_bbox[4,]; dl_var = vars_nirswir
download_study_area <- function(username, password, 
                                dl_product, dl_sensor, dl_correction, dl_processing = NULL,
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
                       dl_product = dl_product, dl_sensor = dl_sensor, 
                       dl_correction = dl_correction, dl_processing = dl_processing,
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

# A wrapper to download all of the EXPERT data for a given sensor
# NB: This is only designed to directly download specific daily EXPERT files
# username = odatis_mr_expert_cred$usrname; password = odatis_mr_expert_cred$psswrd
# dl_sensor = "MODIS"; date_start = "2002-07-04"; date_end = "2024-12-31"
download_expert <- function(username, password, 
                            dl_sensor, date_start, date_end){
  
  # Get a dataframe of download dates by year
  year_range <- year(date_start):year(date_end)
  date_df <- data.frame(date_start = paste0(year_range,"-01-01"),
                        date_end = paste0(year_range,"-12-31"))
  
  # Correct start and end dates accordingly
  date_df$date_start[1] <- ifelse(date_df$date_start[1] < date_start, date_start, date_df$date_start[1])
  date_df$date_end[nrow(date_df)] <- ifelse(date_df$date_end[nrow(date_df)] > date_end, date_end, date_df$date_end[nrow(date_df)])
  
  # Set the full range of possible values
  dl_product <- "ODATIS-MR EXPERT"
  dl_var <- c("CHL", "CHL1", "SPM", "SST", "SST-NIGHT", "T")
  dl_correction <- c("acolite", "nirswir", "polymer")
  dl_processing <- c("OC5", "GONS", "G", "R")

  # Create big grid and clean as necessary
  ply_df <- expand_grid(dl_product, dl_sensor, dl_correction, dl_processing, dl_var) |> 
    mutate(dl_var = case_when(dl_var == "CHL" & dl_processing %in% c("R", "G") ~ NA,
                              dl_var == "SPM" & dl_processing %in% c("OC5", "GONS") ~ NA,
                              # Remove unused corrections by sensor
                              dl_sensor == "MODIS" & dl_correction == "acolite" ~ NA,
                              dl_sensor != "MODIS" & dl_correction == "nirswir" ~ NA,
                              # Only MODIS has SST
                              dl_var %in% c("SST", "SST-NIGHT") & dl_sensor != "MODIS" ~ NA,
                              # No SST-NIGHT for MODIS : polymer
                              dl_var == "SST-NIGHT" & dl_correction == "polymer" ~ NA,
                              # No GONS for MODIS
                              dl_sensor == "MODIS" & dl_processing == "GONS" ~ NA,
                              TRUE ~ dl_var)) |> 
    mutate(dl_processing = case_when(dl_var %in% c("SST", "SST-NIGHT", "CHL1", "T") ~ NA,
                                    TRUE ~ dl_processing)) |> 
    filter(!is.na(dl_var)) |> 
    distinct()
  
  # Bigger grid
  ply_date_df <- expand_grid(ply_df, date_df) |> 
    # Add France bbox, even though it's not used, to keep same structure as other functions
    bind_cols(france_bbox) |> 
    mutate(time_step = "daily",
      username = odatis_mr_expert_cred$usrname,
      password = odatis_mr_expert_cred$psswrd) |> 
    dplyr::select(username, password, dl_product, dl_sensor, dl_correction, dl_processing,
                  time_step, zone, lon_min, lon_max, lat_min, lat_max, dl_var, date_start, date_end)

  # Ply it
  # NB: It appears that the FTP server doesn't accept parallel requests
  plyr::m_ply(.data = ply_date_df, .fun = download_nc_ply, .parallel = FALSE); gc()
}

# Download all SEXTANT ----------------------------------------------------

# Unnecessary as this is all run via the main python script
# Though I would like to remove the file structure complexity in the project generally...


# Download ODATIS-MR ------------------------------------------------------

## MODIS ------------------------------------------------------------------
# 8217 days of data

# Run a loop across all sites
# NB: The multi-core is targeted at running one core for each year of data in the subset by variable
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "MODIS", dl_correction = "nirswir",
                      date_start = "2002-07-04", date_end = "2024-12-31", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_nirswir)
}

# EXPERT products
download_expert("MODIS")


## MERIS -------------------------------------------------------------------
# 3582 days of data

# Run a loop across all sites
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "MERIS", dl_correction = "polymer",
                      date_start = "2002-06-19", date_end = "2012-04-08", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_polymer)
}

# EXPERT products
download_expert("MERIS")


## OLCI-A ------------------------------------------------------------------
# 3172 days of data

# Run a loop across all sites
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "OLCI-A", dl_correction = "polymer",
                      date_start = "2016-04-26", date_end = "2024-12-31", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_polymer)
}

# EXPERT products
download_expert("OLCI-A")


## OLCI-B ------------------------------------------------------------------
# 2423 days of data

# Run a loop across all sites
for(i in 1:nrow(zones_bbox)){
  download_study_area(username = aviso_plus_cred$usrname, password = aviso_plus_cred$psswrd,
                      dl_product = "ODATIS-MR", dl_sensor = "OLCI-B", dl_correction = "polymer",
                      date_start = "2018-05-15", date_end = "2024-12-31", time_step = "daily",
                      zone_info = zones_bbox[i,], dl_var = vars_polymer)
}

# EXPERT products
download_expert("OLCI-B")

