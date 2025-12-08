# func/sentinel.R
# Extract data from sentinel 2/3 NetCDF files and perform basic comparisons


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ncdf4)


# Sentinel2 ---------------------------------------------------------------

# Function that will load one wavelength (nm) from a Sentinel-2 NetCDF file
load_S2_var <- function(file_name, var_choice = NULL, display_vars = FALSE){
  
  # Open NetCDF file connection
  nc_dat <- nc_open(file_name)
  
  # Display variables if desired
  if(display_vars) {
    return(print(names(nc_dat$var)))
  }
  
  # Load the desired variables
  if(!is.null(var_choice)) {
    
    # Extract data
    nc_var <- as.vector(ncvar_get(nc_dat, var_choice))
    nc_lon <- as.vector(ncvar_get(nc_dat, "LON"))
    nc_lat <- as.vector(ncvar_get(nc_dat, "LAT"))
    nc_time <- as.POSIXct(gsub("T", " ", ncatt_get(nc_dat, varid = 0)$ISOTIME), tz = "UTC")
    nc_close(nc_dat)
    
    # Combine, prep, and exit
    nc_df <- data.frame(time = nc_time, lon = nc_lon, lat = nc_lat, var = nc_var)
    colnames(nc_df)[4] <- var_choice
    return(nc_df)
  }
  
  # Else exit
  nc_close(nc_dat)
  message("No variable chosen for extraction. Set a value for `var_choice`.")
  return()
}


## Example use -------------------------------------------------------------

# Check the variables within a S2 NetCDF
load_S2_var("~/Downloads/sentinel_2_3/S2A/S2A_MSIL1C_20170206T102211_N0204_R065_T32TLP_20170206T102733-2_stitch_L2.nc", display_vars = TRUE)

# Extract RHOW 497 data
S2_497 <- load_S2_var("~/Downloads/sentinel_2_3/S2A/S2A_MSIL1C_20170206T102211_N0204_R065_T32TLP_20170206T102733-2_stitch_L2.nc", "RHOW_497")

# Extract RHOW 444 data
S2_444 <- load_S2_var("~/Downloads/sentinel_2_3/S2A/S2A_MSIL1C_20170206T102211_N0204_R065_T32TLP_20170206T102733-2_stitch_L2.nc", "RHOW_444")

# Combine into one dataframe
S2_nm <- left_join(S2_444, S2_497, by = c("time", "lon", "lat"))

# Pivot longer for further analyses and/or plotting
S2_nm_long <- pivot_longer(S2_nm, cols = RHOW_444:RHOW_497, names_to = "nm", values_to = "value")

# Plot values as histogram
ggplot(data = S2_nm_long) +
  geom_histogram(aes(x = value)) +
  facet_wrap(~nm)


# Sentinel3 ---------------------------------------------------------------

# Function that will load Chl a data from Sentinel3 files
## NB: This simply takes the start time of the sample as the single time value
load_S3_var <- function(file_name){
  
  # Search for the lon/lat file
  coord_file <- file.path(dirname(file_name), "geo_coordinates.nc")
  if(file.exists(coord_file)) {
    coord_dat <- nc_open(coord_file)
    nc_lon <- as.vector(ncvar_get(coord_dat, "longitude"))
    nc_lat <- as.vector(ncvar_get(coord_dat, "latitude"))
    nc_close(coord_dat)
  } else {
    stop(paste0("Cannot find coordinates file. Was looking for: ", coord_file))
  }
  
  # Open NetCDF file connection
  nc_dat <- nc_open(file_name)
  
  # Determine which variable to load
  if("CHL_NN" %in% names(nc_dat$var)) {
    nc_var <- as.vector(ncvar_get(nc_dat, "CHL_NN"))
    nc_time <- as.POSIXct(gsub("T", " ", ncatt_get(nc_dat, varid = 0)$start_time), tz = "UTC")
    var_name <- "CHL_NN"
    nc_close(nc_dat)
  } else if("CHL_OC4ME" %in% names(nc_dat$var)) {
    nc_var <- as.vector(ncvar_get(nc_dat, "CHL_OC4ME"))
    nc_time <- as.POSIXct(gsub("T", " ", ncatt_get(nc_dat, varid = 0)$start_time), tz = "UTC")
    var_name <- "CHL_OC4ME"
    nc_close(nc_dat)
  } else {
    stop(paste0("Cannot find a valid Chl a value in the file. Options are 'CHL_NN' or 'CHL_OC4ME'."))
  }
  
  # Combine, prep, and exit
  nc_df <- data.frame(time = nc_time, lon = nc_lon, lat = nc_lat, var = nc_var)
  colnames(nc_df)[4] <- var_name
  return(nc_df)
}


# Example use -------------------------------------------------------------

# Load two different chl files
S3_CHLA_NN <- load_S3_var("~/Downloads/sentinel_2_3/S3B/chl_nn.nc")
S3_CHLA_OC4 <- load_S3_var("~/Downloads/sentinel_2_3/S3B/chl_oc4me.nc")

# Combine into one data frame
## NB: Removing rows with missing values to speed things along
S3_CHLA <- left_join(distinct(filter(S3_CHLA_NN, !is.na(CHL_NN))),
                     distinct(filter(S3_CHLA_OC4, !is.na(CHL_OC4ME))), 
                     by = c("time", "lon", "lat"), 
                     relationship = "many-to-many") |> 
  distinct()

# Calculate difference between Chl a values
S3_CHLA_diff <- S3_CHLA |> 
  mutate(chla_diff = CHL_NN - CHL_OC4ME)

# Plot differences as a histogram
ggplot(data = S3_CHLA_diff) +
  geom_histogram(aes(x = chla_diff))

# NB: These large differences may be caused by the correction not being applied correctly by the NCDF4 package

