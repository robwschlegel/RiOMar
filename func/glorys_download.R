#' Download GLORYS data from Copernicus Marine Data Store
#' 
#' This script uses the CopernicusMarine R package to download GLORYS (Global Ocean Reanalysis and Simulation) data.
#' GLORYS provides global ocean reanalysis products including temperature, salinity, currents, and sea level.
#' 
#' Requirements:
#' - CopernicusMarine R package installed
#' - Valid Copernicus Marine account credentials
#' - Internet connection
#' 
#' Usage:
#' 1. Install the CopernicusMarine package if not already installed
#' 2. Set up your Copernicus Marine credentials
#' 3. Modify the parameters below to specify your data requirements
#' 4. Run the script
#' 
#' NB: This was an initial test of the Mistral vibe coding agent running locally on my machine
#' This script was created in one go via a single command line request written in human language

# Load required library
if (!require("CopernicusMarine")) {
  install.packages("CopernicusMarine")
  library(CopernicusMarine)
}

# Set your Copernicus Marine credentials
# You can get these by registering at https://marine.copernicus.eu/
# It's recommended to store credentials securely, not in the script
Sys.setenv(COPERNICUS_MARINE_USERNAME = "rschlegel1")
Sys.setenv(COPERNICUS_MARINE_PASSWORD = "RobertCMEMS2018")

# GLORYS dataset parameters - modify these according to your needs
# Dataset ID for GLORYS global reanalysis (example for physics)
dataset_id <- "GLOBAL_MULTIYEAR_PHY_001_030"
layer_id <- "cmems_mod_glo_phy_my_0.083deg_P1D-m"

# Time range for data download (format: "YYYY-MM-DD")
start_date <- "2023-01-01"
end_date <- "2023-01-07"

# Geographic bounding box (min_longitude, max_longitude, min_latitude, max_latitude)
# Example: Bay of Biscay region
longitude_range <- c(4, 8)    # Western to Eastern longitude
latitude_range <- c(42, 44)     # Southern to Northern latitude
range_lonlat <- c(longitude_range[1], latitude_range[1], longitude_range[2], latitude_range[2]) # (min_lon, min_lat, max_lon, max_lat1])

# Depth range (in meters)
depth_range <- c(0, 0)        # Surface

# Variables to download (available variables depend on the dataset)
# Common GLORYS variables include:
# "thetao" - Potential temperature
# "so" - Salinity
# "uo" - Eastward sea water velocity
# "vo" - Northward sea water velocity
# "zos" - Sea surface height above geoid
# "mlotst" - Mixed layer depth
variables <- c("thetao", "so", "uo", "vo")

# Output directory and filename
# download_dir <- "/home/calanus/RiOMar/data/GLORYS"
# output_filename <- paste("glorys_", start_date, "_to_", end_date, ".nc", sep = "")

# Create download directory if it doesn't exist
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Download GLORYS data
# Set up the download request
download_request <- CopernicusMarine::cms_download_subset(
  product = dataset_id,
  layer = layer_id,
  variable = variables,
  region = range_lonlat,
  time = c(start_date, end_date),
  verticalrange = depth_range,
  progress = TRUE
)

# Plot
plot(download_request["thetao", drop = TRUE], col = hcl.colors(100), axes = TRUE)

# Additional information about GLORYS datasets:
# GLORYS12V1: global-reanalysis-phys-001-030 (1993-present, 1/12° resolution)
# GLORYS2V4: global-reanalysis-phys-001-025 (1993-2019, 1/4° resolution)
# GLORYS2V3: global-reanalysis-phys-001-012 (1993-2015, 1/4° resolution)

# To find other available datasets, visit:
# https://data.marine.copernicus.eu/products

# Note: Large downloads may take significant time and bandwidth
# Consider downloading smaller time periods or regions if you encounter issues

