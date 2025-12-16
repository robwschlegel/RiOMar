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
Sys.setenv(COPERNICUS_MARINE_USERNAME = "your_username")
Sys.setenv(COPERNICUS_MARINE_PASSWORD = "your_password")

# GLORYS dataset parameters - modify these according to your needs
# Dataset ID for GLORYS global reanalysis (example for physics)
dataset_id <- "global-reanalysis-phys-001-030"

# Time range for data download (format: "YYYY-MM-DD")
start_date <- "2023-01-01"
end_date <- "2023-01-31"

# Geographic bounding box (min_longitude, max_longitude, min_latitude, max_latitude)
# Example: Bay of Biscay region
longitude_range <- c(-10, 0)    # Western to Eastern longitude
latitude_range <- c(40, 50)     # Southern to Northern latitude

# Depth range (in meters)
depth_range <- c(0, 100)        # Surface to 100m depth

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
download_dir <- "/home/calanus/RiOMar/data/GLORYS"
output_filename <- "glorys_data_202301.nc"

# Create download directory if it doesn't exist
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}

# Download GLORYS data
tryCatch({
  # Set up the download request
  download_request <- cm_download(
    dataset_id = dataset_id,
    username = Sys.getenv("COPERNICUS_MARINE_USERNAME"),
    password = Sys.getenv("COPERNICUS_MARINE_PASSWORD"),
    minimum_longitude = longitude_range[1],
    maximum_longitude = longitude_range[2],
    minimum_latitude = latitude_range[1],
    maximum_latitude = latitude_range[2],
    start_datetime = start_date,
    end_datetime = end_date,
    minimum_depth = depth_range[1],
    maximum_depth = depth_range[2],
    force_download = TRUE,
    output_filename = file.path(download_dir, output_filename),
    variables = variables
  )

  # Execute the download
  print(paste("Starting download of GLORYS data for", start_date, "to", end_date))
  print(paste("Region:", longitude_range[1], "to", longitude_range[2], "longitude,", 
              latitude_range[1], "to", latitude_range[2], "latitude"))
  print(paste("Depth range:", depth_range[1], "to", depth_range[2], "meters"))
  print(paste("Variables:", paste(variables, collapse = ", ")))

  # This will download the data and save it to the specified location
  result <- download_request

  print(paste("Download completed successfully!"))
  print(paste("Data saved to:", file.path(download_dir, output_filename)))

}, error = function(e) {
  print(paste("Error downloading GLORYS data:", e))
  # Print more detailed error information
  print("Please check:")
  print("1. Your Copernicus Marine credentials are correct")
  print("2. The dataset ID is valid and available")
  print("3. Your internet connection is working")
  print("4. The date range and geographic area are valid")
})

# Additional information about GLORYS datasets:
# GLORYS12V1: global-reanalysis-phys-001-030 (1993-present, 1/12° resolution)
# GLORYS2V4: global-reanalysis-phys-001-025 (1993-2019, 1/4° resolution)
# GLORYS2V3: global-reanalysis-phys-001-012 (1993-2015, 1/4° resolution)

# To find other available datasets, visit:
# https://data.marine.copernicus.eu/products

# Note: Large downloads may take significant time and bandwidth
# Consider downloading smaller time periods or regions if you encounter issues

