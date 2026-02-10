# earth_data_access
# This script provides the functions, workflow, and examples 
# to download any data product from the earth data server:
# https://www.earthdata.nasa.gov/

# The user guide for MODIS products may be found here :
# https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf


# Libraries ---------------------------------------------------------------

# Luna package is used to access large spatial data products
# NB: It is not available on r-universe not CRAN
# if (!"luna" %in% installed.packages()) {
#   install.packages('luna', repos = 'https://rspatial.r-universe.dev')
# }

# Check for missing libraries and install them if necessary
# We will use th plyr package, but we will not load it explicitly
# if (!all(c("maps", "tidyverse", "ncdf4", "terra", "doParallel", "plyr") %in% installed.packages())) {
#   install.packages(c("maps", "tidyverse", "ncdf4", "terra", "doParallel", "plyr"), repos = "https://cloud.r-project.org/")
# }

library(tidyverse)
library(ncdf4)
library(terra)
library(luna)
library(doParallel); registerDoParallel(cores = detectCores() - 2)


# Setup -------------------------------------------------------------------

# Load username and password
# This is a two column csv file with column names: usrname, psswrd
# It contains one row of data, which is the username and password for an account on earth data
# An account can be created here: https://urs.earthdata.nasa.gov/
# Once you have your account details, create the .csv file shown here and store in a secure location
# earth_up <- read_csv("path/to/file/earthdata_pswd.csv")
earth_up <- read_csv("~/pCloudDrive/Documents/info/earthdata_pswd.csv")


# Functions ---------------------------------------------------------------

# Process MODIS data in a batch
# Doesn't work well with L1 products or very large study areas
MODIS_proc <- function(file_name_df, bbox, band_num, out_dir, land_mask = FALSE){
  
  # Get date from filenames
  files_date <- luna::modisDate(file_name_df$filename)
  
  # Check that only one date of data is being processed
  if(length(unique(files_date$date)) > 1 ) stop("Please ensure only one date of data is being processed.")
  
  # get product ID from file_names
  file_product <- strsplit(file_name_df$filename, "\\/")
  file_product <- sapply(file_product, "[[", length(file_product[[1]]))
  file_product <- unique(sapply(strsplit(file_product, "\\."), "[[", 1))
  
  # Check that only one product type is being processed
  if(length(file_product) > 1 ) stop("Please ensure only one product ID is being processed.")
  
  # Check if file already exists and skip if so
  file_name_out <- file.path(out_dir, paste0("study_area_",file_product,"_",unique(files_date$date),".tif"))
  if(file.exists(file_name_out)){
    message(file_name_out, " already exists. Delete it if you want to reprocess the data. Otherwise, all good.")
    return(NULL)
  }
  
  # Load  and merge the files with desired layers etc.
  if(grepl("MYD01|MYD02", file_product)){
    
    stop("MODIS L1 data are not currently working with this processing workflow. Rather use a different software.")
    
  } else {
    
    # TODO: Get this to detect if it should use 'subds' or 'lyrs'
    data_layer <- rast(file_name_df$filename[1], lyrs = band_num)
    # data_layer
    # plot(data_layer)
    
    # Project to EPSG:4326 and crop
    # data_rectify <- rectify(data_layers[[1]])
    data_proj <- project(data_layer, y = "EPSG:4326")
    # plot(data_proj)
    data_crop <- crop(data_proj, bbox)
    # plot(data_crop)
    
    if(length(file_name_df$filename) > 1){
      # Run and merge each individual file
      data_step <- 2
      while(data_step <= length(file_name_df$filename)){
        data_proj_i <- project(rast(file_name_df$filename[1], lyrs = band_num), y = data_proj)
        data_crop_i <- crop(data_proj_i, bbox)
        data_crop <- terra::merge(data_crop, data_crop_i)
        data_step <- data_step+1
        # plot(data_base)
      }
      rm(data_base_i, data_proj_i, data_crop_i); gc()
    }
    data_base <- data_crop
  }
  
  # Remove unneeded mask layers if compiling the water/land mask
  if(land_mask){
    data_base <- terra::ifel(data_base %in% c(1, 2, 3, 4, 5), NA, data_base)
    # plot(data_base)
  }
  
  # Save and quit
  writeRaster(data_base, file_name_out, overwrite = TRUE)
}


# MODIS data --------------------------------------------------------------

# Lists all products that are currently searchable on Earth Data and related servers
# Or see user guide: https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf
earth_data_catalogue <- luna::getProducts()

# Filter out just the MODIS products
MODIS_catalogue <- earth_data_catalogue[grepl("MOD|MYD", earth_data_catalogue$short_name),]

# MODIS directory that lists all products
# https://nrt3.modaps.eosdis.nasa.gov/archive/allData/61

# Uncomment and run any of these lines to open the product page in your web browser
# MODIS/Aqua Surface Reflectance Daily L2G Global 250m SIN Grid V061
# productInfo("MYD09GQ")

# Level 1
# NB: L1 products are often very difficult to work with in R
# productInfo("MYD01")
# productInfo("MYD021KM")
# productInfo("MYD02HKM")
# productInfo("MYD02QKM")

# Level 2
# productInfo("MYD09GHK") #L2G 500 m
# "MYD09"

# Level 3
# productInfo("MYD09Q1" # MODIS/Aqua Surface Reflectance 8-Day L3 Global 250m SIN Grid V006
# moproductInfo("MYD09A1") # MODIS/Aqua Surface Reflectance 8-Day L3 Global 500m SIN Grid V006

# MODIS/Terra Land Water Mask Derived from MODIS and SRTM L3 Global 250m SIN Grid V061
# https://lpdaac.usgs.gov/documents/1915/MOD44W_User_Guide_ATBD_V61.pdf
# productInfo("MOD44W")


# Workflow ----------------------------------------------------------------

## 1) Setup ---------------------------------------------------------------

# Chose where you would like to save the files
dl_dir <- "~/data/MODIS"

# Chosen start and end dates for downloading
start_date <- "2020-10-01"; end_date <- "2020-10-05"

# Determine what you want your bounding box to be
# NB: The processing functions will fail if too much data are loaded at once
# Here is the area surrounding the Bay of Angels
study_coords <- matrix(c(
  6.9, 43.4,  # Bottom-left corner
  7.7, 43.4, # Bottom-right corner
  7.7, 43.8, # Top-right corner
  6.9, 43.8,  # Top-left corner
  6.9, 43.4   # Close the polygon (same as first point)
), ncol = 2, byrow = TRUE)

# Turn it into the necessary SpatVector object type
study_bbox <- vect(study_coords, crs = "EPSG:4326", type = "polygons")

# Print the object to verify it worked - should be four points that make a box
plot(study_coords)

# Chose the product ID you want to download
product_ID <- "MYD09GA" # L2G 500 m

# Look at the server and version info for the product of choice
earth_data_catalogue[earth_data_catalogue$short_name == product_ID,]

# Then choose the server and version number
# NB: 'LPCLOUD' is the preferred server, but is not always available
# NB: One should generally choose the newest version, i.e. the biggest number
product_server <- "LPCLOUD" # NB: Change this if not shown in the output shown above
product_version <- "061" # NB: Change this if not shown in the output shown above


## 2) Download files -------------------------------------------------------

# First have a peak at what files exist
# NB: If this throws an error, it may be necessary to manually change the download server
# Remove the server and version arguments and run again
# It should show all of the possible servers and versions from which the desired product can be downloaded
luna::getNASA(product_ID, start_date, end_date, aoi = study_bbox, download = FALSE, 
              server = product_server, version = product_version)

# If that looks reasonable, download them
# NB: If this doesn't work, then the product ID, even if it is listed, may not be findable by the luna package
luna::getNASA(product_ID, start_date, end_date, aoi = study_bbox, download = TRUE, overwrite = FALSE, 
              server = product_server, version = product_version,
              path = dl_dir, username = earth_up$usrname, password = earth_up$psswrd)

# To follow the rest of the examples below we also want to download the MODIS mask files
luna::getNASA("MOD44W", start_date, end_date, aoi = study_bbox, download = TRUE, overwrite = FALSE,
              path = dl_dir, username = earth_up$usrname, password = earth_up$psswrd)


## 3) Process files --------------------------------------------------------

# Set file pathways
# NB: Change the directory to where you saved the files if it was changed
# NB: Change the pattern in rast_files to match the product ID you used if it is different
rast_files <- luna::modisDate(list.files(path = "~/data/MODIS", pattern = "MYD09GA\\.", full.names = TRUE)) # Level 2 daily
mask_files <- luna::modisDate(list.files(path = "~/data/MODIS", pattern = "MOD44W\\.", full.names = TRUE))

# Process all of the water mask files
# NB: This requires that this folder exists: ~/data/MODIS
# IF not, create it or change the directories below as desired
plyr::d_ply(.data = mask_files, .variables = c("date"), .fun = MODIS_proc, .parallel = FALSE,
            bbox = study_bbox, out_dir = "~/data/MODIS", band_num = 2, land_mask = TRUE)

# Load the desired mask file
MODIS_mask <- rast("~/data/MODIS/study_area_MOD44W_2020-01-01.tif")

# Check that it looks correct - should show white where land would be
plot(MODIS_mask)

# Prep one day of MODIS data
# NB: This requires that this folder exists: ~/data/MODIS
# IF not, create it or change the directories below to match 
plyr::d_ply(.data = rast_files, .variables = c("date"), .fun = MODIS_proc, .parallel = FALSE,
            bbox = study_bbox, out_dir = "~/data/MODIS", band_num = 1, land_mask = FALSE)

# Load the file
MODIS_rast <- rast("~/data/MODIS/study_area_MYD09GA_2020-10-01.tif")

# Check that it looks correct
plot(MODIS_rast)

# Project the 250 m mask to the same grid as the 500 m raster data
MODIS_mask_proj <- project(MODIS_mask, MODIS_rast)

# Check that it worked
plot(MODIS_mask_proj)

# Mask the raster data
MODIS_water <- mask(MODIS_rast, MODIS_mask_proj)

# Check to see if it looks correct - should show white where land is
plot(MODIS_water)

# Convert to data.frame for easy plotting
MODIS_water_df <- as.data.frame(MODIS_water, xy = TRUE, na.rm = TRUE) |> 
  # NB: Change 'sur_refl_b01_1' to match the correct column name in your data
  dplyr::rename(Rrs = sur_refl_b01_1)


# 4) Plot data ------------------------------------------------------------

# Map
pl_map <- MODIS_water_df |> 
  # Round all surface reflectance values greater than 0.1 down to 0.1 for better plotting
  # mutate(Rrs = case_when(Rrs > 0.1 ~ 0.1, 
  #                        Rrs < 0 ~ 0, TRUE ~ Rrs)) |> 
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = x, y = y, fill = Rrs)) +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
  # NB: Change fill label to correctly indicate which band width was used
  labs(x = "Longitude (°E)", y = "Latitude (°N)", fill = "Surface reflectance (459-479 nm) ") +
  coord_quickmap(xlim = range(MODIS_water_df$x), ylim = range(MODIS_water_df$y)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top", 
        legend.box = "vertical",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

# Save as desired
ggsave("~/data/MODIS/fig_MODIS.png", pl_map, height = 9, width = 14)

