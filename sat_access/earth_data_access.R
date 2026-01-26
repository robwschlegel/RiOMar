# earth_data_access
# This script provides the functions, workflow, and examples 
# to download any data product from the earth data server:
# https://www.earthdata.nasa.gov/

# The user guide for MODIS products may be found here :
# https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ncdf4)
# Luna package is used to access large spatial data products
# install.packages('luna', repos = 'https://rspatial.r-universe.dev')
library(terra)
library(luna)
library(doParallel); registerDoParallel(cores = detectCores() - 2)


# Setup -------------------------------------------------------------------

# Load username and password
earth_up <- read_csv("path/to/file/earthdata_pswd.csv")


# Functions ---------------------------------------------------------------

# Get possible MODIS files
# NB: At the moment this is optimized to work with just one day of data
MODIS_dl <- function(prod_id, dl_start, dl_end, bbox, usrname, psswrd, dl_files = TRUE, dl_dir = "data/MODIS"){
  
  # Check if dl_dir folder exists and create if necessary
  if(!dir.exists(dl_dir)){
    dir.create(dl_dir, recursive = TRUE)
  }
  
  # If download is FALSE, just print possible files
  if(!dl_files){
    message("Data files : ")
    luna::getNASA(prod_id, dl_start, dl_end, aoi = bbox, download = FALSE)
    message("Mask files : ")
    luna::getNASA("MOD44W", dl_start, dl_end, aoi = bbox, download = FALSE)
  } else {
    message("Data files : ")
    luna::getNASA(prod_id, dl_start, dl_end, aoi = bbox, download = TRUE, overwrite = FALSE,
                  path = dl_dir, username = earth_up$usrname, password = earth_up$psswrd)
    message("Mask files : ")
    luna::getNASA("MOD44W", dl_start, dl_end, aoi = bbox, download = TRUE, overwrite = FALSE,
                  path = dl_dir, username = usrname, password = psswrd)
  }
}

# Process MODIS data in a batch
MODIS_proc <- function(file_names, bbox, band_num = NULL, water_mask = FALSE){
  
  # Load files with desired layers etc.
  # NB: If run in parallel, merge() causes a crash to desktop
  if(water_mask){
    data_layers <- lapply(file_names, rast, subds = 2)
    data_merge <- do.call(merge, data_layers)
    # plot(data_merge)
    data_base <- terra::ifel(data_merge %in% c(1, 2, 3, 4, 5), NA, data_merge)
    # plot(data_base)
  } else {
    data_layers <- lapply(file_names, rast, subds = band_num) # Chose the band number (e.g. wavelength range)
    data_base <- do.call(merge, data_layers)
    # plot(data_base)
  }
  
  # Project to EPSG:4326
  data_base_proj <- project(data_base, y = "EPSG:4326")
  # plot(data_base_proj)
  
  # Crop to bbox and exit
  data_crop <- crop(data_base_proj, bbox)
  # plot(data_crop)
  return(data_crop)
}


# MODIS data --------------------------------------------------------------

# Load username and password
earth_up <- read_csv("~/pCloudDrive/Documents/info/earthdata_pswd.csv")

# Lists all products that are currently searchable
# Or see user guide: https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf
MODIS_prod <- luna::getProducts()

# Uncomment and run any of these lines to open the product page in your web browser
# MODIS/Aqua Surface Reflectance Daily L2G Global 250m SIN Grid V061
# productInfo("MYD09GQ")

# Level 1
# productInfo("MYD01")

# MODIS/Aqua Surface Reflectance 8-Day L3 Global 250m SIN Grid V006
# productInfo("MYD09Q1")

# MODIS/Aqua Surface Reflectance 8-Day L3 Global 500m SIN Grid V006
# moproductInfo("MYD09A1")

# MODIS/Terra Land Water Mask Derived from MODIS and SRTM L3 Global 250m SIN Grid V061
# https://lpdaac.usgs.gov/documents/1915/MOD44W_User_Guide_ATBD_V61.pdf
# productInfo("MOD44W")


# Workflow ----------------------------------------------------------------

## 1) Basics ---------------------------------------------------------------

# Determine what you want your bounding box to be
coords <- matrix(c(
  3, 32,  # Bottom-left corner
  27, 32, # Bottom-right corner
  27, 45, # Top-right corner
  3, 45,  # Top-left corner
  3, 32   # Close the polygon (same as first point)
), ncol = 2, byrow = TRUE)

# Turn it into the necessary SpatVector object type
bbox <- vect(coords, crs = "EPSG:4326", type = "polygons")

# Print the object to verify it worked
print(bbox)
plot(bbox)

# Chose the product ID you want to download
product_ID <- "MYD09A1"

# Chosen start and end dates for downloading
start_date <- "2024-08-09"; end_date <- "2024-08-16"


## 2) Download files -------------------------------------------------------

# Download data
# NB: Change the dl_dir if you want to save the files somewhere else
MODIS_dl(prod_id = "MYD09A1", dl_start = start_date, dl_end = end_date, 
         bbox = bbox, usrname = earth_up$usrname, psswrd = earth_up$psswrd,
         dl_files = TRUE, dl_dir = "data/MODIS")


## 3) Process files --------------------------------------------------------

# Set file pathways
# NB: Change the directory to where you saved the files if it was changed
mask_files <- list.files(path = "data/MODIS", pattern = "MOD", full.names = TRUE)
rast_files <- list.files(path = "data/MODIS", pattern = "MYD", full.names = TRUE)

# Process all of the water mask files
# NB: "data/MODIS/study_area_mask.tif" only needs to be created once
if(!file.exists("data/MODIS/study_area_mask.tif")){
  MODIS_mask <- MODIS_proc(file_names = mask_files, bbox = bbox, water_mask = TRUE)
  writeRaster(MODIS_mask, "data/MODIS/study_area_mask.tif", overwrite = TRUE)
}
MODIS_mask <- rast("data/MODIS/study_area_mask.tif")
plot(MODIS_mask)

# Prep one day of MODIS data
# NB: The same for the following two lines
if(!file.exists("data/MODIS/study_area_rast.tif")){
  MODIS_rast <- MODIS_proc(file_names = rast_files, bbox = bbox, band_num = 1) # NB: Change band_num accordingly
  writeRaster(MODIS_rast, "data/MODIS/study_area_rast.tif", overwrite = TRUE)
}
MODIS_rast <- rast("data/MODIS/study_area_rast.tif")
plot(MODIS_rast)

# Projected the 250 m mask to the same grid as the 500 m raster data
MODIS_mask_proj <- project(MODIS_mask, MODIS_rast)

# Mask the raster data
MODIS_water <- mask(MODIS_rast, MODIS_mask_proj)
plot(MODIS_water)

# Convert to data.frame for easy plotting
MODIS_water_df <- as.data.frame(MODIS_water, xy = TRUE, na.rm = TRUE)


# 4) Plot data ------------------------------------------------------------

# Map
pl_map <- ggplot(data = MODIS_water_df) +
  borders(fill = "grey80") +
  # NB: Will need to change 'sur_refl_b01' to match the correct column name in 'MODIS_water_df'
  geom_tile(aes(x = x, y = y, fill = sur_refl_b01)) +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
  # NB: Change fill label to correctly indicate which band width was used
  labs(x = "Longitude (°E)", y = "Latitude (°N)", fill = "Surface reflectance (459-479 nm) ") +
  coord_quickmap(xlim = c(7, 25), ylim = c(35, 42)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top", 
        legend.box = "vertical",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

# Save as desired
ggsave("figures/fig_MODIS.png", pl_map, height = 9, width = 14)

