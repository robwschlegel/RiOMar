# func/VOG.R

# Code used to extract data for studies linked to the VOG


# Setup ------------------------------------------------------------------

# Needed libraries
library(tidyverse)
library(sf)
library(furrr)
# library(doParallel)

# Set cores for use below
# registerDoParallel(cores = parallel::detectCores()-2)


# Functions --------------------------------------------------------------

# Helper to extract data for a given shapefile from a given .csv file
# file_name = panache_files[1]; polygon_sf = VOG_shape
extract_csv <- function(file_name, polygon_sf){

  # Read and convert .csv file
  df_sf <- st_as_sf(read_csv(file_name, show_col_types = FALSE), 
                    coords = c("lon", "lat"), crs = 4326)
  
  # Subset points within the polygon
  df_within <- data.frame()
  for(i in 1:nrow(polygon_sf)){
    df_within_sub <- df_sf[st_within(df_sf, polygon_sf[i,], sparse = FALSE), ]
    df_within_sub$Name <- polygon_sf$Name[i]
    df_within <- rbind(df_within, df_within_sub)
  }

  # Convert back to data.frame
  df_res <- df_within |> 
    mutate(lon = st_coordinates(df_within)[, 1],
           lat = st_coordinates(df_within)[, 2],
           date = as.Date(basename(file_name))) |> 
    st_drop_geometry() |> 
    dplyr::select(Name, date, lon, lat, mask)
  return(df_res)
}


# Load shapes ------------------------------------------------------------

# The overall shape
VOG_shape <- st_polygonize(st_transform(read_sf("data/VOG_shapes/Outline_VOG_2/Outline_VOG.shp"), crs = 4326))
VOG_shape$Name <- "VOG full"
plot(VOG_shape)

# The specific zones of extraction
VOG_zones <- st_transform(read_sf("data/VOG_shapes/Zone_extraction/Zone_Extraction_Sat_VOG.shp"), crs = 4326)
plot(VOG_zones)


# Extract data -----------------------------------------------------------

# List of the Bay of Biscay panache output files
panache_files <- dir("output/REGIONAL_PLUME_DETECTION/BAY_OF_BISCAY/SEXTANT/SPM/merged/Standard/PLUME_DETECTION/DAILY", 
                      recursive = TRUE, pattern = ".csv", full.names = TRUE)

# Prep multicore environment
plan(multisession, workers = parallel::detectCores() - 4)

# Extract VOG shape pixels
pixel_ts_VOG_shape <- future_map_dfr(panache_files, extract_csv, 
  polygon_sf = VOG_shape, .options = furrr_options(seed=TRUE))
write_csv(pixel_ts_VOG_shape, "data/VOG_shapes/pixel_ts_VOG_shape.csv")

# Extract VOG zones pixels
pixel_ts_VOG_zones <- future_map_dfr(panache_files, extract_csv, 
  polygon_sf = VOG_zones, .options = furrr_options(seed=TRUE))
write_csv(pixel_ts_VOG_zones, "data/VOG_shapes/pixel_ts_VOG_zones.csv")

# Close multicore environment
plan(sequential)


# Test visuals -----------------------------------------------------------

ggplot() +
  # Pixels from the VOG shape output
  geom_raster(data = filter(pixel_ts_VOG_shape, date == "1998-01-01"),
              aes(x = lon, y = lat), fill = "blue") +
  # Pixels from the VOG zones output
  geom_raster(data = filter(pixel_ts_VOG_zones, date == "1998-01-01"),
            aes(x = lon, y = lat), fill = "green") +
  # VOG shape
  geom_sf(data = VOG_shape, fill  = NA, colour = "black", linewidth = 0.8) +
  # VOG zones
  geom_sf(data = VOG_zones, fill  = NA, colour = "black", linewidth = 0.8) +
  # Plotting region and sf corrd adjustment
  coord_sf(
    xlim = st_bbox(VOG_shape$geometry)[c("xmin", "xmax")],
    ylim = st_bbox(VOG_shape$geometry)[c("ymin", "ymax")],
    expand = TRUE) +
  # Pretty
  labs(title = "VOG rasters plus extracted data", fill  = "Mask",
       x = NULL, y = NULL) + 
  theme_minimal()

