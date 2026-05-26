# ODATIS-MR_expert_subset.R

# This script subsets ODATIS-MR data directly on the rs-pro2 server
# It will not run locally outside of this environment
# It is designed to streamline the process of subsetting the data for further use

# NB: This script should be run from /data/home/rschlegel/


# Setup ------------------------------------------------------------------

# Activate personal libary
.libPaths("~/R/library")

# Libraries
# library(tidyverse) # Too many dependencies for the version of Debian on Rspro2
library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(furrr)
library(ggplot2)
library(ncdf4)
library(tidync)
library(geosphere)
library(raster)


# Functions --------------------------------------------------------------

# Simple wrapper to extract start and end times of an ODATIS-MR file
# file_name <- "../ODATIS-MR/MODIS/L3m_20020704__FRANCE_03_MOD_CDOM-NS_DAY_00.nc"
get_start_end_time <- function(file_name){
  
  # Get global info
  nc_data <- nc_open(file_name)
  
  # Extract dates accordingly
  df_time <- mutate(start_time = as.POSIXct(gsub("T|Z", " ", ncatt_get(nc_data, varid = 0)[["start_time"]]), format = "%Y%m%d %H%M%S", tz = "GMT"),
                     end_time = as.POSIXct(gsub("T|Z", " ", ncatt_get(nc_data, varid = 0)[["end_time"]]), format = "%Y%m%d %H%M%S", tz = "GMT")) |> 
    dplyr::select(start_time, end_time)

  # Exit
  nc_close(nc_data)
  return(df_time)
}

# Function to extract lon/lat subset from a netCDF file and convert to a dataframe
extract_lon_lat_subset <- function(nc_file, lon_range, lat_range) {

  # Open the netCDF file
  nc_data <- nc_open(nc_file)
  nc_date_val <- as.Date(ncatt_get(nc_data, varid = 0)[["period_start_day"]], format = "%Y%m%d")
  nc_close(nc_data)

  # Get the data within the lon/lat range
  df_sub <- tidync(nc_file) |> 
    hyper_filter(lon = between(lon, lon_range[1], lon_range[2]),
                  lat = between(lat, lat_range[1], lat_range[2])) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
          lat = as.numeric(lat),
          date = as.Date(nc_date_val))
  
  # Exit
  return(df_sub)
}

# Simple wrapper to extract and save full satellite coord grids
get_sat_grid <- function(file_name){
  nc_data <- nc_open(file_name)
  nc_lon <- as.vector(ncvar_get(nc_data, "lon"))
  nc_lat <- as.vector(ncvar_get(nc_data, "lat"))
  nc_close(nc_data)
  coords_sat <- expand.grid(lon = nc_lon, lat = nc_lat, KEEP.OUT.ATTRS = FALSE)
  return(coords_sat)
}

# Create indexes of which pixels match the 1 km range grid around the in situ sites
# target_site = zone_sites[10,]; dist_range = 1; sat_rast = rast_sat; sat_grid = grid_sat
# rm(site_sp, site_buffer, cropped_rast, pixel_coords)
get_pixels <- function(target_site, sat_grid, sat_rast, dist_range = 1){
  
  # Filter the sat_grid by the distance from the target site
  # NB: dist_tange/100*2 is a quick and dirty way to reduce the number of pixels to calculate the distance for
  pixel_coords <- sat_grid |> 
    filter(lon >= target_site$lon - dist_range/100*2, lon <= target_site$lon + dist_range/100*2,
           lat >= target_site$lat - dist_range/100*2, lat <= target_site$lat + dist_range/100*2)
                             
  # Calculate distance of pixels from the target in km
  pixel_coords$dist <- round(distHaversine(pixel_coords, target_site[,c("lon", "lat")])/1000, 2)
  
  # Select the nearest 49 (7x7) pixels
  pixel_coords <- pixel_coords |> arrange(dist) |> slice_head(n = 49)
  
  # Get the pixel IDs and exit
  pixel_coords$cell_numbers <- raster::cellFromXY(sat_rast, xy = pixel_coords)

  # Clean up and exit
  pixel_coords <- pixel_coords |> 
    mutate(zone = target_site$zone, source = target_site$source, site = target_site$site) |> 
    dplyr::select(zone, source, site, lon, lat, dist, cell_numbers)
  return(pixel_coords)
}

# Wrapper for get_pixels() to write the pixels to .csv per sat product
# sat_name <- "OLCI-A"; var_name_full <- "SPM-G-PO"
# sat_name <- "MODIS"; var_name_full <- "CHL-OC5-NS"
write_pixel_coords <- function(sat_name, var_name_full){

  # Load the in situ site metadata
  zone_sites <- read_csv("metadata/in_situ_site_list.csv", show_col_types = FALSE) |> 
    filter(!is.na(zone))

  # Determine file and variable names
  # Logic gate handles the different dates of the fil used for pixels
  if(sat_name == "MODIS"){
    nc_file_base <- paste0("../ODATIS-MR/FRANCE/MODIS/L3m_20200101__FRANCE_03_MOD_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "MERIS"){
    nc_file_base <- paste0("../ODATIS-MR/FRANCE/MERIS/L3m_20090101__FRANCE_03_MER_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "OLCI-A"){
    nc_file_base <- paste0("../ODATIS-MR/FRANCE/OLCI-A/L3m_20200101__FRANCE_03_OLA_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "OLCI-B"){
    nc_file_base <- paste0("../ODATIS-MR/FRANCE/OLCI-B/L3m_20200101__FRANCE_03_OLB_",var_name_full,"_DAY_00.nc")
  }
  sat_var <- paste0(var_name_full,"_mean")
  
  # Get the grid and raster bases
  grid_sat <- get_sat_grid(nc_file_base)
  rast_sat <- raster::raster(nc_file_base, varname = sat_var)
  
  # Extract all pixels
  plan(multisession, workers = parallel::detectCores() - 10)
  zone_pixels <- future_map_dfr(1:nrow(zone_sites),
                                function(i) get_pixels(target_site = zone_sites[i,], 
                                                       sat_grid = grid_sat, sat_rast = rast_sat, 
                                                       dist_range = 1),
                                .options = furrr_options(seed = TRUE))
  plan(sequential)
  
  # Save and exit
  write_csv(zone_pixels, paste0("metadata/zone_pixels_",sat_name,"_",var_name_full,".csv"))
}

# Once the pixels have been determined, use this to extract the data
# file_name <- "../ODATIS-MR/FRANCE/OLCI-B/L3m_20190615__FRANCE_03_OLB_CHL-GONS-AC_DAY_00.nc"
# file_name <- files_var[1]
# df <- zone_pixels
extract_pixels <- function(file_name, df){
  
  # Determine variable names from file pathway
  var_base_name <- str_split(basename(file_name), "_")[[1]][7]
  var_nc_name <- paste0(var_base_name,"_mean")

  # Get date values
  nc_date <- as.Date(str_split(basename(file_name), "_")[[1]][2], format = "%Y%m%d")
  
  # Legacy code kept here for testing purposes
  # Get raster
  # sat_rast <- raster::raster(file_name, varname = var_nc_name)
  
  # Get values
  # sat_vals <- raster::extract(sat_rast, df$cell_numbers)
  
  # Add to data.frame
  # df_res <- df |> 
  #   mutate(date = nc_date,
  #          variable = var_col_name,
  #          value = sat_vals) |> 
  #   dplyr::select(-cell_numbers)
  
  # Extract data and exit
  # NB: A couple of NetCDF files are mysteriously corrupt
  df_res <- tryCatch({
    df |>
      mutate(date = nc_date,
             variable = var_base_name,
             value = raster::extract(raster::raster(file_name, varname = var_nc_name), cell_numbers)) |>
      dplyr::select(-cell_numbers)
  }, error = function(e) {
    message(file_name, " could not be read : ",e$message)
    df |>
      mutate(date = nc_date,
             variable = var_base_name,
             value = NA) |>
      dplyr::select(-cell_numbers)
  })
  return(df_res)
  
  # Load the full nc file to test the raster extraction method
  # df_nc <- tidync::tidync(file_name) |> tidync::hyper_tibble() |>
  #   filter(!is.na(`CHL-OC5-PO_mean`)) |>
  #       # mutate(lon = as.numeric(lon),
  #               # lat = as.numeric(lat)) |>
  #   mutate(lon = round(as.numeric(lon)/0.0005) * 0.0005,
  #          lat = round(as.numeric(lat), 4)) |>
  #   distinct() #|> 
  #   # summarise(SPM_G_PO_mean = mean(`SPM-G-PO_mean`, na.rm = TRUE), .by = c("lon", "lat"))
  
  # Merge and test similarity
  # df_test <- df_res |>
  #   filter(!is.na(value)) |>
  #   mutate(lon = round(as.numeric(lon)/0.0005) * 0.0005,
  #          lat = round(lat, 4)) |>
  #   # summarise(value = mean(value, na.rm = TRUE), .by = c("lon", "lat")) |>
  #   left_join(df_nc, by = c("lon", "lat")) |>
  #       mutate(diff = value-`CHL-OC5-PO_mean`)
  #   # mutate(diff = mean(value-`CHL-OC5-PO_mean`, na.rm = TRUE))
  
  # Test in situ stations with all missing data
  # Eyrac, pk 30, pk 52, Anse du Piquet, Baie d'Yves (a), Cotard, Ile d'Aix, Truscat, Antoine, Luc-sur-Mer
  # is_site <- "Luc-sur-Mer"
  # df_is_test <- df_res |> filter(site == is_site)
  # df_is_test_mean <- summarise(df_is_test, lon = mean(lon), lat = mean(lat), .by = c("zone", "source", "site"))
  # df_nc_test <- df_nc |> filter(lon >= min(df_is_test$lon)-1, lon <= max(df_is_test$lon)+1,
  #                               lat >= min(df_is_test$lat)-1, lat <= max(df_is_test$lat)+1)
  # ggplot() +
  #   annotation_borders() +
  #   geom_raster(data = df_nc_test, aes( x = lon, y = lat, fill = `CHL-OC5-PO_mean`)) +
  #   geom_raster(data = df_is_test, aes( x = lon, y = lat), fill = "red") +
  #   geom_point(data = df_is_test_mean, aes(x = lon, y = lat), colour = "darkred") +
  #   coord_quickmap(xlim = range(df_nc_test$lon), ylim = range(df_nc_test$lat)) +
  #   scale_fill_viridis_c()
  
  # Clean up
  # rm(file_name, df, var_col_name, var_nc_name, df_res, nc_date, df_nc, df_test, df_is_test, df_is_test_mean, df_nc_test, is_site)
}

# Wrapper function to be called per variable
# sat_name <- "OLCI-B"; var_name <- "CHL-OC5-PO"
process_pixels <- function(sat_name, var_name){
  
  # Load the pixel metadata
  file_stub <- paste0(sat_name,"_",var_name)
  
  # Get file pathways
  files_path <- file.path("../ODATIS-MR", "FRANCE", sat_name)
  files_var <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = var_name)

  # Extract data for all files and save
  if(length(files_var) > 0){

    file_name_var <- paste0("output/zone_data_",file_stub,".csv")

    if(!file.exists(file_name_var)){

      # Extract data from all pixels
      message("Started ",var_name," extraction at : ", Sys.time())
      zone_pixels <- read_csv(paste0("metadata/zone_pixels_",file_stub,".csv"), show_col_types = FALSE)

      plan(multisession, workers = parallel::detectCores() - 4)
      zone_data_var <- future_map_dfr(files_var, extract_pixels, df = zone_pixels,
                                      .options = furrr_options(seed = TRUE))
      
      # Save results
      message("Saving ",var_name," extraction at : ", Sys.time())
      data.table::fwrite(zone_data_var, file_name_var)
      plan(sequential)
      gc()
      gc()
    }

    # Create median value time series
    file_name_median_all <- paste0("output/zone_median_",file_stub,"_all.csv")
    file_name_median_small <- paste0("output/zone_median_",file_stub,"_small.csv")
    
    # Create medians etc. from all pixels
    if(!file.exists(file_name_median_all)){
      # Load data from .csv
      if(!exists("zone_data_var")){
        message("Loading ",var_name," extraction at : ", Sys.time())
        zone_data_var <- data.table::fread(file_name_var)
      }
      message("Started median all calculations at : ", Sys.time())
      zone_median_all <- zone_data_var |> 
        filter(value > 0) |>
        summarise(median = median(value, na.rm = TRUE), 
                  mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE),
                  n = n(), .by = c("zone", "source", "site", "date", "variable"))
      data.table::fwrite(zone_median_all, file_name_median_all)
      rm(zone_median_all); gc()
    }

    # Create medians etc. from 'small' pixels
    if(!file.exists(file_name_median_small)){
      # Load data from .csv
      if(!exists("zone_data_var")){
        message("Loading ",var_name," extraction at : ", Sys.time())
        zone_data_var <- data.table::fread(file_name_var)
      }
      message("Started median small calculations at : ", Sys.time())
      slice_n <- 9
      zone_median_small <- zone_data_var |> 
        filter(value > 0) |> # NB: This removes NA pixels before selecting the nearest pixel, which may be incorrect
        group_by(zone, source, site, date, variable) |> 
        arrange(dist) |> 
        slice_head(n = slice_n) |> 
        ungroup() |> 
        summarise(median = median(value, na.rm = TRUE), 
                  mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE),
                  n = n(), .by = c("zone", "source", "site", "date", "variable"))
      data.table::fwrite(zone_median_small, file_name_median_small)
      rm(zone_median_small); gc()
    }

    # Cleanup and exit
    if(exists("zone_data_var")) rm(zone_data_var)
    if(exists("zone_pixels")) rm(zone_pixels)
    rm(file_stub, files_path, files_var, file_name_var, file_name_median_all, file_name_median_small)
    gc()
    gc()
  }
}

# Function that extracts all of the data identified in the product/zone lon/lat metadata files
# sat_name <- "MODIS"
extract_pixels_all <- function(sat_name){

  message("Started run on ", sat_name, " : ", Sys.time())
  
  ## CHL1
  process_pixels(sat_name, "CHL1-AC"); gc(); gc()
  process_pixels(sat_name, "CHL1-NS"); gc(); gc()
  process_pixels(sat_name, "CHL1-PO"); gc(); gc()

  ## CHL-GONS - No NS products
  process_pixels(sat_name, "CHL-GONS-AC"); gc(); gc()
  process_pixels(sat_name, "CHL-GONS-PO"); gc(); gc()
  
  ## CHL-OC5
  process_pixels(sat_name, "CHL-OC5-AC"); gc(); gc()
  process_pixels(sat_name, "CHL-OC5-PO"); gc(); gc()
  process_pixels(sat_name, "CHL-OC5-NS"); gc(); gc()

  ## SPM-G
  process_pixels(sat_name, "SPM-G-AC"); gc(); gc()
  process_pixels(sat_name, "SPM-G-PO"); gc(); gc()
  process_pixels(sat_name, "SPM-G-NS"); gc(); gc()
  
  ## SPM-R
  process_pixels(sat_name, "SPM-R-AC"); gc(); gc()
  process_pixels(sat_name, "SPM-R-PO"); gc(); gc()
  process_pixels(sat_name, "SPM-R-NS"); gc(); gc()

  ## TUR
  process_pixels(sat_name, "T-FNU-AC"); gc(); gc()
  process_pixels(sat_name, "T-FNU-PO"); gc(); gc()
  process_pixels(sat_name, "T-FNU-NS"); gc(); gc()
  
  ## SST - Only NS and PO
  process_pixels(sat_name, "SST-NS"); gc(); gc()
  process_pixels(sat_name, "SST-PO"); gc(); gc()

  ## SST-NIGHT - Only for NS
  process_pixels(sat_name, "SST-NIGHT-NS"); gc(); gc()
  
  # # Clean and exit
  message("Finished ", sat_name," at : ", Sys.time())
  gc(); gc()
}


# File dirs --------------------------------------------------------------

# List all the netCDF files in the given sensor directory
# nc_files_MODIS <- list.files("../ODATIS-MR/FRANCE/MODIS", pattern = "\\.nc$", full.names = TRUE)
# nc_files_MERIS <- list.files("../ODATIS-MR/FRANCE/MERIS", pattern = "\\.nc$", full.names = TRUE)
# nc_files_OLCIA <- list.files("../ODATIS-MR/FRANCE/OLCI-A", pattern = "\\.nc$", full.names = TRUE)
# nc_files_OLCIB <- list.files("../ODATIS-MR/FRANCE/OLCI-B/", pattern = "\\.nc$", full.names = TRUE)

# Filter only certain files based on naming convention (e.g., only OC5 processing)
# In a future version one would not subset by date and rather get everthing
## MODIS
# nc_files_MODIS_SPM_R_NS <- nc_files_MODIS[grepl("(?=.*SPM-R-NS)(?=.*L3m_2019)", nc_files_MODIS, perl = TRUE)]
# nc_files_MODIS_SPM_G_NS <- nc_files_MODIS[grepl("(?=.*SPM-G-NS)(?=.*L3m_2019)", nc_files_MODIS, perl = TRUE)]
# nc_files_MODIS_SPM_R_PO <- nc_files_MODIS[grepl("(?=.*SPM-R-PO)(?=.*L3m_2019)", nc_files_MODIS, perl = TRUE)]
# nc_files_MODIS_SPM_G_PO <- nc_files_MODIS[grepl("(?=.*SPM-G-PO)(?=.*L3m_2019)", nc_files_MODIS, perl = TRUE)]

## MERIS
# nc_files_MERIS_SPM_R_AC <- nc_files_MERIS[grepl("(?=.*SPM-R-AC)(?=.*L3m_2009)", nc_files_MERIS, perl = TRUE)]
# nc_files_MERIS_SPM_G_AC <- nc_files_MERIS[grepl("(?=.*SPM-G-AC)(?=.*L3m_2009)", nc_files_MERIS, perl = TRUE)]
# nc_files_MERIS_SPM_R_PO <- nc_files_MERIS[grepl("(?=.*SPM-R-PO)(?=.*L3m_2009)", nc_files_MERIS, perl = TRUE)]
# nc_files_MERIS_SPM_G_PO <- nc_files_MERIS[grepl("(?=.*SPM-G-PO)(?=.*L3m_2009)", nc_files_MERIS, perl = TRUE)]

## OLCI-A
# nc_files_OLCIA_SPM_R_AC <- nc_files_OLCIA[grepl("(?=.*SPM-R-AC)(?=.*L3m_2019)", nc_files_OLCIA, perl = TRUE)]
# nc_files_OLCIA_SPM_G_AC <- nc_files_OLCIA[grepl("(?=.*SPM-G-AC)(?=.*L3m_2019)", nc_files_OLCIA, perl = TRUE)]
# nc_files_OLCIA_SPM_R_PO <- nc_files_OLCIA[grepl("(?=.*SPM-R-PO)(?=.*L3m_2019)", nc_files_OLCIA, perl = TRUE)]
# nc_files_OLCIA_SPM_G_PO <- nc_files_OLCIA[grepl("(?=.*SPM-G-PO)(?=.*L3m_2019)", nc_files_OLCIA, perl = TRUE)]

## OLCI-B
# nc_files_OLCIB_SPM_R_AC <- nc_files_OLCIB[grepl("(?=.*SPM-R-AC)(?=.*L3m_2019)", nc_files_OLCIB, perl = TRUE)]
# nc_files_OLCIB_SPM_G_AC <- nc_files_OLCIB[grepl("(?=.*SPM-G-AC)(?=.*L3m_2019)", nc_files_OLCIB, perl = TRUE)]
# nc_files_OLCIB_SPM_R_PO <- nc_files_OLCIB[grepl("(?=.*SPM-R-PO)(?=.*L3m_2019)", nc_files_OLCIB, perl = TRUE)]
# nc_files_OLCIB_SPM_G_PO <- nc_files_OLCIB[grepl("(?=.*SPM-G-PO)(?=.*L3m_2019)", nc_files_OLCIB, perl = TRUE)]


# Extract data -----------------------------------------------------------

# Set lon and lat range for data extraction

# lon_range <- c(6.8925000, 7.4200000)
# lat_range <- c(43.2136389, 43.7300000)

# Set core usage
# plan(multisession, workers = parallel::detectCores() - 2)

# Extract the data as a data.frame
## MODIS
# MODIS_SPM_R_NS_sub <- future_map_dfr(nc_files_MODIS_SPM_R_NS, extract_lon_lat_subset, lon_range, lat_range)
# MODIS_SPM_G_NS_sub <- future_map_dfr(nc_files_MODIS_SPM_G_NS, extract_lon_lat_subset, lon_range, lat_range)
# MODIS_SPM_R_PO_sub <- future_map_dfr(nc_files_MODIS_SPM_R_PO, extract_lon_lat_subset, lon_range, lat_range)
# MODIS_SPM_G_PO_sub <- future_map_dfr(nc_files_MODIS_SPM_G_PO, extract_lon_lat_subset, lon_range, lat_range)

## MERIS
# MERIS_SPM_R_AC_sub <- future_map_dfr(nc_files_MERIS_SPM_R_AC, extract_lon_lat_subset, lon_range, lat_range)
# MERIS_SPM_G_AC_sub <- future_map_dfr(nc_files_MERIS_SPM_G_AC, extract_lon_lat_subset, lon_range, lat_range)
# MERIS_SPM_R_PO_sub <- future_map_dfr(nc_files_MERIS_SPM_R_PO, extract_lon_lat_subset, lon_range, lat_range)
# MERIS_SPM_G_PO_sub <- future_map_dfr(nc_files_MERIS_SPM_G_PO, extract_lon_lat_subset, lon_range, lat_range)

## OLCI-A
# OLCIA_SPM_R_AC_sub <- future_map_dfr(nc_files_OLCIA_SPM_R_AC, extract_lon_lat_subset, lon_range, lat_range)
# OLCIA_SPM_G_AC_sub <- future_map_dfr(nc_files_OLCIA_SPM_G_AC, extract_lon_lat_subset, lon_range, lat_range)
# OLCIA_SPM_R_PO_sub <- future_map_dfr(nc_files_OLCIA_SPM_R_PO, extract_lon_lat_subset, lon_range, lat_range)
# OLCIA_SPM_G_PO_sub <- future_map_dfr(nc_files_OLCIA_SPM_G_PO, extract_lon_lat_subset, lon_range, lat_range)

## OLCI-B
# OLCIB_SPM_R_AC_sub <- future_map_dfr(nc_files_OLCIB_SPM_R_AC, extract_lon_lat_subset, lon_range, lat_range)
# OLCIB_SPM_G_AC_sub <- future_map_dfr(nc_files_OLCIB_SPM_G_AC, extract_lon_lat_subset, lon_range, lat_range)
# OLCIB_SPM_R_PO_sub <- future_map_dfr(nc_files_OLCIB_SPM_R_PO, extract_lon_lat_subset, lon_range, lat_range)
# OLCIB_SPM_G_PO_sub <- future_map_dfr(nc_files_OLCIB_SPM_G_PO, extract_lon_lat_subset, lon_range, lat_range)

# Rest core usage to sequential
# plan(sequential)

# Save outputs
## MODIS
# save(MODIS_SPM_R_NS_sub, file = "output/MODIS_SPM_R_NS_sub.RData")
# save(MODIS_SPM_G_NS_sub, file = "output/MODIS_SPM_G_NS_sub.RData")
# save(MODIS_SPM_R_PO_sub, file = "output/MODIS_SPM_R_PO_sub.RData")
# save(MODIS_SPM_G_PO_sub, file = "output/MODIS_SPM_G_PO_sub.RData")

## MERIS
# save(MERIS_SPM_R_AC_sub, file = "output/MERIS_SPM_R_AC_sub.RData")
# save(MERIS_SPM_G_AC_sub, file = "output/MERIS_SPM_G_AC_sub.RData")
# save(MERIS_SPM_R_PO_sub, file = "output/MERIS_SPM_R_PO_sub.RData")
# save(MERIS_SPM_G_PO_sub, file = "output/MERIS_SPM_G_PO_sub.RData")

## OLCI-A
# save(OLCIA_SPM_R_AC_sub, file = "output/OLCIA_SPM_R_AC_sub.RData")
# save(OLCIA_SPM_G_AC_sub, file = "output/OLCIA_SPM_G_AC_sub.RData")
# save(OLCIA_SPM_R_PO_sub, file = "output/OLCIA_SPM_R_PO_sub.RData")
# save(OLCIA_SPM_G_PO_sub, file = "output/OLCIA_SPM_G_PO_sub.RData")

## OLCI-B
# save(OLCIB_SPM_R_AC_sub, file = "output/OLCIB_SPM_R_AC_sub.RData")
# save(OLCIB_SPM_G_AC_sub, file = "output/OLCIB_SPM_G_AC_sub.RData")
# save(OLCIB_SPM_R_PO_sub, file = "output/OLCIB_SPM_R_PO_sub.RData")
# save(OLCIB_SPM_G_PO_sub, file = "output/OLCIB_SPM_G_PO_sub.RData")


# Pixels for validation --------------------------------------------------

# MODIS
# write_pixel_coords("MODIS", "CHL-OC5-NS")
# write_pixel_coords("MODIS", "CHL-OC5-PO")
# write_pixel_coords("MODIS", "CHL1-NS")
# write_pixel_coords("MODIS", "CHL1-PO")
# write_pixel_coords("MODIS", "SPM-R-NS")
# write_pixel_coords("MODIS", "SPM-R-PO")
# write_pixel_coords("MODIS", "SPM-G-NS")
# write_pixel_coords("MODIS", "SPM-G-PO")
# write_pixel_coords("MODIS", "SST-NS")
# write_pixel_coords("MODIS", "SST-PO")
# write_pixel_coords("MODIS", "SST-NIGHT-NS")
# write_pixel_coords("MODIS", "T-FNU-NS")
# write_pixel_coords("MODIS", "T-FNU-PO")

# MERIS
# write_pixel_coords("MERIS", "CHL-GONS-AC")
# write_pixel_coords("MERIS", "CHL-GONS-PO")
# write_pixel_coords("MERIS", "CHL-OC5-AC")
# write_pixel_coords("MERIS", "CHL-OC5-PO")
# write_pixel_coords("MERIS", "CHL1-AC")
# write_pixel_coords("MERIS", "CHL1-PO")
# write_pixel_coords("MERIS", "SPM-R-AC")
# write_pixel_coords("MERIS", "SPM-R-PO")
# write_pixel_coords("MERIS", "SPM-G-AC")
# write_pixel_coords("MERIS", "SPM-G-PO")
# write_pixel_coords("MERIS", "T-FNU-AC")
# write_pixel_coords("MERIS", "T-FNU-PO")

# OLCI-A
# write_pixel_coords("OLCI-A", "CHL-GONS-AC")
# write_pixel_coords("OLCI-A", "CHL-GONS-PO")
# write_pixel_coords("OLCI-A", "CHL-OC5-AC")
# write_pixel_coords("OLCI-A", "CHL-OC5-PO")
# write_pixel_coords("OLCI-A", "CHL1-AC")
# write_pixel_coords("OLCI-A", "CHL1-PO")
# write_pixel_coords("OLCI-A", "SPM-R-AC")
# write_pixel_coords("OLCI-A", "SPM-R-PO")
# write_pixel_coords("OLCI-A", "SPM-G-AC")
# write_pixel_coords("OLCI-A", "SPM-G-PO")
# write_pixel_coords("OLCI-A", "T-FNU-AC")
# write_pixel_coords("OLCI-A", "T-FNU-PO")

# OLCI-B
# write_pixel_coords("OLCI-B", "CHL-GONS-AC")
# write_pixel_coords("OLCI-B", "CHL-GONS-PO")
# write_pixel_coords("OLCI-B", "CHL-OC5-AC")
# write_pixel_coords("OLCI-B", "CHL-OC5-PO")
# write_pixel_coords("OLCI-B", "CHL1-AC")
# write_pixel_coords("OLCI-B", "CHL1-PO")
# write_pixel_coords("OLCI-B", "SPM-R-AC")
# write_pixel_coords("OLCI-B", "SPM-R-PO")
# write_pixel_coords("OLCI-B", "SPM-G-AC")
# write_pixel_coords("OLCI-B", "SPM-G-PO")
# write_pixel_coords("OLCI-B", "T-FNU-AC")
# write_pixel_coords("OLCI-B", "T-FNU-PO")


# Extract data by site pixels --------------------------------------------

extract_pixels_all("MODIS"); gc(); gc()
extract_pixels_all("MERIS"); gc(); gc()
extract_pixels_all("OLCI-A"); gc(); gc()
extract_pixels_all("OLCI-B"); gc(); gc()

