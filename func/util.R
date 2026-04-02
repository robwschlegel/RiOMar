# func/util.R
# The storage point for many functions re-used by other scripts


# Meta-data ---------------------------------------------------------------

# In the future this will be taken from define_parameters() in func/util.py
river_mouths <- data.frame(row_name = 1:4,
                           mouth_name = c("Seine", "Gironde", "Loire", "Grand Rhone"),
                           mouth_lon = c(0.145, -1.05, -2.10, 4.83),
                           mouth_lat = c(49.43, 45.59, 47.29, 43.41))

# The zones en large
zones_list <- c("GULF_OF_LION", "BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY")

# Zone bounding boxes
zones_bbox <- data.frame(zone = zones_list,
                         lon_min = c(3.50, -1.50, -4.00, -5.00),
                         lon_max = c(6.00, 0.50, -0.50, -1.50),
                         lat_min  = c(42.25, 49.25, 44.50, 46.5),
                         lat_max = c(44.00, 50.25, 46.50, 48.00)) |> 
  mutate(zone_pretty = factor(zone, 
                              levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                              labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")), .after = "zone")

# Create FRANCE bounding box with same structure as zones_bbox
france_bbox <- data.frame(zone = "FRANCE",
                         lon_min = c(-7.8),
                         lon_max = c(10.3),
                         lat_min  = c(41.2),
                         lat_max = c(51.5)) 


# Pixels ------------------------------------------------------------------

# Simple wrapper to extract start and end times of an ODATIS-MR file
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_CDOM-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20230424__FRANCE_03_MOD_T-FNU-NS_DAY_00.nc"
# file_name <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/02/01/19980201-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
get_start_end_time <- function(file_name){
  
  # Get global info
  df_info <- ncdump::NetCDF(file_name)$attribute[1]$global
  
  # Extract dates accordingly
  if(grepl("SEXTANT", file_name)){
    df_time <- df_info |> 
      mutate(start_time = ymd_hms(paste(gsub("UTC", "", start_date), gsub("UTC", "", start_time)), tz = "GMT"),
             end_time = ymd_hms(paste(gsub("UTC", "", stop_date), gsub("UTC", "", stop_time)), tz = "GMT")) |>
      dplyr::select(start_time, end_time) 
  } else{
    df_time <- df_info |> 
      mutate(start_time = as.POSIXct(gsub("T|Z", " ", start_time), format = "%Y%m%d %H%M%S", tz = "GMT"),
             end_time = as.POSIXct(gsub("T|Z", " ", end_time), format = "%Y%m%d %H%M%S", tz = "GMT")) |> 
      dplyr::select(start_time, end_time)
  }

  # Exit
  return(df_time)
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
# target_site = zone_site_df[10,]; dist_range = 1; sat_rast = rast_sat
# rm(site_sp, site_buffer, cropped_rast, pixel_coords)
get_pixels <- function(target_site, sat_rast, dist_range = 1){
  
  # Create spatial point
  site_sp <- SpatialPoints(coords = target_site[,c("lon", "lat")], proj4string = CRS(proj4string(sat_rast)))
  
  # Create the chosen buffer 
  site_buffer <- raster::buffer(site_sp, width = dist_range*1000)
  # site_buffer <- st_buffer(site_sp, dist = dist_range*1000)
  
  # Crop and mask the raster
  cropped_rast <- raster::crop(sat_rast, site_buffer)
  # masked_rast <- raster::mask(cropped_rast, site_buffer)
  
  # Get lon/lat coords
  # pixel_values <- getValues(masked_rast)
  pixel_coords <- as.data.frame(raster::xyFromCell(cropped_rast, 1:ncell(cropped_rast)))
  colnames(pixel_coords) <- c("lon", "lat")
                             
  # Calculate distance of pixels from the target in km
  pixel_coords$dist <- round(distHaversine(pixel_coords, target_site[,c("lon", "lat")])/1000, 2)
  
  # Get the pixel IDs and exit
  pixel_coords$cell_numbers <- raster::cellFromXY(sat_rast, xy = pixel_coords)
  return(pixel_coords)
}

# Wrapper for get_pixels() to write the pixels to .csv per sat product
# zone_site_df <- filter(zone_sites, zone == "BAY_OF_BISCAY")
# sat_name <- "MODIS"
# sat_var <- "SPM-G-NS_mean"
# nc_file_base <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc"
write_pixels <- function(zone_site_df, sat_name){
  
  # Check that only one zone is given if not SEXTANT
  if(sat_name != "SEXTANT"){
    if(length(unique(zone_site_df$zone)) > 1) stop("Only pass one zone to this function for non-SEXTANT data structure.")
    zone_unique <- zone_site_df$zone[1]
  }
  
  # Determine file and variable names
  if(sat_name == "SEXTANT"){
    sat_var <- "analysed_spim"
    nc_file_base <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
  } else if(sat_name == "MODIS"){
    sat_var <- "SPM-G-NS_mean"
    nc_file_base <- paste0("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/",zone_unique,"/daily/L3m_20020704__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc")
  } else if(sat_name == "MERIS"){
    sat_var <- "SPM-G-PO_mean"
    nc_file_base <- paste0("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/",zone_unique,"/daily/L3m_20020619__FRANCE_03_MER_SPM-G-PO_DAY_00.nc")
  } else if(sat_name == "OLCI-A"){
    sat_var <- "SPM-G-PO_mean"
    nc_file_base <- paste0("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-A/",zone_unique,"/daily/L3m_20160426__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc")
  } else if(sat_name == "OLCI-B"){
    sat_var <- "SPM-G-PO_mean"
    nc_file_base <- paste0("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-B/",zone_unique,"/daily/L3m_20180515__FRANCE_03_OLB_SPM-G-PO_DAY_00.nc")
  }
  
  # Get the grid and raster bases
  grid_sat <- get_sat_grid(nc_file_base)
  rast_sat <- raster(nc_file_base, varname = sat_var)
  # Extract all pixels
  zone_pixels <- plyr::ddply(.data = zone_site_df, .variables = c("zone", "source", "site"), 
                             .fun = get_pixels, .parallel = TRUE, sat_rast = rast_sat)
  
  # Save and exit
  if(length(unique(zone_pixels$zone)) == 1){
    write_csv(zone_pixels, paste0("metadata/zone_pixels_",sat_name,"_",unique(zone_pixels$zone)[1],".csv"))
  } else {
    write_csv(zone_pixels, paste0("metadata/zone_pixels_",sat_name,".csv"))
  }
}

# Once the pixels have been determined, use this to extract the data
# file_name <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily/L3m_20020704__FRANCE_03_MOD_CHL-OC5-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/BAY_OF_SEINE/daily/L3m_20020619__FRANCE_03_MER_SPM-G-PO_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/SOUTHERN_BRITTANY/daily/L3m_20120908__FRANCE_03_MOD_CDOM-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/SOUTHERN_BRITTANY/daily/L3m_20040825__FRANCE_03_MOD_NRRS555-NS_DAY_00.nc"
# ncdump::NetCDF(file_name)
# df <- zone_pixels
extract_pixels <- function(file_name, df){
  
  # Determine variable names from file pathway
  if(grepl("SEXTANT", file_name)){
    if(grepl("SPM", file_name)){
      var_col_name <- "SPM"
      var_nc_name <- "analysed_spim"
    } else if(grepl("CHL", file_name)){
      var_col_name <- "CHLA"
      var_nc_name <- "analysed_chl_a"
    } else {
      stop("File structure not recognised")
    }
  } else {
    var_base_name <- str_split(basename(file_name), "_")[[1]][7]
    var_col_name <- str_split(var_base_name, "-")[[1]][1]
    var_nc_name <- paste0(var_base_name,"_mean")
  }
  
  # Get date values
  if(grepl("SEXTANT", file_name)){
    nc_date <- as.Date(str_split(basename(file_name), "-")[[1]][1], format = "%Y%m%d")
  } else {
    nc_date <- as.Date(str_split(basename(file_name), "_")[[1]][2], format = "%Y%m%d")
  }
  
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
             variable = var_col_name,
             value = raster::extract(raster(file_name, varname = var_nc_name), cell_numbers)) |>
      dplyr::select(-cell_numbers)
  }, error = function(e) {
    message(file_name, " could not be read : ",e$message)
    df |>
      mutate(date = nc_date,
             variable = var_col_name,
             value = NA) |>
      dplyr::select(-cell_numbers)
  })
  return(df_res)
  
  # Load the full nc file to test the raster extraction method
  # df_nc <- tidync::tidync(file_name) |> tidync::hyper_tibble() |>
  #   mutate(lon = plyr::round_any(as.numeric(lon), 0.005),
  #          lat = round(as.numeric(lat), 2)) |>
  #   dplyr::select(lon, lat, analysed_spim)
  
  # Merge and test similarity
  # df_test <- df_res |>
  #   mutate(lon = plyr::round_any(lon, 0.005),
  #          lat = round(lat, 2)) |>
  #   left_join(df_nc, by = c("lon", "lat")) |>
  #   mutate(diff = mean(value-analysed_spim, na.rm = TRUE))
  
  # Test in situ stations with all missing data
  # Eyrac, pk 30, pk 52, Anse du Piquet, Baie d'Yves (a), Cotard, Ile d'Aix, Truscat, Antoine, Luc-sur-Mer
  # is_site <- "Luc-sur-Mer"
  # df_is_test <- df_res |> filter(site == is_site)
  # df_is_test_mean <- summarise(df_is_test, lon = mean(lon), lat = mean(lat), .by = c("zone", "source", "site"))
  # df_nc_test <- df_nc |> filter(lon >= min(df_is_test$lon)-1, lon <= max(df_is_test$lon)+1,
  #                               lat >= min(df_is_test$lat)-1, lat <= max(df_is_test$lat)+1)
  # ggplot() +
  #   annotation_borders() +
  #   geom_raster(data = df_nc_test, aes( x = lon, y = lat, fill = analysed_spim)) +
  #   geom_raster(data = df_is_test, aes( x = lon, y = lat), fill = "red") +
  #   geom_point(data = df_is_test_mean, aes(x = lon, y = lat), colour = "darkred") +
  #   coord_quickmap(xlim = range(df_nc_test$lon), ylim = range(df_nc_test$lat)) +
  #   scale_fill_viridis_c()
  
  # Clean up
  # rm(file_name, df, var_col_name, var_nc_name, df_res, nc_date, df_nc, df_test, df_is_test, df_is_test_mean, df_nc_test, is_site)
}

# Function that extracts all of the data identified in the product/zone lon/lat metadata files
# sat_name <- "SEXTANT"
# sat_name <- "MODIS"; zone_name <- "SOUTHERN_BRITTANY"
extract_pixels_all <- function(sat_name, zone_name = NULL){#, overwrite = FALSE){
  
  # Load zone pixels
  if(is.null(zone_name)){
    file_stub <- sat_name
  } else {
    file_stub <- paste0(sat_name,"_",zone_name)
  }
  zone_pixels <- read_csv(paste0("metadata/zone_pixels_",file_stub,".csv"))
  message("Started run on ", file_stub, " : ", Sys.time())
  
  # Get file pathways
  if(sat_name == "SEXTANT"){
    files_SPM <- dir("~/pCloudDrive/data/SEXTANT/SPM", pattern = ".nc", full.names = TRUE, recursive = TRUE)
    files_CHL <- dir("~/pCloudDrive/data/SEXTANT/CHLA", pattern = ".nc", full.names = TRUE, recursive = TRUE)
  } else {
    files_path <- file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR",sat_name,zone_name)
    files_CHL <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "CHL")
    files_SPM <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "SPM")
    files_TUR <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "T-FNU")
    files_CDOM <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "CDOM")
    files_RRS <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "RRS")
  } 
  
  # Load SST explicitly for MODIS
  if(sat_name == "MODIS") files_SST <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "SST")
  
  ## CHL
  # TODO: Consider a different way of doing this. I think it could be more streamlined.
  if(exists("files_CHL")){
    message("Started CHL extraction at : ", Sys.time())
    file_name_CHL <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_CHL.csv")
    if(!file.exists(file_name_CHL)){
      zone_data_CHL <- plyr::ldply(.data = files_CHL, .fun = extract_pixels,
                                   .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      data.table::fwrite(zone_data_CHL, file_name_CHL); gc()
    } else {
      if(!exists("zone_data_SEXTANT_CHL")){
        zone_data_CHL <- data.table::fread(file_name_CHL) |> mutate(date = as.Date(date))
      }
    }
  }
  
  # Extract all available data per variable
  ## SPM
  if(exists("files_SPM")){
    message("Started SPM extraction at : ", Sys.time())
    file_name_SPM <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_SPM.csv")
    if(!file.exists(file_name_SPM)){
      # system.time(
      zone_data_SPM <- plyr::ldply(.data = files_SPM, .fun = extract_pixels,
                                   .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      # ) # 57 minutes for all SEXTANT, 3 minutes for all MODIS Bay of Seine
      data.table::fwrite(zone_data_SPM, file_name_SPM); gc()
    } else {
      if(!exists("zone_data_SPM")){
        zone_data_SPM <- data.table::fread(file_name_SPM) |> mutate(date = as.Date(date))
      }
    }
  }
  
  ## TUR
  ### Copy the SPM values as TUR for comparison against in situ values for SEXTANT
  if(sat_name == "SEXTANT"){
    message("Copying TUR from SPM at : ", Sys.time())
    if(!exists("zone_data_TUR")){
      if(!exists("zone_data_SPM")){
        zone_data_SPM <- data.table::fread(file_name_SPM) # Intentionally SPM
      }
      zone_data_TUR <- zone_data_SPM |> mutate(variable = "TUR", date = as.Date(date))
    }
  } else if(exists("files_TUR")){
    message("Started TUR extraction at : ", Sys.time())
    file_name_TUR <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_TUR.csv")
    if(!file.exists(file_name_TUR)){
      # system.time(
      zone_data_TUR <- plyr::ldply(.data = files_TUR, .fun = extract_pixels,
                                   .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      # ) # 5 seconds for 10 turns, 52 minutes for all SEXTANT, 3 minutes for MODIS Bay of Seine
      data.table::fwrite(zone_data_TUR, file_name_TUR); gc()
    } else {
      if(!exists("zone_data_TUR")){
        zone_data_TUR <- data.table::fread(file_name_TUR) |> mutate(date = as.Date(date))
      }
    }
  }

  ## CDOM
  if(exists("files_CDOM")){
    message("Started CDOM extraction at : ", Sys.time())
    file_name_CDOM <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_CDOM.csv")
    if(!file.exists(file_name_CDOM)){
      zone_data_CDOM <- plyr::ldply(.data = files_CDOM, .fun = extract_pixels,
                                    .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      data.table::fwrite(zone_data_CDOM, file_name_CDOM); gc()
    } else {
      if(!exists("zone_data_SEXTANT_CDOM")){
        zone_data_CDOM <- data.table::fread(file_name_CDOM) |> mutate(date = as.Date(date))
      }
    }
  }
  
  ## RRS
  if(exists("files_RRS")){
    message("Started RRS extraction at : ", Sys.time())
    file_name_RRS <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_RRS.csv")
    if(!file.exists(file_name_RRS)){
      zone_data_RRS <- plyr::ldply(.data = files_RRS, .fun = extract_pixels,
                                   .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      data.table::fwrite(zone_data_RRS, file_name_RRS); gc()
    } else {
      if(!exists("zone_data_SEXTANT_RRS")){
        zone_data_RRS <- data.table::fread(file_name_RRS) |> mutate(date = as.Date(date))
      }
    }
  }
  
  ## SST
  if(exists("files_SST")){
    message("Started SST extraction at : ", Sys.time())
    file_name_SST <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_SST.csv")
    if(!file.exists(file_name_SST)){
      zone_data_SST <- plyr::ldply(.data = files_SST, .fun = extract_pixels,
                                   .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      data.table::fwrite(zone_data_SST, file_name_SST); gc()
    } else {
      if(!exists("zone_data_SEXTANT_SST")){
        zone_data_SST <- data.table::fread(file_name_SST) |> mutate(date = as.Date(date))
      }
    }
  }
  
  # Combine all data for median calculations
  if(sat_name == "SEXTANT"){
    zone_median_base <- bind_rows(zone_data_CHL, zone_data_SPM, zone_data_TUR)
    rm(zone_data_SPM, zone_data_TUR, zone_data_CHL); gc()
  } else if(sat_name == "MODIS"){
    zone_median_base <- bind_rows(zone_data_CHL, zone_data_SPM, zone_data_TUR,
                                  zone_data_CDOM, zone_data_RRS, zone_data_SST)
    rm(zone_data_CHL, zone_data_SPM, zone_data_TUR,
       zone_data_CDOM, zone_data_RRS, zone_data_SST); gc()
  } else {
    zone_median_base <- bind_rows(zone_data_CHL, zone_data_SPM, zone_data_TUR,
                                  zone_data_CDOM, zone_data_RRS)
    rm(zone_data_CHL, zone_data_SPM, zone_data_TUR,
       zone_data_CDOM, zone_data_RRS); gc()
  }
  
  # Create median value time series
  file_name_median_all <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,"_all.csv")
  file_name_median_small <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,"_small.csv")
  if(!file.exists(file_name_median_all)){
    message("Started median all calculations at : ", Sys.time())
    
    # Create medians etc. from all pixels
    zone_median_all <- zone_median_base |> 
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
    message("Started median small calculations at : ", Sys.time())
    
    if(sat_name == "SEXTANT"){
      slice_n <- 1
    } else{
      slice_n <- 9
    }
    zone_median_small <- zone_median_base |> 
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
  
  # Clean and exit
  rm(zone_median_base, file_name_median_all, file_name_median_small); gc()
  message("Finished ", file_stub," at : ", Sys.time())
}


# Loading -----------------------------------------------------------------

# Load time series of plume values
load_plume_ts <- function(zone){
  file_name <- paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")
  suppressMessages(
    df_plume <- read_csv(file_name) |> 
      dplyr::select(date:confidence_index_in_perc) |>
      complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
      dplyr::rename(plume_area = area_of_the_plume_mask_in_km2) |> 
      mutate(plume_area = ifelse(plume_area > 20000, NA, plume_area),) |> 
      zoo::na.trim() |> 
      mutate(zone = zone, .before = "date")
  )
  return(df_plume)
}

# Get the date of the plume data from the file name while loading
load_plume_surface <- function(zone){
  # Detect all plume surface csv files
  plume_files <- dir(paste0("output/REGIONAL_PLUME_DETECTION/",zone,"/SEXTANT/SPM/merged/Standard/PLUME_DETECTION/DAILY"), 
                     pattern = ".csv", recursive = TRUE, full.names = TRUE)
  
  # The function to be ply'd across all files
  load_plume_1 <- function(file_name, zone){
    # Get date
    file_date <- as.Date(gsub(".csv", "", basename(file_name)))
    # Load data and switch columns
    df <- read_csv(file_name) |> 
      mutate(zone = zone,
             date = file_date) |> 
      dplyr::select(zone, date, lon, lat) |> 
      filter(!is.na(lon) & !is.na(lat))
    return(df)
  }
  
  # Load all daily maps into one data.frame
  ## NB: There are a lot of files to load, need some heavy lifting to get it done
  df_plume_surface <- plyr::ldply(plume_files, load_plume_1, .parallel = TRUE, zone = zone)
  message(paste0("Loaded ", nrow(df_plume_surface), " plume surface points for zone ", zone))
  return(df_plume_surface)
}

# Load river flow data
# dir_name <- "data/RIVER_FLOW/BAY_OF_SEINE"
load_river_flow <- function(dir_name){
  
  # Get file list
  files_to_load <- list.files(path = dir_name, pattern = "\\.(txt|dat|csv|ascii)$", full.names = TRUE)
  
  # Load each file into a list of dataframes
  data_list <- lapply(files_to_load, 
                      function(file) {
                        ext <- tolower(tools::file_ext(file))
                        if (ext == "csv") {
                          df <- read.csv(file, header = FALSE)
                        } else if (ext == "ascii") {
                          df <- as.data.frame(read.table(file, header = FALSE))
                        } else {
                          df <- read.delim(file, sep = "", header = FALSE)  # Auto-detect separator
                        }
                        return(df)
                      })
  
  # Name the list elements with file names
  names(data_list) <- basename(files_to_load)
  
  # Process each dataframe based on the zone
  for(key in names(data_list)){
    df_river <- data_list[[key]]
    
    if(grepl("GULF_OF_LION", dir_name)){
      colnames(df_river) <- c("flow", "year", "month", "day")
      df_river$date <- as.Date(
        paste0(
          df_river$year, "-",
          str_pad(df_river$month, 2, pad = "0"), "-",
          str_pad(df_river$day, 2, pad = "0")
        )
      )
    }
    
    if (grepl("BAY_OF_BISCAY|SOUTHERN_BRITTANY|BAY_OF_SEINE", dir_name)) {
      colnames(df_river) <- c("date", "time", "flow")
      df_river$date <- as.Date(df_river$date, format = "%d/%m/%Y")
    }
    
    # Select only 'Flow' and 'Date' columns
    data_list[[key]] <- df_river[, c("flow", "date")]
  }
  
  # Combine all dataframes and aggregate by date
  final_df <- bind_rows(data_list) |> 
    summarise(flow = sum(flow, na.rm = TRUE),
              n_rivers = n(), .by = "date") |> 
    filter(n_rivers == length(data_list))
}

# Load tide gauge data
load_tide_gauge <- function(dir_name){
  
  # Tide gauge files
  tide_files <- dir(dir_name, pattern = ".txt", full.names = TRUE)
  
  # Load all files
  suppressMessages(
    df_tide <- map_dfr(tide_files, read_delim, col_names = c("t", "tide", "source"), skip = 14, delim = ";", col_select = c("t", "tide"))
  )
  df_tide_daily <- df_tide |> 
    mutate(t = as.POSIXct(t, format = "%d/%m/%Y %H:%M:%S"),
           date = as.Date(t)) |> 
    summarise(tide_mean = round(mean(tide, na.rm = TRUE), 2),
              tide_range = max(tide, na.rm = TRUE)-min(tide, na.rm = TRUE), .by = "date")
  return(df_tide_daily)
}

# Load wind data
load_wind_sub <- function(file_name, lon_range, lat_range){
  wind_df <- tidync(file_name) |> 
    hyper_filter(longitude = dplyr::between(longitude, lon_range[1], lon_range[2]),
                 latitude = dplyr::between(latitude, lat_range[1], lat_range[2])) |> 
    hyper_tibble() |> 
    dplyr::rename(u = eastward_wind, v = northward_wind, lon = longitude, lat = latitude) |> 
    mutate(date = as.Date(time)) |> 
    dplyr::select(date, lon, lat, u, v) |> 
    summarise(u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE), .by = "date")
  
  # Remove final day of data
  ## it is an artefact from creating daily integrals from hourly data
  final_date <- max(wind_df$date)
  wind_df <- filter(wind_df, date != final_date)
  return(wind_df)
}

# Load ROFI surface NetcDF
load_ROFI <- function(file_name){
  # Get zone from file name
  if(grepl("Seine", file_name)){
    zone <- "BAY_OF_SEINE"
  } else if (grepl("Gironde", file_name)){
    zone <- "BAY_OF_BISCAY"
  } else if (grepl("Loire", file_name)){
    zone <- "SOUTHERN_BRITTANY"
  } else if (grepl("Rhone", file_name)){
    zone <- "GULF_OF_LION"
  } else {
    stop("Zone not recognised from ROFI file name")
  }
  df_ROFI <- tidync(file_name) |> 
    tidync::hyper_tibble() |> 
    mutate(zone = zone,
           date = as.Date(parse_date_time(time_counter, "Ymd HMS"))) |> 
    dplyr::select(zone, date, ROFI_surface) |> 
    summarise(ROFI_surface = mean(ROFI_surface, na.rm = TRUE), .by = c("zone", "date"))
  return(df_ROFI)
}

# Load X11 results
load_X11 <- function(zone_name, plume = TRUE){
  if(plume){
    file_name <- paste0("output/FIXED_THRESHOLD/GULF_OF_LION/X11_ANALYSIS/area_of_the_plume_mask_in_km2/SEXTANT_merged_Standard_WEEKLY.csv")
  } else {
    file_name <- paste0("output/FIXED_THRESHOLD/GULF_OF_LION/X11_ANALYSIS/river_flow/River_flow___WEEKLY.csv")
  }
  suppressMessages(
    df_X11 <- read_csv(file_name) |> 
      dplyr::rename(date = dates) |> 
      complete(date = seq(min(date), max(date), by = "day")) |> 
      mutate(plume_area = ifelse(plume_area > 20000, NA, plume_area)) |> 
      zoo::na.trim() |> 
      mutate(zone = zone, .before = "date")
  )
  return(df_X11)
}


# Statistics --------------------------------------------------------------

# Check for leap year
is_leap_year <- function(year) {
  (year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)
}

# Function to adjust DOY for non-leap years
adjust_doy <- function(year, doy) {
  if (!is_leap_year(year) && doy > 59) {
    doy + 1
  } else {
    doy
  }
}

# Lagged correlations
# NB: The x vector is lagged, meaning it is effectively pushed forward in time
lagged_correlation <- function(x, y, max_lag) {
  df_cor <- tibble(
    lag = 0:max_lag,
    cor = map_dbl(0:max_lag, ~ cor(y, lag(x, .), use = "pairwise.complete.obs"))
  )
  return(df_cor)
}

# Calculate STL
stl_single <- function(x_col, out_col, start_date, ts_freq = 365.25){
  
  # Create ts object and calculate stl
  ts_x <- ts(zoo::na.approx(x_col), frequency = ts_freq, start = c(year(start_date), quarter(start_date)))
  stl_x <- stl(ts_x, s.window = "periodic", t.window = 11, inner = 5, outer = 0)
  #, robust = TRUE) # rather allow outliers to influence results
  
  # Add NA to end if necessary
  if(length(stl_x$time.series[,2]) != length(x_col)){

    # Create a matrix of NA values
    na_matrix <- matrix(NA, nrow = length(x_col)-length(stl_x$time.series[,2]), ncol = ncol(stl_x$time.series))
    
    # Append the NA matrix to the original matrix
    stl_x$time.series <- rbind(stl_x$time.series, na_matrix)
  }
  
  # Decide which column to take
  if(out_col == "seas"){
    return(as.vector(stl_x$time.series[,1]))
  } else if(out_col == "inter"){
    return(as.vector(stl_x$time.series[,2]))
  } else if(out_col == "remain"){
    return(as.vector(stl_x$time.series[,3]))
  } else {
    stop("'out_col' not recognised")
  }
}

# Load panache time series based on river mouth name
plume_clim_calc <- function(zone){
  file_name <- paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")
  suppressMessages(
    df_plume <- read_csv(file_name) |> 
      dplyr::select(date:confidence_index_in_perc) |>
      complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
      dplyr::rename(plume_area = area_of_the_plume_mask_in_km2) |> 
      mutate(plume_area = ifelse(plume_area > 20000, NA, plume_area),) |> 
      zoo::na.trim() |> 
      mutate(zone = zone, .before = "date")
  )
  df_clim <- ts2clm(df_plume, x = date, y = plume_area, climatologyPeriod = c(min(df_plume$date), max(df_plume$date))) |> 
    dplyr::select(zone, date, doy, plume_area, seas, thresh) |> 
    dplyr::rename(plume_seas = seas, plume_thresh = thresh)
}

# The full suite of stats to calculate
compute_stats <- function(x_vec, y_vec){
  
  if(!is.numeric(x_vec)) stop("x_vec is not numeric")
  if(!is.numeric(y_vec)) stop("y_vec is not numeric")
  
  if(length(x_vec) < 3){
    return(data.frame(row.names = NULL,
                      n = length(x_vec),
                      Slope = NA,
                      Slope_log = NA,
                      RMSE = NA,
                      MSA = NA,
                      MAPE = NA,
                      Bias = NA,
                      Error = NA))
  }
  
  # Calculate RMSE (Root Mean Square Error)
  rmse <- sqrt(mean((y_vec - x_vec)^2, na.rm = TRUE))
  
  # Calculate MAPE (Mean Absolute Percentage Error)
  mape <- mean(abs((y_vec - x_vec) / x_vec), na.rm = TRUE) * 100
  
  # Calculate MSA (Mean Squared Adjustment)
  msa <- mean(abs(y_vec - x_vec), na.rm = TRUE)
  
  # Calculate linear slope
  lin_fit <- lm(y_vec ~ x_vec)
  slope <- coef(lin_fit)[2]
  
  # Calculate log-log linear slope
  log_lin_fit <- lm(log10(y_vec) ~ log10(x_vec))
  log_slope <- coef(log_lin_fit)[2]
  
  # Calculate Bias
  log_ratio <- log10(y_vec / x_vec)
  log_ratio_median <- median(log_ratio, na.rm = TRUE)
  bias_perc <- 100 * (sign(log_ratio_median) * (10^abs(log_ratio_median) - 1))
  
  # Calculate error
  log_ratio_median_abs <- median(abs(log_ratio), na.rm = TRUE)
  error_perc <- 100 * (10^log_ratio_median_abs - 1)
  
  # Combine int data.frame and exit
  return(data.frame(row.names = NULL,
                    n = length(x_vec),
                    Slope = round(slope, 2),
                    Slope_log = round(log_slope, 2),
                    RMSE = round(rmse, 6),
                    MSA = round(msa, 6),
                    MAPE = round(mape, 2),
                    Bias = round(bias_perc, 2),
                    Error = round(error_perc, 2)))
}


# Plotting ----------------------------------------------------------------

# Create pretty plot labels from zone values
make_pretty_title <- function(df){
  df <- df |> 
    mutate(plot_title = case_when(zone == "BAY_OF_SEINE" ~ "Bay of Seine",
                                  zone == "SOUTHERN_BRITTANY" ~ "Southern Brittany",
                                  zone == "BAY_OF_BISCAY" ~ "Bay of Biscay",
                                  zone == "GULF_OF_LION" ~ "Gulf of Lion"), .after = "zone") |> 
    mutate(plot_title = factor(plot_title, levels = c("Bay of Seine", "Southern Brittany", "Bay of Biscay", "Gulf of Lion")))
  return(df)
}

# Scale one value to another for tidier double-y-axis plots
sec_axis_adjustement_factors <- function(var_to_scale, var_ref){
  
  index_to_keep <- which(is.finite(var_ref))
  var_ref <- var_ref[index_to_keep]
  
  index_to_keep <- which(is.finite(var_to_scale))
  var_to_scale <- var_to_scale[index_to_keep]
  
  max_var_to_scale <- max(var_to_scale, na.rm = T) 
  min_var_to_scale <- min(var_to_scale, na.rm = T) 
  max_var_ref <- max(var_ref, na.rm = T) 
  min_var_ref <- min(var_ref, na.rm = T) 
  
  diff_to_scale <- max_var_to_scale - min_var_to_scale
  diff_to_scale <- ifelse(diff_to_scale == 0, 1 , diff_to_scale)
  diff_ref <- max_var_ref - min_var_ref
  diff <- diff_ref / diff_to_scale
  
  adjust <- (max_var_ref - max_var_to_scale*diff) 
  
  return(data.frame(diff = diff, adjust = adjust, operation = "scaled var = (var_to_scale * diff) + adjust",
                    trans_axis_operation = "var_to_scale = {scaled_var - adjust} / diff)"))
}

# Consistent theme for project
ggplot_theme <-   function(){
  theme(text = element_text(size = 35, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 55),
        plot.subtitle = element_text(hjust = 0.5, size = 30),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(angle = 0),
        axis.ticks.length = unit(.25, "cm"))
}

# Convenience wrapper for saving png files
save_plot_as_png <- function(plot, name = c(), width = 14, height = 8.27, path, res = 150){
  
  graphics.off()
  
  if (name %>% length() == 1) {
    if (dir.exists(path) == FALSE) {dir.create(path, recursive = TRUE)}
    path <- file.path(path, paste(name, ".png", sep = ""))
  } else {
    path <- paste(path, ".png", sep = "")
  }
  
  if (grepl(pattern = ".png.png", path)) {path <- path %>% gsub(pattern = ".png.png", replacement = ".png", x = .)}
  
  png(path, width = width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
  
}

# Comparison plots
# df <- df_pretty; var_1 <- "plume_inter"; var_2 <- "flow_inter"
# df <- df_seas; var_1 <- "plume_seas"; var_2 <- "flow_seas"
# colour_1 <- "brown"; colour_2 <- "blue"; label_1 <- "Plume area (km^2)"; label_2 <- "River flow (m^3 s-1)"; file_stub <- "comparison_plume_flow_inter"
comparison_plot <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2){
  
  if(grepl("seas", var_1)){
    df_sub <- df[,c("plot_title", "month", var_1, var_2)]
    colnames(df_sub) <- c("plot_title", "month", "var_1", "var_2")
  } else {
    df_sub <- df[,c("plot_title", "date", var_1, var_2)]
    colnames(df_sub) <- c("plot_title", "date", "var_1", "var_2")
  }
  
  # Plot base
  if(grepl("seas", var_1)){
    
    # Scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_sub$var_2, df_sub$var_1)
    df_scaling <- summarise(df_sub, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_sub, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    
    # Get range for ribbon plot
    df_scale_sub <- df_scale |> 
      summarise(var_1_min = min(var_1, na.rm = TRUE),
                var_1_mean = mean(var_1, na.rm = TRUE),
                var_1_max = max(var_1, na.rm = TRUE),
                var_2_min = min(var_2_scaled, na.rm = TRUE),
                var_2_mean = mean(var_2_scaled, na.rm = TRUE),
                var_2_max = max(var_2_scaled, na.rm = TRUE), .by = c("plot_title", "month")) |> 
      mutate(month_int = as.integer(month))
    
    # Plot them
    pl_base <- ggplot(data = df_scale_sub, aes(x = month_int)) + 
      # Var 1
      geom_ribbon(aes(ymin = var_1_min, ymax = var_1_max), fill = colour_1, alpha = 0.3) +
      geom_path(aes(y = var_1_mean), color = colour_1, linewidth = 2) +
      # Var 2
      geom_ribbon(aes(ymin = var_2_min, ymax = var_2_max), fill = colour_2, alpha = 0.3) +
      geom_path(aes(y = var_2_mean), color = colour_2, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_continuous(expand = c(0, 0), breaks = 1:12, labels = month.abb)
  } else {
    
    # Perform rolling mean
    df_roll_mean <- df_sub |>
      mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
      # mutate(date = round_date(date, unit = "months")) |>
      filter(date >= min(df$date)) |>
      group_by(plot_title, date) |>
      summarise(var_1 = mean(var_1, na.rm = TRUE),
                var_2 = mean(var_2, na.rm = TRUE), .groups = "keep") |>
      group_by(plot_title) |>
      mutate(var_1 = roll_mean(var_1, n = 48, fill = NA, align = "center"),
             var_2 = roll_mean(var_2, n = 48, fill = NA, align = "center")) |>
      ungroup()
    
    # Then get the scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_roll_mean$var_2, df_roll_mean$var_1)
    df_scaling <- summarise(df_roll_mean, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_roll_mean, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    unique_years <- df_scale$date |> year() |> unique()
    
    pl_base <- ggplot(data = df_scale) +
      # Var 1 data
      geom_point(aes(x = date, y = var_1), color = colour_1) +
      geom_path(aes(x = date, y = var_1), color = colour_1) +
      # geom_smooth(aes(x = date, y = var_1), method = "lm", se = FALSE, color = colour_1) +
      # Var 2 data
      geom_point(aes(x = date, y = var_2_scaled), color = colour_2) +
      geom_path(aes(x = date, y = var_2_scaled), color = colour_2) +
      # geom_smooth(aes(x = date, y = var_2_scaled), method = "lm", se = FALSE, color = colour_2) +
      # Facet
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      # X-axis labels
      scale_x_date(name = "", expand = c(0, 0),
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) 
  }
  
  # Finish up the comparison plot
  pl_comp <- pl_base +
    # Y-axis labels
    scale_y_continuous(name = label_1,
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
                                           name = label_2)) +
    labs( x = NULL) +
    # Extra bits
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y.left = element_text(color = colour_1),
          axis.ticks.y.left = element_line(color = colour_1),
          axis.line.y.left = element_line(color = colour_1),
          axis.title.y.left = element_text(color = colour_1, margin = unit(c(0, 7.5, 0, 0), "mm")),
          axis.text.y.right = element_text(color = colour_2),
          axis.ticks.y.right = element_line(color = colour_2),
          axis.line.y.right = element_line(color = colour_2),
          axis.title.y.right = element_text(color = colour_2, margin = unit(c(0, 0, 0, 7.5), "mm")),
          panel.border = element_rect(linetype = "solid", fill = NA))
  return(pl_comp)
}

# Convenience wrapper to run and save comparison plots
# NB: this is hard coded to work with four plots
comparison_plot_save <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2, file_stub){
  comp_list <- plyr::dlply(df, c("zone"), comparison_plot, var_1 = var_1, var_2 = var_2, 
                           colour_1 = colour_1, colour_2 = colour_2, label_1 = label_1, label_2 = label_2)
  comp_fig <- comp_list[[2]] + comp_list[[4]] + comp_list[[1]] + comp_list[[3]] + plot_layout(ncol = 1, axes = "collect")
  ggsave(filename = paste0("figures/",file_stub,".png"), plot = comp_fig, width = 24, height = 24, dpi = 300)
}

# It does what it says on the tin
var_labels <- function(var_name){
  
  if(grepl("CHL|RRS", var_name)){
    var_colour <- "green4"
    if(grepl("CHL", var_name)){
      axis_limits <- c(10^-2, 10^2)
      unit <- expression(mg~m^-3)
    } else {
      unit <- expression(sr-1)
      axis_limits <- c(10^-5, 10^-1)
    }
  } else if(grepl("TEMP|SST", var_name)) {
    unit <- expression('"°C"')
    axis_limits <- c(3, 30)
    var_colour <- "orange4"
    # NB: Must be run after SST because the T for turbidity is an issue
  } else if(grepl("SPM|TUR|T|POC|CDOM", var_name)){
    axis_limits <- c(10^-2, 10^3)
    var_colour<- "brown4"
    if(grepl("TUR|T", var_name)){
      unit <- expression(NTU)
    } else if(var_name == "CDOM"){
      unit <- expression(m-1)
    } else if(var_name == "POC"){
      unit <- expression(mg~m^-3)
    } else {
      unit <- expression(g~m^-3)
    }
  } else {
    stop(paste("Could not find ", var_name))
  }
  
  # List and exit
  to_return <- list("unit" = unit,
                    "var_colour" = var_colour, 
                    "axis_limits" = axis_limits)
  return(to_return)
}

# Convenience colour wrapper
colours_of_stations <- function(){
  
  # colour_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
  #                  "Point C" = mako(n = 1,begin = 0.85,end = 0.85), 
  #                  "Luc-sur-Mer" = mako(n = 1,begin = 0.90,end = 0.9), 
  #                  "Smile" = mako(n = 1,begin = 0.95,end = 0.95),
  #                  
  #                  "Bizeux" = viridis(n = 1,begin = 1,end = 1), # Bretagne
  #                  "Le Buron" = viridis(n = 1,begin = 0.925,end = 0.925), 
  #                  "Cézembre" = viridis(n = 1,begin = 0.85,end = 0.85), 
  #                  "Estacade" = viridis(n = 1,begin = 0.775,end = 0.775), 
  #                  "Astan" = viridis(n = 1,begin = 0.7,end = 0.7), 
  #                  "Portzic" = viridis(n = 1,begin = 0.625,end = 0.625), 
  #                  
  #                  "Antioche" = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
  #                  "pk 86" = plasma(n = 1,begin = 0.1,end = 0.1), 
  #                  "pk 52" = plasma(n = 1,begin = 0.15,end = 0.15),
  #                  "pk 30" = plasma(n = 1,begin = 0.2,end = 0.2),
  #                  "Comprian" = plasma(n = 1,begin = 0.25,end = 0.25), 
  #                  "Eyrac" = plasma(n = 1,begin = 0.3,end = 0.3), 
  #                  "Bouee 13" = plasma(n = 1,begin = 0.35,end = 0.35), 
  #                  
  #                  "Sola" = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
  #                  "Sete" = rocket(n = 1,begin = 0.85,end = 0.85), 
  #                  "Frioul" = rocket(n = 1,begin = 0.90,end = 0.90),
  #                  "Point B" = rocket(n = 1,begin = 0.95,end = 0.95)
  # ) 
  
  # colour_values = c('Manche orientale - Mer du Nord' = mako(n = 1,begin = 0.8,end = 0.8), # Manche
  #                  'Baie de Seine' = mako(n = 1,begin = 0.95,end = 0.95),
  # 
  #                  'Manche occidentale' = viridis(n = 1,begin = 1,end = 1), # Bretagne
  #                  'Bretagne Sud' = viridis(n = 1,begin = 0.625,end = 0.625),
  # 
  #                  'Pays de la Loire - Pertuis' = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
  #                  'Sud Golfe de Gascogne' = plasma(n = 1,begin = 0.35,end = 0.35),
  # 
  #                  'Golfe du Lion' = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
  #                  'Mer ligurienne - Corse' = rocket(n = 1,begin = 0.95,end = 0.95)
  # )
  
  # colour_values = c('BAY OF SEINE' = viridis::mako(n = 1, begin = 0.8,end = 0.8),
  #                  'SOUTHERN BRITTANY' = viridis::viridis(n = 1, begin = 0.9,end = 0.9),
  #                  'GULF OF BISCAY' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
  #                  'GULF OF LION' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  colour_values = c('Bay of Seine' = viridis::mako(n = 1, begin = 0.8,end = 0.8),
                    'S. Brittany' = viridis::viridis(n = 1, begin = 0.9,end = 0.9),
                    'Bay of Biscay' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
                    'Gulf of Lion' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  
  return(colour_values)
}

# Plot linear trends and stats for matched data
# var_sub = "TEMP"; df = zone_all_in_situ_monthly; df_stats = zone_all_monthly_lm
# var_sub = "CHLA"
validation_lm_plots <- function(var_sub, sat_name, median_base, df, df_stats){
  
  # Get y-axis label and units
  if(var_sub == "TEMP"){
    y_lab <- "Temperature (°C)"
    y_lim <- c(4, 32)
    var_sat <- "SST"
  } else if(var_sub == "CHLA"){
    y_lab <- "Chlorophyll-a (mg m-3)"
    y_lim <- c(0, 20)
    var_sat <- "CHL"
  } else if(var_sub == "TUR"){
    y_lab <- "Turbidity (NTU)"
    y_lim <- c(0, 20)
  } else if(var_sub == "SPM"){
    y_lab <- "SPM (g m-3)"
    y_lim <- c(0, 20)
  } else {
    # Not worried about other variables for the moment
    y_lab <- NA
    y_lim <- c(0, 20)
  }
  
  # Filter and pivot data
  df_var_sub <- df |> 
    filter(variable == var_sub) |> 
    filter(variable_sat != "NRRS555") |> 
    pivot_longer(cols = value_in_situ:value_satellite) |> 
    mutate(name = gsub("value|_", "", name),
           name = gsub("insitu", "in situ", name)) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  
  # Filter stats labels
  df_stats_var_sub <- df_stats |> 
    filter(variable == var_sub) |> 
    filter(variable_sat != "NRRS555") |> 
    mutate(y = max(y_lim)-2) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  
  # Create plot
  plot_zone_var_TS <- ggplot(data = df_var_sub, aes(x = date, y = value)) +
    geom_point(aes(y = value, colour = name, shape = source), alpha = 0.4) +
    geom_line(aes(y = value, colour = name, linetype = source), alpha = 0.4) +
    geom_smooth(aes(colour = name, linetype = source), method = "lm", se = FALSE) +
    # In situ values
    geom_label(data = filter(df_stats_var_sub, source == "REPHY"), 
               aes(x = as.Date("2005-01-01"), y = y, label = paste0(source," : ", round(slope_is, 2)," / year"), vjust = 1.0)) +
    geom_label(data = filter(df_stats_var_sub, source == "SOMLIT"), 
               aes(x = as.Date("2005-01-01"), y = y, label = paste0(source," : ", round(slope_is, 2)," / year"), vjust = -0.5)) +
    # Satellite values
    geom_label(data = filter(df_stats_var_sub, source == "REPHY"), colour = "red",
               aes(x = as.Date("2020-01-01"), y = y, label = paste0(source," : ", round(slope_sat, 2)," / year"), vjust = 1.0)) +
    geom_label(data = filter(df_stats_var_sub, source == "SOMLIT"),  colour = "red",
               aes(x = as.Date("2020-01-01"), y = y, label = paste0(source," : ", round(slope_sat, 2)," / year"), vjust = -0.5)) +
    facet_wrap(~zone_pretty,nrow = 2) +
    scale_colour_manual(values = c("black", "red")) +
    # scale_y_continuous(limits = c(min(df_var_sub$value, na.rm = TRUE), max(df_var_sub$value, na.rm = TRUE)*1.1)) +
    coord_cartesian(ylim = y_lim) +
    labs(title = paste(sat_name, df_var_sub$variable_sat[1],"vs in situ",df_var_sub$variable[1]),
         y = y_lab, x = NULL, colour = "Type", shape = "Source", linetype = "Source") +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 25, face = "bold"), 
          strip.text = element_text(size = 20),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 20),
          axis.title = element_text(size = 23),
          axis.text = element_text(size = 20),
          legend.position = "bottom")
  # plot_zone_var_TS
  ggsave(paste0("figures/validation/ts/ts_",sat_name,"_",var_sub,"_",median_base,".png"), plot_zone_var_TS, width = 14, height = 8)
  return()
}

# The figure code wrapper
# var_name = "SPM"; match_up_df = zone_in_situ_SEXTANT; match_up_stats = stats_SEXTANT; sat_name = "SEXTANT"
# sat_name = sat_name; match_up_df = zone_all_in_situ; match_up_stats = zone_all_in_situ_stats; var_name = unique(zone_all_in_situ_stats$variabl_combi)[1]
validation_plots <- function(var_name, sat_name, median_base, match_up_df, match_up_stats){
  
  # Split the variable names
  var_name_is <- str_split(var_name, "_")[[1]][1]
  var_name_sat <- str_split(var_name, "_")[[1]][2]
  
  # Get axis labels
  plot_meta_is <- var_labels(var_name_is)
  plot_meta_sat <- var_labels(var_name_sat)

  # Get 1:1 line limits
  identity_line <- data.frame(x = plot_meta_is$axis_limits, 
                              y = plot_meta_sat$axis_limits)

  # Subset datasets for chosen variable
  match_up_df_var <- match_up_df |> 
    filter(variable == var_name_is,
           variable_sat == var_name_sat) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  match_up_stats_var <- match_up_stats |> 
    filter(variable == var_name_is,
           variable_sat == var_name_sat)
  
  # Get colour values
  colour_values <- colours_of_stations()
  
  # Get stats to plot
  if(var_name_sat == 'SST'){
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope
  } else {
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope_log
  }
  
  # Create title and subtitle
  plot_title = paste(sat_name, var_name_sat, "vs.", "in situ", var_name_is)
  
  # Create the plot
  scatterplot <- ggplot(data = match_up_df_var, aes(x = value_in_situ, y = value_satellite)) + 
    
    geom_point(aes(colour = zone_pretty, shape = source), size = 6, show.legend = TRUE) + 
    # TODO: Add 2:1 and 1:2 dashed lines
    geom_line(data = identity_line, aes(x = x, y = y), linetype = "dashed", show.legend = FALSE) +
    
    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('In~situ~measurements~(', plot_meta_is$unit, ')'))) +
    
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Satellite~estimates~(', plot_meta_sat$unit, ')'))) + 
    
    coord_equal(xlim = plot_meta_is$axis_limits, 
                ylim = plot_meta_sat$axis_limits) +
    
    annotate(geom = 'text', x = plot_meta_is$axis_limits[1], y = plot_meta_sat$axis_limits[2], 
             hjust = 0, vjust = 1, color = "black", size = 12.5,
             # label = paste('Error = ', round(ifelse(Error_value |> is.numeric(), Error_value, NA), 1), "%\n",
             #               'Bias = ', round(ifelse(Bias_value |> is.numeric(), Bias_value, NA), 1), " %\n",
             #               # TODO: Change this to show Slope_log or Slope depending on the variable tested (e.g. SST or not)
             #               # 'R²_log = ', round(statistics_values$r2_log, 2),"\n",
             #               'Slope = ', round(ifelse(Slope_value |> is.numeric(), Slope_value, NA), 2),"\n",
             #               'n = ', nrow(match_up_df_var), sep = "")) +
             label = paste("Slope = ", round(ifelse(Slope_value |> is.numeric(), Slope_value, NA), 2),"\n",
                           # 'R² = ', round(statistics_values$r2_log, 2),"\n",
                           "n = ", nrow(match_up_df_var), sep = "")) +
    
    labs(title = plot_title) +
    
    scale_color_manual(name = "zone", values = colour_values, drop = FALSE) +
    
    # scale_linetype_manual(values = c("Identity line" = "dashed",
    #                                  "Linear regression" = "solid"), name = "") +
    
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 10), order = 1),
           shape = guide_legend(ncol = 1, override.aes = list(size = 10), order = 2),
           # linetype = guide_legend(override.aes = list(color = c("black"), shape = c(NA), linetype = c("dashed")), ncol = 2, order = 2)
    ) +
    
    ggplot_theme() + 
    
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 30),
          legend.margin = margin(5, 10, 5, 5),
          plot.subtitle = element_text(hjust = 0.5, color = "black", face = "bold.italic"),
          plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # scatterplot
  
  if(var_name_sat == "SST"){ 
    scatterplot <- scatterplot +  
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      scale_x_continuous(name = parse(text = paste0('In~situ~measurements~(', plot_meta_is$unit, ')'))) +
      scale_y_continuous(name = parse(text = paste0('Satellite~estimates~(', plot_meta_sat$unit, ')')))
  } else {
    scatterplot <- scatterplot + 
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      annotation_logticks()
  }
  # scatterplot
  
  # Add histograms to x and y axes
  scatterplot_with_side_hist <- ggMarginal(scatterplot, type = "histogram", groupFill = TRUE, alpha = 1)
  # scatterplot_with_side_hist
  ggsave(paste0("figures/validation/scatterplot/",sat_name,"_",var_name_sat,"_vs_in_situ_",var_name_is,"_",median_base,".png"), 
         plot = scatterplot_with_side_hist, height = 16, width = 16, bg = "white")
  
  # Barplot of frequency per year
  bar_plot_freq_per_year <- match_up_df_var |> 
    mutate(Year = year(date)) |> 
    dplyr::count(Year) |> 
    ggplot(aes(x = Year)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_is$var_colour, linewidth = 1.5) + 
    scale_x_continuous(breaks = seq(1998, 2025, by = 3), name = "", labels = function(x) substring(x, 3, 4)) +
    scale_y_continuous(expand = c(0,0), name = "n per year") +
    labs(title = plot_title) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_year
  ggsave(paste0("figures/validation/barplot/annual_",sat_name,"_",var_name_sat,"_",median_base,".png"), 
         plot = bar_plot_freq_per_year, height = 10, width = 14, bg = "white")
  
  # Barplot of monthly counts
  bar_plot_freq_per_month <- match_up_df_var |> 
    mutate(Month = month(date)) |> 
    dplyr::count(Month) |> 
    ggplot(aes(x = Month)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_is$var_colour, linewidth = 2) + 
    scale_x_continuous(breaks = 1:12, name = "", labels = function(x) month.abb[x]) +
    scale_y_continuous(expand = c(0,0), name = "n per month") +
    coord_cartesian(xlim = c(1,12)) +
    labs(title = plot_title) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_month
  ggsave(paste0("figures/validation/barplot/monthly_",sat_name,"_",var_name_sat,"_",median_base,".png"), 
         plot = bar_plot_freq_per_month, height = 10, width = 14, bg = "white")
  return()
}

# Run all validation stats and produce the plots
# sat_name = "SEXTANT"; median_base = "small"
# sat_name = "MODIS"; median_base = "all"
# sat_name = "OLCI-A"; median_base = "all"
validate_sensor <- function(sat_name, median_base){
  
  # Get the pixel cutoff based on product
  if(median_base == "small"){
    if(sat_name == "SEXTANT"){
      pixel_n_cut <- 1
    } else {
      pixel_n_cut <- 3
    }
  } else {
    if(sat_name == "SEXTANT"){
      pixel_n_cut <- 3
    } else {
      pixel_n_cut <- 9
    }
  }
  
  # Load prepped data
  sat_files <- dir("output/MATCH_UP_DATA/FRANCE", full.names = TRUE, 
                   pattern = paste0("zone_median_",sat_name))
  sat_files <- sat_files[grepl(paste0("_",median_base), sat_files)]
  zone_median <- map_dfr(sat_files, data.table::fread) |> mutate(date = as.Date(date)) |> 
    # Remove all rows that are below the pixel cutoff and CV cutoff of 20%
    mutate(sd = case_when(n == 1 ~ 0, TRUE ~ sd), # Necessary for CV for 1 pixel count matchups for SEXTANT 'small'
           cv = sd / median) |> 
    filter(n >= pixel_n_cut, cv <= 0.20) |> 
    # Complete all dates
    complete(date = seq(min(date), max(date), by = "day"), fill = list(median = NA), 
             nesting(zone, source, site, variable))
  
  # Make variable name conversions as necessary
  if(sat_name == "SEXTANT"){
    # NB: SEXTANT data had already been tweaked to match the in situ
    # This is rolling that back somewhat and should rather be addressed earlier on in the process
    zone_median <- zone_median |> 
      mutate(variable_is = variable) |> 
      mutate(variable = case_when(variable == "CHLA" ~ "CHL",
                                  variable == "TUR" ~ "SPM", # To see what happens
                                  TRUE ~ variable))

  } else {
    zone_median <- zone_median |> 
      mutate(variable_is = case_when(variable == "CHL" ~ "CHLA",
                                     variable == "SST" ~ "TEMP",
                                     variable == "T" ~ "TUR",
                                     variable == "CDOM" ~ "POC", # Just to see what happens
                                     variable == "NRRS555" ~ "CHLA", # Just to see what happens
                                     variable == "NRRS560" ~ "CHLA", # Just to see what happens
                                     TRUE ~ variable))
  }
  
  # Load in situ data and complete the date column
  zone_data_in_situ <- read_csv("data/INSITU_data/zone_data_in_situ.csv") |> 
    dplyr::select(-lon, -lat) |> 
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA), 
             nesting(zone, zone_pretty, source, site, variable))
  
  # Combine extracted sat data with in situ
  zone_all_in_situ_base <- zone_data_in_situ |> 
    left_join(zone_median, by = c("zone", "source", "site", "date", "variable" = "variable_is")) |>
    dplyr::rename(value_in_situ = value, value_satellite = median, variable_sat = variable.y) |> 
    filter(!is.na(variable_sat)) |> 
    mutate(season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter", month(date) %in% 3:5  ~ "Spring",
      month(date) %in% 6:8  ~ "Summer", month(date) %in% 9:11 ~ "Autumn"), .after = "date") |> 
    dplyr::select(zone, dplyr::everything())

  # Create monthly average TS for lm analysis
  zone_all_in_situ_monthly <- zone_all_in_situ_base |> 
    mutate(date = floor_date(date, unit = "month")) |> 
    summarise(value_in_situ = mean(value_in_situ, na.rm = TRUE),
              value_satellite = mean(value_satellite, na.rm = TRUE),
              .by = c("zone", "zone_pretty", "source", "season", "date", "variable", "variable_sat"))
  
  # Calculate linear model stats to look at change over time
  zone_in_situ_monthly_lm <- zone_all_in_situ_monthly |> 
    filter(value_in_situ > 0) |> 
    group_by(zone, zone_pretty, source, variable, variable_sat) |> 
    do(broom::tidy(lm(value_in_situ ~ date, data = .))) |> 
    filter(term == "date") |> 
    dplyr::rename(slope_is = estimate, p_is = p.value) |> 
    dplyr::select(zone:variable_sat, slope_is, p_is)
  zone_sat_monthly_lm <- zone_all_in_situ_monthly |> 
    filter(value_satellite > 0) |> 
    group_by(zone, zone_pretty, source, variable, variable_sat) |> 
    do(broom::tidy(lm(value_satellite ~ date, data = .))) |> 
    filter(term == "date") |> 
    dplyr::rename(slope_sat = estimate, p_sat = p.value) |> 
    dplyr::select(zone:variable_sat, slope_sat, p_sat)
  
  # Combine results
  zone_all_monthly_lm <- left_join(zone_in_situ_monthly_lm, zone_sat_monthly_lm,
                                   by = join_by(zone, zone_pretty, source, variable, variable_sat)) |> 
    # Convert to values / year
    mutate(slope_is = slope_is*365.25, slope_sat = slope_sat*365.25)
  
  # Save statistics
  write_csv(zone_all_monthly_lm, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/",sat_name,"_lm_stats_",median_base,".csv"))
  
  # Plot the linear TS matchups
  plyr::l_ply(unique(zone_all_in_situ_monthly$variable), validation_lm_plots, 
              sat_name = sat_name, median_base = median_base,
              df = zone_all_in_situ_monthly, df_stats = zone_all_monthly_lm)
  
  # Remove any missing values for stats matchups
  zone_all_in_situ <- zone_all_in_situ_base |> 
    filter(value_in_situ > 0, value_satellite > 0)
  
  # Create big grid of sites to look for obvious outliers
  # ggplot(data = zone_in_situ, aes(x = value_in_situ, y = value_satellite)) +
  #   geom_point(aes(colour = source, shape = variable)) +
  #   facet_wrap(~site, scales = "free")
  
  # TODO: Think of a way to optimise this
  # Stats for all groups together by variable
  zone_in_situ_stats_01 <- zone_all_in_situ |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all groups together by variable and season
  zone_in_situ_stats_02 <- zone_all_in_situ |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones by variable
  zone_in_situ_stats_03 <- zone_all_in_situ |> mutate(source = "ALL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones by variable and season
  zone_in_situ_stats_04 <- zone_all_in_situ |> mutate(source = "ALL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sources by variable
  zone_in_situ_stats_05 <- zone_all_in_situ |> mutate(zone = "GLOBAL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sources by variable and season
  zone_in_situ_stats_06 <- zone_all_in_situ |> mutate(zone = "GLOBAL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones and sources by variable
  zone_in_situ_stats_07 <- zone_all_in_situ |> mutate(site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones and sources by variable and season
  zone_in_situ_stats_08 <- zone_all_in_situ |> mutate(site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sites by variable
  zone_in_situ_stats_09 <- zone_all_in_situ |> mutate(site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sites by variable and season
  zone_in_situ_stats_10 <- zone_all_in_situ |> mutate(site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for each site by variable
  zone_in_situ_stats_11 <- zone_all_in_situ |> mutate(season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for each site by variable and season
  zone_in_situ_stats_12 <- zone_all_in_situ |>
    summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Bind all together
  zone_all_in_situ_stats <- bind_rows(zone_in_situ_stats_01, zone_in_situ_stats_02, zone_in_situ_stats_03,
                                      zone_in_situ_stats_04, zone_in_situ_stats_05, zone_in_situ_stats_06,
                                      zone_in_situ_stats_07, zone_in_situ_stats_08, zone_in_situ_stats_09,
                                      zone_in_situ_stats_10, zone_in_situ_stats_11, zone_in_situ_stats_12) |> 
    mutate(sensor = sat_name, .before = "zone")
  
  # Save results
  write_csv(zone_all_in_situ_stats, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/",sat_name,"_stats_",median_base,".csv"))
  rm(zone_in_situ_stats_01, zone_in_situ_stats_02, zone_in_situ_stats_03,
     zone_in_situ_stats_04, zone_in_situ_stats_05, zone_in_situ_stats_06,
     zone_in_situ_stats_07, zone_in_situ_stats_08, zone_in_situ_stats_09,
     zone_in_situ_stats_10, zone_in_situ_stats_11, zone_in_situ_stats_12); gc()
  
  # Run all of the plots per variable pairing
  zone_all_in_situ_stats$variabl_combi <- paste0(zone_all_in_situ_stats$variable,"_",zone_all_in_situ_stats$variable_sat)
  plyr::l_ply(unique(zone_all_in_situ_stats$variabl_combi), validation_plots, .parallel = TRUE,
              sat_name = sat_name, median_base = median_base,
              match_up_df = zone_all_in_situ, match_up_stats = zone_all_in_situ_stats)
}


# Tables ---------------------------------------------------------------

# Create pretty tables
# file_path = "output/MATCH_UP_DATA/FRANCE/STATISTICS/"; sat_name = "SEXTANT"
# TODO: This will need to be modified to work with other satellite products
validation_tables <- function(file_path, sat_name) {
  
  # Load all output stats
  # TODO: Clean up zone names and order by latitude
  stat_files <- map_dfr(dir(file_path, pattern = sat_name, full.names = TRUE), read_csv) |> 
    dplyr::select(zone, Bias, Error, n, variable) |> 
    mutate_at(vars(-zone, -variable), ~ round(., 1)) |> 
    mutate(var_label = paste("When compared with in situ", ifelse(variable == 'TUR', 'TURB', variable))) |>
    dplyr::select(-variable)
  # dplyr::rename(`Bias (%)` = Bias, `Error (%)` = Error)
  
  desired_colnames <- names(stat_files) |> 
    str_replace_all("Bias", "Bias (%)") |> 
    str_replace_all("Error", "Error (%)") |> 
    str_remove_all("var_label")
  names(desired_colnames) <- names(stat_files)
  
  # Create the SPM/TUR table
  # TODO: Wrap this into a function to be called per variable
  table_SPM_TUR <- stat_files |>
    filter(!grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite SPM', subtitle = 'Compared with SPM and Turbidity in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    # Pretty tabs
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3), # A bit more padding for summaries
                row_group.padding = px(4)) |> # And even more for our groups
    # More styling
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Create the CHLA table
  table_CHLA <- stat_files |>
    filter(grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite CHL', subtitle = 'Compared with CHLA in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3),
                row_group.padding = px(4)) |> 
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Save and exit
  ## SPM/TUR
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.html"), inline_css = TRUE)
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.png"), expand = 10)
  ## CHLA
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.html"), inline_css = TRUE)
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.png"), expand = 10)
}

