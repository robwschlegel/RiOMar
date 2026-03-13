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


# Pixels ------------------------------------------------------------------

# Simple wrapper to extract start and end times of an ODATIS-MR file
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_CDOM-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20230424__FRANCE_03_MOD_T-FNU-NS_DAY_00.nc"
get_start_end_time <- function(file_name){
  df_info <- ncdump::NetCDF(file_name)$attribute[1]$global |> 
    dplyr::select(start_time, end_time) |> 
    mutate(start_time = as.POSIXct(gsub("T|Z", " ", start_time), format = "%Y%m%d %H%M%S"),
           end_time = as.POSIXct(gsub("T|Z", " ", end_time), format = "%Y%m%d %H%M%S"))
  return(df_info)
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
  
  # Extract data and exit
  df_res <- df |> 
    mutate(date = nc_date,
           variable = var_col_name,
           value = raster::extract(raster(file_name, varname = var_nc_name), cell_numbers)) |> 
    dplyr::select(-cell_numbers)
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
    files_SPM <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "SPM")
    files_CHL <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "CHL")
    files_TUR <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "T-FNU")
    files_CDOM <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "CDOM")
    files_RRS <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "RRS")
  } 
  
  # Load SST explicitly for MODIS
  if(sat_name == "MODIS") files_SST <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = "SST")
  
  # Extract all available data per variable
  # TODO: Consider a different way of doing this. I think it could be more streamlined.
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
        zone_data_SPM <- data.table::fread(file_name_SPM)
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
      zone_data_TUR <- zone_data_SPM |> mutate(variable = "TUR")
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
        zone_data_TUR <- data.table::fread(file_name_TUR)
      }
    }
  }
  
  ## CHL
  if(exists("files_CHL")){
    message("Started CHL extraction at : ", Sys.time())
    file_name_CHL <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,"_CHL.csv")
    if(!file.exists(file_name_CHL)){
      zone_data_CHL <- plyr::ldply(.data = files_CHL, .fun = extract_pixels,
                                   .parallel = TRUE, .paropts = list(.inorder = FALSE), df = zone_pixels)
      data.table::fwrite(zone_data_CHL, file_name_CHL); gc()
    } else {
      if(!exists("zone_data_SEXTANT_CHL")){
        zone_data_CHL <- data.table::fread(file_name_CHL)
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
        zone_data_CDOM <- data.table::fread(file_name_CDOM)
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
        zone_data_RRS <- data.table::fread(file_name_RRS)
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
        zone_data_SST <- data.table::fread(file_name_SST)
      }
    }
  }
  
  # Create median value time series
  ## Combine all data as required
  file_name_median <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,".csv")
  message("Started median calculations at : ", Sys.time())
  if(!file.exists(file_name_median)){
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
    
    # Create medians etc.
    zone_median <- zone_median_base |> 
      filter(value > 0) |>
      summarise(median = median(value, na.rm = TRUE), 
                mean = mean(value, na.rm = TRUE),
                sd = sd(value, na.rm = TRUE),
                n = n(), .by = c("zone", "source", "site", "date", "variable"))
    data.table::fwrite(zone_median, file_name_median)
    rm(zone_median_base, zone_median); gc()
  } #else {
  #   if(!exists("zone_median")){
  #     zone_median <- data.table::fread(file_name_median)
  #   }
  # }
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
sec_axis_adjustement_factors <- function(var_to_scale, var_ref) {
  
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
ggplot_theme <-   function() {
  theme(text = element_text(size = 35, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 55),
        plot.subtitle = element_text(hjust = 0.5, size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(angle = 0),
        axis.ticks.length = unit(.25, "cm"))
}

# Convenience wrapper for saving png files
save_plot_as_png <- function(plot, name = c(), width = 14, height = 8.27, path, res = 150) {
  
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

