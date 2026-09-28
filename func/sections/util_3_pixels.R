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
# target_site = zone_site_df[10,]; n_pixels = 9; dist_range = 1; sat_rast = rast_sat
# rm(site_sp, site_buffer, cropped_rast, pixel_coords)
get_pixels <- function(target_site, sat_grid, sat_rast, n_pixels, dist_range = 1){
  
  # Filter the sat_grid by the distance from the target site
  # NB: dist_tange/100*2 is a quick and dirty way to reduce the number of pixels to calculate the distance for
  pixel_coords <- sat_grid |> 
    filter(lon >= target_site$lon - dist_range/100*2, lon <= target_site$lon + dist_range/100*2,
           lat >= target_site$lat - dist_range/100*2, lat <= target_site$lat + dist_range/100*2)
                             
  # Calculate distance of pixels from the target in km
  pixel_coords$dist <- round(distHaversine(pixel_coords, target_site[,c("lon", "lat")])/1000, 2)
  
  # Select the nearest 49 (7x7) or 9 (3x3) pixels depending on product
  pixel_coords <- pixel_coords |> arrange(dist) |> slice_head(n = n_pixels)
  
  # Get the pixel IDs and exit
  pixel_coords$cell_numbers <- as.integer(raster::cellFromXY(sat_rast, xy = pixel_coords))

  # Clean up and exit
  pixel_coords <- pixel_coords |> 
    mutate(zone = target_site$zone, source = target_site$source, site = target_site$site) |> 
    dplyr::select(zone, source, site, lon, lat, dist, cell_numbers)
  return(pixel_coords)
}

# Wrapper for get_pixels() to write the pixels to .csv per sat product
# zone_site_df <- filter(zone_sites, zone == "BAY_OF_BISCAY")
# sat_name <- "MODIS"; var_name_full <- "SPM-G-NS_mean"
# sat_name <- "SEXTANT"; var_name_full <- "SPM"
# nc_file_base <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc"
write_pixel_coords <- function(sat_name, var_name_full){
  
  # Load the in situ site metadata
  zone_sites <- read_csv("metadata/in_situ_site_list.csv", show_col_types = FALSE) |> 
    filter(!is.na(zone))

  # Determine file and variable names
  # Logic gate handles the different dates of the fil used for pixels
  if(sat_name == "SEXTANT"){
    if(var_name_full == "SPM"){
      var_name_stub <- "SPIM"
      sat_var <- "analysed_spim"
    } else if(var_name_full == "CHLA"){
      var_name_stub <- "CHL"
      sat_var <- "analysed_chl_a"
    }
    nc_file_base <- paste0(riomar_data_path("SEXTANT"), "/", var_name_full,
      "/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-",var_name_stub,"-ATL-v01-fv01-OI.nc")
  } else if(sat_name == "MODIS"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/MODIS/L3m_20200101__FRANCE_03_MOD_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "MERIS"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/MERIS/L3m_20090101__FRANCE_03_MER_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "OLCI-A"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/OLCI-A/L3m_20200101__FRANCE_03_OLA_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "OLCI-B"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/OLCI-B/L3m_20200101__FRANCE_03_OLB_",var_name_full,"_DAY_00.nc")
  }
  if(sat_name != "SEXTANT"){
      sat_var <- paste0(var_name_full,"_mean")
  }
  
  # Get number of pixels to extract based on product type
  n_pixels <- ifelse(sat_name == "SEXTANT", 9, 49)

  # Get the grid and raster bases
  grid_sat <- get_sat_grid(nc_file_base)
  rast_sat <- raster::raster(nc_file_base, varname = sat_var)
  
  # Extract all pixels
  plan(multisession, workers = parallel::detectCores() - 4)
  zone_pixels <- future_map_dfr(1:nrow(zone_sites),
                                function(i) get_pixels(target_site = zone_sites[i,], 
                                                       sat_grid = grid_sat, sat_rast = rast_sat, 
                                                       n_pixels = n_pixels, dist_range = 1),
                                .options = furrr_options(seed = TRUE))
  plan(sequential)
  
  # Save and exit
  write_csv(zone_pixels, paste0("metadata/zone_pixels_",sat_name,"_",var_name_full,".csv"))
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
      var_base_name <- "SPM"
      var_nc_name <- "analysed_spim"
    } else if(grepl("CHL", file_name)){
      var_base_name <- "CHLA"
      var_nc_name <- "analysed_chl_a"
    } else {
      stop("File structure not recognised")
    }
  } else {
    var_base_name <- str_split(basename(file_name), "_")[[1]][7]
    # var_col_name <- str_split(var_base_name, "-")[[1]][1]
    var_nc_name <- paste0(var_base_name,"_mean")
  }
  
  # Get date values
  if(grepl("SEXTANT", file_name)){
    nc_date <- as.Date(str_split(basename(file_name), "-")[[1]][1], format = "%Y%m%d")
  } else {
    nc_date <- as.Date(str_split(basename(file_name), "_")[[1]][2], format = "%Y%m%d")
  }
  
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

# Wrapper function to be called per variable
# sat_name <- "SEXTANT"; var_name <- "SPM"
process_pixels <- function(sat_name, var_name){
  
  # Load zone pixels
  file_stub <- paste0(sat_name,"_",var_name)
  
  # Get file pathways
  if(sat_name == "SEXTANT"){
    files_path <- riomar_data_path("SEXTANT", var_name)
    var_name_file <- ifelse(var_name == "SPM", "SPIM", ifelse(var_name == "CHLA", "CHL"))
  } else {
    files_path <- file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR", sat_name, zone_name)
    var_name_file <- var_name
  }
  files_var <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = var_name_file)

  # Filter out .png files from SEXTANT folders
  if(sat_name == "SEXTANT"){
    files_var <- files_var[!grepl(".png", files_var)]
  }

  # Extract data for all files and save
  if(length(files_var) > 0){

    file_name_var <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,".csv")

    # Resume support: pCloud Drive's virtual filesystem periodically detaches
    # under sustained reads, killing long extraction runs partway through.
    # Skip files whose date is already saved in file_name_var, and append in
    # batches so a crash only loses the in-progress batch, not the whole run.
    # NB: extract_pixels() catches a failed file read and fills that file's
    # rows with value = NA rather than raising, so a date where *every* row
    # is NA means the read failed (not just normal cloud/land masking, which
    # only ever hits some of the ~2477 pixel/site rows for a date) and must
    # be re-extracted, not skipped.
    dates_done <- character(0)
    if(file.exists(file_name_var)){
      cached <- data.table::fread(file_name_var, select = c("date", "value"))
      dates_done <- cached[, .(all_na = all(is.na(value))), by = date][all_na == FALSE, as.character(date)]
    }
    file_dates <- if(sat_name == "SEXTANT"){
      purrr::map_chr(files_var, ~ str_split(basename(.x), "-")[[1]][1])
    } else {
      purrr::map_chr(files_var, ~ str_split(basename(.x), "_")[[1]][2])
    }
    file_dates <- as.character(as.Date(file_dates, format = "%Y%m%d"))
    files_remaining <- files_var[!(file_dates %in% dates_done)]

    if(length(files_remaining) > 0){

      message("Started ",var_name," extraction at : ", Sys.time(), " (",
              length(files_remaining)," of ",length(files_var)," files remaining)")
      zone_pixels <- read_csv(paste0("metadata/zone_pixels_",file_stub,".csv"), show_col_types = FALSE)

      # pCloud Drive's virtual filesystem detaches periodically regardless of
      # worker count (observed at 10, 8, and 2) -- resumable batching handles
      # that cheaply now, so favour throughput between crashes instead
      plan(multisession, workers = 8)
      batch_size <- 200
      batches <- split(files_remaining, ceiling(seq_along(files_remaining) / batch_size))
      for(batch in batches){
        zone_data_batch <- future_map_dfr(batch, extract_pixels, df = zone_pixels,
                                          .options = furrr_options(seed = TRUE))
        data.table::fwrite(zone_data_batch, file_name_var, append = file.exists(file_name_var))
        message("  Saved ",var_name," batch of ",length(batch)," files at : ", Sys.time())
      }
      plan(sequential)
    }

    # Create median value time series
    file_name_median_all <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,"_all.csv")
    file_name_median_small <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,"_small.csv")
    
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
                  n = n(), 
                  .by = c("zone", "source", "site", "date", "variable"))
      data.table::fwrite(zone_median_all, file_name_median_all)
      rm(zone_data_var, file_name_median_all, zone_median_all); gc()
    }

    # Create medians etc. from 'small' pixels
    if(!file.exists(file_name_median_small)){
      # Load data from .csv
      if(!exists("zone_data_var")){
        message("Loading ",var_name," extraction at : ", Sys.time())
        zone_data_var <- data.table::fread(file_name_var)
      }
      message("Started median small calculations at : ", Sys.time())
      if(sat_name == "SEXTANT"){
        slice_n <- 1
      } else{
        slice_n <- 9
      }
      zone_median_small <- zone_data_var |> 
        filter(value > 0) |> # NB: This removes NA pixels before selecting the nearest pixel, which may be incorrect
        group_by(zone, source, site, date, variable) |> 
        arrange(dist) |> 
        slice_head(n = slice_n) |> 
        ungroup() |> 
        summarise(median = median(value, na.rm = TRUE), 
                  mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE),
                  n = n(), 
                  .by = c("zone", "source", "site", "date", "variable"))
      data.table::fwrite(zone_median_small, file_name_median_small)
      rm(zone_data_var, file_name_var, file_name_median_small, zone_median_small); gc()
    }

    # Cleanup and exit
    rm(file_stub, zone_pixels)
  }
}

# Function that extracts all of the data identified in the product/zone lon/lat metadata files
# sat_name <- "SEXTANT"
# sat_name <- "MODIS"; zone_name <- "SOUTHERN_BRITTANY"
extract_pixels_all <- function(sat_name, zone_name = NULL){#, overwrite = FALSE){
  
  message("Started run on ", sat_name, " : ", Sys.time())
  
  # Get file pathways
  if(sat_name == "SEXTANT"){
    process_pixels(sat_name, "SPM")
    process_pixels(sat_name, "CHLA")
  } else {

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
  }
  
  # Clean and exit
  message("Finished ", sat_name," at : ", Sys.time())
}


