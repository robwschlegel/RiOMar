# func/util.R
# The storage point for many functions re-used by other scripts
# TODO: test that python calls to other scripts correctly source this one


# Meta-data ---------------------------------------------------------------

# In the future this will be taken from define_parameters() in func/util.py
river_mouths <- data.frame(row_name = 1:4,
                           mouth_name = c("Seine", "Gironde", "Loire", "Grand Rhone"),
                           mouth_lon = c(0.145, -1.05, -2.10, 4.83),
                           mouth_lat = c(49.43, 45.59, 47.29, 43.41))


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
    # hyper_filter(longitude = longitude >= lon_range[1] & longitude <= lon_range[2],
    #              latitude = latitude >= lat_range[1] & latitude <= lat_range[2]) |> 
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
    dplyr::select(zone, date, ROFI_surface)
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
stl_single <- function(x_col, out_col, start_date, ts_freq = 365){
  
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

