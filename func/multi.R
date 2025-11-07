# func/multi.R
# Loads all drivers of plume size, performs stats, plots results

# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(doParallel); doParallel::registerDoParallel(cores = 4)

# Zones
zones <- c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION")


# Functions ---------------------------------------------------------------

# Load all plume and driver data and perform stl
# zone <- zones[1]
multi_stl <- function(zone){
  
  # Determine meta-data based on zone
  if(zone == "BAY_OF_SEINE"){
    gauge = "LE_HAVRE"; mouth_name = "Seine"
    mouth_lon = 0.145; mouth_lat = 49.43
  } else if(zone == "BAY_OF_BISCAY"){
    gauge = "PORT-BLOC"; mouth_name = "Gironde"
    mouth_lon = -1.05; mouth_lat = 45.59
  } else if(zone == "SOUTHERN_BRITTANY"){
    gauge = "SAINT-NAZAIRE"; mouth_name = "Loire"
    mouth_lon = -2.10; mouth_lat = 47.29
  } else if(zone == "GULF_OF_LION"){
    gauge = "MARSEILLE"; mouth_name = "Grand Rhone"
    mouth_lon = 4.83; mouth_lat = 43.41
  } else {
    stop("Zone not recognised.")
  }
  
  # Load panache time series based on river mouth name
  df_plume <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")) |> 
    dplyr::select(date:confidence_index_in_perc) |>
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
    dplyr::rename(plume_area = area_of_the_plume_mask_in_km2) |> 
    zoo::na.trim()
  
  # Standardise last column name so it can be combined across sites
  colnames(df_plume)[ncol(df_plume)] <- "SPM_threshold"
  
  # Load river flow data
  df_river_flow <- load_river_flow(paste0("data/RIVER_FLOW/",zone))

  # Load tide data
  df_tide <- load_tide_gauge(paste0("data/TIDES/",gauge))
  
  # Load and subset wind data to + 0.5 N and E from mouth
  lon_round <- plyr::round_any(mouth_lon, 0.5); lat_round <- plyr::round_any(mouth_lat, 0.5)
  lon_range <- c(lon_round-0.5, lon_round+0.5); lat_range <- c(lat_round-0.5, lat_round+0.5); 
  wind_files <- dir(paste0("~/pCloudDrive/data/WIND/",zone), pattern = "_daily_", full.names = TRUE)
  df_wind <- purrr::map_dfr(wind_files, load_wind_sub, lon_range, lat_range)
  
  # Determine simple upwelling/downwelling index based on coastal direction (based on river mouth name)
  # TODO: Think of a more sophisticated way to do this
  if(zone %in% c("BAY_OF_BISCAY", "SOUTHERN_BRITTANY")){
    df_wind_welling <- df_wind |> 
      mutate(welling = ifelse(u < 0, "up", "down"))
  } else if (zone == "BAY_OF_SEINE"){
    df_wind_welling <- df_wind |> 
      mutate(welling = ifelse(v > 0, "up", "down"))
  } else if (zone == "GULF_OF_LION"){
    df_wind_welling <- df_wind |> 
      mutate(welling = ifelse(v < 0, "up", "down"))
  }
  
  # Calculate wind speed and direction
  # NB: wind_dir is where the wind is coming from, not going to
  df_wind_full <- df_wind_welling |> 
    mutate(wind_spd = round(sqrt(u^2 + v^2), 2),
           wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))
  
  # Combine all dataframes for further stats
  # NB: The trailing NAs are problematic...
  df_all <- left_join(df_plume, df_river_flow, by = "date") |> 
    left_join(df_tide, by = "date") |> 
    left_join(df_wind_full, by = "date") |> 
    mutate(plum_seas = stl_single(plume_area, out_col = "seas", start_date = min(df_plume$date)),
           plum_inter = stl_single(plume_area, out_col = "inter", start_date = min(df_plume$date)),
           flow_seas = stl_single(flow, out_col = "seas", start_date = min(df_plume$date)),
           flow_inter = stl_single(flow, out_col = "inter", start_date = min(df_plume$date)),
           tide_seas = stl_single(tide_range, out_col = "seas", start_date = min(df_plume$date)),
           tide_inter = stl_single(tide_range, out_col = "inter", start_date = min(df_plume$date)),
           wind_seas = stl_single(wind_spd, out_col = "seas", start_date = min(df_plume$date)),
           wind_inter = stl_single(wind_spd, out_col = "inter", start_date = min(df_plume$date))) |> 
    mutate(zone = zone, .before = "date")
  # print(ncol(df_all))

  # Exit
  return(df_all)
}

# Plot the results 
multi_plot <- function(){
  
}


# Run ---------------------------------------------------------------------

# Compute all STL stats and save
stl_all <- plyr::ldply(zones, multi_stl, .parallel = TRUE)
save(stl_all, file = "output/STL/stl_all.RData")

# Create plots

