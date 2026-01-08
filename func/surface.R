# func/surface.R
# Analyses specifically of surface data


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(patchwork)
library(doParallel); doParallel::registerDoParallel(cores = 14)

# Zones
zones <- c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION")

# Location of daily time series per site
# "output/FIXED_THRESHOLD/GULF_OF_LION/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv"

# Location of daily time series per site per year
# "output/FIXED_THRESHOLD/GULF_OF_LION/SEXTANT/SPM/merged/Standard/PLUME_DETECTION/DAILY/1998/Results.csv"

# Location of daily plume maps
# "output/REGIONAL_PLUME_DETECTION/GULF_OF_LION/SEXTANT/SPM/merged/Standard/PLUME_DETECTION/DAILY/1998/MAPS/1998-01-01.csv"


# Functions ---------------------------------------------------------------

# Load all surface data for a RiOMar and plot the surface relation to drivers
surface_plot <- function(zone){
  
  # Determine meta-data based on zone
  if(zone == "BAY_OF_SEINE"){
    gauge = "LE_HAVRE"; mouth_name = "Seine"
    mouth_lon = 0.145; mouth_lat = 49.43
    lon_W = -1.0; lon_E = 0.5; lat_N = 0.5; lat_S = -0.2
  } else if(zone == "BAY_OF_BISCAY"){
    gauge = "PORT-BLOC"; mouth_name = "Gironde"
    mouth_lon = -1.05; mouth_lat = 45.59
    lon_W = -3.0; lon_E = 0.1; lat_N = 1; lat_S = -0.5
  } else if(zone == "SOUTHERN_BRITTANY"){
    gauge = "SAINT-NAZAIRE"; mouth_name = "Loire"
    mouth_lon = -2.10; mouth_lat = 47.29
    lon_W = -3.0; lon_E = 0.2; lat_N = 0.5; lat_S = -1.0
  } else if(zone == "GULF_OF_LION"){
    gauge = "MARSEILLE"; mouth_name = "Grand Rhone"
    mouth_lon = 4.83; mouth_lat = 43.41
    lon_W = -2.0; lon_E = 3.0; lat_N = 0; lat_S = -2.5
  } else {
    stop("Zone not recognised.")
  }
  
  # Set plume dir
  plume_dir <- paste0("output/REGIONAL_PLUME_DETECTION/",zone,"/SEXTANT/SPM/merged/Standard/PLUME_DETECTION/DAILY")
  
  # Detect all csv files
  plume_files <- dir(plume_dir, pattern = ".csv", recursive = TRUE, full.names = TRUE)
  
  # Load all daily maps into one data.frame
  ## NB: There are a lot of files to load, need some heavy lifting to get it done
  df_plume <- plyr::ldply(plume_files, load_plume_surface, .parallel = TRUE)
  
  # Load river flow data
  df_river_flow <- load_river_flow(paste0("data/RIVER_FLOW/",zone))
  
  # Load tide data
  df_tide <- load_tide_gauge(paste0("data/TIDES/",gauge))
  
  # Load and subset wind data to + 0.5 N and E from mouth
  lon_round <- plyr::round_any(mouth_lon, 0.5); lat_round <- plyr::round_any(mouth_lat, 0.5)
  lon_range <- c(lon_round-0.5, lon_round+0.5); lat_range <- c(lat_round-0.5, lat_round+0.5)
  lon_range_wide <- c(lon_round+lon_W, lon_round+lon_E); lat_range_wide <- c(lat_round+lat_S, lat_round+lat_N)
  wind_files <- dir(paste0("~/pCloudDrive/data/WIND/",zone), pattern = "_daily_", full.names = TRUE)
  df_wind <- purrr::map_dfr(wind_files, load_wind_sub, lon_range, lat_range)
  
  # Determine simple on-/off-shore index based on coastal direction (based on river mouth name)
  # TODO: Think of a more sophisticated way to do this
  if(zone %in% c("BAY_OF_BISCAY", "SOUTHERN_BRITTANY")){
    df_wind_direction <- df_wind |> 
      mutate(direction = ifelse(u < 0, "off", "on"))
  } else if (zone == "BAY_OF_SEINE"){
    df_wind_direction <- df_wind |> 
      mutate(direction = ifelse(v > 0, "off", "on"))
  } else if (zone == "GULF_OF_LION"){
    df_wind_direction <- df_wind |> 
      mutate(direction = ifelse(v < 0, "off", "on"))
  }
  
  # Calculate wind speed and direction
  # NB: wind_dir is where the wind is coming from, not going to
  df_wind_full <- df_wind_direction |> 
    mutate(wind_spd = round(sqrt(u^2 + v^2), 2),
           wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))
  
  # Join all by date
  df_full <- left_join(df_plume, df_river_flow, by = "date") |> 
    left_join(df_tide, by = "date") |> 
    left_join(df_wind_full, by = "date")
  
  # Get stats per pixel
  suppressWarnings(
  df_pixel <- summarise(df_full,
                        count = n(),
                        flow_min = min(flow, na.rm = TRUE),
                        flow_mean = mean(flow, na.rm = TRUE),
                        flow_max = max(flow, na.rm = TRUE),
                        tide_range_min = min(tide_range, na.rm = TRUE),
                        tide_range_mean = mean(tide_range, na.rm = TRUE),
                        tide_range_max = max(tide_range, na.rm = TRUE),
                        wind_spd_min = min(wind_spd, na.rm = TRUE),
                        wind_spd_mean = mean(wind_spd, na.rm = TRUE),
                        wind_spd_max = max(wind_spd, na.rm = TRUE),
                        count_on = sum(direction == "on", na.rm = TRUE),
                        count_off = sum(direction == "off", na.rm = TRUE),
                        .by = c("lon", "lat")) |> 
    mutate(prop_n = count/length(unique(df_full$date)),
           prop_on = count_on/count,
           prop_off = count_off/count,
           across(everything(), ~ ifelse(is.finite(.), ., NA)))
  )
  
  # Plot pixel count plus proportion count
  plot_count <- ggplot(df_pixel, aes(x = lon, y = lat)) +
    annotation_borders(regions = "France", fill = "grey70") +
    geom_tile(aes(fill = log10(count))) +
    scale_fill_viridis_c() +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Plume count (n)") +
    theme_bw() +
    theme(legend.position = "bottom")
    # theme(legend.position = "inside",
    #       legend.position.inside = c(0.7, 0.9), # Probably change these per site via first series of logic gates
    #       legend.direction = "horizontal")
  # plot_count
  
  # Plot pixel count plus proportion count
  plot_count_prop <- ggplot(df_pixel, aes(x = lon, y = lat)) +
    annotation_borders(regions = "France", fill = "grey70") +
    geom_tile(aes(fill = prop_n)) +
    scale_fill_viridis_c() +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Plume count proportion (n/all_days)") +
    theme_bw() +
    theme(legend.position = "bottom")
    # theme(legend.position = "inside",
    #       legend.position.inside = c(0.7, 0.9), # Probably change these per site via first series of logic gates
    #       legend.direction = "horizontal")
  # plot_count_prop
  
  # Plot river flow
  plot_flow <- df_pixel |> 
    select(lon, lat, flow_min, flow_mean, flow_max) |> 
    pivot_longer(cols = c(flow_min, flow_mean, flow_max), names_to = "var") |> 
    mutate(var = factor(var, levels = c("flow_min", "flow_mean", "flow_max"))) |> 
    ggplot(aes(x = lon, y = lat)) +
    # annotation_borders() +
    # geom_sf(data = ne_countries(scale = "medium", returnclass = "sf")) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis_c(option = "A") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "River flow range") +
    facet_wrap(~var) +
    theme_bw() +
    theme(legend.position = "bottom")
  # plot_flow
  
  # Plot wind speed
  plot_wind_spd <- df_pixel |> 
    select(lon, lat, wind_spd_min, wind_spd_mean, wind_spd_max) |> 
    pivot_longer(cols = c(wind_spd_min, wind_spd_mean, wind_spd_max), names_to = "var") |> 
    mutate(var = factor(var, levels = c("wind_spd_min", "wind_spd_mean", "wind_spd_max"))) |> 
    ggplot(aes(x = lon, y = lat)) +
    # annotation_borders() +
    # geom_sf(data = ne_countries(scale = "medium", returnclass = "sf")) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis_c(option = "C") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, , title = "Wind speed range") +
    facet_wrap(~var) +
    theme_bw() +
    theme(legend.position = "bottom")
  # plot_wind_spd
  
  # Plot wind direction
  plot_wind_dir <- df_pixel |> 
    select(lon, lat, prop_on, prop_off) |> 
    pivot_longer(cols = c(prop_on, prop_off), names_to = "var") |> 
    mutate(var = factor(var, levels = c("prop_on", "prop_off"))) |> 
    ggplot(aes(x = lon, y = lat)) +
    # annotation_borders() +
    # geom_sf(data = ne_countries(scale = "medium", returnclass = "sf")) +
    geom_tile(aes(fill = value)) +
    # scale_fill_viridis_c(option = "C") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Proportion of on- or off-shore winds") +
    facet_wrap(~var) +
    theme_bw() +
    theme(legend.position = "bottom")
  # plot_wind_dir
  
  # Plot tide range
  plot_tide_range <- df_pixel |> 
    select(lon, lat, tide_range_min, tide_range_mean, tide_range_max) |> 
    pivot_longer(cols = c(tide_range_min, tide_range_mean, tide_range_max), names_to = "var") |> 
    mutate(var = factor(var, levels = c("tide_range_min", "tide_range_mean", "tide_range_max"))) |> 
    ggplot(aes(x = lon, y = lat)) +
    # annotation_borders() +
    # geom_sf(data = ne_countries(scale = "medium", returnclass = "sf")) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis_c(option = "E") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Tidal range range") +
    facet_wrap(~var) +
    theme_bw() +
    theme(legend.position = "bottom")
  # plot_tide_range
  
  # Put it all together
  plot_multi <- (plot_count | plot_count_prop) /
    (plot_flow | plot_tide_range) / 
    (plot_wind_spd | plot_wind_dir)
  ggsave(filename = paste0("figures/surface_stats_",zone,".png"), plot = plot_multi, height = 14, width = 26)
  return()
}


# Run ---------------------------------------------------------------------

plyr::l_ply(zones, surface_plot)

