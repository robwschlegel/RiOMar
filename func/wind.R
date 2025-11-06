# func/wind.R
# Comparisons of wind forcing against plume size


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(doParallel); doParallel::registerDoParallel(cores = 4)


# Meta-data ---------------------------------------------------------------

# In the future this will be taken from define_parameters() in func/util.py
river_mouths <- data.frame(row_name = 1:4,
                           mouth_name = c("Seine", "Gironde", "Loire", "Grand Rhone"),
                           mouth_lon = c(0.145, -1.05, -2.10, 4.83),
                           mouth_lat = c(49.43, 45.59, 47.29, 43.41))


# Functions ---------------------------------------------------------------

# Function for loading subsets of wind NetCDF files
load_wind_sub <- function(file_name, lon_range, lat_range){
  wind_df <- tidync(file_name) |> 
    hyper_filter(longitude = dplyr::between(longitude, lon_range[1], lon_range[2]),
                 latitude = dplyr::between(latitude, lat_range[1], lat_range[2])) |> 
    # hyper_filter(longitude = longitude >= lon_range[1] & longitude <= lon_range[2],
    #              latitude = latitude >= lat_range[1] & latitude <= lat_range[2]) |> 
    hyper_tibble() |> 
    dplyr::rename(u = eastward_wind, v = northward_wind, lon = longitude, lat = latitude) |> 
    mutate(t = as.Date(time)) |> 
    dplyr::select(t, lon, lat, u, v)
  
  # Remove final day of data
  ## it is an artefact from creating daily integrals from hourly data
  final_date <- max(wind_df$t)
  wind_df <- filter(wind_df, t != final_date)
  return(wind_df)
}

# Caluclate wind stats and create plots
# mouth_info <- river_mouths[1,]
spatial_wind_calc <- function(mouth_info){
  
  # Get zone name from river mouth
  if(mouth_info$mouth_name == "Seine"){
    zone <- "BAY_OF_SEINE"
  } else if(mouth_info$mouth_name == "Gironde"){
    zone <- "BAY_OF_BISCAY"
  } else if(mouth_info$mouth_name == "Loire"){
    zone <- "SOUTHERN_BRITTANY"
  } else if(mouth_info$mouth_name == "Grand Rhone"){
    zone <- "GULF_OF_LION"
  } else {
    stop("River mouth not recognised.")
  }
  
  # Load and subset wind data to + 0.5 N and E from mouth
  lon_round <- plyr::round_any(mouth_info$mouth_lon, 0.5); lat_round <- plyr::round_any(mouth_info$mouth_lat, 0.5)
  lon_range <- c(lon_round-0.5, lon_round+0.5); lat_range <- c(lat_round-0.5, lat_round+0.5); 
  wind_files <- dir(paste0("~/pCloudDrive/data/WIND/",zone), pattern = "_daily_", full.names = TRUE)
  wind_df <- purrr::map_dfr(wind_files, load_wind_sub, lon_range, lat_range)
  
  # Calculate spatial average of wind vectors by day
  wind_df_mean <- wind_df |> 
    summarise(u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE), .by = "t")
  
  # Determine simple upwelling/downwelling index based on coastal direction (based on river mouth name)
  # TODO: Think of a more sophisticated way to do this
  if(zone %in% c("BAY_OF_BISCAY", "SOUTHERN_BRITTANY")){
    wind_df_updown <- wind_df_mean |> 
      mutate(welling = ifelse(u < 0, "up", "down"))
  } else if (zone == "BAY_OF_SEINE"){
    wind_df_updown <- wind_df_mean |> 
      mutate(welling = ifelse(v > 0, "up", "down"))
  } else if (zone == "GULF_OF_LION"){
    wind_df_updown <- wind_df_mean |> 
      mutate(welling = ifelse(v < 0, "up", "down"))
  }
    
  # Calculate wind speed and direction
  # NB: wind_dir is where the wind is coming from, not going to
  wind_df_full <- wind_df_updown |> 
    mutate(wind_spd = round(sqrt(u^2 + v^2), 2),
           wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))
  
  # Load panache time series based on river mouth name
  plume_daily <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")) |> 
    dplyr::select(date:path_to_file) |> dplyr::select(-path_to_file) |> 
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA))
  
  # Eventually it would be ideal to load river flow
  
  # Combine dataframes for further analyses
  wind_plume_df <- left_join(plume_daily, wind_df_full, join_by(date == t)) |> 
    zoo::na.trim()
  # wind_plume_df_month <- wind_plume_df |> 
  #   mutate()
  
  # Create time series objects for stl
  ts_plume <- ts(zoo::na.approx(wind_plume_df$area_of_the_plume_mask_in_km2), frequency = 365, start = c(year(min(wind_plume_df$date)), 1))
  stl_plume <- stl(ts_plume, s.window = "periodic")
  ts_wind <- ts(zoo::na.approx(wind_plume_df$wind_spd), frequency = 365, start = c(year(min(wind_plume_df$date)), 1))
  stl_wind <- stl(ts_wind, s.window = "periodic")
  
  # Add trend elements back into dataframe for further stats
  wind_plume_df$wind_stl <- as.vector(stl_wind$time.series[,2])
  wind_plume_df$plume_stl <- as.vector(stl_plume$time.series[,2])
  
  # Compare panache size against wind speed
  # Two separate comparisons based on upwelling or downwelling times
  wind_plume_stats_all <- wind_plume_df |> 
    mutate(welling = "all") |> 
    summarise(r = cor(wind_spd, area_of_the_plume_mask_in_km2, use = "pairwise.complete.obs"), .by = "welling")
  wind_plume_stats <- wind_plume_df |> 
    summarise(r = cor(wind_spd, area_of_the_plume_mask_in_km2, use = "pairwise.complete.obs"), .by = "welling") |> 
    rbind(wind_plume_stats_all)
  
  # Lagged correlations
  wind_plume_lag_cor <- tibble(
    lag = 0:30,
    cor = map_dbl(0:30, ~ cor(wind_plume_df$wind_spd, lag(wind_plume_df$area_of_the_plume_mask_in_km2, .), use = "complete.obs"))
  )
  
  # Plot wind speed and up/downwelling times
  wind_plot <- ggplot(wind_plume_df, aes(x = date, y = wind_spd)) +
    geom_vline(aes(xintercept = date, colour = welling), show.legend = FALSE) +
    geom_line() +
    labs(y = "wind speed (m s-1)", x = NULL) +
    scale_x_date(expand = 0) +
    # guides(colour = guide_legend(override.aes = list(shape = 15))) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
   # wind_plot
  
  # Panache size
  panache_plot <- ggplot(wind_plume_df, aes(x = date, y = area_of_the_plume_mask_in_km2)) +
    geom_line() +
    labs(y = "plume area (km^2)", x = NULL) +
    scale_x_date(expand = 0) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # panache_plot
  
  # Plot wind speed and panache size correlation
  wind_plume_cor_plot <- ggplot(wind_plume_df, aes(x = wind_spd, y = area_of_the_plume_mask_in_km2)) + 
    geom_point(aes(colour = welling), alpha = 0.7) +
    geom_smooth(method = "lm") +
    geom_smooth(method = "lm", aes(colour = welling)) +
    labs(y = "plume area (km^2)", x = "wind speed (m s-1)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  # wind_plume_cor_plot
  
  # Plot lag results
  wind_plume_cor_lag_plot <- ggplot(wind_plume_lag_cor, aes(x = lag, y = cor)) +
    geom_point() +
    labs(x = "lag plume after wind (days)", y = "correlation (r)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # wind_plume_cor_lag_plot
  
  # Combine plots and save
  wind_plume_title <- grid::textGrob(paste0(mouth_info$mouth_name," : wind vs plume size"), gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(wind_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
  cor_plot <- ggpubr::ggarrange(wind_plume_cor_plot, wind_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
  full_plot_title <- ggpubr::ggarrange(wind_plume_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
  ggsave(filename = paste0("figures/wind_plume_cor_plot_",mouth_info$mouth_name,".png"), width = 12, height = 6, dpi = 600)

  # Get scaling factor for plotting
  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = wind_plume_df$wind_stl, 
                                                 var_ref = wind_plume_df$plume_stl)
  wind_plume_df <- wind_plume_df |> 
    mutate(wind_scaled = wind_stl * scaling_factor$diff + scaling_factor$adjust)
  unique_years <- wind_plume_df$date |> year() |> unique()
  
  # Plots for STL decomposition
  ggplot(data = wind_plume_df) + 
    # Plume data
    geom_point(aes(x = date, y = plume_stl), color = "brown") + 
    geom_path(aes(x = date, y = plume_stl), color = "brown") + 
    # Wind data
    geom_point(aes(x = date, y = wind_scaled), color = "blue") + 
    geom_path(aes(x = date, y = wind_scaled), color = "blue") + 
    # X-axis labels
    scale_x_date(name = "", 
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "Plume area (km²)",
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
                                           name = "Wind speed (m/s)")) +
    # Extra bits
    labs(title = zone) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y.left = element_text(color = "brown"), 
          axis.ticks.y.left = element_line(color = "brown"),
          axis.line.y.left = element_line(color = "brown"),
          axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
          axis.text.y.right = element_text(color = "blue"), 
          axis.ticks.y.right = element_line(color = "blue"),
          axis.line.y.right = element_line(color = "blue"),
          axis.title.y.right = element_text(color = "blue", margin = unit(c(0, 0, 0, 7.5), "mm")),
          panel.border = element_rect(linetype = "solid", fill = NA))
}


# Run ---------------------------------------------------------------------

plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = spatial_wind_calc)

