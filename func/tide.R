# func/tide.R
# Comparisons of tides against plume size


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(tidync)
library(doParallel); doParallel::registerDoParallel(cores = 4)


# Meta-data ---------------------------------------------------------------

# In the future this will be taken from define_parameters() in func/util.py
river_mouths <- data.frame(row_name = 1:4,
                           mouth_name = c("Seine", "Gironde", "Loire", "Grand Rhone"),
                           mouth_lon = c(0.145, -1.05, -2.10, 4.83),
                           mouth_lat = c(49.43, 45.59, 47.29, 43.41))


# Functions ---------------------------------------------------------------

# Function for loading annual tidal range data files
# file_name <- "data/TIDES/LE_HAVRE/4_1938.txt"
# dir_name <- "data/TIDES/LE_HAVRE"
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

# Calculate relationship between tides and panache size
# mouth_info <- river_mouths[1,]
tide_calc <- function(mouth_info){
  
  # Get zone name from river mouth
  if(mouth_info$mouth_name == "Seine"){
    gauge <- "LE_HAVRE"
    zone <- "BAY_OF_SEINE"
  } else if(mouth_info$mouth_name == "Gironde"){
    gauge <- "PORT-BLOC"
    zone <- "BAY_OF_BISCAY"
  } else if(mouth_info$mouth_name == "Loire"){
    gauge <- "SAINT-NAZAIRE"
    zone <- "SOUTHERN_BRITTANY"
  } else if(mouth_info$mouth_name == "Grand Rhone"){
    gauge <- "MARSEILLE"
    zone <- "GULF_OF_LION"
  } else {
    stop("River mouth not recognised.")
  }
  
  # Load tide data
  tide_df <- load_tide_gauge(paste0("data/TIDES/",gauge))
  
  # Load panache time series based on river mouth name
  plume_daily <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")) |> 
    dplyr::select(date:path_to_file) |> dplyr::select(-path_to_file)
  
  # Combine
  tide_plume_df <- left_join(plume_daily, tide_df, join_by(date))
  
  # Compare panache size against wind speed
  # Two separate comparisons based on upwelling or downwelling times
  tide_plume_stats_all <- tide_plume_df |> 
    summarise(r = cor(tide_range, area_of_the_plume_mask_in_km2, use = "pairwise.complete.obs"))
  
  # Lagged correlations
  tide_plume_lag_cor <- tibble(
    lag = 0:30,
    cor = map_dbl(0:30, ~ cor(tide_plume_df$tide_range, lag(tide_plume_df$area_of_the_plume_mask_in_km2, .), use = "complete.obs"))
  )
  
  # Plot wind speed and up/downwelling times
  tide_plot <- ggplot(tide_plume_df, aes(x = date, y = tide_mean)) +
    geom_ribbon(aes(ymin = (tide_mean-(tide_range/2)), ymax = tide_mean+tide_range/2)) +
    geom_line() +
    labs(y = "tidal range (m)", x = NULL) +
    scale_x_date(expand = 0) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  tide_plot
  
  # Panache size
  panache_plot <- ggplot(tide_plume_df, aes(x = date, y = area_of_the_plume_mask_in_km2)) +
    geom_line() +
    labs(y = "plume area (km^2)", x = NULL) +
    scale_x_date(expand = 0) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # panache_plot
  
  # Plot wind speed and panache size correlation
  tide_plume_cor_plot <- ggplot(tide_plume_df, aes(x = tide_range, y = area_of_the_plume_mask_in_km2)) + 
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm") +
    labs(y = "plume area (km^2)", x = "tidal range (m)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  # tide_plume_cor_plot
  
  # Plot lag results
  tide_plume_cor_lag_plot <- ggplot(tide_plume_lag_cor, aes(x = lag, y = cor)) +
    geom_point() +
    labs(x = "lag plume after tidal range (days)", y = "correlation (r)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # tide_plume_cor_lag_plot
  
  # Combine plots and save
  tide_plume_title <- grid::textGrob(paste0(mouth_info$mouth_name," : tide vs plume size"), gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(tide_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
  cor_plot <- ggpubr::ggarrange(tide_plume_cor_plot, tide_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
  full_plot_title <- ggpubr::ggarrange(tide_plume_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
  ggsave(filename = paste0("figures/tide_plume_cor_plot_",mouth_info$mouth_name,".png"), width = 12, height = 6, dpi = 600)
}


# Run ---------------------------------------------------------------------

plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = tide_calc)

