# func/flow.R
# Comparisons of river flow against plume size


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(doParallel); doParallel::registerDoParallel(cores = 4)


# Functions ---------------------------------------------------------------

# Calculate relationship between tides and panache size
# mouth_info <- river_mouths[1,]
flow_calc <- function(mouth_info){
  
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
  
  # Load river flow data
  flow_df <- load_river_flow(paste0("data/RIVER_FLOW/",zone))
  
  # Load panache time series based on river mouth name
  plume_daily <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")) |> 
    dplyr::select(date:path_to_file) |> dplyr::select(-path_to_file) |> 
    mutate(area_of_the_plume_mask_in_km2 = ifelse(area_of_the_plume_mask_in_km2 > 20000, NA, area_of_the_plume_mask_in_km2))
  
  # Combine
  flow_plume_df <- left_join(plume_daily, flow_df, join_by(date)) |> 
    zoo::na.trim()
  
  # Create time series objects for stl
  ts_plume <- ts(zoo::na.approx(flow_plume_df$area_of_the_plume_mask_in_km2), frequency = 365, start = c(year(min(flow_plume_df$date)), 1))
  stl_plume <- stl(ts_plume, s.window = "periodic")
  ts_flow <- ts(zoo::na.approx(flow_plume_df$flow), frequency = 365, start = c(year(min(flow_plume_df$date)), 1))
  stl_flow <- stl(ts_flow, s.window = "periodic")
  
  # Add trend elements back into dataframe for further stats
  flow_plume_df$flow_stl <- as.vector(stl_flow$time.series[,2])
  flow_plume_df$plume_stl <- as.vector(stl_plume$time.series[,2])
  # cor(flow_plume_df$flow_stl, flow_plume_df$plume_stl)
  
  # Compare panache size against wind speed
  # Two separate comparisons based on upwelling or downwelling times
  flow_plume_stats_all <- flow_plume_df |> 
    summarise(r = cor(flow, area_of_the_plume_mask_in_km2, use = "pairwise.complete.obs"))
  
  # Lagged correlations
  flow_plume_lag_cor <- tibble(
    lag = 0:30,
    cor = map_dbl(0:30, ~ cor(flow_plume_df$flow, lag(flow_plume_df$area_of_the_plume_mask_in_km2, .), use = "pairwise.complete.obs"))
  )
  
  # Plot wind speed and up/downwelling times
  flow_plot <- ggplot(flow_plume_df, aes(x = date, y = flow)) +
    # geom_ribbon(aes(ymin = (flow-(flow/2)), ymax = flow+flow/2)) +
    geom_line() +
    labs(y = "River flow (m^3 s-1)", x = NULL) +
    scale_x_date(expand = 0) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # flow_plot
  
  # Panache size
  panache_plot <- ggplot(flow_plume_df, aes(x = date, y = area_of_the_plume_mask_in_km2)) +
    geom_line() +
    labs(y = "plume area (km^2)", x = NULL) +
    scale_x_date(expand = 0) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # panache_plot
  
  # Plot wind speed and panache size correlation
  flow_plume_cor_plot <- ggplot(flow_plume_df, aes(x = flow, y = area_of_the_plume_mask_in_km2)) + 
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm") +
    labs(y = "plume area (km^2)", x = "River flow (m^3 s-1)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  # flow_plume_cor_plot
  
  # Plot lag results
  flow_plume_cor_lag_plot <- ggplot(flow_plume_lag_cor, aes(x = lag, y = cor)) +
    geom_point() +
    labs(x = "lag plume after tidal range (days)", y = "correlation (r)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # flow_plume_cor_lag_plot
  
  # Combine plots and save
  flow_plume_title <- grid::textGrob(paste0(mouth_info$mouth_name," : river flow vs plume size"), gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(flow_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
  cor_plot <- ggpubr::ggarrange(flow_plume_cor_plot, flow_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
  full_plot_title <- ggpubr::ggarrange(flow_plume_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
  ggsave(filename = paste0("figures/cor_plot_flow_plume_",mouth_info$mouth_name,".png"), width = 12, height = 6, dpi = 600)
  
  # Get scaling factor for plotting
  # scaling_factor <- sec_axis_adjustement_factors(var_to_scale = flow_plume_df$flow_stl, 
  #                                                var_ref = flow_plume_df$plume_stl)
  # flow_plume_df <- flow_plume_df |> 
  #   mutate(flow_scaled = flow_stl * scaling_factor$diff + scaling_factor$adjust)
  # unique_years <- flow_plume_df$date |> year() |> unique()
  # 
  # # Plots for STL decomposition
  # ggplot(data = flow_plume_df) + 
  #   # Plume data
  #   geom_point(aes(x = date, y = plume_stl), color = "brown") + 
  #   geom_path(aes(x = date, y = plume_stl), color = "brown") + 
  #   # Wind data
  #   geom_point(aes(x = date, y = flow_scaled), color = "blue") + 
  #   geom_path(aes(x = date, y = flow_scaled), color = "blue") + 
  #   # X-axis labels
  #   scale_x_date(name = "", 
  #                breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
  #                labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
  #   # Y-axis labels
  #   scale_y_continuous(name = "Plume area (km²)",
  #                      sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
  #                                          name = "Wind speed (m/s)")) +
  #   # Extra bits
  #   labs(title = zone) +
  #   ggplot_theme() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  #         plot.subtitle = element_text(hjust = 0.5),
  #         axis.text.y.left = element_text(color = "brown"), 
  #         axis.ticks.y.left = element_line(color = "brown"),
  #         axis.line.y.left = element_line(color = "brown"),
  #         axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
  #         axis.text.y.right = element_text(color = "blue"), 
  #         axis.ticks.y.right = element_line(color = "blue"),
  #         axis.line.y.right = element_line(color = "blue"),
  #         axis.title.y.right = element_text(color = "blue", margin = unit(c(0, 0, 0, 7.5), "mm")),
  #         panel.border = element_rect(linetype = "solid", fill = NA))
}


# Run ---------------------------------------------------------------------

plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_calc)

