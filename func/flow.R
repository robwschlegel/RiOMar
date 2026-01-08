# func/flow.R
# Comparisons of river flow against plume size


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(patchwork)
library(doParallel); doParallel::registerDoParallel(cores = 4)


# Functions ---------------------------------------------------------------

# Calculate relationship between river flow and panache size
# mouth_info <- river_mouths[1,]
flow_comp <- function(mouth_info){
  
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
  
  # Compare panache size against river flow
  flow_plume_stats_all <- flow_plume_df |> 
    summarise(r = cor(flow, area_of_the_plume_mask_in_km2, use = "pairwise.complete.obs"))
  
  # Lagged correlations
  flow_plume_lag_cor <- tibble(
    lag = 0:30,
    cor = map_dbl(0:30, ~ cor(flow_plume_df$flow, lag(flow_plume_df$area_of_the_plume_mask_in_km2, .), use = "pairwise.complete.obs"))
  )
  
  # Plot river flow
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
  
  # Plot river flow and panache size correlation
  flow_plume_cor_plot <- ggplot(flow_plume_df, aes(x = flow, y = area_of_the_plume_mask_in_km2)) + 
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linewidth = 2, linetype = "dashed", color = "black") +
    geom_smooth(method = "lm",  se = FALSE, colour = "black", linewidth = 2) +
    labs(y = "plume area (km^2)", x = "River flow (m^3 s-1)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  # flow_plume_cor_plot
  
  # Plot lag results
  flow_plume_cor_lag_plot <- ggplot(flow_plume_lag_cor, aes(x = lag, y = cor)) +
    geom_point() +
    labs(x = "lag plume after river flow (days)", y = "correlation (r)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # flow_plume_cor_lag_plot
  
  # Combine plots and save
  flow_plume_title <- grid::textGrob(paste0(mouth_info$mouth_name," : river flow vs plume size"), gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(flow_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
  cor_plot <- ggpubr::ggarrange(flow_plume_cor_plot, flow_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
  full_plot_title <- ggpubr::ggarrange(flow_plume_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
  ggsave(filename = paste0("figures/cor_plot_flow_plume_",mouth_info$mouth_name,".png"), width = 12, height = 6, dpi = 600)
}

# Calculate the linear trends for river flow and panache size
# mouth_info <- river_mouths[1,]
flow_trend <- function(mouth_info){
  
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
  flow_plume_daily_all <- left_join(plume_daily, flow_df, join_by(date)) |> 
    zoo::na.trim() |> 
    dplyr::select(date, area_of_the_plume_mask_in_km2, flow) |> 
    dplyr::rename(plume_area = area_of_the_plume_mask_in_km2)
  ## Monthly
  flow_plume_monthly_all <- flow_plume_daily_all |> 
    mutate(date = round_date(date, "month") + days(14)) |> 
    summarise(plume_area = mean(plume_area, na.rm = TRUE), 
              flow = mean(flow, na.rm = TRUE), .by = c("date"))
  ## Annual
  flow_plume_annual_all <- flow_plume_daily_all |> 
    mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
    summarise(plume_area = mean(plume_area, na.rm = TRUE), 
              flow = mean(flow, na.rm = TRUE), .by = c("date"))
  
  # Calculate trend 
  trend_daily_all <- flow_plume_daily_all |> 
    summarise(plume_slope = coef(lm(plume_area ~ date))["date"] * 365.25,
              flow_slope = coef(lm(flow ~ date))["date"] * 365.25)
  trend_monthly_all <- flow_plume_monthly_all |> 
    summarise(plume_slope = coef(lm(plume_area ~ date))["date"] * 365.25,
              flow_slope = coef(lm(flow ~ date))["date"] * 365.25)
  trend_annual_all <- flow_plume_annual_all |> 
    summarise(plume_slope = coef(lm(plume_area ~ date))["date"] * 365.25,
              flow_slope = coef(lm(flow ~ date))["date"] * 365.25)
  trend_labels_all <- trend_daily_all |> 
    reshape2::melt() |> 
    mutate(x = c(as.Date("2000-01-01"), as.Date("2015-01-01")),
           y = round(max(flow_plume_daily_all$plume_area, na.rm = TRUE), -2))
  
  # Scale river flow to plume size
  scaling_factor <- sec_axis_adjustement_factors(flow_plume_daily_all$flow, flow_plume_daily_all$plume_area)
  flow_plume_daily_all <- flow_plume_daily_all |> 
    mutate(flow_scaled = flow * scaling_factor$diff + scaling_factor$adjust, .after = "flow")
  
  # Calculate 90th percentile for each column
  plume_90 <- quantile(flow_plume_daily_all$plume_area, probs = 0.9, na.rm = TRUE)
  flow_90 <- quantile(flow_plume_daily_all$flow, probs = 0.9, na.rm = TRUE)
  flow_scaled_90 <- quantile(flow_plume_daily_all$flow_scaled, probs = 0.9, na.rm = TRUE)
  
  # Plot daily data for both variables
  line_trend_base_all <- ggplot(flow_plume_daily_all, aes(x = date, y = plume_area)) +
    geom_point(alpha = 0.3, aes(colour = "plume")) +
    geom_point(alpha = 0.3, aes(y = flow_scaled, colour = "flow")) +
    geom_smooth(method = "lm", colour = "black", linewidth = 2) +
    geom_smooth(method = "lm", colour = "brown", linewidth = 1.8) +
    geom_smooth(aes(y = flow_scaled), method = "lm", colour = "black", linewidth = 2) +
    geom_smooth(aes(y = flow_scaled), method = "lm", colour = "blue", linewidth = 1.8) +
    geom_hline(yintercept = plume_90, linetype = "dashed", colour = "brown", linewidth = 1) +
    geom_hline(yintercept = flow_scaled_90, linetype = "dashed", colour = "blue", linewidth = 1) +
    geom_label(data = filter(trend_labels_all, variable == "plume_slope"), size = 5, hjust = 0,
               aes(x = x, y = y, label = paste0("Plume area slope = ", round(value[1], 2), " km^2 yr-1", sep = ""))) +
    geom_label(data = filter(trend_labels_all, variable == "flow_slope"), size = 5, hjust = 0,
                 aes(x = x, y = y, label = paste0("River flow slope = ", round(value[1], 2), " m^3 s-1 y-1", sep = ""))) +
    scale_y_continuous(name = "Plume area (km^2)",
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = "River flow (m^3 s-1)")) +
    labs(x = NULL, title = paste(mouth_info$mouth_name,": Trends for daily plume area and river flow"),
         subtitle = "Solid line = linear trend; dashed line = 90th percentile") +
    scale_color_manual(name = "Variable",
                       values = c("plume" = "brown", "flow" = "blue"),
                       breaks = c("plume", "flow")) +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA, colour = "black"))
  # line_trend_base_all
  
  # Get the count of the 90th percentile days changing over time
  flow_plume_90th_all <- flow_plume_daily_all |> 
    mutate(plume_90th = ifelse(plume_area >= plume_90, 1, 0),
           flow_90th = ifelse(flow >= flow_90, 1, 0),
           year = year(date),
           date = round_date(date, "month") + days(14)) |> 
    summarise(plume_90th_count = sum(plume_90th, na.rm = TRUE),
              flow_90th_count = sum(flow_90th, na.rm = TRUE), .by = c("date"))
  
  # The trends in 90th percentile exceedance
  trend_90th_all <- flow_plume_90th_all |> 
    summarise(plume_90th_slope = coef(lm(plume_90th_count ~ date))["date"] * 365.25,
              flow_90th_slope = coef(lm(flow_90th_count ~ date))["date"] * 365.25)
  trend_90th_labels_all <- trend_90th_all |> 
    reshape2::melt() |> 
    mutate(x = c(as.Date("2000-01-01"), as.Date("2015-01-01")),
           y = max(flow_plume_90th_all$flow_90th_count, na.rm = TRUE))
  
  # Plot the 90th perc. exceedances as barplots
  line_trend_90th_all <-  flow_plume_90th_all |>
    dplyr::rename(plume = plume_90th_count, flow = flow_90th_count) |> 
    pivot_longer(plume:flow) |> 
    ggplot(aes(x = date, y = value)) +
    geom_col(alpha = 0.5, aes(fill = name), position = "dodge") +
    geom_smooth(method = "lm", linewidth = 2, se = FALSE, colour = "black", aes(group = name)) +
    geom_smooth(method = "lm", linewidth = 1.8, se = FALSE, aes(colour = name)) +
    # geom_hline(yintercept = mean(flow_plume_90th_all$plume_90th_count), linetype = "dashed", colour = "brown", linewidth = 1) +
    # geom_hline(yintercept = mean(flow_plume_90th_all$flow_90th_count), linetype = "dashed", colour = "blue", linewidth = 1) +
    geom_label(data = filter(trend_90th_labels_all, variable == "plume_90th_slope"), size = 5, hjust = 0,
               aes(x = x, y = y, label = paste0("Plume 90th slope = ", round(value[1], 2), " days yr-1", sep = ""))) +
    geom_label(data = filter(trend_90th_labels_all, variable == "flow_90th_slope"), size = 5, hjust = 0,
               aes(x = x, y = y, label = paste0("River flow 90th slope = ", round(value[1], 2), " days yr-1", sep = ""))) +
    labs(x = NULL, y = "Days per month above 90th perc. thresh. (n)",
         title = paste(mouth_info$mouth_name,": Trends for daily exceedance of 90th percentile"),
         subtitle = "Solid line = linear trend") +
    scale_color_manual(name = "Variable",
                       values = c("plume" = "brown", "flow" = "blue"),
                       breaks = c("plume", "flow"), aesthetics = c("colour", "fill")) +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA, colour = "black"))
  
  # Combine plots and save
  flow_plume_combi <- line_trend_base_all / line_trend_90th_all
  ggsave(filename = paste0("figures/trends_plume_flow_",mouth_info$mouth_name,".png"), width = 12, height = 8, dpi = 600)
}


# Run ---------------------------------------------------------------------

# Compare panache size and river flow data
plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_comp)

# Calculate the linear trends and 90th percentile days for panache and river flow
plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_trend, .parallel = TRUE)

