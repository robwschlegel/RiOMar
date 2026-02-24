# func/flow.R
# Comparisons of river flow against plume size


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(sandwich) # For HAC covariance tests
library(lmtest) # For more detailed linear model tests
library(patchwork)
library(doParallel); doParallel::registerDoParallel(cores = 14)


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
              plume_perc = round((plume_slope/mean(flow_plume_daily_all$plume_area, na.rm = TRUE))*100, 2),
              plume_p = round(summary(lm(plume_area ~ date))[["coefficients"]][2,4], 4),
              flow_slope = coef(lm(flow ~ date))["date"] * 365.25,
              flow_perc = round((flow_slope/mean(flow_plume_daily_all$flow, na.rm = TRUE))*100, 2),
              flow_p = round(summary(lm(flow ~ date))[["coefficients"]][2,4], 4))
  trend_monthly_all <- flow_plume_monthly_all |> 
    summarise(plume_slope = coef(lm(plume_area ~ date))["date"] * 365.25,
              flow_slope = coef(lm(flow ~ date))["date"] * 365.25)
  trend_annual_all <- flow_plume_annual_all |> 
    summarise(plume_slope = coef(lm(plume_area ~ date))["date"] * 365.25,
              flow_slope = coef(lm(flow ~ date))["date"] * 365.25)
  trend_labels_all <- trend_daily_all |> 
    reshape2::melt() |> 
    mutate(x = c(rep(as.Date("1999-01-01"), 3), rep(as.Date("2012-01-01"), 3)),
           y = round(max(flow_plume_daily_all$plume_area, na.rm = TRUE)-quantile(flow_plume_daily_all$plume_area, 0.1, na.rm = TRUE), -2))
  
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
    # geom_hline(yintercept = plume_90, linetype = "dashed", colour = "brown", linewidth = 1) +
    # geom_hline(yintercept = flow_scaled_90, linetype = "dashed", colour = "blue", linewidth = 1) +
    geom_label(data = filter(trend_labels_all, x == "1999-01-01"), size = 5, hjust = 0,
               aes(x = x[1], y = y[1], label = paste0("Plume area slope = ", round(value[variable == "plume_slope"], 2), " km^2 yr-1\n",
                                                "Change per year = ",round(value[variable == "plume_perc"], 2), "% ; ",
                                                "p-value = ",round(value[variable == "plume_p"], 2), "", sep = ""))) +
    geom_label(data = filter(trend_labels_all, x == "2012-01-01"), size = 5, hjust = 0,
                 aes(x = x, y = y, label = paste0("River flow slope = ", round(value[variable == "flow_slope"], 2), " m^3 s-1 y-1\n", 
                                                  "Change per year = ",round(value[variable == "flow_perc"], 2), "% ; ",
                                                  "p-value = ",round(value[variable == "flow_p"], 2), "", sep = ""))) +
    scale_y_continuous(name = "Plume area (km^2)",
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = "River flow (m^3 s-1)")) +
    labs(x = NULL, title = paste(mouth_info$mouth_name,": Trends for daily plume area and river flow"),
         subtitle = "Solid line = linear trend") +#; dashed line = 90th percentile") +
    scale_color_manual(name = "Variable",
                       values = c("plume" = "brown", "flow" = "blue"),
                       breaks = c("plume", "flow")) +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = NA, colour = "black"))
  # line_trend_base_all
  
  # Get the count of the 90th percentile days changing over time
  # flow_plume_90th_all <- flow_plume_daily_all |> 
  #   mutate(plume_90th = ifelse(plume_area >= plume_90, 1, 0),
  #          flow_90th = ifelse(flow >= flow_90, 1, 0),
  #          year = year(date),
  #          date = round_date(date, "month") + days(14)) |> 
  #   summarise(plume_90th_count = sum(plume_90th, na.rm = TRUE),
  #             flow_90th_count = sum(flow_90th, na.rm = TRUE), .by = c("date"))
  
  # The trends in 90th percentile exceedance
  # trend_90th_all <- flow_plume_90th_all |> 
  #   summarise(plume_90th_slope = coef(lm(plume_90th_count ~ date))["date"] * 365.25,
  #             flow_90th_slope = coef(lm(flow_90th_count ~ date))["date"] * 365.25)
  # trend_90th_labels_all <- trend_90th_all |> 
  #   reshape2::melt() |> 
  #   mutate(x = c(as.Date("2000-01-01"), as.Date("2015-01-01")),
  #          y = max(flow_plume_90th_all$flow_90th_count, na.rm = TRUE))
  
  # Plot the 90th perc. exceedances as barplots
  # line_trend_90th_all <-  flow_plume_90th_all |>
  #   dplyr::rename(plume = plume_90th_count, flow = flow_90th_count) |> 
  #   pivot_longer(plume:flow) |> 
  #   ggplot(aes(x = date, y = value)) +
  #   geom_col(alpha = 0.5, aes(fill = name), position = "dodge") +
  #   geom_smooth(method = "lm", linewidth = 2, se = FALSE, colour = "black", aes(group = name)) +
  #   geom_smooth(method = "lm", linewidth = 1.8, se = FALSE, aes(colour = name)) +
  #   # geom_hline(yintercept = mean(flow_plume_90th_all$plume_90th_count), linetype = "dashed", colour = "brown", linewidth = 1) +
  #   # geom_hline(yintercept = mean(flow_plume_90th_all$flow_90th_count), linetype = "dashed", colour = "blue", linewidth = 1) +
  #   geom_label(data = filter(trend_90th_labels_all, variable == "plume_90th_slope"), size = 5, hjust = 0,
  #              aes(x = x, y = y, label = paste0("Plume 90th slope = ", round(value[1], 2), " days yr-1", sep = ""))) +
  #   geom_label(data = filter(trend_90th_labels_all, variable == "flow_90th_slope"), size = 5, hjust = 0,
  #              aes(x = x, y = y, label = paste0("River flow 90th slope = ", round(value[1], 2), " days yr-1", sep = ""))) +
  #   labs(x = NULL, y = "Days per month above 90th perc. thresh. (n)",
  #        title = paste(mouth_info$mouth_name,": Trends for daily exceedance of 90th percentile"),
  #        subtitle = "Solid line = linear trend") +
  #   scale_color_manual(name = "Variable",
  #                      values = c("plume" = "brown", "flow" = "blue"),
  #                      breaks = c("plume", "flow"), aesthetics = c("colour", "fill")) +
  #   theme(legend.position = "bottom",
  #         panel.border = element_rect(fill = NA, colour = "black"))
  
  # Combine plots and save
  # flow_plume_combi <- line_trend_base_all / line_trend_90th_all
  ggsave(filename = paste0("figures/trends_plume_flow_",mouth_info$mouth_name,".png"), 
         line_trend_base_all, width = 12, height = 6, dpi = 600)
}



# Calculate flow analyses -------------------------------------------------

# Compare panache size and river flow data
plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_comp)

# Calculate the linear trends for panache and river flow
plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_trend, .parallel = TRUE)


# Trends and components ---------------------------------------------------

# TODO: Decompose flow and panache time series with STL and X11 and then fit linear trends

# Linear and seasonal decomposition matching Toulouse workflow
# https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2022.1045667/full
# https://github.com/NOAA-PMEL/TOATS/blob/main/TOATS.ipynb

# Get the trends for river and panache using a more rigorous approach
# mouth_info <- river_mouths[4,]
flow_plume_trend_plus <- function(mouth_info){
  
  # Get zone info
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
  plume_df <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")) |> 
    dplyr::select(date:path_to_file) |> dplyr::select(-path_to_file) |> 
    mutate(area_of_the_plume_mask_in_km2 = ifelse(area_of_the_plume_mask_in_km2 > 20000, NA, area_of_the_plume_mask_in_km2))
  
  # Combine
  flow_plume_df <- left_join(plume_df, flow_df, join_by(date)) |> 
    zoo::na.trim() |> 
    dplyr::select(date, area_of_the_plume_mask_in_km2, flow) |> 
    dplyr::rename(plume_area = area_of_the_plume_mask_in_km2) |> 
    mutate(doy = yday(date),
           month = month(date),
           year = year(date))
  
  # Calculate residuals
  flow_plume_resid_df <- flow_plume_df |> 
    mutate(flow_resid = residuals(lm(flow ~ date, na.action = na.exclude)),
           plume_area_resid = residuals(lm(plume_area ~ date, na.action = na.exclude)))
  
  # Test visual
  # ggplot(data = flow_plume_resids_df, aes(x = date, y = plume_area)) +
  #   geom_point() +
  #   # geom_point(aes(y = flow_resid), colour = "red") +
  #   geom_point(aes(y = plume_area_flat), colour = "blue") +
  #   # geom_smooth(method = "lm") #+
  #   geom_smooth(aes(y = plume_area_flat), method = "lm", colour = "purple")
  
  # Calculate daily, monthly and annual clims
  flow_plume_daily_clim_df <- flow_plume_resid_df |> 
    summarise(flow_resid_doy = mean(flow_resid, na.rm = TRUE),
              plume_area_resid_doy = mean(plume_area_resid, na.rm = TRUE), .by = "doy") |> 
    mutate(flow_resid_doy_clim = flow_resid_doy-mean(flow_resid_doy, na.rm = TRUE),
           plume_area_resid_doy_clim = plume_area_resid_doy-mean(plume_area_resid_doy, na.rm = TRUE))
  flow_plume_month_clim_df <- flow_plume_resid_df |> 
    summarise(flow_resid_monthly = mean(flow_resid, na.rm = TRUE),
              plume_area_resid_monthly = mean(plume_area_resid, na.rm = TRUE), .by = "month") |> 
    mutate(flow_resid_monthly_clim = flow_resid_monthly-mean(flow_resid_monthly, na.rm = TRUE),
           plume_area_resid_monthly_clim = plume_area_resid_monthly-mean(plume_area_resid_monthly, na.rm = TRUE))
  
  # Calculate seasonal amplitudes
  # NB: Not used for final stats
  flow_plume_seas_amp <- flow_plume_month_clim_df |> 
    summarise(flow_resid_min = min(flow_resid_monthly_clim, na.rm = TRUE),
              flow_resid_max = max(flow_resid_monthly_clim, na.rm = TRUE),
              plume_area_resid_min = min(plume_area_resid_monthly_clim, na.rm = TRUE),
              plume_area_resid_max = max(plume_area_resid_monthly_clim, na.rm = TRUE))
  
  # Get annual averages
  # NB: Not used for final stats
  flow_plume_annual_df <- flow_plume_month_clim_df |> 
    summarise(flow_resid_annual = mean(flow_resid_monthly, na.rm = TRUE),
              plume_area_resid_annual = mean(plume_area_resid_monthly, na.rm = TRUE))
  
  # Calculate non-de-trended doy and monthly means and apply monthly clim correction
  flow_plume_daily_df <- flow_plume_df |> 
    left_join(flow_plume_daily_clim_df, by = join_by(doy)) |> 
    mutate(flow_doy_adj = flow-flow_resid_doy_clim,
           plume_area_doy_adj = plume_area-plume_area_resid_doy_clim) |> 
    # TODO: Make this dynamic based on input data
    # To update once the newer flow data have been added
    filter(date <= as.Date("2023-12-31")) |> 
    mutate(date_int = seq(1:n()), .after = "date") # Necessary for linear model
  flow_plume_monthly_df <- flow_plume_df |> 
    mutate(date = floor_date(date, "month")) |> 
    summarise(flow_monthly = mean(flow, na.rm = TRUE),
              plume_area_monthly = mean(plume_area, na.rm = TRUE), .by = c("year", "month", "date")) |> 
    left_join(flow_plume_month_clim_df, by = join_by(month)) |> 
    mutate(flow_monthly_adj = flow_monthly-flow_resid_monthly_clim,
           plume_area_monthly_adj = plume_area_monthly-plume_area_resid_monthly_clim) |> 
    # TODO: Make this dynamic based on input data
    # To update once the newer flow data have been added
    filter(date <= as.Date("2023-12-31")) |> 
    mutate(date_int = seq(1:n()), .after = "date") # Necessary for linear model
  
  # Fit Ar(1) model to residuals to get autocorrelation structure
  ar_weights_func <- function(val_col, start_year, time_step){
    
    # Create ts object
    ts_obj <- ts(zoo::na.approx(val_col), frequency = time_step, start = c(start_year, 1))
    
    # Run AR model
    ar_model <- ar(ts_obj, order.max = 1)
    
    # Estimated AR(1) coefficient
    phi_est <- ar_model$ar
    
    # Estimated noise variance
    sigma_est <- sqrt(phi_est)
    error_variance <- sigma_est^2 / (1 - phi_est^2)
    ar_weights <- rep((1 / (error_variance^2)), length(val_col))
    return(ar_weights)
  }
  # ar_weights_func(flow_plume_monthly_df$plume_area_monthly_adj, 1998, 12)
  
  # Fit STL to get the error variance for weighting
  stl_weights_func <- function(val_col, start_year, time_step){
    
    # Create ts object
    ts_obj <- ts(zoo::na.approx(val_col), frequency = time_step, start = c(start_year, 1))
    
    # Run STL
    stl_ts <- stl(ts_obj, s.window = "periodic")
    
    # get weights from residuals as an expression of variance around trend and seasonal cycle
    stl_var <- as.vector(stl_ts$time.series[,"remainder"])
    stl_weights <- 1 / (stl_var^2)
    return(stl_weights)
  }
  # stl_weights_func(flow_plume_monthly_df$plume_area_monthly_adj, 1998, 12)

  # Linear model function that runs the desired weights and then runs an HAC correction on the output
  lm_HAC_weights_func <- function(weight_choice, val_col, date_col){
    
    # Determine start_year from date_col
    start_year <- year(min(date_col))
    
    # Determine time_step from date_col
    if(length(val_col) < 1000){
      time_step <- 12
    } else {
      time_step <- 365
    }
    
    # First calculate weights
    if(weight_choice == "ar"){
      weights <- ar_weights_func(val_col, start_year, time_step)
    } else if(weight_choice == "stl"){
      weights <- stl_weights_func(val_col, start_year, time_step)
    } else {
      weights <- rep(1, length(val_col))
    }
    
    lm_model <- lm(val_col ~ date_col, weights = weights)
    lm_model_HAC <- coeftest(lm_model, vcov = vcovHAC(lm_model))
    
    # Create dataframe of results
    lm_res <- tibble(
      n = length(val_col),
      time_step = time_step,
      start_year = start_year,
      weight_choice = weight_choice,
      intercept = lm_model_HAC[1,1],
      slope = lm_model_HAC[2,1],
      slope_se = lm_model_HAC[2,2],
      slope_t = lm_model_HAC[2,3],
      slope_p = lm_model_HAC[2,4]
    )
    return(lm_res)
  }
  
  # Calculate linear models for daily and monthly flows
  wls_flow_daily <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, 
                                  val_col = flow_plume_daily_df$flow_doy_adj, 
                                  date_col = flow_plume_daily_df$date)
  wls_flow_monthly <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, 
                                 val_col = flow_plume_monthly_df$flow_monthly_adj, 
                                 date_col = flow_plume_monthly_df$date)
  
  # Calculate linear models for daily and monthly plumes
  wls_plume_area_daily <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, 
                                val_col = flow_plume_daily_df$plume_area_doy_adj, 
                                date_col = flow_plume_daily_df$date)
  wls_plume_area_monthly <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, 
                                  val_col = flow_plume_monthly_df$plume_area_monthly_adj, 
                                  date_col = flow_plume_monthly_df$date)
  
  # Get labels for plotting
  trend_labels_plume <- rbind(wls_plume_area_daily, wls_plume_area_monthly) |> 
    filter(weight_choice == "ar") |> # Just working with error estimates that match the paper
    reshape2::melt(id.vars = c("weight_choice", "time_step")) |> 
    filter(variable %in% c("slope", "slope_p")) |> 
    arrange(time_step) |>
    mutate(x = c(rep(as.Date("1999-01-01"), 2), rep(as.Date("2012-01-01"), 2)),
           y = round(max(flow_plume_daily_df$plume_area, na.rm = TRUE)-quantile(flow_plume_daily_df$plume_area, 0.2, na.rm = TRUE), -2))
  trend_labels_flow <- rbind(wls_flow_daily, wls_flow_monthly) |> 
    filter(weight_choice == "ar") |> # Just working with error estimates that match the paper
    reshape2::melt(id.vars = c("weight_choice", "time_step")) |> 
    filter(variable %in% c("slope", "slope_p")) |> 
    arrange(time_step) |>
    mutate(x = c(rep(as.Date("1999-01-01"), 2), rep(as.Date("2012-01-01"), 2)),
           y = round(max(flow_plume_daily_df$flow, na.rm = TRUE)-quantile(flow_plume_daily_df$flow, 0.05, na.rm = TRUE), -2))
  
  # Plot the plume results with labels
  pl_plume <- ggplot(data = flow_plume_daily_df, aes(x = date, y = plume_area)) +
    # Add daily+monthly points
    geom_point(colour = "sienna", alpha = 0.1) +
    geom_point(aes(y = plume_area_doy_adj), colour = "darkblue", alpha = 0.1) +
    geom_point(data = flow_plume_monthly_df, aes(y = plume_area_monthly), colour = "sienna", alpha = 0.6, size = 3) +
    geom_point(data = flow_plume_monthly_df, aes(y = plume_area_monthly_adj), colour = "darkred", size = 3) +
    # Add slopes
    geom_abline(intercept = wls_plume_area_daily$intercept[1],
                slope = wls_plume_area_daily$slope[1], linewidth = 2, colour = "darkblue") +
    geom_abline(intercept = wls_plume_area_monthly$intercept[1],
                slope = wls_plume_area_monthly$slope[1], linewidth = 2, colour = "darkred") +
    # Add labels
    geom_label(data = filter(trend_labels_plume, time_step == 365), size = 5, hjust = 0, colour = "darkblue",
               aes(x = x[1], y = y[1], label = paste0("Plume area slope = ", round(value[variable == "slope"]*365.25, 2), " km^2 yr-1\n",
                                                      "p-value = ",round(value[variable == "slope_p"], 2), "", sep = ""))) +
    geom_label(data = filter(trend_labels_plume, time_step == 12), size = 5, hjust = 0, colour = "darkred",
               aes(x = x[1], y = y[1], label = paste0("Plume area slope = ", round(value[variable == "slope"]*365.25, 2), " km^2 yr-1\n",
                                                "p-value = ",round(value[variable == "slope_p"], 2), "", sep = ""))) +
    labs(x = NULL, y = "Plume area [km^2]", 
         title = paste0(mouth_info$mouth_name," : Plume area after statistical treatment"),
         subtitle = "Red values show adjusted monthly values; blue shows adjusted daily values; brown are original data") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # pl_plume
    
    # Plot the river flow results with labels
  pl_flow <- ggplot(data = flow_plume_daily_df, aes(x = date, y = flow)) +
      # Add daily+monthly points
      geom_point(colour = "purple", alpha = 0.1) +
      geom_point(aes(y = flow_doy_adj), colour = "darkblue", alpha = 0.1) +
      geom_point(data = flow_plume_monthly_df, aes(y = flow_monthly), colour = "purple", alpha = 0.6, size = 3) +
      geom_point(data = flow_plume_monthly_df, aes(y = flow_monthly_adj), colour = "darkred", size = 3) +
      # Add slopes
      geom_abline(intercept = wls_flow_daily$intercept[1],
                  slope = wls_flow_daily$slope[1], linewidth = 2, colour = "darkblue") +
      geom_abline(intercept = wls_flow_monthly$intercept[1],
                  slope = wls_flow_monthly$slope[1], linewidth = 2, colour = "darkred") +
      # Add labels
      geom_label(data = filter(trend_labels_flow, time_step == 365), size = 5, hjust = 0, colour = "darkblue",
                 aes(x = x[1], y = y[1], label = paste0("River flow slope = ", round(value[variable == "slope"]*365.25, 2), " m^3 s-1 y-1\n", 
                                                        "p-value = ",round(value[variable == "slope_p"], 2), "", sep = ""))) +
      geom_label(data = filter(trend_labels_flow, time_step == 12), size = 5, hjust = 0, colour = "darkred",
                 aes(x = x[1], y = y[1], label = paste0("River flow slope = ", round(value[variable == "slope"]*365.25, 2), " m^3 s-1 y-1\n", 
                                                        "p-value = ",round(value[variable == "slope_p"], 2), "", sep = ""))) +
      labs(x = NULL, y = "River flow [m^3 s-1]", 
           title = paste0(mouth_info$mouth_name," : River flow after statistical treatment"),
           subtitle = "Red values show adjusted monthly values; blue shows adjusted daily values; purple are original data") +
      theme(panel.border = element_rect(fill = NA, colour = "black"))
  # pl_flow
  
  # Combine, save, exit
  pl_combi <- ggpubr::ggarrange(pl_plume, pl_flow, ncol = 1, nrow = 2)
  ggsave(filename = paste0("figures/trends_plume_flow_adj_", mouth_info$mouth_name,".png"), 
         pl_combi, width = 12, height = 10)
}

# Run for all zones
plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_plume_trend_plus, .parallel = TRUE)

