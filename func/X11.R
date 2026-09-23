# func/X11.R
# This script provides graphical creation for the X11.py script


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(viridis)
library(ggpubr)
library(scales)
library(ggnewscale)
library(zoo)

# The central functions -- multi.R (which itself sources util.R) rather than
# util.R alone, needed by plot_driver_x11_trend_comparison() below for
# driver_display and zone_title().
source("func/multi.R")

# where_are_saved_X11_results = "~/RiOMar/output/FIXED_THRESHOLD/BAY_OF_BISCAY"
# Zone= "BAY_OF_SEINE"
# Data_source = "SEXTANT"
# sensor_name = "merged"
# atmospheric_correction = "Standard"
# Temporal_resolution = "WEEKLY"
### Load the files 


# Utils -------------------------------------------------------------------

get_X11_data <- function(where_are_saved_X11_results, Zone, Data_source, sensor_name, atmospheric_correction, Temporal_resolution) {
  
  path_to_X11_data <- where_are_saved_X11_results %>% file.path(Zone, "X11_ANALYSIS")
  
  X11_plume_area_ts <- path_to_X11_data %>% 
    file.path("area_of_the_plume_mask_in_km2", 
              paste(Data_source, "_", sensor_name, "_", atmospheric_correction, "_", Temporal_resolution, ".csv", sep = "")) %>% 
    read_csv()
  
  X11_river_flow_ts <- path_to_X11_data %>% 
    file.path("river_flow", paste("River_flow___", Temporal_resolution, ".csv", sep = "")) %>% 
    read_csv()
  
  X11_ts <- X11_plume_area_ts %>% inner_join(X11_river_flow_ts, by = "dates", suffix = c("_plume_area", "_river_flow"))
  
  return( X11_ts )
  
}


make_the_plot <- function(X11_data, type_of_signal) {
  
  unique_years <- X11_data$dates %>% year() %>% unique()
  
  X11_data_for_plot <- X11_data %>% 
    rename(river_flow = !!sym(paste(type_of_signal, "signal_river_flow", sep = "_")),
           plume_area = !!sym(paste(type_of_signal, "signal_plume_area", sep = "_"))) %>% 
    select(dates, river_flow, plume_area)
  
  if (type_of_signal %in% c("Seasonal", "Residual")) {
    X11_data_for_plot <- X11_data_for_plot %>% 
      mutate(river_flow = river_flow + mean(X11_data$Raw_signal_river_flow, na.rm = T),
             plume_area = plume_area + mean(X11_data$Raw_signal_plume_area, na.rm = T))
  }
  
  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = X11_data_for_plot$river_flow, 
                                                 var_ref = X11_data_for_plot$plume_area)
  
  X11_data_for_plot <- X11_data_for_plot %>% mutate(river_flow_scaled = river_flow * scaling_factor$diff + scaling_factor$adjust)

  r_value <- cor(X11_data_for_plot$plume_area, X11_data_for_plot$river_flow, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  the_plot <- ggplot() +

    geom_point(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown") +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown") +

    geom_point(data = X11_data_for_plot, aes(x = dates, y = river_flow_scaled), color = "blue") +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = river_flow_scaled), color = "blue") +

    annotate("text", x = min(X11_data_for_plot$dates), y = Inf, label = r_label,
            hjust = 0, vjust = 1.5, size = 6, colour = "black") +

    scale_x_date(name = "",
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    
    scale_y_continuous(name = "Plume area (km²)",
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
                                           name = "River flow (m³ s⁻¹)")) +
    
    labs(title = paste(type_of_signal, "signal")) +
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
  
  return(the_plot)

}


# Driver comparison (generalised, wind/tide/wave/current) -------------------
# Rollback (2026-09-23) of the Aug-11 driver-vs-plume X11 migration's Python
# plotting (func/X11.py::plot_driver_x11_dual_axis(), called from
# func/compute_driver_x11_figures.py) -- Python stays limited to the X11
# *calculation* (still func/compute_x11_driver_signals.py, which persists
# each driver's weekly Seasonal_signal/Interannual_signal to CSV); this is
# the R replacement for the *plotting* half, matching make_the_plot()'s
# house style (dual-axis, scaled second axis, r annotation) rather than
# matplotlib's. Uses the Interannual_signal component only, matching what
# the Python version compared (see decompose_driver_series() upstream).

plot_driver_x11_trend_comparison <- function(where_are_saved_X11_results, zone_name, driver_name){

  plume_path <- file.path(where_are_saved_X11_results, zone_name, "X11_ANALYSIS",
                          "area_of_the_plume_mask_in_km2", "SEXTANT_merged_Standard_WEEKLY.csv")
  driver_path <- file.path(where_are_saved_X11_results, zone_name, "X11_ANALYSIS",
                           driver_name, paste0(driver_name, "_WEEKLY.csv"))

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)

  if(!file.exists(plume_path) || !file.exists(driver_path)){
    return(ggplot() + labs(title = paste0(zone_title(zone_name), " (insufficient data)")) + ggplot_theme())
  }

  plume_ts <- readr::read_csv(plume_path, show_col_types = FALSE) |>
    dplyr::transmute(dates = date, plume_area = Interannual_signal)
  driver_ts <- readr::read_csv(driver_path, show_col_types = FALSE) |>
    dplyr::transmute(dates = date, driver_value = Interannual_signal)

  X11_data_for_plot <- dplyr::inner_join(plume_ts, driver_ts, by = "dates")
  unique_years <- X11_data_for_plot$dates |> year() |> unique()

  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = X11_data_for_plot$driver_value,
                                                 var_ref = X11_data_for_plot$plume_area)
  X11_data_for_plot <- X11_data_for_plot |>
    mutate(driver_scaled = driver_value * scaling_factor$diff + scaling_factor$adjust)

  r_value <- cor(X11_data_for_plot$plume_area, X11_data_for_plot$driver_value, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  ggplot() +
    geom_point(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown") +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown") +
    geom_point(data = X11_data_for_plot, aes(x = dates, y = driver_scaled), color = disp$driver_colour) +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = driver_scaled), color = disp$driver_colour) +
    annotate("text", x = min(X11_data_for_plot$dates), y = Inf, label = r_label,
            hjust = 0, vjust = 1.5, size = 6, colour = "black") +
    scale_x_date(name = "",
                breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist()) +
    scale_y_continuous(name = "Plume area (km²)",
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = disp$driver_label)) +
    labs(title = zone_title(zone_name)) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
         axis.text.y.left = element_text(color = "brown"),
         axis.ticks.y.left = element_line(color = "brown"),
         axis.line.y.left = element_line(color = "brown"),
         axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
         axis.text.y.right = element_text(color = disp$driver_colour),
         axis.ticks.y.right = element_line(color = disp$driver_colour),
         axis.line.y.right = element_line(color = disp$driver_colour),
         axis.title.y.right = element_text(color = disp$driver_colour, margin = unit(c(0, 0, 0, 7.5), "mm")),
         panel.border = element_rect(linetype = "solid", fill = NA))
}


# Main --------------------------------------------------------------------

plot_time_series_of_plume_area_and_river_flow <- function(where_are_saved_X11_results,
                                                          Zone, Data_source, sensor_name,
                                                          atmospheric_correction, Temporal_resolution) {
  
  X11_data <- get_X11_data(where_are_saved_X11_results, Zone, Data_source, sensor_name, atmospheric_correction, Temporal_resolution)
  
  Raw_plot <- make_the_plot(X11_data, type_of_signal = 'Raw')
  
  Interannual_plot <- make_the_plot(X11_data, type_of_signal = 'Interannual')
  
  Seasonal_plot <- make_the_plot(X11_data, type_of_signal = 'Seasonal')
  
  Residual_plot <- make_the_plot(X11_data, type_of_signal = 'Residual')
  
  final_plot <- ggarrange(Raw_plot, Interannual_plot, Seasonal_plot, Residual_plot, ncol = 1, nrow = 4, align = "v")
  
  final_plot <- annotate_figure(final_plot, top=text_grob(Zone %>% str_replace_all("_", " "), face = "bold", size = 60, color = "black"))
  
  save_plot_as_png(plot = final_plot, width = 40, height = 25, # = 35
                   path = file.path(where_are_saved_X11_results, Zone, "X11_ANALYSIS", "plume_area_vs_river_flow"),
                   name = paste(Data_source, sensor_name, atmospheric_correction, Temporal_resolution, sep = "_"))
  
}

