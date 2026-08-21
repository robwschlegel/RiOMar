# func/X11.R
# This script provides graphical creation for the X11.py script


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(viridis)
library(ggpubr)
library(scales)
library(ggnewscale)
library(zoo)

# The central functions
source("func/util.R")

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
                                           name = "River flow (m³/s)")) +
    
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

