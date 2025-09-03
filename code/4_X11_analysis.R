list_of_packages <- c("plyr", "tidyverse", "ggpubr", "viridis", "doParallel", "zoo", "ggnewscale", "ggpubr", "scales")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, require, character.only = TRUE)

# where_are_saved_X11_results = '/home/terrats/Desktop/RIOMAR/TEST/RESULTS/TEST4'
# Zone= "BAY_OF_SEINE"
# Data_source = "SEXTANT"
# sensor_name = "merged"
# atmospheric_correction = "Standard"
# Temporal_resolution = "WEEKLY"
### Load the files 



### Plot them based on the script used for plotting plume area ts


#### Main ####

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


#### Utils ####



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

ggplot_theme <-   function() {
  theme(text = element_text(size=35, colour = "black"), #25
        plot.title = element_text(hjust = 0.5, size = 55),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x=element_text(angle=0),
        axis.ticks.length=unit(.25, "cm"))}

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
  
  the_plot <- ggplot() + 
    
    geom_point(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown") + 
    geom_path(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown") + 
    
    geom_point(data = X11_data_for_plot, aes(x = dates, y = river_flow_scaled), color = "blue") + 
    geom_path(data = X11_data_for_plot, aes(x = dates, y = river_flow_scaled), color = "blue") + 
    
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