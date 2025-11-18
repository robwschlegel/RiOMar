# func/plume.R
# This script visualises the results generated from plume.py
# It is called via code/3_plumes.py


# Libraries ---------------------------------------------------------------

# list_of_packages <- c("plyr", "tidyverse", "ggpubr", "viridis", "doParallel", "zoo")
# new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
# lapply(list_of_packages, require, character.only = TRUE)

library(plyr)
library(tidyverse)
library(viridis)
library(ggpubr)
library(zoo)
library(doParallel)


# Ideas for analysing plume dynamics
# Loire + Gironde plumes
# plumes generally turn to the right as they exit the estuary and follow the coastline northwards
## the Gironde plume may merge with the Loire plume
# A south-easterly wind, which favours downwellings, will have the effect of constraining the plume to the coast and reinforcing the circulation of the coastal current
# A north-westerly wind, favouring upwellings, will stop the propagation of the plumes and disperse them offshore
# winter winds tend to favour plume confinement at the coast
## while a change in wind direction at the end of March or early April favours southwards propagation offshore of the plumes
# Seine plume
# The dynamics of the Seine plume in the bay of Seine depend mainly on bathymetry and tides, and to a lesser extent on wind
# The annual cycle is the largest.
# - Perhaps it would be interesting to measure the distance from the mouth to the edge of the plume
# - Or how much plume can be attributed to SPM along the coastal zone / within a shallow depth area
# [bottom] SPM concentrations are on average twice the surface SPM concentrations [Seine]

# Tides
# For each ROFI, one tide gauge is processed for: Havre for the Seine river, Saint-Nazaire for the Loire river, and Port-Bloc for the Gironde river
# At the diurnal timescale, the three ROFIs studied are therefore mainly influenced by the semi-diurnal tidal cycle, 
## although there are slight differences in their response to this forcing.
# - I will need to think of a way to search forward and backward through matched time series to find neap and spring tides
## - to then associate these to times when there are high runoff, that then quickly shifts to smaller/larger plume areas
# - Create an index value based on the build-up/-down from neap to spring tide
## - Meaning, don't measure the tidal range on a given day, but rather where on the progression along the larger tidal scale that day is
# Gironde
# The bulge is created in moderate to high discharge conditions and relatively weak wind conditions
# (1) formation of a bulge during the transition from spring to neap tides, 
# (2) reinforcement of the bulge in terms of stratification during neap tides, 
# (3) destabilization of the bulge during the neap to spring tides transition, with increased mixing that breaks stratification. 
# When the bulge is relaxed, upwelling favorable winds can export the freshwater offshore
# Seine
# For a given tidal range, increasing river discharge leads to increase SPM concentrations by a factor of 2 due to the shift of the 
# estuarine turbidity maximum toward the mouth, and high waves again multiply by 2 SPM concentrations for a given tidal range and river discharge.

# Wind
# The zonal and meridional wind speeds were averaged spatially for each region of interest
# The favourable downwelling and upwelling winds were defined as the quarter sector of wind parallel to the coast, 
## directed with the coast to the right and to the left, respectively
# As ROFI is stretched into a thin layer on the surface, it loses stability and becomes much more sensitive to mixing.
# These wind events on the scale of a few days can shape the ROFI and have an impact on the whole month ahead.
# However, for the ROFI of the Seine, the wind rarely has an impact except during neap tides. 
## During these periods, only strong upwelling favourable wind can enhance the ROFI extent.
# Weekly to monthly variability of ROFIs is dominantly driven by tidal dynamics for the Seine region and by wind and runoff dynamics in Loire and Gironde regions.
# The dynamics of the Loire and Gironde ROFIs are significantly impacted by wind between January and June, during the periods of highest runoff
# However, wind exerts minimal influence when river runoff is insufficient to generate a substantial ROFI. 
# - Think of a way to show how wind may proportionately increase plume size, based on the pre-existing plume size
## - Meaning, the plume will get relatively bigger due to wind forcing if it is already beefy
## - A small plume will not be magnified much by wind

# Discharge
# Flood events trigger a more or less rapid response from the ROFI, from few days to a month
# The extent of the ROFI is amplified when neap tides occur after a flood event
# At the onset of the spring tide, the freshwater volume decreases in all salinity classes as the surface area decreases. [Seine]
# The ROFIs for the Loire and Gironde regions do not show such a marked response to the spring and neap tide cycles.
# Notably larger surface extensions during the neap periods following a flood


# Utils -------------------------------------------------------------------

ggplot_theme <-   function() {
  theme(text = element_text(size = 35, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 55),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(angle = 0),
        axis.ticks.length = unit(.25, "cm"))}

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

save_file_as_csv <- function(data, file_name) {
  path <- file_name %>% str_split(pattern = "/") %>% .[[1]] %>% head(-1) %>% paste(collapse = "/")
  if (dir.exists(path) == FALSE) {dir.create(path = path, recursive = TRUE)}
  if (grepl(".csv", file_name) == FALSE) {file_name <- paste(file_name, ".csv", sep = "")}
  data.table::fwrite(data %>% as.data.frame(), file = file_name, row.names = TRUE, col.names = TRUE)
}

plot_function <- function(data, y_variable, Temporal_resolution) {
  
  if ("Satellite_sensor" %in% names(data)) {
    ggplot_base <- ggplot(data, aes(x = date, y = .data[[y_variable]], group = Satellite_sensor, colour = Satellite_sensor))  
  } else {
    ggplot_base <- ggplot(data, aes(x = date, y = .data[[y_variable]])) #+ labs(colour = "")
  }
  
  the_plot <- ggplot_base +
    geom_line(linewidth = 1.5, na.rm = TRUE) + geom_point(size = 2, na.rm = TRUE) + 
    labs(y = case_when(grepl("^SPM_threshold", y_variable) ~ y_variable %>% str_replace_all("_", " ") %>% str_replace_all("threshold ", "threshold\n"),
                       y_variable == "area_of_the_plume_mask_in_km2" ~ "Plume area (km²)")) + 
    ggplot_theme() 
  
  if (str_detect(y_variable, "threshold")) {
    
    window_size <- case_when(Temporal_resolution == "WEEKLY" ~ 8,
                             Temporal_resolution == "MONTHLY" ~ 2,
                             Temporal_resolution == "DAILY" ~ 60)+1
    
    smoothed_threshold_data <- rollmean(data[,y_variable], k = window_size, fill = "NA")
    # smoothed_threshold_data <- rollmean(smoothed_threshold_data, k = window_size, fill = "NA")
    data['smoothed_threshold_data'] <- smoothed_threshold_data
    
    the_plot <- the_plot + geom_line(aes(y = smoothed_threshold_data), 
                                     linewidth = 1.5, na.rm = TRUE, colour = "black")
    
  }
  
  return(the_plot)
  
}

make_the_plot_from_the_df <- function(data, row, where_to_save_the_file, plot_name, Temporal_resolution) {
  
  SPM_thresholds_columns <- data %>% names() %>% grep(pattern = "^SPM_threshold", value = TRUE)
  
  time_series_plot <- ggarrange(plotlist = c("area_of_the_plume_mask_in_km2", SPM_thresholds_columns) %>% 
                                  llply(plot_function, data = data, Temporal_resolution = Temporal_resolution),
                                ncol = 1, align = 'v', common.legend = TRUE)
  
  where_to_save_the_file %>% 
    l_ply(save_plot_as_png, name = plot_name, width = 28, height = 22, plot = time_series_plot)
  
}


# Main --------------------------------------------------------------------

# where_are_saved_plume_results = "output/FIXED_THRESHOLD"
# where_to_save_plume_time_series = "output/FIXED_THRESHOLD"
# # Zone = c('BAY_OF_SEINE', 'GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY')
# Zone = c('BAY_OF_SEINE')
# Data_source = 'SEXTANT'
# Satellite_sensor = c('merged')
# atmospheric_correction = 'Standard'
# Temporal_resolution = 'DAILY'
# Years = 1998:2024
# # Plumes = list('BAY_OF_BISCAY' = list('Gironde', 'Charente', 'Sevre'), 'SOUTHERN_BRITTANY' = list('Loire', 'Vilaine'),
# #               'GULF_OF_LION' = list('Grand Rhone', 'Petit Rhone'), 'BAY_OF_SEINE' = list('Seine'))
# Plumes = list('BAY_OF_SEINE' = list('Seine'))
# nb_cores = 6

plot_time_series_of_plume_area_and_thresholds <- function(where_are_saved_plume_results, where_to_save_plume_time_series, 
                                                          Zone, Data_source, Satellite_sensor, 
                                                          atmospheric_correction, Temporal_resolution, Years, Plumes, nb_cores) {
  
  registerDoParallel(cores = nb_cores)
  
  cases_to_process <- expand.grid(list("where_are_saved_plume_results" = where_are_saved_plume_results, 
                                       "Zone" = Zone, 
                                       "Data_source" = Data_source, 
                                       "Variable" = "SPM", 
                                       "Satellite_sensor" = Satellite_sensor, 
                                       "atmospheric_correction" = atmospheric_correction, 
                                       'PLUME_DETECTION' = 'PLUME_DETECTION', 
                                       "Temporal_resolution" = Temporal_resolution, 
                                       "Years" = Years %>% as.character()) 
                                  )
  
  ts_data <- cases_to_process %>%
    
    adply(1, function(row) {
      
      file_to_load <- list.files(row %>% as_vector() %>% paste(collapse = "/"), pattern = ".csv", full.names = TRUE, recursive = TRUE )
      
      if (file_to_load %>% is_empty()) {return()}
      
      data <- read.csv(file_to_load, stringsAsFactors = FALSE, colClasses = c("date" = "character")) %>%
        mutate(date = date %>% as.Date(format = "%Y-%m-%d"))
      
      where_to_save_the_file = file_to_load %>% str_replace(where_are_saved_plume_results, where_to_save_plume_time_series) %>% str_remove('/[A-Za-z]*.csv')
      
      make_the_plot_from_the_df(data, row, where_to_save_the_file, "Time_series_of_plume_area_and_SPM_threshold", Temporal_resolution)
      
      return(data)
      
    }, .parallel = FALSE, .inform = TRUE)
  
  # if (nrow(ts_data) == 0) {print("No file found"); return()}
  if (nrow(ts_data) == 0) {print(cases_to_process); return()}
  
  global_cases_to_process <- expand.grid(list("where_are_saved_plume_results" = where_are_saved_plume_results, 
                                              "Zone" = Zone, 
                                              "Data_source" = "*", 
                                              "Variable" = "SPM", 
                                              "Satellite_sensor" = "*", 
                                              "atmospheric_correction" = "*",
                                              'PLUME_DETECTION' = 'PLUME_DETECTION', 
                                              "Temporal_resolution" = Temporal_resolution)
                                         )
  
  plume_names <- ts_data %>% select_at(vars(starts_with('SPM_threshold_'))) %>% names() %>% str_replace_all('SPM_threshold_', '')
  
  global_cases_to_process %>%
    
    a_ply(1, function(row) {
      
      file_path <- row %>% as_vector() %>% paste(collapse = "/")
      
      ts_data_of_the_case <- ts_data %>%
        filter(Zone == row$Zone, Temporal_resolution == row$Temporal_resolution) %>%
        select(-matches( paste0("SPM_threshold_", setdiff(Plumes %>% unlist() %>% str_replace(" ", ".") %>% append(" "), plume_names) ) )) %>% 
        select_if(~sum(!is.na(.)) > 0) %>% 
        # filter_at(vars(starts_with('SPM_threshold')), ~ is.finite(.)) %>% 
        mutate(path_to_file = file_path %>% 
                 str_replace(as.character(where_are_saved_plume_results), where_to_save_plume_time_series) %>% 
                 str_remove_all('[*][A-Za-z0-9/_-]*')) 
      
      sub_cases <- expand.grid( list("Data_source" = unique(ts_data_of_the_case$Data_source), 
                                     "Satellite_sensor" = unique(ts_data_of_the_case$Satellite_sensor), 
                                     "atmospheric_correction" = unique(ts_data_of_the_case$atmospheric_correction)) )
      
      sub_cases %>% a_ply(1, function(sub_case) {
        
        ts_data_of_the_sub_case <- ts_data_of_the_case %>% filter(Data_source == sub_case$Data_source, 
                                                                  Satellite_sensor == sub_case$Satellite_sensor, 
                                                                  atmospheric_correction == sub_case$atmospheric_correction)
        
        where_to_save_sub_case_data <- file.path(ts_data_of_the_sub_case$path_to_file[1], sub_case$Data_source, 'SPM', 
                                                 sub_case$Satellite_sensor, sub_case$atmospheric_correction, 
                                                 row$PLUME_DETECTION, row$Temporal_resolution)
        
        make_the_plot_from_the_df(data = ts_data_of_the_sub_case, row, 
                                  where_to_save_the_file = where_to_save_sub_case_data,
                                  plot_name = 'Time_series_of_plume_area_and_SPM_threshold',
                                  Temporal_resolution)
        
        save_file_as_csv(ts_data_of_the_sub_case, file.path(where_to_save_sub_case_data, 'Time_series_of_plume_area_and_SPM_threshold'))    
        
      })
      
      save_file_as_csv(ts_data_of_the_case,
                       file.path(ts_data_of_the_case$path_to_file[1], row$PLUME_DETECTION,
                                 paste('Time_series_of', row$Temporal_resolution ,'plume_area_and_SPM_threshold', sep = "_")))
      
      make_the_plot_from_the_df(data = ts_data_of_the_case, row, 
                                where_to_save_the_file = file.path(ts_data_of_the_case$path_to_file[1], row$PLUME_DETECTION),
                                plot_name = paste('Time_series_of', row$Temporal_resolution, 'plume_area_and_SPM_threshold', sep = "_"),
                                Temporal_resolution)
      
    }, .parallel = FALSE, .inform = TRUE)
  
}

