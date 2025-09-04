list_of_packages <- c("plyr", "tidyverse", "maps", "scales", "ggpubr")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, require, character.only = TRUE)

#### Main ####
# where_to_save_the_figure <- "/home/terrats/Desktop/RIOMAR/TEST"

Figure_1 <- function(where_to_save_the_figure) {
  
  main_folder_of_Figure_1 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_1")

  SPM_map <- file.path( main_folder_of_Figure_1, "DATA", "SPM_map.csv" ) %>% read_csv()
  insitu_stations <- file.path( main_folder_of_Figure_1, "DATA", "Stations_position.csv" ) %>% read_csv()
  RIOMAR_limits <- file.path( main_folder_of_Figure_1, "DATA", "RIOMAR_limits.csv" ) %>% read_csv() %>% 
                    pivot_longer(-...1, names_to = "Zone", values_to = "Value") %>%
                    pivot_wider(names_from = ...1, values_from = Value)

  basic_map <- create_the_basic_map(map_df = SPM_map, var_name = 'SPM', in_situ_fixed_station = insitu_stations)
  
  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'),
                                      longitude = c(0,0),
                                      latitude = c(0,0))
  
  Figure_1 <- basic_map + 
      geom_point(data = insitu_stations %>% filter(SOURCE == 'REPHY'), 
                 aes(x = LONGITUDE, y = LATITUDE), 
                 fill = "red", color = "black", size = 4, shape = 24, stroke = 1) +
      geom_point(data = insitu_stations %>% filter(SOURCE == 'SOMLIT'), 
                 aes(x = LONGITUDE, y = LATITUDE), 
                 fill = "red", color = "black", size = 10, shape = 21, stroke = 2) + 
      geom_rect(data = RIOMAR_limits, aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max),
                fill = "transparent", color = "red", linetype = "dashed", size = 2) +
      
      geom_point(data = points_for_the_legend, aes(x = longitude, y = latitude, shape = SOURCE), size = 0.1) +
      
      scale_shape_manual(values = c('SOMLIT' = 21, "REPHY" = 24), breaks=c('SOMLIT', 'REPHY'), 
                         labels = c(paste('SOMLIT (n=', length(which(insitu_stations$SOURCE == "SOMLIT")), ")", sep = ""),
                                    paste('REPHY (n=', length(which(insitu_stations$SOURCE == "REPHY")), ")", sep = ""))) +
      guides(
        shape = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                             override.aes = list(size = c(10, 4),
                                                 shape = c(21,24),
                                                 fill = c("red", "red"),
                                                 color = c('black', 'black'),
                                                 stroke = c(2, 1)),
                             order = 1),
        fill = guide_colorbar(barwidth = 30, barheight = 2)) +
      labs(shape = "In-situ stations") + 
      theme(legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(angle = 0, hjust = 0.5),
            legend.spacing.x = unit(5, "cm"))
    
  save_plot_as_png(Figure_1, "Figure_1", width = 18, height = 16, path = main_folder_of_Figure_1)
  
}

# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/'
Figure_2 <- function(where_to_save_the_figure, include_station_points) {
  
  main_folder_of_Figure_2 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_2")
  
  SPM_map_data <- file.path( main_folder_of_Figure_2, "DATA" ) %>% 
                list.files(pattern = "*.csv", full.names = TRUE) %>% 
                llply(read_csv) %>% 
                keep(~ 'analysed_spim' %in% names(.))
  
  insitu_stations <- file.path( main_folder_of_Figure_2, "DATA", "Stations_position.csv" ) %>% read_csv()

  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'), longitude = c(0,0), latitude = c(0,0))
  
  SPM_maps <- SPM_map_data %>% 
                llply(function(x) {
                  
                  insitu_stations_of_the_map <- insitu_stations %>% 
                                                  filter((LATITUDE %>% between(min(x$lat), max(x$lat))) & 
                                                         (LONGITUDE %>% between(min(x$lon), max(x$lon))))
                  
                  the_map <- create_the_basic_map(x, 'SPM', legend_limits = c(4,10))
                  
                  if (include_station_points) {
                  
                    the_map <- the_map + 
                      
                      geom_point(data = insitu_stations_of_the_map %>% filter(SOURCE == 'REPHY'), 
                                 aes(x = LONGITUDE, y = LATITUDE), 
                                 fill = "red", color = "black", size = 6, shape = 24, stroke = 1) +
                      geom_point(data = insitu_stations_of_the_map %>% filter(SOURCE == 'SOMLIT'), 
                                 aes(x = LONGITUDE, y = LATITUDE), 
                                 fill = "red", color = "black", size = 14, shape = 21, stroke = 2) + 
                      
                      geom_point(data = points_for_the_legend, aes(x = longitude, y = latitude, shape = SOURCE), size = 0.1) +
                      
                      scale_shape_manual(values = c('SOMLIT' = 21, "REPHY" = 24), breaks=c('SOMLIT', 'REPHY'), 
                                         labels = c('SOMLIT', 'REPHY')) +
                      guides(
                        shape = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                                             override.aes = list(size = c(14, 6),
                                                                 shape = c(21,24),
                                                                 fill = c("red", "red"),
                                                                 color = c('black', 'black'),
                                                                 stroke = c(2, 1)),
                                             order = 1)) +
                      labs(shape = "In-situ stations") 
                    
                  }
                  
                  the_map <- the_map + 
                                guides(fill = guide_colorbar(barwidth = 45, barheight = 2)) +
                                theme(legend.position = "bottom",
                                      legend.title.position = "top",
                                      legend.title = element_text(angle = 0, hjust = 0.5),
                                      legend.spacing.x = unit(5, "cm"),
                                      axis.text = element_text(size=25, colour = "black"))

                    return(the_map)
                  }) %>% 
                ggarrange(plotlist = ., common.legend = TRUE)
  
  save_plot_as_png(SPM_maps, paste("Figure_2", ifelse(include_station_points, "with_stations", "wo_stations"), sep = "_"), 
                   width = 28, height = 16, path = main_folder_of_Figure_2)
  
}


# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURES/FIGURE_4'
# name_of_the_plot <- "C"
Figure_4 <- function(where_to_save_the_figure, name_of_the_plot) {

  SPM_map_data <- read_csv(file.path( where_to_save_the_figure, "DATA", paste(name_of_the_plot, ".csv", sep = "") ))
  
  if (name_of_the_plot %in% c("A", "B")) {
    the_map <- create_the_basic_map(SPM_map_data, 'SPM', legend_limits = c(4,10))
  } else {
    the_map <- create_the_basic_map(SPM_map_data, 'plume', legend_limits = c(4,10))
  }
    
  the_map <- the_map + 
              guides(fill = guide_colorbar(barwidth = 60, barheight = 2, title.position = "top")) +
              theme(legend.position = "top",
                    legend.title = element_text(angle = 0, hjust = 0.5),
                    axis.text = element_text(size=25, colour = "black"))
  
  if (name_of_the_plot == "B") {
    points_used_for_finding_SPM_threshold <- read_csv(file.path( where_to_save_the_figure, "DATA", "B_points_used_for_finding_SPM_threshold.csv" ))
    all_points_used_for_finding_SPM_threshold <- read_csv(file.path( where_to_save_the_figure, "DATA", "B_all_points_used_for_finding_SPM_threshold.csv" ))
    the_map <- the_map +
      geom_point(data = all_points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "grey50", size = 5) +
      geom_point(data = points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "red", size = 5)
  }
  
  save_plot_as_png(the_map, name_of_the_plot, width = 28, height = 16, path = where_to_save_the_figure)
  
}


# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURES/FIGURE_5'
Figure_5 <- function(where_to_save_the_figure) {
  
  SPM_map_data <- where_to_save_the_figure %>% file.path('DATA') %>% 
                    list.files(pattern = "*.csv", full.names = TRUE) %>% 
                    llply(read_csv)
  
  SPM_maps <- SPM_map_data %>% llply(function(SPM_map) {
    
    create_the_basic_map(SPM_map, 'plume', legend_limits = c(4,10)) + 
      guides(fill = guide_colorbar(barwidth = 60, barheight = 2, title.position = "top")) +
      theme(legend.position = "top",
            legend.title = element_text(angle = 0, hjust = 0.5),
            axis.text = element_text(size=25, colour = "black"))
    
  })
  
  save_plot_as_png(ggarrange(plotlist = SPM_maps, common.legend = TRUE), 
                   'Figure_5', width = 28, height = 16, path = where_to_save_the_figure)
  
}


# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURES/FIGURES_6_7'
Figures_6_7 <- function(where_to_save_the_figure) {
  
  SPM_map_data <- where_to_save_the_figure %>% file.path('DATA', 'ts_data.csv') %>% read_csv()
  SPM_map_data$Dynamic_threshold <- ifelse(SPM_map_data$Dynamic_threshold, 'Dynamic threshold', 'Fixed threshold')
  
  SPM_map_ts <- SPM_map_data %>% filter(Dynamic_threshold == 'Dynamic threshold') %>% dlply(.(Zone), function(df_zone) {
    
    unique_years <- df_zone$Years %>% unique()
    
    points_for_the_legend <- data.frame(Satellite_sensor = c('merged', 'modis'),
                                        date = c('2020-01-01','2020-01-01') %>% as.Date(),
                                        area_of_the_plume_mask_in_km2 = c(-9999,-9999))
    
    index_to_remove <- which((df_zone$Satellite_sensor == "modis") & 
                             (df_zone$area_of_the_plume_mask_in_km2 > quantile(df_zone$area_of_the_plume_mask_in_km2, probs = 0.999, na.rm = TRUE)))
    
    if (index_to_remove %>% length() > 0) {df_zone <- df_zone[-index_to_remove,]}
    
    the_ts_plot_wo_modis <- ggplot() + 
      geom_point(data = df_zone %>% filter(Satellite_sensor == "merged"), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") + 
      geom_path(data = df_zone %>% filter(Satellite_sensor == "merged"), 
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") + 
      
      scale_x_date(name = "", 
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist(),
                   expand = c(0.01,0.01)) +
      
      coord_cartesian(ylim = c(0, max(df_zone$area_of_the_plume_mask_in_km2, na.rm = TRUE))) +
      labs(y = "Plume area (km²)", x = "", title = df_zone$Zone[1] %>% str_replace_all("_", " ")) + 
      ggplot_theme() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = c(.9,.9),
            legend.background = element_rect(fill = "transparent"),
            plot.title = element_text(size=30, colour = "black"),
            text = element_text(size=25, colour = "black"),
            axis.text = element_text(size=20, colour = "black"),
            axis.title = element_text(size=30, colour = "black")) 
    
    
    the_ts_plot_with_modis <- the_ts_plot_wo_modis + 
      geom_point(data = df_zone %>% filter(Satellite_sensor == "modis"), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "blue3", alpha = 0.5) + 
      geom_path(data = df_zone %>% filter(Satellite_sensor == "modis"), 
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "blue3", alpha = 0.5) + 
      geom_point(data = points_for_the_legend, aes(x = date, y = area_of_the_plume_mask_in_km2, color = Satellite_sensor), size = 0.1) +
        
      scale_color_manual(values = c('merged' = "red3", "modis" = "blue3"), name = "") +
        
      guides(
        color = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                             override.aes = list(size = c(5, 5),
                                                 alpha = c(1, 0.5)))) 
        
    return(list("wo_modis" = the_ts_plot_wo_modis, "w_modis" = the_ts_plot_with_modis))
    
  })
    
  save_plot_as_png(ggarrange(plotlist = SPM_map_ts %>% llply(function(x) {x$wo_modis}), common.legend = FALSE, ncol = 1, nrow = 4, align = "v"), 
                   'Figure_6', width = 20, height = 16, path = where_to_save_the_figure)
  
  save_plot_as_png(ggarrange(plotlist = SPM_map_ts %>% llply(function(x) {x$w_modis}), common.legend = FALSE, ncol = 1, nrow = 4, align = "v"), 
                   'Figure_7_merged_modis', width = 20, height = 16, path = where_to_save_the_figure)
  
  
  SPM_map_ts <- SPM_map_data %>% filter(Satellite_sensor == "merged") %>% dlply(.(Zone), function(df_zone) {
    
    unique_years <- df_zone$Years %>% unique()
    
    points_for_the_legend <- data.frame(Dynamic_threshold = c('Dynamic threshold', 'Fixed threshold'),
                                        date = c('2020-01-01','2020-01-01') %>% as.Date(),
                                        area_of_the_plume_mask_in_km2 = c(-9999,-9999))
    
    # index_to_remove <- which((df_zone$Satellite_sensor == "modis") & 
    #                            (df_zone$area_of_the_plume_mask_in_km2 > quantile(df_zone$area_of_the_plume_mask_in_km2, probs = 0.999, na.rm = TRUE)))
    # 
    # if (index_to_remove %>% length() > 0) {df_zone <- df_zone[-index_to_remove,]}
    
    the_ts_plot <- ggplot() + 
      
      geom_point(data = df_zone %>% filter(Dynamic_threshold == 'Dynamic threshold'), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") + 
      geom_path(data = df_zone %>% filter(Dynamic_threshold == 'Dynamic threshold'), 
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") + 
      
      scale_x_date(name = "", 
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist(),
                   expand = c(0.01,0.01)) +
      
      coord_cartesian(ylim = c(0, max(df_zone$area_of_the_plume_mask_in_km2, na.rm = TRUE))) +
      labs(y = "Plume area (km²)", x = "", title = df_zone$Zone[1] %>% str_replace_all("_", " ")) + 
      ggplot_theme() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = c(.9,.9),
            legend.background = element_rect(fill = "transparent"),
            plot.title = element_text(size=30, colour = "black"),
            text = element_text(size=25, colour = "black"),
            axis.text = element_text(size=20, colour = "black"),
            axis.title = element_text(size=30, colour = "black"))  +
      
      geom_point(data = df_zone %>% filter(Dynamic_threshold == 'Fixed threshold'), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "chartreuse4", alpha = 0.5) + 
      geom_path(data = df_zone %>% filter(Dynamic_threshold == 'Fixed threshold'), 
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "chartreuse4", alpha = 0.5) + 
      geom_point(data = points_for_the_legend, aes(x = date, y = area_of_the_plume_mask_in_km2, color = Dynamic_threshold), size = 0.1) +
      
      scale_color_manual(values = c('Dynamic threshold'= "red3", 'Fixed threshold' = "chartreuse4"), name = "") +
      
      guides(
        color = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                             override.aes = list(size = c(5, 5),
                                                 alpha = c(1, 0.5)))) 
    
    return(the_ts_plot)
    
  })
  
  save_plot_as_png(ggarrange(plotlist = SPM_map_ts, common.legend = FALSE, ncol = 1, nrow = 4, align = "v"), 
                   'Figure_7_threshold', width = 20, height = 16, path = where_to_save_the_figure)
  
}



# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURES/FIGURES_8_9_10'
Figures_8_9_10 <- function(where_to_save_the_figure) {
  
  plume_data <- where_to_save_the_figure %>% file.path('DATA', 'ts_plume_data.csv') %>% read_csv()
  river_data <- where_to_save_the_figure %>% file.path('DATA', 'ts_river_data.csv') %>% read_csv()

  regions <- unique(plume_data$Zone) %>% sort()
  
  X11_all_ts <- regions %>% llply(function(region) {
    
    plume_data_region <- plume_data %>% filter(Zone == region)
    river_data_region <- river_data %>% filter(Zone == region)
    
    X11_ts <- plume_data_region %>% inner_join(river_data_region, by = "dates", suffix = c("_plume_area", "_river_flow"))
    
    Raw_plot <- make_the_X11_plot_of_river_and_plume(X11_ts, type_of_signal = 'Raw')
    
    Interannual_plot <- make_the_X11_plot_of_river_and_plume(X11_ts, type_of_signal = 'Interannual')
    
    Seasonal_plot <- make_the_X11_plot_of_river_and_plume(X11_ts, type_of_signal = 'Seasonal')
    
    Residual_plot <- make_the_X11_plot_of_river_and_plume(X11_ts, type_of_signal = 'Residual')
    
    return(list("Interannual" = Interannual_plot + labs(title = region %>% str_replace_all("_", " ")),
                "Seasonal" = Seasonal_plot + labs(title = region %>% str_replace_all("_", " ")),
                "Residual" = Residual_plot + labs(title = region %>% str_replace_all("_", " "))))
    
  })
  
  save_plot_as_png(ggarrange(plotlist = X11_all_ts %>% llply(function(x) {x$Seasonal}), ncol = 1, nrow = 4, align = "v"), 
                   'Figure_8', width = 20, height = 16, path = where_to_save_the_figure)
  
  save_plot_as_png(ggarrange(plotlist = X11_all_ts %>% llply(function(x) {x$Interannual}), ncol = 1, nrow = 4, align = "v"), 
                   'Figure_9', width = 20, height = 16, path = where_to_save_the_figure)
  
  save_plot_as_png(ggarrange(plotlist = X11_all_ts %>% llply(function(x) {x$Residual}), ncol = 1, nrow = 4, align = "v"), 
                   'Figure_10', width = 20, height = 16, path = where_to_save_the_figure)
    
}


#### Utils ####

create_the_basic_map <- function(map_df, var_name, 
                           in_situ_fixed_station = NULL, cruise_stations = NULL, 
                           glider_stations = NULL,
                           legend_limits = NULL) {
  
  if (str_detect(var_name, 'chl|CHL')) {
    title = "[Chl-a]"
    unit = "mg m³"
    if (legend_limits %>% is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  if (str_detect(var_name, 'tsm|SPM|TSM|plume')) {
    title = "[SPM]"
    unit = "g m³"
    # legend_limits <- map_df$analysed_spim[which(map_df$plume)] %>% quantile(probs = c(0.1, 0.9), na.rm = TRUE)
    if (legend_limits %>% is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  FRANCE_shapefile <- map_data('world')[map_data('world')$region == "France",]

  the_base_map <- ggplot() + 
    geom_raster(data = map_df, aes(x = lon, y = lat, fill = analysed_spim), interpolate = FALSE) + 
    scale_fill_viridis_c(na.value = "transparent", option = "viridis", trans = "log10", 
                         limits = c(legend_limits[1], legend_limits[2]), oob = scales::squish, 
                         n.breaks = 5, name = paste(title, " (", unit, ")", sep = "")) +
    guides(fill = guide_colourbar(title.position = "right"))
    
  if (var_name == 'plume') {
    the_base_map <- the_base_map + geom_raster(data = map_df[which(map_df$plume),], aes(x = lon, y = lat), fill = "red", interpolate = FALSE) 
  }
  
  the_map <- the_base_map + 
    
    ## First layer: worldwide map
    geom_polygon(data = map_data("world"), aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black') +
    ## Second layer: Country map
    geom_polygon(data = FRANCE_shapefile, aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black') +
    coord_cartesian(xlim = range(map_df$lon), ylim = range(map_df$lat), expand = FALSE) +
    
    scale_x_continuous(name = "", labels = function(x) paste(x, "°E", sep = "")) +
    scale_y_continuous(name = "", labels = function(x) paste(x, "°N", sep = "")) +
    ggplot_theme() + 
    
    theme(plot.title = element_text(size = 45),
          legend.position = "right",
          legend.title = element_text(angle = -90, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key.height = unit(6, "lines"),
          legend.key.width = unit(3, "lines")) 
  
  return(the_map)
  
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

make_the_X11_plot_of_river_and_plume <- function(X11_data, type_of_signal) {
  
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
          plot.title = element_text(size=30, colour = "black"),
          text = element_text(size=25, colour = "black"),
          axis.text = element_text(size=20, colour = "black"),
          axis.title = element_text(size=30, colour = "black")) +
    
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

