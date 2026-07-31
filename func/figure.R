# func/figure.R
# This is called by func/figure.py to provide graphical outputs


# Libraries ---------------------------------------------------------------

# list_of_packages <- c("plyr", "tidyverse", "maps", "scales", "ggpubr")
# new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
# lapply(list_of_packages, require, character.only = TRUE)

library(plyr)
library(tidyverse)
library(scales)
library(maps)
library(ggpubr)
library(cowplot)
library(magick)

# For zones_bbox -- the authoritative zone boundaries (also used by
# func/validate.R for actual zone classification), as opposed to util.py's
# separate, wider lat/lon_range_of_the_map_to_plot (map-cropping extent only).
source("func/util.R")


# Tests -------------------------------------------------------------------

# test1 <- tidync::tidync("~/pCloudDrive/data/WIND/BAY_OF_BISCAY/wind_202501_202509.nc") |> tidync::hyper_tibble()
# test2 <- tidync::tidync("~/pCloudDrive/data/WIND/BAY_OF_BISCAY/wind_daily_202501_202509.nc") |> tidync::hyper_tibble()
# test3 <- tidync::tidync("~/pCloudDrive/data/GLORYS/BAY_OF_BISCAY/glorys_202107_202507.nc") |> tidync::hyper_tibble()


# Utils -------------------------------------------------------------------

create_the_basic_map <- function(map_df, var_name, 
                                 in_situ_fixed_station = NULL, 
                                 cruise_stations = NULL, 
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


# Main --------------------------------------------------------------------

# where_to_save_the_figure <- "~/RiOMar/figures/"

Figure_1 <- function(where_to_save_the_figure) {

  main_folder_of_Figure_1 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_1")

  SPM_map <- file.path( main_folder_of_Figure_1, "DATA", "SPM_map.csv" ) %>% read_csv()
  insitu_stations <- file.path( main_folder_of_Figure_1, "DATA", "Stations_position.csv" ) %>% read_csv()

  # Zone boundaries: use util.R's zones_bbox (the boxes func/validate.R
  # actually classifies stations/match-ups against), NOT the separate,
  # wider lat_range_of_the_map_to_plot/lon_range_of_the_map_to_plot in
  # util.py (used only for cropping the regional map extent to load/plot,
  # not the zone definition itself) -- the old RIOMAR_limits.csv was built
  # from the latter and had drifted out of sync with the former.
  RIOMAR_limits <- zones_bbox %>% dplyr::rename(Zone = zone)

  basic_map <- create_the_basic_map(map_df = SPM_map, var_name = 'SPM', in_situ_fixed_station = insitu_stations)

  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'),
                                      longitude = c(0,0),
                                      latitude = c(0,0))

  national_map <- basic_map +
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

  # Zoomed regional insets, floating near each river mouth -------------------
  # Manuscript Figure 2 is the satellite-vs-in-situ validation scatterplot
  # (produced by func/validate.R during 1_validate.py) -- these insets are
  # the former standalone Figure_2(), now folded into Figure_1() per Robert's
  # request. Each inset crops the SAME national SPM_map grid used for the
  # main map (not a separate per-zone snapshot date), so the insets and the
  # backdrop are guaranteed to show the same underlying data.
  # River name shown in the corner of each inset, so it's identifiable
  # without having to trace the dashed line back to the national map.
  zone_river_labels <- c(BAY_OF_SEINE = "Seine", SOUTHERN_BRITTANY = "Loire",
                         BAY_OF_BISCAY = "Gironde", GULF_OF_LION = "Rhône")

  build_zone_inset <- function(zone_name) {
    limits <- RIOMAR_limits %>% dplyr::filter(Zone == zone_name)

    zone_SPM <- SPM_map %>%
      dplyr::filter(lon >= limits$lon_min, lon <= limits$lon_max,
                    lat >= limits$lat_min, lat <= limits$lat_max)

    zone_stations <- insitu_stations %>%
      dplyr::filter(LONGITUDE >= limits$lon_min, LONGITUDE <= limits$lon_max,
                    LATITUDE >= limits$lat_min, LATITUDE <= limits$lat_max)

    create_the_basic_map(zone_SPM, 'SPM') +
      geom_point(data = zone_stations %>% dplyr::filter(SOURCE == 'REPHY'),
                 aes(x = LONGITUDE, y = LATITUDE),
                 fill = "red", color = "black", size = 2, shape = 24, stroke = 0.6) +
      geom_point(data = zone_stations %>% dplyr::filter(SOURCE == 'SOMLIT'),
                 aes(x = LONGITUDE, y = LATITUDE),
                 fill = "red", color = "black", size = 3.5, shape = 21, stroke = 0.8) +
      ggtitle(zone_river_labels[[zone_name]]) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(size = 28, face = "bold", hjust = 0.5,
                                      colour = "black", margin = margin(b = 4)),
            plot.margin = margin(t = 8, r = 6, b = 6, l = 6),
            plot.background = element_rect(fill = "white", colour = "red", linewidth = 1.8))
  }

  # x/y = bottom-left corner of each inset, w/h = width/height, all as
  # fractions of the whole canvas (cowplot::draw_plot() convention), over the
  # empty land in the middle of the map. Deliberately irregular sizing/
  # spacing rather than an even grid, per Robert: Seine upper right, Rhone
  # underneath it, Loire where Seine used to sit, Gironde roughly in place --
  # shifted further right than that as a whole group so no inset covers any
  # coloured SPM pixels (checked against zones_bbox's true-box fractions).
  inset_layout <- tibble::tribble(
    ~Zone,               ~x,    ~y,    ~w,    ~h,
    "BAY_OF_SEINE",       0.66,  0.62,  0.26,  0.22,
    "SOUTHERN_BRITTANY",  0.43,  0.49,  0.21,  0.22,
    "BAY_OF_BISCAY",      0.46,  0.24,  0.23,  0.21,
    "GULF_OF_LION",       0.72,  0.34,  0.24,  0.20
  )

  Figure_1 <- ggdraw() + draw_plot(national_map, x = 0, y = 0, width = 1, height = 1)

  for (i in seq_len(nrow(inset_layout))) {
    Figure_1 <- Figure_1 +
      draw_plot(build_zone_inset(inset_layout$Zone[i]),
                x = inset_layout$x[i], y = inset_layout$y[i],
                width = inset_layout$w[i], height = inset_layout$h[i])
  }

  save_plot_as_png(Figure_1, "Figure_1", width = 28, height = 22, path = main_folder_of_Figure_1)

}

# Builds the 4-zone regional SPM map grid (optionally with SOMLIT/REPHY
# station overlay), shared by Figure_1() (with stations, embedded as insets)
# and regional_zone_maps() below (without stations, a standalone reference
# figure used by Figure_5()).
regional_zone_maps_panels <- function(data_folder, include_station_points) {

  SPM_map_data <- file.path(data_folder, "DATA") %>%
    list.files(pattern = "*.csv", full.names = TRUE) %>%
    llply(read_csv) %>%
    keep(~ 'analysed_spim' %in% names(.))

  insitu_stations <- file.path( data_folder, "DATA", "Stations_position.csv" ) %>% read_csv()

  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'), longitude = c(0,0), latitude = c(0,0))

  SPM_map_data %>%
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

          scale_shape_manual(values = c('SOMLIT' = 21, "REPHY" = 24),
                             breaks = c('SOMLIT', 'REPHY'),
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
}

# Standalone regional-zone-maps figure, kept for Figure_5()'s
# without-stations reference-map byproduct. The with-stations version is no
# longer produced standalone -- it is embedded directly in Figure_1().
regional_zone_maps <- function(where_to_save_the_figure, include_station_points) {

  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_1")

  SPM_maps <- regional_zone_maps_panels(main_folder, include_station_points)

  save_plot_as_png(SPM_maps, paste("regional_zone_maps", ifelse(include_station_points, "with_stations", "wo_stations"), sep = "_"),
                   width = 28, height = 16, path = main_folder)

}


# Manuscript Figure 2: satellite-vs-in-situ validation scatterplots, panel
# (a) SPM and panel (b) Chl-a. Combines the two already-rendered, 3x3
# grid-size (spatially averaged over the 3x3 pixel window around each in
# situ station -- not the 1x1 single-native-pixel version) SEXTANT
# scatterplots from func/validate.py's match-up pipeline (run during
# 1_validate.py) -- does not regenerate the scatterplots themselves, so
# re-tweaking them there and re-running 1_validate.py is picked up
# automatically next time this runs.
Figure_2 <- function(spm_scatterplot_path, chla_scatterplot_path, where_to_save_the_figure) {

  main_folder_of_Figure_2 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_2")
  if (!dir.exists(main_folder_of_Figure_2)) dir.create(main_folder_of_Figure_2, recursive = TRUE)

  label_panel <- function(path, label) {
    image_annotate(image_read(path), label, size = 150, weight = 700,
                   gravity = "northwest", location = "+40+20", color = "black")
  }

  panel_a <- label_panel(spm_scatterplot_path, "a)")
  panel_b <- label_panel(chla_scatterplot_path, "b)")

  # Side-by-side (not stacked): each panel is a square 4800x4800 scatterplot
  # grid, so stacking them made a 1:2 aspect-ratio image too tall to fit on
  # one manuscript page alongside its caption. Side-by-side gives a 2:1
  # image instead, at the same per-panel resolution.
  Figure_2 <- image_append(c(panel_a, panel_b), stack = FALSE)

  image_write(Figure_2, file.path(main_folder_of_Figure_2, "Figure_2.png"))

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

  # Plume-area trend line (manuscript Figure 4) uses the same AR(1)/HAC-weighted
  # fit as Table 5's "Surface area" row -- see func/compute_area_trend.R.
  area_trend <- read_csv("output/STATS/area_trend_summary.csv", show_col_types = FALSE)

  SPM_map_ts <- SPM_map_data %>% filter(Dynamic_threshold == 'Dynamic threshold') %>% dlply(.(Zone), function(df_zone) {

    unique_years <- df_zone$Years %>% unique()

    points_for_the_legend <- data.frame(Satellite_sensor = c('merged', 'modis'),
                                        date = c('2020-01-01','2020-01-01') %>% as.Date(),
                                        area_of_the_plume_mask_in_km2 = c(-9999,-9999))

    index_to_remove <- which((df_zone$Satellite_sensor == "modis") &
                               (df_zone$area_of_the_plume_mask_in_km2 > quantile(df_zone$area_of_the_plume_mask_in_km2, probs = 0.999, na.rm = TRUE)))

    if (index_to_remove %>% length() > 0) {df_zone <- df_zone[-index_to_remove,]}

    zone_trend <- area_trend %>% filter(zone == df_zone$Zone[1]) %>%
      mutate(trend_label = paste0("Trend: ", sprintf("%+.2f", slope_annualised), " km² yr⁻¹ (",
                                  ifelse(slope_p < 0.001, "p < 0.001", paste0("p = ", signif(slope_p, 2))), ")"))

    df_merged <- df_zone %>% filter(Satellite_sensor == "merged")
    spm_scaling <- sec_axis_adjustement_factors(var_to_scale = df_merged$mean_SPM_in_the_plume_area,
                                                var_ref = df_merged$area_of_the_plume_mask_in_km2)
    df_merged <- df_merged %>% mutate(SPM_scaled = mean_SPM_in_the_plume_area * spm_scaling$diff + spm_scaling$adjust)

    the_ts_plot_wo_modis <- ggplot() +
      geom_point(data = df_merged,
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") +
      geom_path(data = df_merged,
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") +
      geom_point(data = df_merged, aes(x = date, y = SPM_scaled), color = "steelblue", alpha = 0.4) +
      geom_path(data = df_merged, aes(x = date, y = SPM_scaled), color = "steelblue", alpha = 0.4) +
      geom_abline(data = zone_trend, aes(intercept = intercept, slope = slope),
                  colour = "black", linewidth = 1.2, linetype = "dashed") +
      geom_text(data = zone_trend, aes(x = -Inf, y = Inf, label = trend_label), inherit.aes = FALSE,
               hjust = -0.03, vjust = 1.5, size = 6, colour = "black") +

      scale_x_date(name = "",
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist(),
                   expand = c(0.01,0.01)) +

      coord_cartesian(ylim = c(0, max(df_zone$area_of_the_plume_mask_in_km2, na.rm = TRUE))) +
      scale_y_continuous(name = "Plume area (km²)",
                         sec.axis = sec_axis(transform = ~ {. - spm_scaling$adjust} / spm_scaling$diff,
                                            name = "Mean SPM in plume (g m⁻³)")) +
      labs(x = "", title = df_zone$Zone[1] %>% str_replace_all("_", " ")) +
      ggplot_theme() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.subtitle = element_text(hjust = 0.5),
            plot.margin = margin(t = 20, r = 10, b = 5, l = 5),
            legend.position = c(.9,.9),
            legend.background = element_rect(fill = "transparent"),
            plot.title = element_text(size=30, colour = "black"),
            text = element_text(size=25, colour = "black"),
            axis.text = element_text(size=20, colour = "black"),
            axis.title.y.right = element_text(size = 22, colour = "steelblue"),
            axis.text.y.right = element_text(colour = "steelblue"),
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


# manuscript Figure 5: daily plume area vs. river flow, one row per zone --
# a) scatter with a linear trend line, b) lagged correlation (plume lagged
# behind flow, max_lag_daily days). Reuses func/multi.R's
# combine_plume_driver()/driver_plume_correlation() (same machinery behind
# run_driver_suite("flow")'s cor_plot_flow_plume_*.png files), but keeps only
# the scatter + lag panels -- the raw driver/plume time series (panels a/b of
# plot_driver_plume_comparison()) would duplicate manuscript Figure 4 -- and
# caps the lag search at 14 days (vs. that function's default 30) per
# Robert's request. Distinct from the Figure_5() function above (regional
# zone maps feeding the Figure 3 composite): that name was already taken
# before this figure existed, and manuscript.tex's own "Figure 5" caption
# had no generator at all until now (see [[project_manuscript_pipeline_gaps]]
# item 4).
Figure_5_driver_comparison <- function(where_to_save_the_figure, max_lag_daily = 14){

  # Sourced here rather than at file scope: multi.R pulls in tidync/ncdf4
  # (for load_plume_surface(), unused by this function) which fails to load
  # under rpy2's embedded R in this environment, and would otherwise break
  # every other figure.R function called from Python via source(figure.R).
  source("func/multi.R")

  main_folder_of_Figure_5_driver <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_5")
  if (!dir.exists(main_folder_of_Figure_5_driver)) dir.create(main_folder_of_Figure_5_driver, recursive = TRUE)

  panels <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    df <- combine_plume_driver("flow", meta)

    cor_df <- driver_plume_correlation(df, max_lag_daily = max_lag_daily) %>%
      dplyr::filter(timestep == "daily")
    peak <- cor_df %>% dplyr::slice_max(cor, n = 1)

    # ggplot_theme()'s font sizes are tuned for the single-column, 4-row
    # figures elsewhere (e.g. Figures_6_7); an 8-panel 2x4 grid needs smaller
    # text and explicit margins so axis titles don't clip against the plot
    # edge or the panel above.
    panel_theme <- theme(plot.title = element_text(hjust = 0.5, size = 20),
                         axis.title = element_text(size = 15, colour = "black"),
                         axis.text = element_text(size = 12, colour = "black"),
                         plot.margin = margin(t = 8, r = 12, b = 5, l = 5),
                         panel.background = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(linetype = "solid", fill = NA))

    scatter_plot <- ggplot(df, aes(x = value, y = plume_area)) +
      geom_point(alpha = 0.3, colour = "grey30", size = 0.8) +
      geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1.2) +
      labs(x = "River flow (m³ s⁻¹)", y = "Plume area (km²)", title = meta$mouth_name) +
      panel_theme

    lag_plot <- ggplot(cor_df, aes(x = lag, y = cor)) +
      geom_line(colour = "grey30") +
      geom_point(colour = "grey30") +
      geom_point(data = peak, colour = "firebrick", size = 3) +
      geom_text(data = peak, aes(label = paste0(lag, "d")), vjust = -1.2, colour = "firebrick", size = 4) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(x = "Lag, plume after flow (days)", y = "Correlation (r)", title = meta$mouth_name) +
      panel_theme

    list(scatter = scatter_plot, lag = lag_plot)
  })

  plotlist <- panels %>% purrr::map(function(p) list(p$scatter, p$lag)) %>% purrr::flatten()

  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(panels), align = "v")

  save_plot_as_png(full_plot, "Figure_5", width = 16, height = 18, path = main_folder_of_Figure_5_driver)
}


# manuscript Figure 7: wind and wave direction/magnitude roses, one row per
# zone, coloured by the flow-controlled plume-area response
# (multi.R::plot_driver_rose()). Replaces the original Figure 7/8 concept
# (X11-decomposed wind/wave magnitude time series) -- Robert's call, since
# wind and wave direction matter as much or more than magnitude, which a
# magnitude-only time series can't show.
Figure_7_driver_rose <- function(where_to_save_the_figure, n_sectors = 16){

  source("func/multi.R")  # see Figure_5_driver_comparison() for why this is lazy, not file-scope

  main_folder_of_Figure_7 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_7")
  if (!dir.exists(main_folder_of_Figure_7)) dir.create(main_folder_of_Figure_7, recursive = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    list(wind = plot_driver_rose("wind", meta, n_sectors), wave = plot_driver_rose("wave", meta, n_sectors))
  }) %>% purrr::map(function(p) list(p$wind, p$wave)) %>% purrr::flatten()

  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = nrow(zone_meta), align = "v")

  save_plot_as_png(full_plot, "Figure_7", width = 16, height = 20, path = main_folder_of_Figure_7)
}


# manuscript Figure 8: flow-controlled plume-area residual vs. wave height,
# coloured by on/off-shore wind category, one panel per zone
# (multi.R::plot_driver_category_scatter()). Replaces the original Figure 8
# concept (X11-decomposed wave-height magnitude time series).
Figure_8_driver_category <- function(where_to_save_the_figure){

  source("func/multi.R")  # see Figure_5_driver_comparison() for why this is lazy, not file-scope

  main_folder_of_Figure_8 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_8")
  if (!dir.exists(main_folder_of_Figure_8)) dir.create(main_folder_of_Figure_8, recursive = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    plot_driver_category_scatter(meta)
  })

  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

  save_plot_as_png(full_plot, "Figure_8", width = 14, height = 12, path = main_folder_of_Figure_8)
}


# manuscript Figure 9: GAM partial-dependence curves for flow, wind, wave,
# and current, one row per zone (driver_interactions.R::fit_gam()/
# gam_partial_effect()). Tide is intentionally excluded from the plot --
# Robert's call, since tidal range is an essentially fixed astronomical
# property of each site rather than something worth a dedicated panel; it
# stays in the underlying GAM/Table 6 statistics, just not visualised here.
# Reuses the same corrected (ROFI-free) GAM already needed for Fig. 10/
# Table 6 rather than fitting a separate model -- refit here from the
# already-saved daily_driver_matrix_<zone>.csv (Stage 4 output, see
# func/driver_interactions.R::run_full_analysis()) instead of rerunning the
# full driver-interactions pipeline (GLM/regime/RF steps not needed here).
Figure_9_gam_partial <- function(where_to_save_the_figure, stats_dir = "output/STATS"){

  source("func/driver_interactions.R")  # see Figure_5_driver_comparison() for why this is lazy, not file-scope

  main_folder_of_Figure_9 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_9")
  if (!dir.exists(main_folder_of_Figure_9)) dir.create(main_folder_of_Figure_9, recursive = TRUE)

  driver_labels <- c(flow = "River flow (m³ s⁻¹)", wind_spd = "Wind speed (m s⁻¹)",
                     wave_height = "Wave height (m)", current = "Current speed (m s⁻¹)")

  plotlist <- purrr::map(zones, function(zone_name){
    df <- readr::read_csv(file.path(stats_dir, paste0("daily_driver_matrix_", zone_name, ".csv")), show_col_types = FALSE)
    gam_model <- fit_gam(df)
    drivers_to_show <- intersect(names(driver_labels), available_drivers(df))

    purrr::imap(drivers_to_show, function(d, i){
      curve <- gam_partial_effect(gam_model, d, df)
      ggplot(curve, aes(x = x, y = fit)) +
        geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), fill = "grey80", alpha = 0.5) +
        geom_line(colour = "black", linewidth = 1) +
        labs(x = driver_labels[[d]], y = "Partial effect on plume area (km²)",
            title = if(i == 1) stringr::str_replace_all(zone_name, "_", " ") else NULL) +
        theme(panel.border = element_rect(fill = NA, colour = "black"))
    })
  }) %>% purrr::flatten()

  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = length(driver_labels), nrow = length(zones), align = "v")

  save_plot_as_png(full_plot, "Figure_9", width = 20, height = 16, path = main_folder_of_Figure_9)
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

