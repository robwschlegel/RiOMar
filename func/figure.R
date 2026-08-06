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

# For zone_meta/get_zone_meta/combine_plume_driver/plot_driver_rose/etc,
# used by Figure_5_driver_comparison()/Figure_7_driver_rose()/
# Figure_8_driver_category() below. (func/driver_interactions.R, needed only
# by Figure_9_gam_partial(), is sourced lazily inside that one function
# instead -- it pulls in several heavyweight modelling packages (mgcv,
# gratia, ranger, iml) not worth loading for every other figure in this file.)
source("func/multi.R")


# Tests -------------------------------------------------------------------

# test1 <- tidync::tidync("~/pCloudDrive/data/WIND/BAY_OF_BISCAY/wind_202501_202509.nc") |> tidync::hyper_tibble()
# test2 <- tidync::tidync("~/pCloudDrive/data/WIND/BAY_OF_BISCAY/wind_daily_202501_202509.nc") |> tidync::hyper_tibble()
# test3 <- tidync::tidync("~/pCloudDrive/data/GLORYS/BAY_OF_BISCAY/glorys_202107_202507.nc") |> tidync::hyper_tibble()


# Utils -------------------------------------------------------------------

# High-resolution France coastline (GADM level-0 boundary, ~216,000 vertices)
# for Figure 1's zoomed insets -- Robert asked for this in place of
# maps::map_data("world")'s coastline, which is visibly blocky once cropped
# to a single zone. Same file panache's own zone configs already point to as
# coast_shapefile (metadata/zone_config_dynamic_*.json), so this isn't a new
# input, just reusing what the real detection pipeline already trusts for
# coastline masking.
#
# Loaded lazily (library(sf) is NOT called at file scope, unlike the other
# packages above) and cached in .high_res_france_coastline_cache, rather than
# read once eagerly when this file is sourced: figure.R is source()'d for
# every manuscript figure via func/figure.py's rpy2 calls, and loading `sf`
# unconditionally broke every one of them, not just Figure 1 -- the RiOMar
# conda env bundles an older GEOS (3.10.6) than the system R/GDAL stack
# expects (3.12.1), and rpy2 embeds R inside the conda Python process, so R's
# dynamic linker picks up the conda env's older GEOS and sf's dyn.load() then
# fails with "undefined symbol: GEOSConcaveHull_r" against the system
# libgdal.so.34 -- confirmed 2026-08-02 by reproducing the failure via
# Figure_1() called through func/figure.py (rpy2), while a plain
# `Rscript -e 'source("func/figure.R"); Figure_1(...)'` (no conda Python
# involved) works fine. This lazy load means only Figure_1()'s insets are
# still affected when run through the real python code/5_figures.py
# pipeline, until that conda/GEOS conflict is resolved -- see the pipeline
# map's "Still open" list.
.high_res_france_coastline_cache <- NULL
high_res_coastline <- function(){
  if(is.null(.high_res_france_coastline_cache)){
    library(sf)
    .high_res_france_coastline_cache <<- sf::st_read("data/FRANCE_shapefile/gadm41_FRA_0.shp", quiet = TRUE) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      dplyr::transmute(long = X, lat = Y, group = interaction(L1, L2, L3, drop = TRUE))
  }
  .high_res_france_coastline_cache
}

create_the_basic_map <- function(map_df, var_name,
                                 in_situ_fixed_station = NULL,
                                 cruise_stations = NULL,
                                 glider_stations = NULL,
                                 legend_limits = NULL,
                                 log_scale = TRUE,
                                 high_res_coast = FALSE) {
  
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
    scale_fill_viridis_c(na.value = "transparent", option = "viridis", trans = if(log_scale) "log10" else "identity",
                         limits = c(legend_limits[1], legend_limits[2]), oob = scales::squish,
                         n.breaks = 5, name = paste(title, " (", unit, ")", sep = "")) +
    guides(fill = guide_colourbar(title.position = "right"))
  
  if (var_name == 'plume') {
    the_base_map <- the_base_map + geom_raster(data = map_df[which(map_df$plume),], aes(x = lon, y = lat), fill = "red", interpolate = FALSE) 
  }
  
  the_map <- the_base_map +

    (if(high_res_coast){
      ## High-resolution France coastline (GADM), for zoomed insets where
      ## the crude "world" database coastline is visibly blocky -- other
      ## countries aren't relevant at inset zoom levels, so no world layer.
      geom_polygon(data = high_res_coastline(), aes(x = long, y = lat, group = group), color = 'grey60', fill = 'black')
    } else {
      list(
        ## First layer: worldwide map
        geom_polygon(data = map_data("world"), aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black'),
        ## Second layer: Country map
        geom_polygon(data = FRANCE_shapefile, aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black')
      )
    }) +
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


# ggplot_theme(), save_plot_as_png(), and sec_axis_adjustement_factors()
# used to be redefined here identically (ggplot_theme differed slightly --
# see func/util.R's copy for the note); now sourced from func/util.R above.

plot_x11_river_and_plume <- function(X11_data, type_of_signal) {
  
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

  main_folder_of_Figure_1 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_1")

  SPM_map <- file.path( main_folder_of_Figure_1, "DATA", "SPM_map.csv" ) %>% read_csv()
  insitu_stations <- file.path( main_folder_of_Figure_1, "DATA", "Stations_position.csv" ) %>% read_csv()

  # Zone boundaries: use util.R's zones_bbox (the boxes func/validate.R
  # actually classifies stations/match-ups against), NOT the separate,
  # wider lat_range_of_the_map_to_plot/lon_range_of_the_map_to_plot in
  # util.py (used only for cropping the regional map extent to load/plot,
  # not the zone definition itself) -- the old RIOMAR_limits.csv was built
  # from the latter and had drifted out of sync with the former.
  RIOMAR_limits <- zones_bbox %>% dplyr::rename(Zone = zone)

  basic_map <- create_the_basic_map(map_df = SPM_map, var_name = 'SPM', in_situ_fixed_station = insitu_stations, log_scale = FALSE)

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
  # request. Changed 2026-08-01 (per Robert): each inset now reads its own
  # zone-specific date-of-maximum-plume-area CSV (figure.py::
  # dates_for_each_zone(), written by load_the_regional_maps_and_save_them_
  # for_plotting() into DATA/<zone>.csv) rather than cropping the national
  # SPM_map (now the temporal-mean field used only for the backdrop) -- the
  # point of using each zone's biggest-plume day is to actually show that,
  # which the shared national map's date/mean can no longer do.
  # Zone name shown in the corner of each inset (previously the single
  # representative river's name -- changed 2026-08-01 per Robert, since a
  # river name mislabels zones with more than one contributing river).
  # River-mouth coordinates below, marked with a red x on each inset per
  # Robert's request, are the exact starting_points panache.utils.define_parameters()
  # uses for plume detection (verified directly, not RiOMar's older,
  # now-unmaintained local copy of that function) -- several zones have more
  # than one contributing river (Southern Brittany: Loire, Vilaine; Bay of
  # Biscay: Gironde, Charente, Sevre; Gulf of Lion: Grand Rhone, Petit Rhone).
  zone_river_mouths <- tibble::tribble(
    ~zone,               ~river,         ~lat,   ~lon,
    "BAY_OF_SEINE",       "Seine",        49.43,  0.145,
    "SOUTHERN_BRITTANY",  "Loire",        47.24, -2.2,
    "SOUTHERN_BRITTANY",  "Vilaine",      47.48, -2.55,
    "BAY_OF_BISCAY",      "Gironde",      45.61, -1.14,
    "BAY_OF_BISCAY",      "Charente",     45.98, -1.15,
    "BAY_OF_BISCAY",      "Sevre",        46.26, -1.2,
    "GULF_OF_LION",       "Grand Rhone",  43.32,  4.85,
    "GULF_OF_LION",       "Petit Rhone",  43.45,  4.39
  )

  build_zone_inset <- function(zone_name) {
    zone_SPM <- file.path(main_folder_of_Figure_1, "DATA", paste0(zone_name, ".csv")) %>% read_csv()
    mouths <- zone_river_mouths %>% dplyr::filter(zone == zone_name)

    create_the_basic_map(zone_SPM, 'SPM', log_scale = FALSE, high_res_coast = TRUE) +
      geom_point(data = mouths, aes(x = lon, y = lat), shape = 4, colour = "red", size = 4, stroke = 2) +
      ggtitle(zone_title(zone_name)) +
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
zone_maps_panels <- function(data_folder, include_station_points) {

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

  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_1")

  SPM_maps <- zone_maps_panels(main_folder, include_station_points)

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
Figure_2 <- function(spm_scatterplot_path, turb_scatterplot_path, where_to_save_the_figure) {

  main_folder_of_Figure_2 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_2")
  if (!dir.exists(main_folder_of_Figure_2)) dir.create(main_folder_of_Figure_2, recursive = TRUE)

  label_panel <- function(path, label) {
    image_annotate(image_read(path), label, size = 150, weight = 700,
                   gravity = "northwest", location = "+40+20", color = "black")
  }

  panel_a <- label_panel(spm_scatterplot_path, "a)")
  panel_b <- label_panel(turb_scatterplot_path, "b)")

  # Side-by-side (not stacked): each panel is a square 4800x4800 scatterplot
  # grid, so stacking them made a 1:2 aspect-ratio image too tall to fit on
  # one manuscript page alongside its caption. Side-by-side gives a 2:1
  # image instead, at the same per-panel resolution.
  Figure_2 <- image_append(c(panel_a, panel_b), stack = FALSE)

  image_write(Figure_2, file.path(main_folder_of_Figure_2, "Figure_2.png"))

}


# Renders one methodology panel (A-E) for manuscript Figure 3 -- renamed
# from Figure_4 (2026-08-01): it has never produced manuscript Figure 4,
# that name was left over from before the manuscript was renumbered.
# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURE_3'
# name_of_the_plot <- "C"
Figure_3_panel <- function(where_to_save_the_figure, name_of_the_plot) {
  
  SPM_map_data <- read_csv(file.path( where_to_save_the_figure, "DATA", paste(name_of_the_plot, ".csv", sep = "") ))
  
  # legend_limits matches Figure_3_zone_maps()'s panels E-H (0.1-10 g/m3,
  # displayed/labelled as "0-10") -- this row's own colour bar was removed
  # 2026-08-05 (per Robert) since it duplicated across all four panels, so
  # A-D now rely on the shared E-H legend below; that only works if both
  # rows are on the same colour scale, hence matching the range here too
  # (was c(4,10), inconsistent with E-H's 0-10 even before the legend was
  # hidden). Lower bound is 0.1, not 0: create_the_basic_map()'s fill scale
  # is log10-transformed, on which log10(0) = -Inf breaks ggplot's break
  # computation -- 0.1 is the same near-zero floor already used as this
  # function's own default legend_limits elsewhere in this file.
  if (name_of_the_plot %in% c("A", "B")) {
    the_map <- create_the_basic_map(SPM_map_data, 'SPM', legend_limits = c(0.1,10))
  } else {
    the_map <- create_the_basic_map(SPM_map_data, 'plume', legend_limits = c(0.1,10))
  }
  
  # No per-panel colour bar or lon/lat axis text on this top methodology
  # row (2026-08-05, per Robert): each of A-D previously carried its own
  # copy of the same SPM legend, which visually duplicated across the row;
  # the zone maps panel below already carries the shared legend, and the
  # lon/lat axis ticks are already shown there too.
  the_map <- the_map +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.tag = element_text(size = 35, face = "bold"),
          plot.tag.position = c(0.02, 0.98)) +
    labs(tag = name_of_the_plot)

  if (name_of_the_plot == "B") {
    points_used_for_finding_SPM_threshold <- read_csv(file.path( where_to_save_the_figure, "DATA", "B_points_used_for_finding_SPM_threshold.csv" ))
    all_points_used_for_finding_SPM_threshold <- read_csv(file.path( where_to_save_the_figure, "DATA", "B_all_points_used_for_finding_SPM_threshold.csv" ))
    the_map <- the_map +
      geom_point(data = all_points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "grey50", size = 5) +
      geom_point(data = points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "red", size = 5)
  }
  
  save_plot_as_png(the_map, name_of_the_plot, width = 28, height = 16, path = where_to_save_the_figure)
  
}


# Renders the per-zone plume-maps panel feeding manuscript Figure 3 --
# renamed from Figure_5 (2026-08-01): it produces zone_maps_panel.png, not
# manuscript Figure 5 (that's Figure_5_driver_comparison()'s Figure_5.png).
# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURE_3'
Figure_3_zone_maps <- function(where_to_save_the_figure) {

  # Read only the four per-zone SPM-map CSVs figure.py's Figure_3_zone_maps()
  # writes here (Zone.csv, via zone_meta$zone for canonical zone naming/order)
  # -- a plain "*.csv" glob on this shared FIGURE_3/DATA/ folder also picks up
  # Figure_3_panels()' A-E.csv (which lack a `plume` column entirely) and its
  # *_threshold.csv debug files, crashing create_the_basic_map()'s
  # which(map_df$plume) on the first file missing that column.
  SPM_map_data <- where_to_save_the_figure %>% file.path('DATA', paste0(zone_meta$zone, ".csv")) %>%
    llply(read_csv)

  # Continues the lettering from Figure_3_panel()'s methodology row (A-D)
  # -- zone_meta$zone is already arranged by ZONE_ORDER (north to south:
  # Seine, Southern Brittany, Bay of Biscay, Gulf of Lion), matching the
  # order these zones are listed in the Figure 3 caption.
  panel_letters <- c("E", "F", "G", "H")
  SPM_maps <- purrr::map2(SPM_map_data, panel_letters, function(SPM_map, letter) {

    create_the_basic_map(SPM_map, 'plume', legend_limits = c(0.1,10)) +
      guides(fill = guide_colorbar(barwidth = 60, barheight = 2, title.position = "top")) +
      theme(legend.position = "top",
            legend.title = element_text(angle = 0, hjust = 0.5),
            axis.text = element_text(size=25, colour = "black"),
            plot.tag = element_text(size = 35, face = "bold"),
            plot.tag.position = c(0.02, 0.98)) +
      labs(tag = letter)

  })
  
  save_plot_as_png(ggarrange(plotlist = SPM_maps, common.legend = TRUE),
                   'zone_maps_panel', width = 28, height = 16, path = where_to_save_the_figure)

}


# manuscript Figure 4: daily plume area + mean SPM concentration time
# series (dynamic threshold, merged sensor), with an AR(1)/HAC-weighted
# trend line, one panel per zone. Renamed from Figures_6_7 and split from
# the threshold-comparison plot below (2026-08-01) -- that single function
# used to produce two different manuscript figures (4 and S1) under one
# name matching neither cleanly. Both still read the same
# FIGURE_4/DATA/ts_data.csv (prepped once by figure.py's
# Figure_4_S1_timeseries(), shared between this function and
# Figure_S1_thresholds() below).
# where_to_save_the_figure <- 'figures'
Figure_4_timeseries <- function(where_to_save_the_figure){
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_4")

  SPM_map_data <- main_folder %>% file.path('DATA', 'ts_data.csv') %>% read_csv()
  SPM_map_data$Dynamic_threshold <- ifelse(SPM_map_data$Dynamic_threshold, 'Dynamic threshold', 'Fixed threshold')

  # Plume-area trend line uses the same AR(1)/HAC-weighted fit as Table 5's
  # "Surface area" row -- see func/compute_area_trend.R.
  area_trend <- read_csv("output/STATS/area_trend_summary.csv", show_col_types = FALSE)

  # Panel order follows ZONE_ORDER (north to south), matching every other
  # multi-zone figure/table in the project -- dlply()'s default grouping
  # would otherwise sort panels alphabetically.
  SPM_map_data$Zone <- factor(SPM_map_data$Zone, levels = ZONE_ORDER)

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

    the_ts_plot_wo_modis <- ggplot() +
      geom_point(data = df_merged,
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") +
      geom_path(data = df_merged,
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") +
      geom_abline(data = zone_trend, aes(intercept = intercept, slope = slope),
                  colour = "black", linewidth = 1.2, linetype = "dashed") +
      geom_text(data = zone_trend, aes(x = -Inf, y = Inf, label = trend_label), inherit.aes = FALSE,
               hjust = -0.03, vjust = 1.5, size = 6, colour = "black") +

      scale_x_date(name = "",
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist(),
                   expand = c(0.01,0.01)) +

      coord_cartesian(ylim = c(0, max(df_zone$area_of_the_plume_mask_in_km2, na.rm = TRUE))) +
      # No per-panel y-axis title -- a single shared label is added once via
      # annotate_figure() on the assembled composite below, rather than
      # repeating "Plume area (km²)" on all four stacked panels.
      scale_y_continuous(name = NULL) +
      labs(x = "", title = zone_title(df_zone$Zone[1])) +
      ggplot_theme() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.subtitle = element_text(hjust = 0.5),
            plot.margin = margin(t = 20, r = 10, b = 5, l = 5),
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
  
  save_plot_as_png(annotate_figure(
                     ggarrange(plotlist = SPM_map_ts %>% llply(function(x) {x$wo_modis}), common.legend = FALSE, ncol = 1, nrow = 4, align = "v"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30)),
                   'Figure_4', width = 20, height = 16, path = main_folder)

  save_plot_as_png(annotate_figure(
                     ggarrange(plotlist = SPM_map_ts %>% llply(function(x) {x$w_modis}), common.legend = FALSE, ncol = 1, nrow = 4, align = "v"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30)),
                   'Figure_7_merged_modis', width = 20, height = 16, path = main_folder)
}

# manuscript Figure S1: plume-area time series, fixed vs. dynamic threshold
# comparison, one panel per zone. Reads the same ts_data.csv Figure_4_timeseries()
# does (both share one Python-side data prep), but writes to its own
# FIGURE_S1/ folder rather than FIGURE_4/ -- every manuscript figure has its
# own FIGURE_* folder (2026-08-01).
# where_to_save_the_figure <- 'figures'
Figure_S1_thresholds <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_4")
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_S1")

  SPM_map_data <- data_dir %>% file.path('DATA', 'ts_data.csv') %>% read_csv()
  SPM_map_data$Dynamic_threshold <- ifelse(SPM_map_data$Dynamic_threshold, 'Dynamic threshold', 'Fixed threshold')

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
      labs(y = "Plume area (km²)", x = "", title = zone_title(df_zone$Zone[1])) +
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
                   'Figure_S1', width = 20, height = 16, path = main_folder)

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
# zone maps, now folded into the Figure 3 composite rather than a standalone
# figure): manuscript.tex's own "Figure 5" caption had no generator at all
# until now (see [[project_manuscript_pipeline_gaps]] item 4).
Figure_5_driver_comparison <- function(where_to_save_the_figure, max_lag_daily = 14){

  # FIGURE_5 (not FIGURE_5_DRIVER): Figure_5() above no longer writes here --
  # its regional-zone-maps panel now feeds the Figure 3 composite directly
  # (see Figure_3() / figure.py::Figure_5()), freeing this folder for the
  # actual manuscript Figure 5.
  main_folder_of_Figure_5_driver <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_5")
  if (!dir.exists(main_folder_of_Figure_5_driver)) dir.create(main_folder_of_Figure_5_driver, recursive = TRUE)

  # ggplot_theme()'s font sizes are tuned for the single-column, 4-row
  # figures elsewhere (e.g. Figures_6_7); an 8-panel 2x4 grid needs smaller
  # text and explicit margins so axis titles don't clip against the plot
  # edge or the panel above. Same theme for every zone, so built once here
  # rather than inside the per-zone loop below.
  panel_theme <- theme(plot.title = element_text(hjust = 0.5, size = 20),
                       axis.title = element_text(size = 15, colour = "black"),
                       axis.text = element_text(size = 12, colour = "black"),
                       plot.margin = margin(t = 8, r = 12, b = 5, l = 5),
                       panel.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_rect(linetype = "solid", fill = NA))

  panels <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    df <- combine_plume_driver("flow", meta)

    cor_df <- driver_plume_correlation(df, max_lag_daily = max_lag_daily) %>%
      dplyr::filter(timestep == "daily")
    peak <- cor_df %>% dplyr::slice_max(cor, n = 1)

    # Panel title is the zone (the flow series compared here is already
    # zone-summed across every contributing river -- see load_river_flow()),
    # not the representative river/mouth name.
    zone_label <- zone_title(meta$zone)

    scatter_plot <- ggplot(df, aes(x = value, y = plume_area)) +
      geom_point(alpha = 0.3, colour = "grey30", size = 0.8) +
      geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1.2) +
      labs(x = "River flow (m³ s⁻¹)", y = "Plume area (km²)", title = zone_label) +
      panel_theme

    lag_plot <- ggplot(cor_df, aes(x = lag, y = cor)) +
      geom_line(colour = "grey30") +
      geom_point(colour = "grey30") +
      geom_point(data = peak, colour = "firebrick", size = 3) +
      geom_text(data = peak, aes(label = paste0(lag, "d")), vjust = -1.2, colour = "firebrick", size = 4) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(x = "Lag, plume after flow (days)", y = "Correlation (r)", title = zone_label) +
      panel_theme

    list(scatter = scatter_plot, lag = lag_plot)
  })

  plotlist <- panels %>% purrr::map(function(p) list(p$scatter, p$lag)) %>% purrr::flatten()

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(panels), align = "v",
                                 labels = panel_labels, font.label = list(size = 18, face = "bold"),
                                 hjust = -0.3, vjust = 1.3)

  save_plot_as_png(full_plot, "Figure_5", width = 16, height = 18, path = main_folder_of_Figure_5_driver)
}


# manuscript Figure 7: wind and wave direction/magnitude roses, one row per
# zone, coloured by the flow-controlled plume-area response
# (multi.R::plot_driver_rose()). Replaces the original Figure 7/8 concept
# (X11-decomposed wind/wave magnitude time series) -- Robert's call, since
# wind and wave direction matter as much or more than magnitude, which a
# magnitude-only time series can't show.
Figure_7_driver_rose <- function(where_to_save_the_figure, n_sectors = 8){

  main_folder_of_Figure_7 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_7")
  if (!dir.exists(main_folder_of_Figure_7)) dir.create(main_folder_of_Figure_7, recursive = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    # One flow join/STL per zone, shared across the wind and wave roses below
    # (plot_driver_rose() would otherwise recompute the identical df_flow twice).
    df_flow <- combine_plume_driver("flow", meta) %>% dplyr::select(date, plume_area, flow = value)
    list(wind = plot_driver_rose("wind", meta, n_sectors, df_flow), wave = plot_driver_rose("wave", meta, n_sectors, df_flow))
  }) %>% purrr::map(function(p) list(p$wind, p$wave)) %>% purrr::flatten()

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = nrow(zone_meta), align = "v",
                                 labels = panel_labels, font.label = list(size = 18, face = "bold"),
                                 hjust = -0.3, vjust = 1.3)

  save_plot_as_png(full_plot, "Figure_7", width = 16, height = 20, path = main_folder_of_Figure_7)
}


# manuscript Figure 8: flow-controlled plume-area residual vs. wave height,
# coloured by on/off-shore wind category, one panel per zone
# (multi.R::plot_category_scatter()). Replaces the original Figure 8
# concept (X11-decomposed wave-height magnitude time series).
Figure_8_driver_category <- function(where_to_save_the_figure){

  main_folder_of_Figure_8 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_8")
  if (!dir.exists(main_folder_of_Figure_8)) dir.create(main_folder_of_Figure_8, recursive = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    plot_category_scatter(meta)
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

  # Sourced here rather than at file scope (unlike multi.R above): this pulls
  # in several heavyweight modelling packages (mgcv, gratia, ranger, iml)
  # that only this one figure function needs -- not worth loading for every
  # other figure in this file.
  source("func/driver_interactions.R")

  main_folder_of_Figure_9 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_9")
  if (!dir.exists(main_folder_of_Figure_9)) dir.create(main_folder_of_Figure_9, recursive = TRUE)

  driver_labels <- c(flow = "River flow (m³ s⁻¹)", wind_spd = "Wind speed (m s⁻¹)",
                     wave_height = "Wave height (m)", current = "Current speed (m s⁻¹)")

  plotlist <- purrr::map(zones, function(zone_name){
    df <- readr::read_csv(file.path(stats_dir, paste0("daily_driver_matrix_", zone_name, ".csv")), show_col_types = FALSE)
    gam_model <- fit_gam(df)
    drivers_to_show <- intersect(names(driver_labels), available_drivers(df))

    zone_plots <- purrr::map(drivers_to_show, function(d){
      curve <- gam_partial_effect(gam_model, d, df)
      ggplot(curve, aes(x = x, y = fit)) +
        geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), fill = "grey80", alpha = 0.5) +
        geom_line(colour = "black", linewidth = 1) +
        labs(x = driver_labels[[d]], y = "Partial effect on plume area (km²)") +
        theme(panel.border = element_rect(fill = NA, colour = "black"))
    })
    # Only the row's first panel gets the zone name as a title.
    zone_plots[[1]] <- zone_plots[[1]] + labs(title = zone_title(zone_name))
    zone_plots
  }) %>% purrr::flatten()

  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = length(driver_labels), nrow = length(zones), align = "v")

  save_plot_as_png(full_plot, "Figure_9", width = 20, height = 16, path = main_folder_of_Figure_9)
}


# Computes the per-zone X11 Interannual/Seasonal/Residual plots for one
# threshold (dynamic or static) from the weekly plume/river time series
# prepped by figure.py's Figure_X11_weekly_results(). Pure computation, no
# saving -- shared by the four Figure_*_x11_*() functions below so the
# per-zone join/plot logic isn't duplicated per manuscript figure.
# Replaces the old Figures_8_9_10() (2026-08-01, deprecated): that function
# also computed a "Raw" signal plot per zone that was never saved anywhere
# (dead computation), dropped here.
compute_x11_zone_plots <- function(data_dir){
  plume_data <- data_dir %>% file.path('DATA', 'ts_plume_data.csv') %>% read_csv()
  river_data <- data_dir %>% file.path('DATA', 'ts_river_data.csv') %>% read_csv()

  regions <- unique(plume_data$Zone) %>% order_zones()

  regions %>% llply(function(region) {

    plume_data_region <- plume_data %>% filter(Zone == region)
    river_data_region <- river_data %>% filter(Zone == region)

    X11_ts <- plume_data_region %>% inner_join(river_data_region, by = "dates", suffix = c("_plume_area", "_river_flow"))

    zone_label <- zone_title(region)
    list("Interannual" = plot_x11_river_and_plume(X11_ts, type_of_signal = 'Interannual') + labs(title = zone_label),
        "Seasonal" = plot_x11_river_and_plume(X11_ts, type_of_signal = 'Seasonal') + labs(title = zone_label),
        "Residual" = plot_x11_river_and_plume(X11_ts, type_of_signal = 'Residual') + labs(title = zone_label))

  })
}

# Stacks one X11 component (Interannual/Seasonal/Residual) across all 4
# zones into a single column -- the layout shared by every figure below.
stack_x11_component <- function(zone_plots, component){
  ggarrange(plotlist = zone_plots %>% llply(function(x) x[[component]]), ncol = 1, nrow = 4, align = "v")
}

# manuscript Figure 6: X11 interannual (long-term) signal of plume area vs.
# river flow, dynamic threshold (main results), all four zones.
Figure_6_x11_interannual <- function(where_to_save_the_figure){
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_6")
  zone_plots <- compute_x11_zone_plots(main_folder)
  save_plot_as_png(stack_x11_component(zone_plots, "Interannual"), "Figure_6", width = 20, height = 16, path = main_folder)
}

# Renders each of two X11-component stacks (e.g. Seasonal + Residual) to its
# own PNG first, then composites them with magick::image_append(stack=TRUE)
# -- nesting a second ggarrange() around two already-4-row ggarrange() stacks
# instead corrupts the rotated secondary-axis text (checked visually: the
# y-axis labels overlap into an unreadable smear). Compositing at the image
# level, like the pre-2026-08-01 Figure_8/Figure_10 assembly did, avoids this.
save_x11_component_composite <- function(zone_plots, components, name, path){
  panel_files <- purrr::map_chr(components, function(component){
    save_plot_as_png(stack_x11_component(zone_plots, component), paste0(name, "_", tolower(component)),
                     width = 20, height = 16, path = path)
    file.path(path, paste0(name, "_", tolower(component), ".png"))
  })
  composite <- magick::image_append(magick::image_read(panel_files), stack = TRUE)
  magick::image_write(composite, file.path(path, paste0(name, ".png")))
}

# manuscript Figure S2: X11 seasonal + residual components, dynamic
# threshold, all four zones -- the two components not shown in main Figure 6.
# Shares Figure 6's DATA/ prep (data_dir points at FIGURE_6/, not FIGURE_S2/).
# Renamed 2026-07-31 from Figure_S3_x11_components() when the daily-vs-weekly
# Figure S2 was removed from the manuscript and every later supplementary
# figure renumbered down by one (S3->S2, S4->S3, S5->S4, S6->S5, S7->S6).
Figure_S2_x11_components <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_6")
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_S2")
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_zone_plots(data_dir)
  save_x11_component_composite(zone_plots, c("Seasonal", "Residual"), "Figure_S2", main_folder)
}

# manuscript Figure S5: X11 interannual signal, static threshold
# (supplementary robustness check) -- mirrors Figure 6. Renamed 2026-07-31
# from Figure_S6_x11_interannual_static(), see Figure_S2_x11_components() note.
Figure_S5_x11_interannual_static <- function(where_to_save_the_figure){
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_S5")
  zone_plots <- compute_x11_zone_plots(main_folder)
  save_plot_as_png(stack_x11_component(zone_plots, "Interannual"), "Figure_S5", width = 20, height = 16, path = main_folder)
}

# manuscript Figure S6: X11 seasonal + residual components, static threshold
# -- mirrors Figure S2. Shares Figure S5's DATA/ prep. Renamed 2026-07-31
# from Figure_S7_x11_components_static(), see Figure_S2_x11_components() note.
Figure_S6_x11_components_static <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_S5")
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_S6")
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_zone_plots(data_dir)
  save_x11_component_composite(zone_plots, c("Seasonal", "Residual"), "Figure_S6", main_folder)
}

# manuscript Figure S3: seasonal (JJA vs. NDJ) boxplots of plume metrics,
# dynamic vs. static threshold, all four zones. Migrated 2026-08-01 from
# manuscript/make_figures_tables.R::generate_figure_s4_seasonal_thresholds()
# into the real pipeline (called from code/5_figures.py) so it writes
# straight to figures/ARTICLE/FIGURE_S3/ instead of via the
# manuscript/figures/ copy step. Renamed 2026-07-31 from
# Figure_S4_seasonal_boxplots(), see Figure_S2_x11_components() note.
Figure_S3_seasonal_boxplots <- function(where_to_save_the_figure){
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_S3")

  metrics <- c(
    area_of_the_plume_mask_in_km2           = "Plume area (km²)",
    mean_SPM_in_the_plume_area              = "Mean SPM (g m⁻³)",
    lon_weighted_centroid_of_the_plume_area = "Centroid longitude (°E)",
    lat_weighted_centroid_of_the_plume_area = "Centroid latitude (°N)"
  )

  all_data <- bind_rows(lapply(c("dynamic", "static"), function(thresh){
    bind_rows(lapply(zones, function(zone){
      csv_path <- file.path("output", "panache", thresh, zone, "Results.csv")
      if (!file.exists(csv_path)) {
        message("Figure S3: skipping (not found): ", csv_path)
        return(NULL)
      }
      read_csv(csv_path, show_col_types = FALSE) %>%
        mutate(date = as.Date(date), month = as.integer(format(date, "%m")), zone = zone, threshold = thresh)
    }))
  }))

  if (nrow(all_data) == 0) {
    message("Figure S3: no Results.csv files found -- skipping.")
    return(invisible(FALSE))
  }

  seasonal_data <- all_data %>%
    filter(month %in% c(6L, 7L, 8L, 11L, 12L, 1L)) %>%
    mutate(season = case_when(month %in% c(6L, 7L, 8L) ~ "Summer (JJA)",
                              month %in% c(11L, 12L, 1L) ~ "Winter (NDJ)"),
           zone = factor(zone, levels = zones, labels = zone_title(zones)),
           threshold = factor(threshold, levels = c("dynamic", "static"), labels = c("Dynamic", "Static"))) %>%
    select(date, zone, threshold, season, all_of(names(metrics))) %>%
    pivot_longer(cols = names(metrics), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = names(metrics), labels = unname(metrics)))

  p <- ggplot(seasonal_data, aes(x = season, y = value, fill = threshold)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, position = position_dodge(0.75)) +
    facet_grid(metric ~ zone, scales = "free_y") +
    scale_fill_manual(values = c("Dynamic" = "#2166ac", "Static" = "#d6604d"), name = "Threshold") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 8), legend.position = "bottom",
         axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())

  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  ggsave(file.path(main_folder, "Figure_S3.png"), p, width = 14, height = 12, dpi = 150)
  message("Wrote Figure_S3.png (", nrow(seasonal_data), " observations across available zones and thresholds)")
  invisible(TRUE)
}

