# Main --------------------------------------------------------------------

# where_to_save_the_figure <- "~/RiOMar/figures/"

Figure_1 <- function(where_to_save_the_figure) {

  main_folder_of_Figure_1 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_1")

  SPM_map <- file.path(main_folder_of_Figure_1, "DATA", "SPM_map.csv") |> read_csv()
  insitu_stations <- file.path(main_folder_of_Figure_1, "DATA", "Stations_position.csv") |> read_csv()

  RIOMAR_limits <- zones_bbox |> dplyr::rename(Zone = zone)

  basic_map <- create_the_basic_map(map_df = SPM_map, var_name = 'SPM', in_situ_fixed_station = insitu_stations,
                                    log_scale = FALSE, legend_limits = c(0.1, 10), high_res_coast = TRUE)

  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'),
                                      longitude = c(0,0),
                                      latitude = c(0,0))

  national_map <- basic_map +
    geom_point(data = insitu_stations |> filter(SOURCE == 'REPHY'),
               aes(x = LONGITUDE, y = LATITUDE),
               fill = "red", color = "black", size = 4, shape = 24, stroke = 1) +
    geom_point(data = insitu_stations |> filter(SOURCE == 'SOMLIT'),
               aes(x = LONGITUDE, y = LATITUDE),
               fill = "red", color = "black", size = 10, shape = 21, stroke = 2) +
    geom_rect(data = RIOMAR_limits, aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max),
              fill = "transparent", color = "red", linetype = "dashed", size = 2) +

    # Open-water sea/ocean labels, positioned clear of the zone bboxes and
    # zoomed insets above (Atlantic: open water west of Southern Brittany/
    # Gironde shelf; Mediterranean: open water south of the Rhône shelf coast).
    annotate("text", x = -4.6, y = 45.5, label = "Atlantic\nOcean",
             colour = "white", fontface = "italic", size = 15, hjust = 0.5) +
    annotate("text", x = 6.6, y = 42.6, label = "Mediterranean\nSea",
             colour = "white", fontface = "italic", size = 15, hjust = 0.5) +

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
  # `primary` marks the river carrying the bulk of each zone's discharge
  # (Seine: only river in its zone; Loire vs. Vilaine and Gironde vs.
  # Charente/Sevre: the named river the manuscript treats as the zone's main
  # discharge series throughout.
  zone_river_mouths <- tibble::tribble(
    ~zone,               ~river,         ~lat,   ~lon,    ~primary,
    "BAY_OF_SEINE",       "Seine",        49.43,  0.145,  TRUE,
    "SOUTHERN_BRITTANY",  "Loire",        47.24, -2.2,    TRUE,
    "SOUTHERN_BRITTANY",  "Vilaine",      47.48, -2.55,   FALSE,
    "BAY_OF_BISCAY",      "Gironde",      45.61, -1.14,   TRUE,
    "BAY_OF_BISCAY",      "Charente",     45.98, -1.15,   FALSE,
    "BAY_OF_BISCAY",      "Sevre",        46.26, -1.2,    FALSE,
    "GULF_OF_LION",       "Grand\nRhône",  43.32,  4.85,  TRUE,
    "GULF_OF_LION",       "Petit\nRhône",  43.45,  4.39,  FALSE
  )

  # One geom_label() call per river (not vectorised per zone) so each
  # label's position can be tuned individually against the actual coastline
  # geometry, rather than sharing a single per-zone offset.
  mouth <- function(river_name) dplyr::filter(zone_river_mouths, river == river_name)
  river_label_style <- function(river_name, primary, colour, ...){
    geom_label(data = mouth(river_name), aes(x = lon, y = lat, label = river),
              fontface = if(primary) "bold" else "plain", size = if(primary) 11 else 9,
              colour = colour, fill = if(colour == "white") "black" else "white", alpha = 0.4, ...)
  }

  river_labels_by_zone <- list(
    BAY_OF_SEINE = list(
      river_label_style("Seine", primary = TRUE, colour = "white", hjust = 0, nudge_x = 0.0, nudge_y = 0.12)
    ),
    SOUTHERN_BRITTANY = list(
      river_label_style("Loire", primary = TRUE, colour = "white", hjust = 0, nudge_x = 0.08, nudge_y = 0.00),
      river_label_style("Vilaine", primary = FALSE, colour = "white", hjust = 0, nudge_x = 0.08, nudge_y = 0.05)
    ),
    BAY_OF_BISCAY = list(
      river_label_style("Gironde", primary = TRUE, colour = "black", hjust = 1, nudge_x = -0.10, nudge_y = -0.02),
      river_label_style("Charente", primary = FALSE, colour = "black", hjust = 1, nudge_x = -0.10, nudge_y = 0.11),
      river_label_style("Sevre", primary = FALSE, colour = "black", hjust = 1, nudge_x = -0.10, nudge_y = 0.05)
    ),
    GULF_OF_LION = list(
      river_label_style("Grand\nRhône", primary = TRUE, colour = "black", hjust = 0, nudge_x = -0.08, nudge_y = -0.22),
      river_label_style("Petit\nRhône", primary = FALSE, colour = "black", hjust = 0, nudge_x = -0.14, nudge_y = -0.20)
    )
  )

  build_zone_inset <- function(zone_name) {
    zone_SPM <- file.path(main_folder_of_Figure_1, "DATA", paste0(zone_name, ".csv")) |> read_csv()
    mouths <- zone_river_mouths |> dplyr::filter(zone == zone_name)

    create_the_basic_map(zone_SPM, 'SPM', log_scale = FALSE, legend_limits = c(0.1, 10), high_res_coast = TRUE) +
      geom_point(data = mouths, aes(x = lon, y = lat), shape = 4, colour = "red", size = 4, stroke = 2) +
      river_labels_by_zone[[zone_name]] +
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
  # spacing rather than an even grid: Seine upper right, Rhone
  # underneath it, Loire where Seine used to sit, Gironde roughly in place --
  # shifted further right than that as a whole group so no inset covers any
  # coloured SPM pixels (checked against zones_bbox's true-box fractions).
  inset_layout <- tibble::tribble(
    ~Zone,               ~x,    ~y,    ~w,    ~h,
    "BAY_OF_SEINE",       0.68,  0.64,  0.23,  0.20,
    "SOUTHERN_BRITTANY",  0.41,  0.50,  0.25,  0.21,
    "BAY_OF_BISCAY",      0.45,  0.26,  0.24,  0.21,
    "GULF_OF_LION",       0.72,  0.34,  0.23,  0.20
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

# Builds the 4-zone regional SPM map grid
zone_maps_panels <- function(data_folder, include_station_points) {

  SPM_map_data <- file.path(data_folder, "DATA") |>
    list.files(pattern = "*.csv", full.names = TRUE) |>
    plyr::llply(read_csv) |>
    keep(~ 'analysed_spim' %in% names(.))

  insitu_stations <- file.path( data_folder, "DATA", "Stations_position.csv" ) |> read_csv()

  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'), longitude = c(0,0), latitude = c(0,0))

  SPM_maps <- SPM_map_data |>
    plyr::llply(function(x) {
      insitu_stations_of_the_map <- insitu_stations |>
        filter((LATITUDE |> between(min(x$lat), max(x$lat))) &
                 (LONGITUDE |> between(min(x$lon), max(x$lon))))

      the_map <- create_the_basic_map(x, 'SPM', legend_limits = c(4,10))

      if (include_station_points) {

        the_map <- the_map +

          geom_point(data = insitu_stations_of_the_map |> filter(SOURCE == 'REPHY'),
                     aes(x = LONGITUDE, y = LATITUDE),
                     fill = "red", color = "black", size = 6, shape = 24, stroke = 1) +
          geom_point(data = insitu_stations_of_the_map |> filter(SOURCE == 'SOMLIT'),
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
    })

  ggarrange(plotlist = SPM_maps, common.legend = TRUE)
}

# Standalone regional-zone-maps figure
regional_zone_maps <- function(where_to_save_the_figure, include_station_points) {

  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURE_1")

  SPM_maps <- zone_maps_panels(main_folder, include_station_points)

  save_plot_as_png(SPM_maps, paste("regional_zone_maps", ifelse(include_station_points, "with_stations", "wo_stations"), sep = "_"),
                   width = 28, height = 16, path = main_folder)

}


# Satellite-vs-in-situ validation scatterplots, panel (a) SPM and panel (b)
# Turbidity. Manuscript slot "validation_scatterplot_panel" -- see
# metadata/figure_table_registry.csv for its current figure number and
# output folder (get_registry_row(), func/util.R).
plot_validation_scatterplot_panel <- function(spm_scatterplot_path, turb_scatterplot_path, where_to_save_the_figure) {

  output_subdir <- get_registry_row("validation_scatterplot_panel")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)

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
  combined <- image_append(c(panel_a, panel_b), stack = FALSE)

  image_write(combined, file.path(main_folder, registry_filename(output_subdir)))

}


# Renders one methodology panel (A-D) for the plume_methodology_panel figure
# (metadata/figure_table_registry.csv)
# where_to_save_the_figure <- '/figures/ARTICLE/' + that slot's output_subdir
# name_of_the_plot <- "C"
plot_methodology_worked_example_panel <- function(where_to_save_the_figure, name_of_the_plot) {
  
  SPM_map_data <- read_csv(file.path(where_to_save_the_figure, "DATA", paste(name_of_the_plot, ".csv", sep = "")))
  
  # legend_limits matches plot_methodology_zone_maps_panel()'s panels E-H
  if (name_of_the_plot %in% c("A", "B")) {
    the_map <- create_the_basic_map(SPM_map_data, 'SPM', legend_limits = c(0.1,10), high_res_coast = TRUE)
  } else {
    the_map <- create_the_basic_map(SPM_map_data, 'plume', legend_limits = c(0.1,10), high_res_coast = TRUE)
  }

  # Convert name_of_plot to pretty labels
  tag_label <- paste0(tolower(name_of_the_plot),")")
  
  # No per-panel colour bar or lon/lat axis text on the top row
  the_map <- the_map +
    labs(tag = tag_label) +
    theme(legend.position = "none",
          # NB: For some reason it is necessary to call x and y explicitly
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.tag = element_text(size = 60, face = "bold"),
          plot.tag.position = c(0.02, 0.98),
          plot.margin = margin(t = 40, r = 10, b = 10, l = 10))

  if (name_of_the_plot == "B") {
    points_used_for_finding_SPM_threshold <- read_csv(file.path(where_to_save_the_figure,
                                                                "DATA", "B_points_used_for_finding_SPM_threshold.csv"))
    all_points_used_for_finding_SPM_threshold <- read_csv(file.path(where_to_save_the_figure,
                                                                    "DATA", "B_all_points_used_for_finding_SPM_threshold.csv"))
    the_map <- the_map +
      geom_point(data = all_points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "grey50", size = 3) +
      geom_point(data = points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "red", size = 3)
  }

  if (name_of_the_plot == "A") {
    # Marks the Grand Rhone river mouth (panel a's "daily map of SPM for the
    # study region"), per Robert's request -- river_mouths comes from
    # func/util.R (panache_zone_metadata), already loaded.
    grand_rhone_mouth <- dplyr::filter(river_mouths, mouth_name == "Grand Rhone")
    the_map <- the_map +
      geom_point(data = grand_rhone_mouth, aes(x = mouth_lon, y = mouth_lat),
                colour = "red", shape = 4, size = 7, stroke = 3)

    # 15-pixel-radius circle showing panache's near-mouth sampling window
    # (near_mouth_radius_pixels, panache/src/panache/config.py -- the
    # circular window used to derive the SPM_threshold's minimal/maximal
    # bounds; same 15-pixel default already used by figure.py's panel-B
    # data export). Radius converted from pixels to degrees using this
    # scene's own grid spacing (SPM_map_data's lon/lat are on a regular
    # grid), not a fixed km distance -- an ellipse rather than a true
    # circle wherever the grid's lon/lat pixel spacing differ, matching
    # what a 15-pixel *index* radius (panache samples in row/column space)
    # actually looks like in lon/lat map space.
    lon_grid <- sort(unique(SPM_map_data$lon))
    lat_grid <- sort(unique(SPM_map_data$lat))
    near_mouth_radius_pixels <- 15
    radius_lon <- near_mouth_radius_pixels * mean(diff(lon_grid))
    radius_lat <- near_mouth_radius_pixels * mean(diff(lat_grid))
    theta <- seq(0, 2 * pi, length.out = 100)
    near_mouth_circle <- data.frame(
      lon = grand_rhone_mouth$mouth_lon + radius_lon * cos(theta),
      lat = grand_rhone_mouth$mouth_lat + radius_lat * sin(theta))
    the_map <- the_map +
      geom_path(data = near_mouth_circle, aes(x = lon, y = lat),
               colour = "red", linewidth = 1.8, inherit.aes = FALSE)
  }

  # 20 m bathymetric exclusion boundary:
  # `shallow` is all-False outside the Gulf of Lion (the only zone using
  # this general exclusion), so geom_contour() draws nothing there.
  if (any(SPM_map_data$shallow)) {
    the_map <- the_map +
      geom_contour(data = SPM_map_data, aes(x = lon, y = lat, z = as.numeric(shallow)),
                  breaks = 0.5, colour = "white", linewidth = 1, linetype = "dashed")
  }

  save_plot_as_png(the_map, name_of_the_plot, width = 12, height = 8, path = where_to_save_the_figure)

}


# New methodology panel, inserted between the
# A-D worked-example row and the per-zone f)-i) grid: shows the actual SPM
# value found at every point tested along panel B's transects, plotted
# against distance from the river mouth, coloured by whether the point
# survived the gradient-cutoff filter (an edge candidate) or was rejected,
# together with the near-mouth minimal/maximal_threshold bounds and the
# final derived SPM_threshold. Makes explicit what panel B's grey/red points
# only show spatially: how the gradient cutoff and quantile bounds actually
# combine to pick the scene-specific plume-edge threshold.
plot_methodology_transect_panel <- function(where_to_save_the_figure) {

  transect_values <- read_csv(file.path(where_to_save_the_figure, "DATA", "B_transect_values.csv"))
  threshold_values <- read_csv(file.path(where_to_save_the_figure, "DATA", "B_threshold_values.csv"))

  ref_lines <- tibble::tribble(
    ~label,               ~value,
    "maximal_threshold",  threshold_values$maximal_threshold[1],
    "minimal_threshold",  threshold_values$minimal_threshold[1],
    "SPM_threshold",      threshold_values$SPM_threshold[1]
  )
  ref_line_styles <- c(maximal_threshold = "dotted", minimal_threshold = "dashed", SPM_threshold = "solid")
  ref_line_labels <- c(maximal_threshold = "maximal_threshold", minimal_threshold = "minimal_threshold",
                       SPM_threshold = "SPM_threshold (final)")

  the_plot <- ggplot(transect_values, aes(x = distance_km, y = analysed_spim)) +
    geom_point(aes(colour = kept), size = 2, alpha = 0.7) +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "kept (edge candidate)", `FALSE` = "rejected"), name = NULL) +
    geom_hline(data = ref_lines, aes(yintercept = value, linetype = label), colour = "black", linewidth = 1) +
    scale_linetype_manual(values = ref_line_styles, labels = ref_line_labels, name = NULL) +
    labs(x = "Distance from river mouth (km)", y = expression(SPM~(g~m^{-3})), tag = "e)") +
    ggplot_theme() +
    # Legend moved inside the plot area 2026-08-11 (was legend.position =
    # "right", eating into the panel's plotting width) and its text sized up
    # for legibility, per metadata/TODO.md. Anchored top-right: rendered
    # against the real transect data, SPM decays sharply with distance from
    # the mouth, so nothing is ever plotted in the high-distance/high-SPM
    # corner -- confirmed empty, not assumed.
    theme(plot.tag = element_text(size = 35, face = "bold"), plot.tag.position = c(0.01, 0.98),
          legend.position = "inside", legend.position.inside = c(0.88, 0.8),
          legend.background = element_rect(fill = "white", colour = "grey70"),
          legend.text = element_text(size = 24), legend.title = element_text(size = 24),
          text = element_text(size = 20, colour = "black"),
          axis.text = element_text(size = 18, colour = "black"))

  save_plot_as_png(the_plot, "transect_panel", width = 24, height = 7.2, path = where_to_save_the_figure)

}


# Renders the per-zone plume-maps panel feeding the plume_methodology_panel
# figure (metadata/figure_table_registry.csv)
plot_methodology_zone_maps_panel <- function(where_to_save_the_figure) {

  # Read only the four per-zone SPM-map CSVs figure.py's plot_methodology_zone_maps_panel()
  # writes here (Zone.csv, via zone_meta$zone for canonical zone naming/order)
  # -- a plain "*.csv" glob on this shared DATA/ folder also picks up
  # Figure_2_methodology_panels()' A-E.csv (which lack a `plume` column entirely) and its
  # *_threshold.csv debug files, crashing create_the_basic_map()'s
  # which(map_df$plume) on the first file missing that column.
  SPM_map_data <- where_to_save_the_figure |>
    file.path('DATA', paste0(zone_meta$zone, ".csv")) |>
    plyr::llply(read_csv)

  # Continues the lettering from plot_methodology_worked_example_panel()'s methodology row (A-D)
  # and plot_methodology_transect_panel()'s e) zone_meta$zone is already arranged by
  # ZONE_ORDER (north to south) Seine, Southern Brittany, Bay of Biscay, Gulf of Lion),
  # matching the order these zones are listed in the figure's caption.
  panel_letters <- c("f)", "g)", "h)", "i)")
  SPM_maps <- purrr::map2(SPM_map_data, panel_letters, function(SPM_map, letter) {

    the_map <- create_the_basic_map(SPM_map, 'plume', legend_limits = c(0.1,10), high_res_coast = TRUE) +
      guides(fill = guide_colorbar(barwidth = 60, barheight = 2, title.position = "top")) +
      theme(legend.position = "top",
            legend.title = element_text(angle = 0, hjust = 0.5),
            axis.text = element_text(size=25, colour = "black"),
            plot.tag = element_text(size = 35, face = "bold"),
            plot.tag.position = c(0.02, 0.98)) +
      labs(tag = letter)

    # Bathymetric exclusion boundary matching panache's actual per-zone
    # resuspension-removal depth (`shallow` computed in figure.py from
    # maximal_bathymetric_for_zone_with_resuspension: Seine 20 m,
    # Gironde/Charente/Sevre 10 m, Loire/Vilaine 10 m, Grand/Petit Rhone
    # 20 m -- see figure.py for detail).
    if (any(SPM_map$shallow)) {
      the_map <- the_map +
        geom_contour(data = SPM_map, aes(x = lon, y = lat, z = as.numeric(shallow)),
                    breaks = 0.5, colour = "white", linewidth = 1, linetype = "dashed")
    }

    the_map

  })
  
  save_plot_as_png(ggarrange(plotlist = SPM_maps, common.legend = TRUE),
                   'zone_maps_panel', width = 28, height = 16, path = where_to_save_the_figure)

}


