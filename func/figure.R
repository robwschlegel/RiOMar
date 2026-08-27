# func/figure.R
# This is called by func/figure.py to provide graphical outputs


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(scales)
library(maps)
library(ggpubr)
library(cowplot)
library(magick)

# For zones_bbox
source("func/util.R")

# For zone_meta/get_zone_meta/combine_plume_driver/plot_driver_rose/etc
source("func/multi.R")

# For the tide QC diagnostics (tide_qc_all(), save_tide_qc_examples(), etc.)
# -- .load_tide_raw() and the QC logic load_tide_gauge() calls both live in
# func/util.R itself, not here; nothing in the pipeline path requires this
# file, it is sourced only for its standalone diagnostic/reporting functions.
source("func/tide.R")


# Utils -------------------------------------------------------------------

# High-resolution France coastline (GADM level-0 boundary, ~216,000 vertices)
# for the study-zone-map figure's zoomed insets.
# Loaded lazily (library(sf) is NOT called at file scope

# Step size for four_ticks_from_zero() below, factored out so a dual-axis
# figure (e.g. plot_plume_area_timeseries()) can also compute the step needed to
# reach an arbitrary target -- not just a series' own max(x) -- letting two
# independently-scaled axes share the same top tick instead of each rounding
# up separately and leaving whichever axis rounds up less stranded well
# below the panel top (fixed 2026-08-11, the plume-area-timeseries figure's
# "dead space" tweak).
# Unlike base pretty()/scales::breaks_pretty(), which restrict the step to a
# 1/2/5 x 10^k multiple, a step of 2.5 x 10^k is allowed too -- needed so 3
# equal steps can land close to the target (e.g. target ~7500 -> steps of
# 2500, giving 0, 2500, 5000, 7500) rather than overshooting it by a lot.
# Picks the *smallest* nice step whose 3 steps still reach the target --
# picking whichever nice step is merely closest to target/3 (tried first)
# can round down and clip real data above the top tick, found via the Gulf
# of Lion's actual max area (8987 km^2) rounding down to a 7500 top tick.
nice_step_for_target <- function(target){
  raw_step <- target / 3
  magnitude <- 10 ^ floor(log10(raw_step))
  candidate_steps <- c(1, 2, 2.5, 5, 10) * magnitude
  min(candidate_steps[candidate_steps * 3 >= target])
}

# Four evenly spaced y-axis ticks from 0 up to a "nice" round number at or
# above max(x). See nice_step_for_target() above for the step-selection logic.
four_ticks_from_zero <- function(x){
  step <- nice_step_for_target(max(x, na.rm = TRUE))
  seq(0, step * 3, by = step)
}

# Natural Earth 10m "Land" (naturalearthdata.com, public domain), vendored
# into data/EUROPE_shapefile/ 2026-08-11: a single global land-polygon layer
# at the same 10m-class resolution as the GADM France file this replaced,
# covering all of Europe and the UK (and everywhere else) in one file. Was
# previously France-only (GADM ADM0), overplotted on the crude `maps`-
# package "world" coastline for every other country -- visibly mismatched
# resolution at the UK/Belgium/Spain edges of every map. No per-country
# subsetting needed: "Land" is one continuous landmass polygon, matching
# how it's used below (a single filled black layer, not coloured by
# country), so create_the_basic_map() no longer needs the low-res
# map_data("world") layer at all when high_res_coast = TRUE.
.high_res_coastline_cache <- NULL
high_res_coastline <- function(){
  if(is.null(.high_res_coastline_cache)){
    library(sf)
    .high_res_coastline_cache <<- sf::st_read("data/EUROPE_shapefile/ne_10m_land.shp", quiet = TRUE) |>
      sf::st_coordinates() |>
      as.data.frame() |>
      dplyr::transmute(long = X, lat = Y, group = interaction(L1, L2, L3, drop = TRUE))
  }
  .high_res_coastline_cache
}

create_the_basic_map <- function(map_df, var_name,
                                 in_situ_fixed_station = NULL,
                                 cruise_stations = NULL,
                                 glider_stations = NULL,
                                 legend_limits = NULL,
                                 log_scale = TRUE,
                                 high_res_coast = FALSE) {
  
  if (str_detect(var_name, 'chl|CHL')) {
    title = "Chl-a"
    unit = "mg m³"
    if (legend_limits |> is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  if (str_detect(var_name, 'tsm|SPM|TSM|plume')) {
    title = "SPM"
    unit = "g m³"
    # legend_limits <- map_df$analysed_spim[which(map_df$plume)] |> quantile(probs = c(0.1, 0.9), na.rm = TRUE)
    if (legend_limits |> is.null()) {legend_limits <- c(1e-1, 5e0)} 
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
      # High-resolution coastline (Natural Earth 10m land, all of Europe +
      # UK and beyond -- see high_res_coastline() above) as the only layer.
      # Fixed 2026-08-11: previously this was France-only (GADM), layered
      # OVER the crude low-res map_data("world") coastline underneath, so
      # every neighbouring country's coast (UK, Belgium, Spain) stayed at
      # visibly cruder resolution than France's own -- both on the main
      # national panel and every zoomed inset. The new file covers the
      # whole extent any of these maps needs, so the low-res layer is no
      # longer drawn at all when high_res_coast = TRUE.
      list(
        geom_polygon(data = high_res_coastline(), aes(x = long, y = lat, group = group), color = 'grey60', fill = 'black')
      )
    } else {
      list(
        # First layer: worldwide map
        geom_polygon(data = map_data("world"), aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black'),
        # Second layer: Country map
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

plot_x11_river_and_plume <- function(X11_data, type_of_signal, show_axis_titles = TRUE) {

  unique_years <- X11_data$dates |> year() |> unique()

  X11_data_for_plot <- X11_data |>
    rename(river_flow = !!sym(paste(type_of_signal, "signal_river_flow", sep = "_")),
           plume_area = !!sym(paste(type_of_signal, "signal_plume_area", sep = "_"))) |>
    select(dates, river_flow, plume_area)

  if (type_of_signal %in% c("Seasonal", "Residual")) {
    X11_data_for_plot <- X11_data_for_plot |>
      mutate(river_flow = river_flow + mean(X11_data$Raw_signal_river_flow, na.rm = T),
             plume_area = plume_area + mean(X11_data$Raw_signal_plume_area, na.rm = T))
  }

  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = X11_data_for_plot$river_flow,
                                                 var_ref = X11_data_for_plot$plume_area)

  X11_data_for_plot <- X11_data_for_plot |> mutate(river_flow_scaled = river_flow * scaling_factor$diff + scaling_factor$adjust)

  r_value <- cor(X11_data_for_plot$plume_area, X11_data_for_plot$river_flow, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  the_plot <- ggplot() +

    geom_path(data = X11_data_for_plot, aes(x = dates, y = plume_area), color = "brown", linewidth = 1) +

    geom_path(data = X11_data_for_plot, aes(x = dates, y = river_flow_scaled), color = "blue", linewidth = 1) +

    annotate("text", x = min(X11_data_for_plot$dates), y = Inf, label = r_label,
            hjust = 0, vjust = 1.5, size = 6, colour = "black") +

    scale_x_date(name = "",
                 breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                 labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist()) +

    scale_y_continuous(name = if(show_axis_titles) "Plume area (km²)" else NULL,
                       breaks = scales::breaks_pretty(n = 3),
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = if(show_axis_titles) "River flow (m³ s⁻¹)" else NULL,
                                           breaks = scales::breaks_pretty(n = 3))) +

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
  # discharge series throughout; Grand vs. Petit Rhone: Grand Rhone carries
  # ~90% of the Rhone's combined discharge, see the Rhone-apportionment note
  # in CLAUDE.md) -- used below to bold/enlarge the primary river's label
  # and de-emphasise the others.
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
  # geometry, rather than sharing a single per-zone offset (2026-08-10; the
  # old shared-offset approach clipped the Seine label against its inset's
  # right edge and crowded the two Rhone labels together).
  mouth <- function(river_name) dplyr::filter(zone_river_mouths, river == river_name)
  # fill/alpha added 2026-08-11 (per TODO.md): the labels had no background
  # fill at all (geom_label()'s box border was visible but its interior was
  # effectively transparent), making text hard to read against the busy SPM
  # colour raster underneath. fill is the opposite of the label's own text
  # colour (dark fill behind white text, light fill behind black text) so
  # each stays legible against its own box, at a shared alpha = 0.2 so the
  # coastline/SPM data is still visible through it.
  river_label_style <- function(river_name, primary, colour, ...){
    geom_label(data = mouth(river_name), aes(x = lon, y = lat, label = river),
              fontface = if(primary) "bold" else "plain", size = if(primary) 8 else 6,
              colour = colour, fill = if(colour == "white") "black" else "white", alpha = 0.2, ...)
  }

  river_labels_by_zone <- list(
    BAY_OF_SEINE = list(
      river_label_style("Seine", primary = TRUE, colour = "white", hjust = 1, nudge_x = -0.05, nudge_y = 0.08)
    ),
    SOUTHERN_BRITTANY = list(
      river_label_style("Loire", primary = TRUE, colour = "white", hjust = 0, nudge_x = 0.08, nudge_y = 0.00),
      river_label_style("Vilaine", primary = FALSE, colour = "white", hjust = 0, nudge_x = 0.08, nudge_y = 0.05)
    ),
    BAY_OF_BISCAY = list(
      river_label_style("Gironde", primary = TRUE, colour = "black", hjust = 1, nudge_x = -0.30, nudge_y = 0.05),
      river_label_style("Charente", primary = FALSE, colour = "black", hjust = 1, nudge_x = -0.30, nudge_y = 0.05),
      river_label_style("Sevre", primary = FALSE, colour = "black", hjust = 1, nudge_x = -0.30, nudge_y = 0.05)
    ),
    GULF_OF_LION = list(
      river_label_style("Grand\nRhône", primary = TRUE, colour = "black", hjust = 0, nudge_x = -0.08, nudge_y = -0.22),
      # Nudged down/left rather than up: the Gulf of Lion inset's data only
      # extends to lat 43.6 (0.15 above this river's own mouth), so an
      # upward nudge clips the label against the panel's top edge.
      river_label_style("Petit\nRhône", primary = FALSE, colour = "black", hjust = 1, nudge_x = -0.05, nudge_y = -0.10)
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
    "BAY_OF_SEINE",       0.66,  0.62,  0.24,  0.20,
    "SOUTHERN_BRITTANY",  0.43,  0.50,  0.22,  0.23,
    "BAY_OF_BISCAY",      0.46,  0.25,  0.23,  0.21,
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
# manuscript/figure_table_registry.csv for its current figure number and
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
# (manuscript/figure_table_registry.csv)
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
          plot.tag = element_text(size = 50, face = "bold"),
          plot.tag.position = c(0.02, 0.98))

  if (name_of_the_plot == "B") {
    points_used_for_finding_SPM_threshold <- read_csv(file.path(where_to_save_the_figure,
                                                                "DATA", "B_points_used_for_finding_SPM_threshold.csv"))
    all_points_used_for_finding_SPM_threshold <- read_csv(file.path(where_to_save_the_figure,
                                                                    "DATA", "B_all_points_used_for_finding_SPM_threshold.csv"))
    the_map <- the_map +
      geom_point(data = all_points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "grey50", size = 3) +
      geom_point(data = points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "red", size = 3)
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
    labs(x = "Distance from river mouth (km)", y = "SPM (g m⁻³)", tag = "e)") +
    ggplot_theme() +
    # Legend moved inside the plot area 2026-08-11 (was legend.position =
    # "right", eating into the panel's plotting width) and its text sized up
    # for legibility, per manuscript/TODO.md. Anchored top-right: rendered
    # against the real transect data, SPM decays sharply with distance from
    # the mouth, so nothing is ever plotted in the high-distance/high-SPM
    # corner -- confirmed empty, not assumed.
    theme(plot.tag = element_text(size = 35, face = "bold"), plot.tag.position = c(0.01, 0.98),
          legend.position = "inside", legend.position.inside = c(0.88, 0.8),
          legend.background = element_rect(fill = "white", colour = "grey70"),
          legend.text = element_text(size = 18), legend.title = element_text(size = 18),
          text = element_text(size = 20, colour = "black"),
          axis.text = element_text(size = 18, colour = "black"))

  save_plot_as_png(the_plot, "transect_panel", width = 24, height = 8, path = where_to_save_the_figure)

}


# Renders the per-zone plume-maps panel feeding the plume_methodology_panel
# figure (manuscript/figure_table_registry.csv)
plot_methodology_zone_maps_panel <- function(where_to_save_the_figure) {

  # Read only the four per-zone SPM-map CSVs figure.py's plot_methodology_zone_maps_panel()
  # writes here (Zone.csv, via zone_meta$zone for canonical zone naming/order)
  # -- a plain "*.csv" glob on this shared DATA/ folder also picks up
  # Figure_3_panels()' A-E.csv (which lack a `plume` column entirely) and its
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

    # 20 m bathymetric exclusion boundary:
    # `shallow` is all-False outside the Gulf of Lion (the only zone using
    # this general exclusion), so geom_contour() draws nothing at e)-g).
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


# Daily plume area + SPM mass time series (dynamic threshold, merged
# sensor), with an AR(1)/HAC-weighted trend line, one panel per zone.
# Manuscript slot "plume_area_timeseries" -- see
# manuscript/figure_table_registry.csv for its current figure number.
# where_to_save_the_figure <- 'figures'
plot_plume_area_timeseries <- function(where_to_save_the_figure){
  output_subdir <- get_registry_row("plume_area_timeseries")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)

  # Manual top of each panel's left-hand (plume area) y-axis, set by hand
  # per zone rather than derived from that zone's own data max -- see
  # comment at its call site below.
  manual_area_ylim_top <- c(
    "BAY_OF_SEINE" = 3000,
    "SOUTHERN_BRITTANY" = 7000,
    "BAY_OF_BISCAY" = 11000,
    "GULF_OF_LION" = 11000
  )

  # Static left-hand (plume area) axis tick breaks per zone, set by hand
  # independently of manual_area_ylim_top above -- these don't necessarily
  # reach the axis top (e.g. Southern Brittany's top tick is 6000 against
  # a 7000 limit), they're just the labelled ticks.
  manual_area_breaks <- list(
    "BAY_OF_SEINE" = c(0, 1000, 2000, 3000),
    "SOUTHERN_BRITTANY" = c(0, 2000, 4000, 6000),
    "BAY_OF_BISCAY" = c(0, 3000, 6000, 9000),
    "GULF_OF_LION" = c(0, 3000, 6000, 9000)
  )

  SPM_map_data <- main_folder |> file.path('DATA', 'ts_data.csv') |> read_csv()
  SPM_map_data$Dynamic_threshold <- ifelse(SPM_map_data$Dynamic_threshold, 'Dynamic threshold', 'Fixed threshold')

  # Panel order follows ZONE_ORDER (north to south), matching every other
  # multi-zone figure/table in the project
  SPM_map_data$Zone <- factor(SPM_map_data$Zone, levels = ZONE_ORDER)

  SPM_map_ts <- SPM_map_data |> filter(Dynamic_threshold == 'Dynamic threshold') |> plyr::dlply(c("Zone"), function(df_zone) {

    unique_years <- df_zone$Years |> unique()

    points_for_the_legend <- data.frame(Satellite_sensor = c('merged', 'modis'),
                                        date = c('2020-01-01','2020-01-01') |> as.Date(),
                                        area_of_the_plume_mask_in_km2 = c(-9999,-9999))

    index_to_remove <- which((df_zone$Satellite_sensor == "modis") &
                               (df_zone$area_of_the_plume_mask_in_km2 > quantile(df_zone$area_of_the_plume_mask_in_km2, probs = 0.999, na.rm = TRUE)))

    if (index_to_remove |> length() > 0) {df_zone <- df_zone[-index_to_remove,]}

    df_merged <- df_zone |> filter(Satellite_sensor == "merged") |>
      mutate(mass_t = mass_SPM_in_the_plume_area_in_t)  # already tonnes

    # Mass plotted on a secondary axis, scaled into area's own range --
    # same dual-axis pattern as multi.R::plot_x11_river_and_plume().
    scaling_factor <- sec_axis_adjustement_factors(var_to_scale = df_merged$mass_t,
                                                    var_ref = df_merged$area_of_the_plume_mask_in_km2)
    df_merged <- df_merged |> mutate(mass_scaled = mass_t * scaling_factor$diff + scaling_factor$adjust)

    r_value <- cor(df_merged$area_of_the_plume_mask_in_km2, df_merged$mass_t, use = "complete.obs")
    r_label <- paste0("r (area, mass) = ", sprintf("%.2f", r_value))

    # Left-hand (area) axis top and tick breaks are both set manually per
    # zone, rather than via nice_step_for_target()/four_ticks_from_zero()
    # on the data max, per Robert's request 2026-08-21/2026-08-22. Labels
    # are still rounded to the nearest integer below in case that ever
    # changes.
    zone_key <- as.character(df_zone$Zone[1])
    ylim_top <- manual_area_ylim_top[[zone_key]]
    area_breaks <- manual_area_breaks[[zone_key]]

    mass_breaks <- four_ticks_from_zero(df_merged$mass_t)
    # mass_breaks' own top tick, transformed into area's scale -- rescale
    # mass_breaks to span the same manual ylim_top whenever its natural
    # top tick doesn't already land there, so the right-hand axis' last
    # tick/label isn't stranded below the panel top.
    mass_breaks_top_scaled <- max(mass_breaks) * scaling_factor$diff + scaling_factor$adjust
    if (mass_breaks_top_scaled != ylim_top) {
      mass_target_native <- (ylim_top - scaling_factor$adjust) / scaling_factor$diff
      mass_breaks <- seq(0, mass_target_native, length.out = 4)
    }

    the_ts_plot_wo_modis <- ggplot() +
      geom_point(data = df_merged,
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3", alpha = 0.6) +
      geom_path(data = df_merged,
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") +

      geom_point(data = df_merged, aes(x = date, y = mass_scaled), color = "steelblue4", alpha = 0.6) +
      geom_path(data = df_merged, aes(x = date, y = mass_scaled), color = "steelblue4") +

      annotate("text", x = mean(range(df_merged$date, na.rm = TRUE)), y = Inf, label = r_label,
              hjust = 0.5, vjust = 1.5, size = 6, colour = "black") +

      scale_x_date(name = "",
                   breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                   labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist(),
                   expand = c(0.01,0.01)) +

      coord_cartesian(ylim = c(0, ylim_top)) +
      # No per-panel y-axis title
      scale_y_continuous(name = NULL, breaks = area_breaks, labels = function(b) round(b),
                         sec.axis = sec_axis(transform = ~ (. - scaling_factor$adjust) / scaling_factor$diff,
                                             name = NULL, breaks = mass_breaks,
                                             labels = function(b) format(round(b / 1e5, 1), trim = TRUE))) +
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
            axis.title = element_text(size=30, colour = "black"),
            axis.text.y.left = element_text(color = "red3"),
            axis.ticks.y.left = element_line(color = "red3"),
            axis.text.y.right = element_text(color = "steelblue4"),
            axis.ticks.y.right = element_line(color = "steelblue4"))
    
    
    the_ts_plot_with_modis <- the_ts_plot_wo_modis + 
      geom_point(data = df_zone |> filter(Satellite_sensor == "modis"), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "blue3", alpha = 0.5) + 
      geom_path(data = df_zone |> filter(Satellite_sensor == "modis"), 
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
                     ggarrange(plotlist = SPM_map_ts |> plyr::llply(function(x) {x$wo_modis}), common.legend = FALSE, ncol = 1, nrow = 4, align = "v"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30, color = "red3"),
                     right = text_grob("SPM mass (t x 10⁵)", rot = -90, size = 30, color = "steelblue4")),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}


# Monthly heatmap (sec:results_seasonal) of all four plume properties and
# five drivers, per zone, dynamic threshold. Each tile is that month's median
# value divided by the zone's own all-time median (dynamic threshold only) --
# i.e. how a typical day in that month compares to a typical day across the
# full record for that zone x variable -- so zones of very different raw
# magnitude (see the panache_stats_table slot) are comparable in one figure
# on a scale centred on 1 (= typical). Real (unscaled) interquartile values
# are annotated as text. Also writes the shared long-format data (both thresholds) that
# plot_seasonal_boxplots_dynamic_vs_static() below reads back in, so the
# static-threshold pass is computed once rather than twice. Manuscript slot
# "seasonal_boxplot_heatmap" -- see manuscript/figure_table_registry.csv for
# its current figure number.
plot_seasonal_boxplot_heatmap <- function(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                                       where_are_saved_plume_results_with_static_threshold = "output/panache/static",
                                       where_to_save_the_figure){

  figure_5_output_subdir <- get_registry_row("seasonal_boxplot_heatmap")$output_subdir
  figure_5_dir <- file.path(where_to_save_the_figure, "ARTICLE", figure_5_output_subdir)
  data_dir <- file.path(figure_5_dir, "DATA")
  if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  mass_col <- "mass_SPM_in_the_plume_area_in_t"  # tonnes, see compute_mass_spm_trend.R
  drivers <- c("flow", "wind", "tide", "wave", "current")
  variable_display <- c(
    plume_area    = "Plume area (km²)",
    SPM_mass      = "SPM mass (t)",
    compactness   = "Compactness",
    alongcoast_km = "Along-coast drift (km)",
    flow          = "River flow (m³ s⁻¹)",
    wind          = "Wind speed (m s⁻¹)",
    tide          = "Tidal range (m)",
    wave          = "Wave height (m)",
    current       = "Current speed (m s⁻¹)"
  )

  thresholds <- c(dynamic = where_are_saved_plume_results_with_dynamic_threshold,
                  static = where_are_saved_plume_results_with_static_threshold)

  long_data <- purrr::map_dfr(names(thresholds), function(threshold_label){
    plume_dir <- thresholds[[threshold_label]]

    purrr::pmap_dfr(zone_meta, function(...){
      meta <- tibble::tibble(...)

      df_area <- load_plume_ts(meta$zone, plume_dir = plume_dir, outlier_max = 20000) |>
        dplyr::transmute(date, variable = "plume_area", value = plume_area)

      df_mass <- load_plume_ts(meta$zone, plume_dir = plume_dir, metric_col = mass_col, outlier_max = NULL) |>
        dplyr::transmute(date, variable = "SPM_mass", value = plume_area)  # already tonnes

      # func/compute_plume_shape.py must be run (now wired into
      # code/4_time_series.py) before this figure -- compactness is a
      # required panel, not an optional one, so a missing file is a hard
      # error rather than a silently dropped panel.
      shape_path <- paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv")
      if(!file.exists(shape_path)){
        stop("plot_seasonal_boxplot_heatmap: missing ", shape_path,
             " -- run func/compute_plume_shape.py before regenerating this figure.")
      }
      df_shape <- read_csv(shape_path, show_col_types = FALSE) |>
        dplyr::mutate(date = as.Date(date)) |>
        dplyr::transmute(date, variable = "compactness", value = compactness)

      df_coast <- compute_alongcoast_ts(meta$zone, meta, plume_dir) |>
        dplyr::transmute(date, variable = "alongcoast_km", value = value)

      df_drivers <- purrr::map_dfr(drivers, function(driver_name){
        load_driver(driver_name, meta) |>
          dplyr::transmute(date, variable = driver_name, value = value)
      })

      dplyr::bind_rows(df_area, df_mass, df_shape, df_coast, df_drivers) |>
        dplyr::mutate(threshold = threshold_label, zone = meta$zone, .before = 1)
    })
  }) |>
    dplyr::mutate(category = ifelse(variable %in% drivers, "driver", "property"),
                  month = lubridate::month(date))

  readr::write_csv(long_data, file.path(data_dir, "monthly_boxplot_data.csv"))

  # Zone x variable's own all-time median (dynamic threshold only) -- the
  # "typical day" denominator each month's median is compared against below.
  # A ratio, unlike the old 2nd/98th-percentile-of-range scaling this
  # replaced, is naturally robust to heavy-tailed variables (SPM mass, river
  # flow) and rare extreme-event runs (e.g. Bay of Seine along-coast drift,
  # dominated by the Feb-2014 storm cluster) without needing any winsorising:
  # the median denominator itself already ignores those days.
  overall_median <- long_data |>
    dplyr::filter(threshold == "dynamic") |>
    dplyr::summarise(overall_median = stats::median(value, na.rm = TRUE),
                     .by = c(zone, variable))

  # geom_tile's discrete y-axis places the first factor level at the bottom,
  # so levels are reversed from ZONE_ORDER here to read north (Bay of Seine)
  # at top -> south (Gulf of Lion) at bottom, top-to-bottom -- matching the
  # project's north-to-south panel convention (see ZONE_ORDER, func/util.R)
  # used by every facet_wrap multi-zone figure elsewhere. "Southern Brittany"
  # is abbreviated to "S. Brittany" for this axis label only (not zone_title()
  # itself, which other figures/tables still use unabbreviated) to save
  # left-margin white space in this 3x3 panel grid.
  zone_labels <- zone_title(rev(zones))
  zone_labels[zone_labels == "Southern Brittany"] <- "S. Brittany"

  # Heatmap of monthly medians. One small zone x
  # month heatmap per variable, all 9 (4 properties + 5 drivers) in a single
  # 3x3 panel grid ; the full distributional detail (this exact
  # median plus IQR/range) is still in monthly_boxplot_data.csv (written
  # above), and the per-month linear trend (as opposed to the median shown
  # here) is in the monthly_trend_pct_heatmap slot's figure
  # (func/generate_monthly_trend_pct_heatmap.R) for anyone who needs it.
  heat_stats <- long_data |>
    dplyr::filter(threshold == "dynamic") |>
    dplyr::summarise(month_median = stats::median(value, na.rm = TRUE), .by = c(zone, variable, month)) |>
    dplyr::left_join(overall_median, by = c("zone", "variable")) |>
    dplyr::mutate(ratio = month_median / overall_median,
                  zone = factor(zone, levels = rev(zones), labels = zone_labels),
                  month = factor(month, levels = 1:12, labels = month.abb),
                  variable = factor(variable, levels = names(variable_display), labels = unname(variable_display)))

  # Diverging scale centred on 1 (= that month matches the zone's own
  # all-time typical day); purple/orange rather than the blue/red diverging
  # scale used by the monthly_trend_pct_heatmap slot's %-change-per-year
  # figure (func/generate_monthly_trend_pct_heatmap.R), so the two are not
  # visually conflated.
  p_heatmap <- ggplot(heat_stats, aes(x = month, y = zone, fill = ratio)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_wrap(~variable, ncol = 3) +
    scale_fill_gradient2(name = "Month / overall\nmedian", low = "#7b3294", mid = "white", high = "#e66101", midpoint = 1) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
         axis.text.y = element_text(size = 10), panel.grid = element_blank())

  save_plot_as_png(p_heatmap, registry_basename(figure_5_output_subdir), width = 12, height = 10, path = figure_5_dir)
  message("Wrote ", registry_filename(figure_5_output_subdir))
  invisible(TRUE)
}


# Supplementary "Sx. Lagged daily correlations" figure (fig:daily_flow):
# daily plume area vs. river flow scatter + lagged correlation, per zone.
# Misplaced under a "Deprecated" heading by the 2026-08-11 figure.R cleanup
# pass (only its old main-text role, superseded by plot_seasonal_boxplot_heatmap()
# above, was ever deprecated -- this Supplementary figure itself is still
# live, called from code/5_figures.py); moved back out here. Manuscript slot
# "daily_flow_lagged_correlation" -- see manuscript/figure_table_registry.csv
# for its current figure number.
plot_daily_flow_lagged_correlation <- function(where_to_save_the_figure, max_lag_daily = 14){

  output_subdir <- get_registry_row("daily_flow_lagged_correlation")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)

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

    cor_df <- driver_plume_correlation(df, max_lag_daily = max_lag_daily) |>
      dplyr::filter(timestep == "daily")
    peak <- cor_df |> dplyr::slice_max(cor, n = 1)

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

  plotlist <- panels |> purrr::map(function(p) list(p$scatter, p$lag)) |> purrr::flatten()

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(panels), align = "v",
                                 labels = panel_labels, font.label = list(size = 18, face = "bold"),
                                 hjust = -0.3, vjust = 1.3)

  save_plot_as_png(full_plot, registry_basename(output_subdir), width = 16, height = 18, path = main_folder)
}


# X11 interannual (long-term) signal of plume area vs. river flow, dynamic
# threshold (main results), all four zones. Per-panel axis titles are
# suppressed (show_axis_titles = FALSE) in favour of one shared left/right
# label on the assembled composite, matching the plume_methodology_panel
# figure's convention (annotate_figure(), not a title repeated on all four
# panels). Manuscript slot "x11_interannual_river_flow" -- see
# manuscript/figure_table_registry.csv for its current figure number.
plot_x11_interannual_river_flow <- function(where_to_save_the_figure){
  output_subdir <- get_registry_row("x11_interannual_river_flow")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  zone_plots <- compute_x11_zone_plots(main_folder, show_axis_titles = FALSE)
  save_plot_as_png(annotate_figure(
                     stack_x11_component(zone_plots, "Interannual"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30, color = "brown"),
                     right = text_grob("River flow (m³ s⁻¹)", rot = -90, size = 30, color = "blue")),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}


# Wind, wave, and current direction/magnitude roses, one row per zone,
# coloured by the flow-controlled plume-area response
# (multi.R::plot_driver_rose()). Current column added 2026-08-11 (per
# Robert), a plotting-only addition -- current speed/direction was already
# loaded elsewhere in the pipeline (e.g. the driver_stats_table's driver set)
# under the same column-naming convention plot_driver_rose() already expects.
# Manuscript slot "driver_rose_diagram" -- see
# manuscript/figure_table_registry.csv for its current figure number.
plot_driver_rose_diagram <- function(where_to_save_the_figure, n_sectors = 8){

  output_subdir <- get_registry_row("driver_rose_diagram")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)

  zone_results <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    # One flow join/STL per zone, shared across the wind/wave/current roses
    # below (plot_driver_rose() would otherwise recompute the identical
    # df_flow three times).
    df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)

    summaries <- purrr::map(c("wind", "wave", "current"), compute_driver_rose_summary,
                            meta = meta, n_sectors = n_sectors, df_flow = df_flow) |>
      purrr::keep(~ nrow(.x) > 0)

    # Shared, symmetric colour scale across this zone's wind/wave/current
    # roses, so the same colour means the same residual magnitude in every
    # panel of the row -- each panel would otherwise auto-scale to its own
    # data independently. Falls back to a fixed range only if all three
    # drivers are missing direction data for this zone.
    lim <- if (length(summaries) > 0) max(abs(unlist(purrr::map(summaries, ~ range(.x$mean_area_resid_plot, na.rm = TRUE))))) else 1
    fill_limits <- c(-lim, lim)

    panels <- list(wind = plot_driver_rose("wind", meta, n_sectors, df_flow, fill_limits = fill_limits, show_legend = FALSE),
                  wave = plot_driver_rose("wave", meta, n_sectors, df_flow, fill_limits = fill_limits, show_legend = FALSE),
                  current = plot_driver_rose("current", meta, n_sectors, df_flow, fill_limits = fill_limits, show_legend = FALSE))

    # A legend-only grob for this zone's shared scale, extracted from a
    # throwaway plot built with the same scale -- becomes the row's single
    # colourbar, replacing the 3 per-panel legends suppressed above.
    legend_plot <- ggplot(data.frame(x = 0, y = fill_limits), aes(x = x, y = y, fill = y)) +
      geom_point() +
      scale_fill_gradient2(low = "steelblue", mid = "grey90", high = "firebrick", midpoint = 0,
                          name = "Plume-area\nresidual (km²)", limits = fill_limits) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14),
            legend.key.size = unit(1.1, "cm"))
    colorbar <- ggpubr::as_ggplot(cowplot::get_legend(legend_plot))

    list(panels = panels, colorbar = colorbar)
  })

  plotlist <- purrr::map(zone_results, ~ list(.x$panels$wind, .x$panels$wave, .x$panels$current)) |> purrr::flatten()

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  panel_grid <- ggpubr::ggarrange(plotlist = plotlist, ncol = 3, nrow = nrow(zone_meta), align = "v",
                                 labels = panel_labels, font.label = list(size = 18, face = "bold"),
                                 hjust = -0.3, vjust = 1.3)

  row_labels <- ggpubr::ggarrange(plotlist = purrr::map(zone_meta$zone, ~ ggpubr::text_grob(zone_title(.x), face = "bold", size = 18, rot = 90)),
                                  ncol = 1, nrow = nrow(zone_meta))
  colorbars <- ggpubr::ggarrange(plotlist = purrr::map(zone_results, "colorbar"), ncol = 1, nrow = nrow(zone_meta))
  row_and_panel_grid <- ggpubr::ggarrange(row_labels, colorbars, panel_grid, ncol = 3, widths = c(0.04, 0.13, 1))

  col_labels <- ggpubr::ggarrange(plotlist = purrr::map(c("wind", "wave", "current"),
                                                        ~ ggpubr::text_grob(dplyr::filter(driver_display, driver_name == .x)$driver_label, face = "bold", size = 18)),
                                  ncol = 3, nrow = 1)
  col_labels_row <- ggpubr::ggarrange(ggpubr::text_grob(""), col_labels, ncol = 2, widths = c(0.17, 1))

  full_plot <- ggpubr::ggarrange(col_labels_row, row_and_panel_grid, nrow = 2, heights = c(0.04, 1))

  save_plot_as_png(full_plot, registry_basename(output_subdir), width = 24, height = 20, path = main_folder)
}


# GAM partial-dependence curves for flow, wind, wave, and current, one row
# per zone (driver_interactions.R::fit_gam()/gam_partial_effect()). Tide is
# intentionally excluded from the plot -- since tidal range is an
# essentially fixed astronomical property of each site rather than something
# worth a dedicated panel; it stays in the underlying GAM/driver_stats_table
# statistics, just not visualised here. Manuscript slot "gam_partial_effects"
# -- see manuscript/figure_table_registry.csv for its current figure number.
plot_gam_partial_effects <- function(where_to_save_the_figure, stats_dir = "output/STATS"){

  # Sourced here rather than at file scope (unlike multi.R above): this pulls
  # in several heavyweight modelling packages (mgcv, gratia, ranger, iml)
  # that only this one figure function needs -- not worth loading for every
  # other figure in this file.
  source("func/driver_interactions.R")

  output_subdir <- get_registry_row("gam_partial_effects")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)

  driver_labels <- c(flow = "River flow (m³ s⁻¹)", wind_spd = "Wind speed (m s⁻¹)",
                     wave_height = "Wave height (m)", current = "Current speed (m s⁻¹)")
  n_drivers <- length(driver_labels)

  # No per-panel titles/axis titles here -- a title on only the first column
  # of each row (the old approach) gives that panel a shorter plot area than
  # its title-less row-mates, so the panels no longer line up vertically
  # despite ggarrange's align = "v". Row/column labelling is added once,
  # externally, after the grid is built (zone name per row below; x-axis
  # title on the bottom row only, since every row in a column shares that
  # driver's x-axis; single shared y-axis title via annotate_figure()).
  zone_rows <- purrr::map(zones, function(zone_name){
    df <- readr::read_csv(file.path(stats_dir, paste0("daily_driver_matrix_", zone_name, ".csv")), show_col_types = FALSE)
    gam_model <- fit_gam(df)
    drivers_to_show <- intersect(names(driver_labels), available_drivers(df))
    if(!setequal(drivers_to_show, names(driver_labels))) stop("plot_gam_partial_effects(): zone ", zone_name,
      " does not have the full driver set (", paste(names(driver_labels), collapse = ", "),
      "); the fixed 4-column grid assumes every zone does.")

    purrr::map(drivers_to_show, function(d){
      curve <- gam_partial_effect(gam_model, d, df)
      ggplot(curve, aes(x = x, y = fit)) +
        geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), fill = "grey80", alpha = 0.5) +
        geom_line(colour = "black", linewidth = 1) +
        labs(x = NULL, y = NULL) +
        theme(panel.border = element_rect(fill = NA, colour = "black"))
    })
  })

  # Bottom row only: one x-axis title per column.
  zone_rows[[length(zones)]] <- purrr::map2(zone_rows[[length(zones)]], names(driver_labels),
    function(p, d) p + labs(x = driver_labels[[d]]))

  panel_labels <- paste0(letters[seq_len(length(zones) * n_drivers)], ")")
  # hjust/vjust pushed further in than the driver_rose_diagram/daily_flow_lagged_correlation figures' convention
  # (hjust=-0.3, vjust=1.3): this grid's panels are much smaller (4x4 vs.
  # 2-column), so the same absolute offset landed the tag on top of the
  # topmost y-axis tick label instead of clear of it.
  panel_grid <- ggpubr::ggarrange(plotlist = purrr::flatten(zone_rows), ncol = n_drivers, nrow = length(zones),
                                  align = "hv", labels = panel_labels, font.label = list(size = 14, face = "bold"),
                                  hjust = -0.8, vjust = 1.8)

  row_labels <- ggpubr::ggarrange(plotlist = purrr::map(zones, ~ ggpubr::text_grob(zone_title(.x), face = "bold", size = 16, rot = 90)),
                                  ncol = 1, nrow = length(zones))

  full_plot <- ggpubr::annotate_figure(
    ggpubr::ggarrange(row_labels, panel_grid, ncol = 2, widths = c(0.05, 1)),
    left = ggpubr::text_grob("Partial effect on plume area (km²)", rot = 90, size = 20))

  save_plot_as_png(full_plot, registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}


# Computes the per-zone X11 Interannual/Seasonal/Residual plots for one
# threshold (dynamic or static) from the weekly plume/river time series
# prepped by figure.py's Figure_X11_weekly_results(). Pure computation, no
# saving -- shared by the four Figure_*_x11_*() functions below so the
# per-zone join/plot logic isn't duplicated per manuscript figure.
compute_x11_zone_plots <- function(data_dir, show_axis_titles = TRUE){
  plume_data <- data_dir |> file.path('DATA', 'ts_plume_data.csv') |> read_csv()
  river_data <- data_dir |> file.path('DATA', 'ts_river_data.csv') |> read_csv()

  regions <- unique(plume_data$Zone) |> order_zones()

  regions |> plyr::llply(function(region) {

    plume_data_region <- plume_data |> filter(Zone == region)
    river_data_region <- river_data |> filter(Zone == region)

    X11_ts <- plume_data_region |> inner_join(river_data_region, by = "dates", suffix = c("_plume_area", "_river_flow"))

    zone_label <- zone_title(region)
    list("Interannual" = plot_x11_river_and_plume(X11_ts, type_of_signal = 'Interannual', show_axis_titles = show_axis_titles) + labs(title = zone_label),
        "Seasonal" = plot_x11_river_and_plume(X11_ts, type_of_signal = 'Seasonal', show_axis_titles = show_axis_titles) + labs(title = zone_label),
        "Residual" = plot_x11_river_and_plume(X11_ts, type_of_signal = 'Residual', show_axis_titles = show_axis_titles) + labs(title = zone_label))

  })
}

# Stacks one X11 component (Interannual/Seasonal/Residual) across all 4
# zones into a single column -- the layout shared by every figure below.
stack_x11_component <- function(zone_plots, component){
  ggarrange(plotlist = zone_plots |> plyr::llply(function(x) x[[component]]), ncol = 1, nrow = 4, align = "v")
}

# Plume-area time series, fixed vs. dynamic threshold comparison, one panel
# per zone. Reads the same ts_data.csv plot_plume_area_timeseries() does
# (both share one Python-side data prep), but writes to its own folder.
# Manuscript slot "thresholds_comparison" -- see
# manuscript/figure_table_registry.csv for its current figure number.
# where_to_save_the_figure <- 'figures'
plot_threshold_comparison <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("plume_area_timeseries")$output_subdir)
  output_subdir <- get_registry_row("thresholds_comparison")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)

  SPM_map_data <- data_dir |> file.path('DATA', 'ts_data.csv') |> read_csv()
  SPM_map_data$Dynamic_threshold <- ifelse(SPM_map_data$Dynamic_threshold, 'Dynamic threshold', 'Fixed threshold')

  # Panel order follows ZONE_ORDER (north to south), matching the
  # plume_area_timeseries figure and every other multi-zone figure/table --
  # dlply()'s default grouping would otherwise sort panels alphabetically.
  SPM_map_data$Zone <- factor(SPM_map_data$Zone, levels = ZONE_ORDER)

  SPM_map_ts <- SPM_map_data |> filter(Satellite_sensor == "merged") |> plyr::dlply(c("Zone"), function(df_zone) {
    
    unique_years <- df_zone$Years |> unique()
    
    points_for_the_legend <- data.frame(Dynamic_threshold = c('Dynamic threshold', 'Fixed threshold'),
                                        date = c('2020-01-01','2020-01-01') |> as.Date(),
                                        area_of_the_plume_mask_in_km2 = c(-9999,-9999))
    
    # index_to_remove <- which((df_zone$Satellite_sensor == "modis") & 
    #                            (df_zone$area_of_the_plume_mask_in_km2 > quantile(df_zone$area_of_the_plume_mask_in_km2, probs = 0.999, na.rm = TRUE)))
    # 
    # if (index_to_remove |> length() > 0) {df_zone <- df_zone[-index_to_remove,]}
    
    the_ts_plot <- ggplot() + 
      
      geom_point(data = df_zone |> filter(Dynamic_threshold == 'Dynamic threshold'), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") + 
      geom_path(data = df_zone |> filter(Dynamic_threshold == 'Dynamic threshold'), 
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "red3") + 
      
      scale_x_date(name = "", 
                   breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(), 
                   labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist(),
                   expand = c(0.01,0.01)) +
      
      coord_cartesian(ylim = c(0, max(df_zone$area_of_the_plume_mask_in_km2, na.rm = TRUE))) +
      # No per-panel y-axis title -- a single shared label is added once via
      # annotate_figure() on the assembled composite below, matching the
      # plume_area_timeseries figure.
      labs(y = NULL, x = "", title = zone_title(df_zone$Zone[1])) +
      ggplot_theme() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = c(.9,.9),
            legend.background = element_rect(fill = "transparent"),
            plot.title = element_text(size=30, colour = "black"),
            text = element_text(size=25, colour = "black"),
            axis.text = element_text(size=20, colour = "black"),
            axis.title = element_text(size=30, colour = "black"))  +
      
      geom_point(data = df_zone |> filter(Dynamic_threshold == 'Fixed threshold'), 
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "chartreuse4", alpha = 0.5) + 
      geom_path(data = df_zone |> filter(Dynamic_threshold == 'Fixed threshold'), 
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "chartreuse4", alpha = 0.5) + 
      geom_point(data = points_for_the_legend, aes(x = date, y = area_of_the_plume_mask_in_km2, color = Dynamic_threshold), size = 0.1) +
      
      scale_color_manual(values = c('Dynamic threshold'= "red3", 'Fixed threshold' = "chartreuse4"), name = "") +
      
      guides(
        color = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                             override.aes = list(size = c(5, 5),
                                                 alpha = c(1, 0.5)))) 
    
    return(the_ts_plot)
    
  })
  
  save_plot_as_png(annotate_figure(
                     ggarrange(plotlist = SPM_map_ts, common.legend = FALSE, ncol = 1, nrow = 4, align = "v"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30)),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)

}

# Plots one X11 component of plume area under the dynamic vs. static thresholds
plot_x11_dynamic_vs_static <- function(X11_data, type_of_signal) {

  unique_years <- X11_data$dates |> year() |> unique()

  X11_data_for_plot <- X11_data |>
    rename(dynamic = !!sym(paste(type_of_signal, "signal_dynamic", sep = "_")),
           static = !!sym(paste(type_of_signal, "signal_static", sep = "_"))) |>
    select(dates, dynamic, static)

  if (type_of_signal %in% c("Seasonal", "Residual")) {
    X11_data_for_plot <- X11_data_for_plot |>
      mutate(dynamic = dynamic + mean(X11_data$Raw_signal_dynamic, na.rm = T),
             static = static + mean(X11_data$Raw_signal_static, na.rm = T))
  }

  r_value <- cor(X11_data_for_plot$dynamic, X11_data_for_plot$static, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  ggplot() +

    geom_point(data = X11_data_for_plot, aes(x = dates, y = static), color = "chartreuse4", alpha = 0.5) +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = static), color = "chartreuse4", alpha = 0.5) +

    geom_point(data = X11_data_for_plot, aes(x = dates, y = dynamic), color = "red3") +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = dynamic), color = "red3") +

    annotate("text", x = min(X11_data_for_plot$dates), y = Inf, label = r_label,
            hjust = 0, vjust = 1.5, size = 6, colour = "black") +

    scale_x_date(name = "",
                 breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                 labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist()) +

    scale_y_continuous(name = "Plume area (km²)") +

    labs(title = paste(type_of_signal, "signal")) +
    ggplot_theme() +

    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(size=30, colour = "black"),
          text = element_text(size=25, colour = "black"),
          axis.text = element_text(size=20, colour = "black"),
          axis.title = element_text(size=30, colour = "black"),
          panel.border = element_rect(linetype = "solid", fill = NA))
}

# Computes the per-zone dynamic-vs-static plume-area comparison plots for
# one X11 component -- mirrors compute_x11_zone_plots() but joins the
# dynamic- and static-threshold plume-area series (by date) instead of
# plume area and river flow.
compute_x11_dynamic_vs_static_plots <- function(data_dir){
  ts_data <- data_dir |> file.path('DATA', 'ts_plume_dynamic_vs_static.csv') |> read_csv()

  regions <- unique(ts_data$Zone) |> order_zones()

  regions |> plyr::llply(function(region) {
    ts_dynamic <- ts_data |> filter(Zone == region, threshold == "dynamic") |> select(-Zone, -threshold)
    ts_static  <- ts_data |> filter(Zone == region, threshold == "static") |> select(-Zone, -threshold)

    X11_ts <- ts_dynamic |> inner_join(ts_static, by = "dates", suffix = c("_dynamic", "_static"))

    zone_label <- zone_title(region)
    list("Interannual" = plot_x11_dynamic_vs_static(X11_ts, type_of_signal = 'Interannual') + labs(title = zone_label),
        "Seasonal" = plot_x11_dynamic_vs_static(X11_ts, type_of_signal = 'Seasonal') + labs(title = zone_label),
        "Residual" = plot_x11_dynamic_vs_static(X11_ts, type_of_signal = 'Residual') + labs(title = zone_label))
  })
}

# X11 seasonal component of plume area vs. river flow, dynamic threshold,
# all four zones. Shares x11_interannual_river_flow's DATA/ prep. Manuscript
# slot "x11_seasonal_river_flow" -- see manuscript/figure_table_registry.csv
# for its current figure number. Split 2026-08-26 out of the former
# plot_x11_components_dynamic(), which rendered this and the residual
# component as one stacked seasonal-on-top-of-residual composite image via
# save_x11_component_composite() (now deleted) -- that stacked two already
# 4-zone-panel figures into one 8-panel image, causing page overflow in the
# compiled PDF. Each component now gets its own standalone 4-panel figure,
# same as every other X11 figure in this file.
plot_x11_seasonal_river_flow <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("x11_interannual_river_flow")$output_subdir)
  output_subdir <- get_registry_row("x11_seasonal_river_flow")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_zone_plots(data_dir)
  save_plot_as_png(stack_x11_component(zone_plots, "Seasonal"), registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# X11 residual (short-term) component of plume area vs. river flow, dynamic
# threshold, all four zones. Shares x11_interannual_river_flow's DATA/ prep.
# Manuscript slot "x11_residual_river_flow" -- see
# manuscript/figure_table_registry.csv for its current figure number.
# Renamed/split 2026-08-26 from plot_x11_components_dynamic() (see
# plot_x11_seasonal_river_flow() above for why); this function renders
# residual only, matching the slot's narrowed role now that seasonal has
# moved to the new main-text x11_seasonal_river_flow slot.
plot_x11_residual_river_flow <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("x11_interannual_river_flow")$output_subdir)
  output_subdir <- get_registry_row("x11_residual_river_flow")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_zone_plots(data_dir)
  save_plot_as_png(stack_x11_component(zone_plots, "Residual"), registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# X11 interannual signal of plume area, dynamic vs. static threshold, all
# four zones. Manuscript slot "x11_interannual_dynamic_vs_static" -- see
# manuscript/figure_table_registry.csv for its current figure number.
plot_x11_interannual_dynamic_vs_static <- function(where_to_save_the_figure){
  output_subdir <- get_registry_row("x11_interannual_dynamic_vs_static")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  zone_plots <- compute_x11_dynamic_vs_static_plots(main_folder)
  save_plot_as_png(stack_x11_component(zone_plots, "Interannual"), registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# X11 seasonal signal of plume area, dynamic vs. static threshold. Shares
# x11_interannual_dynamic_vs_static's DATA/ prep. One of three separate
# figures (interannual/seasonal/residual) requested in manuscript/TODO.md,
# rather than a seasonal+residual composite. Manuscript slot
# "x11_seasonal_dynamic_vs_static" -- see manuscript/figure_table_registry.csv
# for its current figure number.
plot_x11_seasonal_dynamic_vs_static <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("x11_interannual_dynamic_vs_static")$output_subdir)
  output_subdir <- get_registry_row("x11_seasonal_dynamic_vs_static")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_dynamic_vs_static_plots(data_dir)
  save_plot_as_png(stack_x11_component(zone_plots, "Seasonal"), registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# X11 residual variance of plume area, dynamic vs. static threshold. Shares
# x11_interannual_dynamic_vs_static's DATA/ prep. Third of the three separate
# figures requested in manuscript/TODO.md. Manuscript slot
# "x11_residual_dynamic_vs_static" -- see manuscript/figure_table_registry.csv
# for its current figure number.
plot_x11_residual_dynamic_vs_static <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("x11_interannual_dynamic_vs_static")$output_subdir)
  output_subdir <- get_registry_row("x11_residual_dynamic_vs_static")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_dynamic_vs_static_plots(data_dir)
  save_plot_as_png(stack_x11_component(zone_plots, "Residual"), registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# Monthly (previously JJA vs. NDJ) dynamic-vs-static threshold comparison of
# the four plume properties, all four zones. Drivers
# are deliberately not shown here: driver values
# don't depend on the plume-detection threshold at all, so their dynamic and
# static boxes are only ever near-identical up to sampling noise -- not an
# informative threshold comparison the way the plume properties are.
# Reads the shared long-format data plot_seasonal_boxplot_heatmap() above
# already wrote to seasonal_boxplot_heatmap's DATA/monthly_boxplot_data.csv,
# rather than re-reading Results.csv independently, so the two-threshold
# computation only happens once -- this function must therefore be called
# after plot_seasonal_boxplot_heatmap() (Figure_5_seasonal_analysis()) in the
# code/5_figures.py pipeline. Values are shown on the same
# 0-100%-of-dynamic-range scale as that figure, so a static-threshold box
# sitting outside 0-100% is a direct visual signal that the two thresholds
# disagree, not scaling noise. Manuscript slot
# "seasonal_boxplots_dynamic_vs_static" -- see
# manuscript/figure_table_registry.csv for its current figure number.
plot_seasonal_boxplots_dynamic_vs_static <- function(where_to_save_the_figure){
  output_subdir <- get_registry_row("seasonal_boxplots_dynamic_vs_static")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  data_path <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("seasonal_boxplot_heatmap")$output_subdir, "DATA", "monthly_boxplot_data.csv")

  if(!file.exists(data_path)){
    message("seasonal_boxplots_dynamic_vs_static: shared data file not found (", data_path,
           ") -- run plot_seasonal_boxplot_heatmap() first. Skipping.")
    return(invisible(FALSE))
  }

  variable_display <- c(
    plume_area    = "Plume area (km²)",
    SPM_mass      = "SPM mass (t)",
    compactness   = "Compactness",
    alongcoast_km = "Along-coast drift (km)",
    flow          = "River flow (m³ s⁻¹)",
    wind          = "Wind speed (m s⁻¹)",
    tide          = "Tidal range (m)",
    wave          = "Wave height (m)",
    current       = "Current speed (m s⁻¹)"
  )

  long_data <- readr::read_csv(data_path, show_col_types = FALSE)

  # Own dynamic-threshold-only 2nd/98th-percentile-of-range scale (fixed
  # 2026-08-11): boxplot whiskers below use per-day pct extremes directly, so
  # a robust-but-fixed range matters here to keep rare extreme-event days
  # from setting the whole 28-year scale. plot_seasonal_boxplot_heatmap()
  # above now uses a different scale (month median / zone's own all-time
  # median, since it only plots that single more-robust summary per tile),
  # so the two figures' scales are no longer directly comparable -- this one
  # stands alone.
  scale_range <- long_data |>
    dplyr::filter(threshold == "dynamic") |>
    dplyr::summarise(range_min = stats::quantile(value, 0.02, na.rm = TRUE),
                     range_max = stats::quantile(value, 0.98, na.rm = TRUE),
                     .by = c(zone, variable))

  df <- long_data |>
    dplyr::left_join(scale_range, by = c("zone", "variable")) |>
    dplyr::mutate(pct = 100 * (value - range_min) / (range_max - range_min),
                  zone = factor(zone, levels = zones, labels = zone_title(zones)),
                  month = factor(month, levels = 1:12, labels = month.abb),
                  variable = factor(variable, levels = names(variable_display), labels = unname(variable_display)),
                  threshold = factor(threshold, levels = c("dynamic", "static"), labels = c("Dynamic", "Static")))

  box_stats <- df |>
    dplyr::filter(category == "property") |>
    dplyr::summarise(
      ymin = min(pct, na.rm = TRUE), lower = stats::quantile(pct, 0.25, na.rm = TRUE),
      middle = stats::median(pct, na.rm = TRUE), upper = stats::quantile(pct, 0.75, na.rm = TRUE),
      ymax = max(pct, na.rm = TRUE),
      .by = c(zone, variable, month, threshold)
    )

  p_properties <- ggplot(box_stats, aes(x = month, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
                                        fill = threshold)) +
    geom_boxplot(stat = "identity", position = position_dodge(0.75), width = 0.65, linewidth = 0.3) +
    facet_grid(variable ~ zone, scales = "free_x") +
    scale_fill_manual(values = c("Dynamic" = "#2166ac", "Static" = "#d6604d"), name = "Threshold") +
    labs(x = NULL, y = "% of zone's own observed dynamic-threshold range") +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0, size = 7), strip.text.x = element_text(size = 9),
         axis.text.x = element_text(angle = 45, hjust = 1, size = 6), legend.position = "bottom",
         panel.grid.minor = element_blank())

  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  save_plot_as_png(p_properties, registry_basename(output_subdir), width = 12, height = 9, path = main_folder)

  message("Wrote Figure_S3.png (plume properties only, monthly dynamic-vs-static comparison)")
  invisible(TRUE)
}


# Deprecated -------------------------------------------------------------

# Figure XXX: flow-controlled plume-area residual vs. wave height,
# coloured by on/off-shore wind category, one panel per zone
# (multi.R::plot_category_scatter()).
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

