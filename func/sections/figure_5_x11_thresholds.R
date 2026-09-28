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
# common_legend/legend_position are only used by the dynamic-vs-static
# comparison plots (plot_x11_component_dynamic_vs_static()), whose panels
# carry a Dynamic/Fixed threshold colour legend; the other callers' panels
# have no legend, so the defaults leave their output unchanged.
stack_x11_component <- function(zone_plots, component, common_legend = FALSE, legend_position = "none"){
  ggarrange(plotlist = zone_plots |> plyr::llply(function(x) x[[component]]), ncol = 1, nrow = 4, align = "v",
           common.legend = common_legend, legend = legend_position)
}

# Plume-area time series, fixed vs. dynamic threshold comparison, one panel
# per zone. Reads the same ts_data.csv plot_plume_area_timeseries() does
# (both share one Python-side data prep), but writes to its own folder.
# Manuscript slot "thresholds_comparison" -- see
# metadata/figure_table_registry.csv for its current figure number.
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

    # Per-zone Pearson correlation between the dynamic- and fixed-threshold
    # area series (paired by date), annotated the same way as the X11
    # dynamic-vs-static comparison plots (plot_x11_dynamic_vs_static() below).
    joined_vals <- df_zone |> filter(Dynamic_threshold == 'Dynamic threshold') |> select(date, area_of_the_plume_mask_in_km2) |>
      inner_join(df_zone |> filter(Dynamic_threshold == 'Fixed threshold') |> select(date, area_of_the_plume_mask_in_km2),
                 by = "date", suffix = c("_dynamic", "_fixed"))
    r_value <- cor(joined_vals$area_of_the_plume_mask_in_km2_dynamic, joined_vals$area_of_the_plume_mask_in_km2_fixed, use = "complete.obs")
    r_label <- paste0("r = ", sprintf("%.2f", r_value))

    the_ts_plot <- ggplot() +

      geom_point(data = df_zone |> filter(Dynamic_threshold == 'Dynamic threshold'),
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "#E69F00") +
      geom_path(data = df_zone |> filter(Dynamic_threshold == 'Dynamic threshold'),
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "#E69F00") +

      annotate("text", x = min(df_zone$date), y = Inf, label = r_label,
              hjust = 0, vjust = 1.5, size = 6, colour = "black") +

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
            legend.position = "bottom",
            legend.text = element_text(size = 30, colour = "black"),
            plot.title = element_text(size=30, colour = "black"),
            text = element_text(size=25, colour = "black"),
            axis.text = element_text(size=20, colour = "black"),
            axis.title = element_text(size=30, colour = "black"))  +

      geom_point(data = df_zone |> filter(Dynamic_threshold == 'Fixed threshold'),
                 aes(x = date, y = area_of_the_plume_mask_in_km2), color = "#CC79A7", alpha = 0.5) +
      geom_path(data = df_zone |> filter(Dynamic_threshold == 'Fixed threshold'),
                aes(x = date, y = area_of_the_plume_mask_in_km2), color = "#CC79A7", alpha = 0.5) +
      geom_point(data = points_for_the_legend, aes(x = date, y = area_of_the_plume_mask_in_km2, color = Dynamic_threshold), size = 0.1) +

      scale_color_manual(values = c('Dynamic threshold'= "#E69F00", 'Fixed threshold' = "#CC79A7"), name = "") +

      guides(
        color = guide_legend(keyheight = unit(1, "cm"), keywidth = unit(1.5, "cm"), byrow = TRUE,
                             override.aes = list(size = c(5, 5),
                                                 alpha = c(1, 0.5))))

    return(the_ts_plot)

  })

  save_plot_as_png(annotate_figure(
                     ggarrange(plotlist = SPM_map_ts, common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 4, align = "v"),
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

    geom_point(data = X11_data_for_plot, aes(x = dates, y = static, color = "Fixed threshold"), alpha = 0.5) +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = static, color = "Fixed threshold"), alpha = 0.5) +

    geom_point(data = X11_data_for_plot, aes(x = dates, y = dynamic, color = "Dynamic threshold")) +
    geom_path(data = X11_data_for_plot, aes(x = dates, y = dynamic, color = "Dynamic threshold")) +

    annotate("text", x = min(X11_data_for_plot$dates), y = Inf, label = r_label,
            hjust = 0, vjust = 1.5, size = 6, colour = "black") +

    scale_color_manual(values = c('Dynamic threshold' = "#E69F00", 'Fixed threshold' = "#CC79A7"), name = "") +

    # Dot-only legend keys (no line segment), matching the
    # plot_threshold_comparison() (Fig. S2) legend style above.
    guides(
      color = guide_legend(keyheight = unit(1, "cm"), keywidth = unit(1.5, "cm"), byrow = TRUE,
                           override.aes = list(linetype = 0, shape = 16, size = c(5, 5),
                                               alpha = c(1, 0.5)))) +

    scale_x_date(name = "",
                 breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                 labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist()) +

    scale_y_continuous(name = "Plume area (km²)") +

    labs(title = paste(type_of_signal, "signal")) +
    ggplot_theme() +

    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size = 30, colour = "black"),
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
# slot "x11_seasonal_river_flow" -- see metadata/figure_table_registry.csv
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
  zone_plots <- compute_x11_zone_plots(data_dir, show_axis_titles = FALSE)
  save_plot_as_png(annotate_figure(
                     stack_x11_component(zone_plots, "Seasonal"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30, color = "brown"),
                     right = text_grob("River flow (m³ s⁻¹)", rot = -90, size = 30, color = "blue")),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# X11 residual (short-term) component of plume area vs. river flow, dynamic
# threshold, all four zones. Shares x11_interannual_river_flow's DATA/ prep.
# Manuscript slot "x11_residual_river_flow" -- see
# metadata/figure_table_registry.csv for its current figure number.
# Renamed/split 2026-08-26 from plot_x11_components_dynamic() (see
# plot_x11_seasonal_river_flow() above for why); this function renders
# residual only, matching the slot's narrowed role now that seasonal has
# moved to the new main-text x11_seasonal_river_flow slot.
plot_x11_residual_river_flow <- function(where_to_save_the_figure){
  data_dir <- file.path(where_to_save_the_figure, "ARTICLE", get_registry_row("x11_interannual_river_flow")$output_subdir)
  output_subdir <- get_registry_row("x11_residual_river_flow")$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)
  zone_plots <- compute_x11_zone_plots(data_dir, show_axis_titles = FALSE)
  save_plot_as_png(annotate_figure(
                     stack_x11_component(zone_plots, "Residual"),
                     left = text_grob("Plume area (km²)", rot = 90, size = 30, color = "brown"),
                     right = text_grob("River flow (m³ s⁻¹)", rot = -90, size = 30, color = "blue")),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# X11 signal of plume area, dynamic vs. static threshold, all four zones --
# one component (interannual/seasonal/residual) per call. Refactored
# 2026-09-23 (Robert's call) from three near-identical functions
# (plot_x11_interannual_dynamic_vs_static/plot_x11_seasonal_dynamic_vs_static/
# plot_x11_residual_dynamic_vs_static, ~10 duplicated lines each differing
# only in registry slot_key and stack_x11_component()'s component argument)
# down to this one shared implementation. The three registry-facing entry
# points below are kept as one-line dispatchers since
# func/figure.py::Figure_X11_weekly_results() looks up and calls each slot's
# r_function by name generically (no type_of_signal argument passed through)
# -- collapsing further would mean special-casing that dispatch loop, a
# larger change for no real benefit. All three share the same DATA/ prep
# (data_dir always resolves to the interannual slot's own folder, since
# func/figure.py::_prep_x11_dynamic_vs_static_data() only ever writes
# ts_plume_dynamic_vs_static.csv there -- true for the interannual call too,
# since its own output_subdir *is* that slot).
plot_x11_component_dynamic_vs_static <- function(where_to_save_the_figure, type_of_signal){
  slot_key <- switch(type_of_signal,
                     "Interannual" = "x11_interannual_dynamic_vs_static",
                     "Seasonal"    = "x11_seasonal_dynamic_vs_static",
                     "Residual"    = "x11_residual_dynamic_vs_static")
  output_subdir <- get_registry_row(slot_key)$output_subdir
  main_folder <- file.path(where_to_save_the_figure, "ARTICLE", output_subdir)
  if (!dir.exists(main_folder)) dir.create(main_folder, recursive = TRUE)

  data_dir <- file.path(where_to_save_the_figure, "ARTICLE",
                        get_registry_row("x11_interannual_dynamic_vs_static")$output_subdir)
  zone_plots <- compute_x11_dynamic_vs_static_plots(data_dir)
  save_plot_as_png(stack_x11_component(zone_plots, type_of_signal, common_legend = TRUE, legend_position = "bottom"),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}

# Manuscript slot "x11_interannual_dynamic_vs_static" -- see
# metadata/figure_table_registry.csv for its current figure number.
plot_x11_interannual_dynamic_vs_static <- function(where_to_save_the_figure){
  plot_x11_component_dynamic_vs_static(where_to_save_the_figure, "Interannual")
}

# Manuscript slot "x11_seasonal_dynamic_vs_static" -- see
# metadata/figure_table_registry.csv for its current figure number.
plot_x11_seasonal_dynamic_vs_static <- function(where_to_save_the_figure){
  plot_x11_component_dynamic_vs_static(where_to_save_the_figure, "Seasonal")
}

# Manuscript slot "x11_residual_dynamic_vs_static" -- see
# metadata/figure_table_registry.csv for its current figure number.
plot_x11_residual_dynamic_vs_static <- function(where_to_save_the_figure){
  plot_x11_component_dynamic_vs_static(where_to_save_the_figure, "Residual")
}
