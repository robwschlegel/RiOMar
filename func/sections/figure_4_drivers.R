# Wind, wave, and current direction/magnitude roses, one row per zone,
# coloured by the flow-controlled plume-area response
# (multi.R::plot_driver_rose()). Current column added 2026-08-11 (per
# Robert), a plotting-only addition -- current speed/direction was already
# loaded elsewhere in the pipeline (e.g. the driver_stats_table's driver set)
# under the same column-naming convention plot_driver_rose() already expects.
# Manuscript slot "driver_rose_diagram" -- see
# metadata/figure_table_registry.csv for its current figure number.
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

    summaries <- purrr::map(c("wind", "wave", "current"), function(d)
                              compute_driver_rose_summary(d, meta = meta, n_sectors = n_sectors, df_flow = df_flow) |>
                                dplyr::mutate(zone = meta$zone, driver = d, .before = 1)) |>
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

    list(panels = panels, colorbar = colorbar, summaries = dplyr::bind_rows(summaries))
  })

  # Per-sector day share and mean flow-controlled residual behind every rose,
  # so the text can cite them (written alongside the figure, as other figures' DATA/)
  if (!dir.exists(file.path(main_folder, "DATA"))) dir.create(file.path(main_folder, "DATA"))
  readr::write_csv(dplyr::bind_rows(purrr::map(zone_results, "summaries")),
                   file.path(main_folder, "DATA", "rose_sector_summary.csv"))

  plotlist <- purrr::map(zone_results, ~ list(.x$panels$wind, .x$panels$wave, .x$panels$current)) |> purrr::flatten()

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  panel_grid <- ggpubr::ggarrange(plotlist = plotlist, ncol = 3, nrow = nrow(zone_meta), align = "v",
                                 labels = panel_labels, font.label = list(size = 18, face = "bold"),
                                 hjust = -0.3, vjust = 1.3)

  row_labels <- ggpubr::ggarrange(plotlist = purrr::map(zone_meta$zone, ~ ggpubr::text_grob(zone_title(.x), face = "bold", size = 18, rot = 90)),
                                  ncol = 1, nrow = nrow(zone_meta))
  colorbars <- ggpubr::ggarrange(plotlist = purrr::map(zone_results, "colorbar"), ncol = 1, nrow = nrow(zone_meta))
  row_and_panel_grid <- ggpubr::ggarrange(row_labels, colorbars, panel_grid, ncol = 3, widths = c(0.04, 0.13, 1))

  # The roses show direction, not magnitude, so the columns are titled by
  # direction (driver_display's labels are the speed/height ones, and their
  # superscript minus also rendered as a missing glyph here).
  col_labels <- ggpubr::ggarrange(plotlist = purrr::map(c("Wind direction", "Wave direction", "Current direction"),
                                                        ~ ggpubr::text_grob(.x, face = "bold", size = 18)),
                                  ncol = 3, nrow = 1)
  col_labels_row <- ggpubr::ggarrange(ggpubr::text_grob(""), col_labels, ncol = 2, widths = c(0.17, 1))

  full_plot <- ggpubr::ggarrange(col_labels_row, row_and_panel_grid, nrow = 2, heights = c(0.04, 1))

  save_plot_as_png(full_plot, registry_basename(output_subdir), width = 20, height = 22, path = main_folder)
}


# GAM partial-dependence curves for flow, wind, wave, and current, one row
# per zone (driver_interactions.R::fit_gam()/gam_partial_effect()). Tide is
# intentionally excluded from the plot -- since tidal range is an
# essentially fixed astronomical property of each site rather than something
# worth a dedicated panel; it stays in the underlying GAM/driver_stats_table
# statistics, just not visualised here. Manuscript slot "gam_partial_effects"
# -- see metadata/figure_table_registry.csv for its current figure number.
plot_gam_partial_effects <- function(where_to_save_the_figure, stats_dir = "output/STATS"){

  # Sourced here rather than at file scope (unlike multi.R above): this pulls
  # in heavyweight modelling packages (mgcv, ranger)
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
      # Clipped to the driver's 2nd-98th percentile of observed values: the
      # sparse tails beyond it carry SE bands several times the effect itself
      # (co-author meeting 2026-10-01), which set the y-range for nothing.
      x_lim <- stats::quantile(df[[d]], c(0.02, 0.98), na.rm = TRUE)
      curve <- gam_partial_effect(gam_model, d, df) |>
        dplyr::filter(x >= x_lim[[1]], x <= x_lim[[2]])
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


