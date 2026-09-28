# Daily plume area + SPM mass time series (dynamic threshold, merged
# sensor), with an AR(1)/HAC-weighted trend line, one panel per zone.
# Manuscript slot "plume_area_timeseries" -- see
# metadata/figure_table_registry.csv for its current figure number.
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
                     right = text_grob(expression(SPM~mass~(t~x~10^{5})), rot = -90, size = 30, color = "steelblue4")),
                   registry_basename(output_subdir), width = 20, height = 16, path = main_folder)
}


# Monthly heatmap (sec:results_seasonal) of all four plume properties and
# five drivers, per zone, dynamic threshold. Each tile is that month's median
# value divided by the zone's own all-time median (dynamic threshold only) --
# i.e. how a typical day in that month compares to a typical day across the
# full record for that zone x variable -- so zones of very different raw
# magnitude (see the panache_stats_table slot) are comparable in one figure
# on a scale centred on 1 (= typical). Real (unscaled) interquartile values
# are annotated as text. Also writes the long-format data (both thresholds)
# behind it to DATA/monthly_boxplot_data.csv. Manuscript slot
# "seasonal_boxplot_heatmap" -- see metadata/figure_table_registry.csv for
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
  # Plotmath source strings (parsed via facet_wrap(labeller = label_parsed)
  # below), not literal Unicode superscripts -- U+207B/U+00B3 intermittently
  # render as missing-glyph boxes in this figure's strip text depending on
  # the R session's font/graphics-device state (reproduced 2026-09-24: same
  # code, same machine, correct in one render and broken in the next two),
  # so plotmath avoids depending on that glyph being available at all.
  variable_display <- c(
    plume_area    = '"Plume area ("*km^2*")"',
    SPM_mass      = '"SPM mass (t)"',
    compactness   = '"Compactness"',
    alongcoast_km = '"Along-coast drift (km)"',
    flow          = '"River flow ("*m^3~s^{-1}*")"',
    wind          = '"Wind speed ("*m~s^{-1}*")"',
    tide          = '"Tidal range (m)"',
    wave          = '"Wave height (m)"',
    current       = '"Current speed ("*m~s^{-1}*")"'
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

      # func/analysis/compute_plume_shape.py must be run (now wired into
      # code/4_time_series.py) before this figure -- compactness is a
      # required panel, not an optional one, so a missing file is a hard
      # error rather than a silently dropped panel.
      shape_path <- paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv")
      if(!file.exists(shape_path)){
        stop("plot_seasonal_boxplot_heatmap: missing ", shape_path,
             " -- run func/analysis/compute_plume_shape.py before regenerating this figure.")
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
  # (func/analysis/generate_monthly_trend_pct_heatmap.R) for anyone who needs it.
  heat_stats <- long_data |>
    dplyr::filter(threshold == "dynamic") |>
    dplyr::summarise(month_median = stats::median(value, na.rm = TRUE), .by = c(zone, variable, month)) |>
    dplyr::left_join(overall_median, by = c("zone", "variable")) |>
    dplyr::mutate(ratio = month_median / overall_median,
                  zone = factor(zone, levels = rev(zones), labels = zone_labels),
                  month = factor(month, levels = 1:12, labels = month.abb),
                  variable = factor(variable, levels = names(variable_display), labels = unname(variable_display)))

  # Diverging scale centred on 1 (= that month matches the zone's own
  # all-time typical day). Co-authors flagged the original purple/orange
  # scale as still reading too close to the blue/red diverging scale used by
  # the monthly_trend_pct_heatmap slot's %-change-per-year figure
  # (func/analysis/generate_monthly_trend_pct_heatmap.R) -- both are a cool colour
  # against a warm colour, so the two heatmaps were visually conflated at a
  # glance even though the hues differ. Switched to ColorBrewer's PRGn
  # (purple-green), a colourblind-safe diverging palette (verified
  # deuteranopia/protanopia/tritanopia-safe at colorbrewer2.org) that reads
  # as a genuinely different colour axis -- purple vs. green, not
  # warm vs. cool -- rather than another warm/cool pair.
  p_heatmap <- ggplot(heat_stats, aes(x = month, y = zone, fill = ratio)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_wrap(~variable, ncol = 3, labeller = label_parsed) +
    scale_fill_gradient2(name = "Month / overall\nmedian", low = "#762A83", mid = "white", high = "#1B7837", midpoint = 1) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
         axis.text.y = element_text(size = 10), panel.grid = element_blank())

  save_plot_as_png(p_heatmap, registry_basename(figure_5_output_subdir), width = 12, height = 8, path = figure_5_dir)
  message("Wrote ", registry_filename(figure_5_output_subdir))
  invisible(TRUE)
}


# X11 interannual (long-term) signal of plume area vs. river flow, dynamic
# threshold (main results), all four zones. Per-panel axis titles are
# suppressed (show_axis_titles = FALSE) in favour of one shared left/right
# label on the assembled composite, matching the plume_methodology_panel
# figure's convention (annotate_figure(), not a title repeated on all four
# panels). Manuscript slot "x11_interannual_river_flow" -- see
# metadata/figure_table_registry.csv for its current figure number.
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


