# func/compute_driver_x11_figures.R
#
# NOTE (2026-08-11): Section 1 (the x11_interannual_river_flow-style stacked
# X11 trend-cycle comparison, one PNG per driver) moved to func/compute_driver_x11_figures.py,
# which uses func/X11.py's weekly Census-Pezzulli decomposition instead of
# the monthly-only seasonal::seas()/X-13ARIMA-SEATS this file used to call
# directly (X-13 rejects weekly-frequency input outright -- confirmed by
# testing a weekly/frequency=48 series, rejected with "Seasonal period too
# large"). This file now contains only Section 2 below (the raw daily
# wind/wave category-split plots), which are not X11 decomposition and were
# always independent of the X-13 constraint above.
#
# Run from repo root: Rscript func/compute_driver_x11_figures.R

source("func/multi.R")

OUT_DIR <- "figures/driver_x11_comparison"
CATEGORY_DIR <- file.path(OUT_DIR, "by_category")
dir.create(CATEGORY_DIR, recursive = TRUE, showWarnings = FALSE)


# 2. Category-split versions for wind and wave height -----------------------
# "The categories we have established" = the calm (<3 m/s) / onshore /
# offshore wind classification already used for the gam_partial_effects
# figure (plot_category_scatter(), func/multi.R) -- applied here to both the
# wind and the wave-height comparison, since that's the only categorical
# scheme in the pipeline that applies uniformly to all four zones (the
# Rhone-only calm/onshore-easterly/Mistral/other scheme in
# rhone_wind_wave_effect() is a one-off side-study, not the
# manuscript's own categorisation, and isn't defined for the other 3 zones).
#
# Unlike Section 1 above, this stays at *daily*, raw (not STL-smoothed)
# resolution: a wind category is inherently a per-day label (a given day is
# either calm, onshore, or offshore), and the monthly STL Interannual signal
# moves far too slowly for a day-scale category split to show anything
# against it. A thin grey line shows the full daily record for context, and
# coloured points mark just the days in that category, so it's visible
# whether the driver still tracks plume area when only e.g. calm days are
# considered.

CATEGORY_LEVELS <- c("calm (<3 m/s)", "onshore", "offshore")

plot_by_category <- function(zone_name, driver_name, category){
  meta <- get_zone_meta(zone_name = zone_name)
  df <- combine_plume_driver(driver_name, meta)  # daily, raw plume_area + value

  # Distinct column names here (not wind_spd/direction) because combine_plume_driver()
  # already returns a same-named "direction" column of its own when driver_name
  # is "wind" (load_driver()'s wind branch includes it) -- joining a second
  # "direction" column in that case silently produces direction.x/direction.y
  # instead of erroring at the join, which then breaks the drop_na() below.
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, cat_wind_spd = value, cat_wind_direction = direction)
  df <- df |> dplyr::left_join(df_wind, by = "date") |> tidyr::drop_na(cat_wind_spd, cat_wind_direction)

  df$wind_category <- dplyr::case_when(
    df$cat_wind_spd < 3            ~ "calm (<3 m/s)",
    df$cat_wind_direction == "off" ~ "offshore",
    TRUE                           ~ "onshore"
  )
  df$wind_category <- factor(df$wind_category, levels = CATEGORY_LEVELS)

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)
  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = df$value, var_ref = df$plume_area)
  df <- df |> dplyr::mutate(driver_scaled = value * scaling_factor$diff + scaling_factor$adjust)
  df_cat <- dplyr::filter(df, wind_category == category)
  unique_years <- df$date |> lubridate::year() |> unique()

  # r computed within the category's own days only (not the full grey
  # background series) -- this is the number that answers "does the
  # driver-plume relationship hold up within just this wind regime",
  # which is the whole point of the category split.
  r_value <- cor(df_cat$plume_area, df_cat$value, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  ggplot2::ggplot() +
    ggplot2::geom_line(data = df, ggplot2::aes(x = date, y = plume_area), colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_line(data = df, ggplot2::aes(x = date, y = driver_scaled), colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_point(data = df_cat, ggplot2::aes(x = date, y = plume_area), colour = "brown", size = 0.6, alpha = 0.6) +
    ggplot2::geom_point(data = df_cat, ggplot2::aes(x = date, y = driver_scaled), colour = disp$driver_colour, size = 0.6, alpha = 0.6) +
    ggplot2::annotate("text", x = min(df$date), y = Inf, label = r_label,
                      hjust = 0, vjust = 1.5, size = 6, colour = "black") +
    ggplot2::scale_x_date(name = "",
                          breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                          labels = unique_years |> stringr::str_extract_all('[0-9][0-9]$') |> unlist()) +
    ggplot2::scale_y_continuous(name = "Plume area (km²)",
                                sec.axis = ggplot2::sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                                             name = disp$driver_label)) +
    ggplot2::labs(title = paste0(zone_title(zone_name), ": ", category),
                 subtitle = paste0("n = ", nrow(df_cat), " / ", nrow(df), " days (raw daily values, not STL-smoothed)")) +
    ggplot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                  plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, colour = "grey40"),
                  axis.text.y.left = ggplot2::element_text(color = "brown"),
                  axis.title.y.left = ggplot2::element_text(color = "brown"),
                  axis.text.y.right = ggplot2::element_text(color = disp$driver_colour),
                  axis.title.y.right = ggplot2::element_text(color = disp$driver_colour),
                  panel.border = ggplot2::element_rect(linetype = "solid", fill = NA))
}

message("Building category-split comparisons for wind and wave height...")
for(driver_name in c("wind", "wave")){
  for(category in CATEGORY_LEVELS){
    zone_plots <- purrr::map(zones, plot_by_category, driver_name = driver_name, category = category)
    composite <- ggpubr::ggarrange(plotlist = zone_plots, ncol = 1, nrow = length(zones), align = "v")
    category_slug <- stringr::str_extract(category, "^[a-z]+")
    out_file <- file.path(CATEGORY_DIR, paste0(driver_name, "_", category_slug, ".png"))
    ggplot2::ggsave(out_file, composite, width = 20, height = 16, dpi = 200)
    message("Wrote ", out_file)
  }
}
