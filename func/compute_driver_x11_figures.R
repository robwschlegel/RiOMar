# func/compute_driver_x11_figures.R
#
# Manuscript Figure 6 (func/figure.R::Figure_6_x11_interannual(), driven by
# the bespoke Python Census-X11 port in func/X11.py) shows the long-term
# ("Interannual") signal of plume area against river discharge, one panel
# per zone, stacked top-to-bottom. Robert asked for the same figure for the
# other drivers (wind, tide, wave height, current speed), plus a version of
# the wind and wave-height comparisons split by the on-/off-shore/calm wind
# categories already established for manuscript Figure 8
# (func/multi.R::plot_category_scatter()).
#
# This reuses the R-side driver-comparison machinery already built in
# func/multi.R (combine_plume_driver()), rather than extending the Python
# Census-X11 pipeline itself: that pipeline
# (func/X11.py::Apply_X11_method_on_time_series()) only has a path for
# plume area + river flow, and porting it to four more drivers would be a
# large undertaking for what's meant to be an explanatory figure set.
#
# Section 1 originally used STL (util.R::stl_single()) in place of a true
# X11 decomposition. Robert asked for real X11 instead, and pointed at
# something already in the codebase for this -- that's
# `library(seasonal)` (func/multi.R, loaded but, per its own comment there,
# "currently not used"). The `seasonal` package is an R wrapper around the
# U.S. Census Bureau's official X-13ARIMA-SEATS binary (via `x13binary`,
# confirmed installed and working with seasonal::checkX13()), which can be
# forced into its X-11 method specifically with seas(x, x11 = list())
# rather than the default SEATS method -- a genuine Census-X11 decomposition,
# not a STL stand-in, and the same family of method (Cleveland 1976) the
# manuscript cites for Figure 6.
#
# One constraint this ran into: X-13 only accepts monthly (frequency = 12)
# or quarterly (frequency = 4) input -- confirmed by testing a weekly
# (frequency = 48) series, which X-13 rejects with "Seasonal period too
# large" -- so the weekly binning Figure 6's own Python pipeline uses isn't
# reproducible via this route. Monthly aggregation (driver_plume_timesteps(),
# already built for this) is used instead, which is also coarser than
# Figure 6's weekly bins but is what let the STL version's smoothness match
# Figure 6's own look reasonably well (checked visually) before this switch.
# X-13 also rejects series with internal NAs, hence the na.approx() below
# (the same gap-filling stl_single() already does for STL).
#
# Run from repo root: Rscript func/compute_driver_x11_figures.R

source("func/multi.R")

OUT_DIR <- "figures/driver_x11_comparison"
CATEGORY_DIR <- file.path(OUT_DIR, "by_category")
dir.create(CATEGORY_DIR, recursive = TRUE, showWarnings = FALSE)

# Monthly-mean plume/driver series with a monthly-resolution Census-X11
# trend-cycle component (the "Interannual" signal, in the manuscript's
# terminology for Figure 6) -- see the file header for the X-13 constraints
# behind the monthly aggregation and na.approx() calls here.
x11_dual_axis_df <- function(driver_name, meta){
  df_daily <- combine_plume_driver(driver_name, meta)
  df_month <- driver_plume_timesteps(df_daily)$monthly
  start_year  <- lubridate::year(min(df_month$date))
  start_month <- lubridate::month(min(df_month$date))

  ts_plume  <- ts(zoo::na.approx(df_month$plume_area), frequency = 12, start = c(start_year, start_month))
  ts_driver <- ts(zoo::na.approx(df_month$value),      frequency = 12, start = c(start_year, start_month))

  # arima.model/regression.aictest/outlier fix X-13's internal pre-adjustment
  # to the standard "airline" model instead of letting it auto-identify one:
  # automatic identification failed to converge for at least one series
  # (current speed, Gulf of Lion) with "Estimation failed to converge during
  # the automatic model identification procedure" -- the fixed airline model
  # is the standard fallback for exactly this and converged on every series
  # tested here.
  plume_seas  <- seasonal::seas(ts_plume,  x11 = list(), arima.model = "(0 1 1)(0 1 1)",
                                regression.aictest = NULL, outlier = NULL)
  driver_seas <- seasonal::seas(ts_driver, x11 = list(), arima.model = "(0 1 1)(0 1 1)",
                                regression.aictest = NULL, outlier = NULL)
  df_month$plume_trend  <- as.numeric(seasonal::trend(plume_seas))
  df_month$driver_trend <- as.numeric(seasonal::trend(driver_seas))
  df_month
}

# Dual-axis X11 trend-cycle plot with an r annotation -- the same visual
# convention as Figure 6 itself (figure.R::plot_x11_river_and_plume()),
# generalised to any driver via driver_display, and adding the r-value
# label that function has but plot_driver_plume_dual_axis() (used for the
# STL version previously) does not.
plot_x11_dual_axis <- function(df, driver_name, zone_name){
  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)
  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = df$driver_trend, var_ref = df$plume_trend)
  df <- df |> dplyr::mutate(driver_scaled = driver_trend * scaling_factor$diff + scaling_factor$adjust)
  unique_years <- df$date |> lubridate::year() |> unique()

  r_value <- cor(df$plume_trend, df$driver_trend, use = "complete.obs")
  r_label <- paste0("r = ", sprintf("%.2f", r_value))

  ggplot2::ggplot(data = df) +
    ggplot2::geom_point(ggplot2::aes(x = date, y = plume_trend), colour = "brown") +
    ggplot2::geom_path(ggplot2::aes(x = date, y = plume_trend), colour = "brown") +
    ggplot2::geom_point(ggplot2::aes(x = date, y = driver_scaled), colour = disp$driver_colour) +
    ggplot2::geom_path(ggplot2::aes(x = date, y = driver_scaled), colour = disp$driver_colour) +
    ggplot2::annotate("text", x = min(df$date), y = Inf, label = r_label,
                      hjust = 0, vjust = 1.5, size = 6, colour = "black") +
    ggplot2::scale_x_date(name = "",
                          breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                          labels = unique_years |> stringr::str_extract_all('[0-9][0-9]$') |> unlist()) +
    ggplot2::scale_y_continuous(name = "Plume area (km²)",
                                sec.axis = ggplot2::sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                                             name = disp$driver_label)) +
    ggplot2::labs(title = zone_name) +
    ggplot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y.left = ggplot2::element_text(color = "brown"),
                  axis.title.y.left = ggplot2::element_text(color = "brown"),
                  axis.text.y.right = ggplot2::element_text(color = disp$driver_colour),
                  axis.title.y.right = ggplot2::element_text(color = disp$driver_colour),
                  panel.border = ggplot2::element_rect(linetype = "solid", fill = NA))
}


# 1. Figure 6-style stacked comparison, one new figure per driver -----------
# River flow is Figure 6 itself, so it's not repeated here. ROFI is excluded
# per the standing project rule that it's downstream of discharge, not an
# independent driver (RiOMar memory feedback_rofi_not_a_driver.md).

other_drivers <- c("wind", "tide", "wave", "current")

message("Building Figure-6-style stacked comparison for each driver (Census X11)...")
for(driver_name in other_drivers){
  zone_plots <- purrr::map(zones, function(zone_name){
    meta <- get_zone_meta(zone_name = zone_name)
    df <- x11_dual_axis_df(driver_name, meta)
    plot_x11_dual_axis(df, driver_name, zone_title(zone_name))
  })
  composite <- ggpubr::ggarrange(plotlist = zone_plots, ncol = 1, nrow = length(zones), align = "v")
  out_file <- file.path(OUT_DIR, paste0("Figure6_style_", driver_name, ".png"))
  ggplot2::ggsave(out_file, composite, width = 20, height = 16, dpi = 200)
  message("Wrote ", out_file)
}


# 2. Category-split versions for wind and wave height -----------------------
# "The categories we have established" = the calm (<3 m/s) / onshore /
# offshore wind classification already used for manuscript Figure 8
# (plot_category_scatter(), func/multi.R) -- applied here to both the
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
