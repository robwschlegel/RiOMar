# One-off: compass-octant occupancy trends (compute_octant_trend() /
# compute_monthly_octant_trend(), func/multi.R) for the three directional
# drivers (wind, wave, current), annual and per-calendar-month, for every
# zone. Dynamic-threshold-independent -- direction drivers don't depend on
# the plume detection threshold, unlike compute_seasonal_trend.R's plume
# properties, so there is no threshold loop here.
#
# Feeds a new Supplementary table (\label{tab:s_octant_trends}): compact
# summary combining each octant's annual occupancy trend with a roll-up of
# its 12 within-month trends (months significant out of 12, majority
# direction, range of monthly trend magnitude). Full annual and monthly
# detail is written to CSV only, not printed in the manuscript.
#
# Also writes the underlying annual per-octant day-counts
# (output/STATS/driver_octant_annual_counts.csv) and, from those, one
# diagnostic barplot PNG per driver (figures/driver_octant_trend/) -- bars =
# days per year in that octant, line = simple OLS trend fit for visual
# inspection only (the manuscript-facing trend is the AR/HAC fit above).
#
# Run from repo root: Rscript func/compute_direction_octant_trend.R
source("func/multi.R")
# util.R::load_tide_gauge() needs tide.R::.load_tide_raw() -- see
# compute_seasonal_trend.R's identical comment (load_driver() relies on
# util.R already having tide.R's QC functions in scope even for non-tide
# drivers).
source("func/tide.R")

direction_drivers <- c("wind", "wave", "current")

# Annual per-octant day-counts. compute_octant_trend() (multi.R) bins the
# same way internally but only returns the fitted trend summary -- this
# keeps the counts themselves, for the barplots below.
count_annual_octant_days <- function(degrees, date){
  octant <- compass_octant(degrees)
  labels <- levels(octant)
  daily <- tibble::tibble(year = lubridate::year(date), octant = octant) |>
    dplyr::filter(!is.na(octant))
  years <- unique(daily$year)

  daily |>
    dplyr::summarise(n_days = dplyr::n(), .by = c("year", "octant")) |>
    tidyr::complete(year = years, octant = labels, fill = list(n_days = 0))
}

# Single load_driver() pass per zone/driver, shared by the trend fit and the
# counts table below it (load_driver() re-reads the raw NetCDFs each call --
# reusing df_driver here avoids doing that twice).
octant_annual_results <- purrr::map(direction_drivers, function(driver_name){
  dir_col <- paste0(driver_name, "_dir")

  purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    message("annual / ", driver_name, " / ", meta$zone)

    df_driver <- load_driver(driver_name, meta)
    trend <- compute_octant_trend(df_driver[[dir_col]], df_driver$date) |>
      dplyr::mutate(driver = driver_name, zone = meta$zone, mouth_name = meta$mouth_name, .before = 1)
    counts <- count_annual_octant_days(df_driver[[dir_col]], df_driver$date) |>
      dplyr::mutate(driver = driver_name, zone = meta$zone, mouth_name = meta$mouth_name, .before = 1)
    list(trend = trend, counts = counts)
  })
}) |> purrr::flatten()

octant_trend_annual_detail <- purrr::map_dfr(octant_annual_results, "trend")
driver_octant_annual_counts <- purrr::map_dfr(octant_annual_results, "counts")
readr::write_csv(driver_octant_annual_counts, "output/STATS/driver_octant_annual_counts.csv")

octant_trend_monthly_detail <- purrr::map_dfr(direction_drivers, function(driver_name){
  dir_col <- paste0(driver_name, "_dir")

  purrr::pmap_dfr(zone_meta, function(...){
    meta <- tibble::tibble(...)
    message("monthly / ", driver_name, " / ", meta$zone)

    df_driver <- load_driver(driver_name, meta)
    compute_monthly_octant_trend(df_driver[[dir_col]], df_driver$date) |>
      dplyr::mutate(driver = driver_name, zone = meta$zone, mouth_name = meta$mouth_name, .before = 1)
  })
})

readr::write_csv(octant_trend_annual_detail, "output/STATS/octant_trend_annual_summary.csv")
readr::write_csv(octant_trend_monthly_detail, "output/STATS/octant_trend_monthly_summary.csv")

# Compact summary (new Supplementary table): one row per zone x driver x
# octant, combining the annual-level trend with a roll-up of the 12 monthly
# trends -- same aggregation pattern as compute_seasonal_trend.R's
# monthly_trend_compact.
monthly_rollup <- octant_trend_monthly_detail |>
  dplyr::summarise(
    n_months_fit = sum(!is.na(slope_p)),
    n_significant_months = sum(slope_p < 0.05, na.rm = TRUE),
    n_positive_months = sum(slope > 0, na.rm = TRUE),
    n_negative_months = sum(slope < 0, na.rm = TRUE),
    monthly_slope_min = if(all(is.na(slope))) NA_real_ else min(slope, na.rm = TRUE),
    monthly_slope_max = if(all(is.na(slope))) NA_real_ else max(slope, na.rm = TRUE),
    .by = c(driver, zone, mouth_name, octant)
  ) |>
  dplyr::mutate(monthly_majority_direction = dplyr::case_when(
    n_positive_months > n_negative_months ~ "+",
    n_negative_months > n_positive_months ~ "-",
    TRUE ~ "mixed"
  ))

octant_trend_compact_summary <- octant_trend_annual_detail |>
  dplyr::select(driver, zone, mouth_name, octant, n_years, mean_occurrence,
                annual_slope = slope, annual_slope_p = slope_p) |>
  dplyr::left_join(monthly_rollup, by = c("driver", "zone", "mouth_name", "octant")) |>
  dplyr::select(driver, zone, mouth_name, octant, n_years, mean_occurrence,
                annual_slope, annual_slope_p, n_months_fit, n_significant_months,
                monthly_majority_direction, monthly_slope_min, monthly_slope_max)

readr::write_csv(octant_trend_compact_summary, "output/STATS/octant_trend_compact_summary.csv")
print(octant_trend_compact_summary, n = Inf)


# Plotting --------------------------------------------------------------

# ggplot_theme()'s font sizes are tuned for single-panel, article-scale
# figures; an 8-facet-per-zone grid needs its own much smaller theme (same
# reasoning as figure.R::Figure_S_daily_flow's panel_theme).
.octant_barplot_theme <- theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
       strip.background = element_rect(fill = "grey90"),
       strip.text = element_text(size = 7),
       axis.title = element_text(size = 8, colour = "black"),
       axis.text = element_text(size = 6, colour = "black"),
       panel.grid.minor = element_blank(),
       panel.border = element_rect(fill = NA, colour = "black"),
       plot.margin = margin(t = 4, r = 6, b = 2, l = 2))

# Bars = count of days per year the zone-day-averaged direction fell in that
# octant (driver_octant_annual_counts, above). Black line = simple OLS trend
# (geom_smooth(method = "lm")), fit purely for visual inspection -- not the
# AR/HAC trend already fit by compute_octant_trend() and written to
# octant_trend_annual_summary.csv/octant_trend_compact_summary.csv.
.octant_barplot_panel <- function(df_counts, driver_name, zone_name){
  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)
  ggplot(df_counts, aes(x = year, y = n_days)) +
    geom_col(fill = disp$driver_colour) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
    facet_wrap(~octant, nrow = 2) +
    labs(x = NULL, y = "Days per year", title = zone_title(zone_name)) +
    .octant_barplot_theme
}

# One PNG per driver: 8-octant barplot of annual day-counts, one row per zone.
plot_octant_barplots <- function(driver_name, output_dir = "figures/driver_octant_trend"){
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    df_counts <- dplyr::filter(driver_octant_annual_counts, driver == driver_name, zone == meta$zone)
    .octant_barplot_panel(df_counts, driver_name, meta$zone)
  })

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 1, nrow = length(plotlist),
                                 labels = panel_labels, font.label = list(size = 14, face = "bold"))

  save_plot_as_png(full_plot, paste0("octant_barplot_", driver_name), width = 14, height = 16, path = output_dir)
}

purrr::walk(direction_drivers, plot_octant_barplots)
