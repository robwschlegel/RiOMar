# One-off: compute per-calendar-month linear trends (compute_monthly_trend(),
# func/multi.R) for all four plume properties (area, SPM mass, compactness,
# along-coast centroid drift) and five known drivers (discharge, wind speed,
# tidal range, wave height, current speed), at both the dynamic and static
# detection thresholds, for every zone. Threshold is looped rather than
# duplicated as separate code (func/multi.R::combine_plume_driver()/
# load_plume_ts() both take a plume_dir argument for this).
#
# Feeds manuscript Table 7 (compact summary: months significant out of 12,
# majority direction, range of monthly trend magnitude, dynamic threshold
# only) and the companion Supplementary table (full month-by-month detail,
# both thresholds).
#
# Run from repo root: Rscript func/compute_seasonal_trend.R
source("func/multi.R")

thresholds <- c(dynamic = "output/panache/dynamic", static = "output/panache/static")
mass_col <- "mass_SPM_in_the_plume_area_in_g_m"  # kg, see compute_mass_spm_trend.R
drivers <- c("flow", "wind", "tide", "wave", "current")

monthly_trend_detail <- purrr::map_dfr(names(thresholds), function(threshold_label){
  plume_dir <- thresholds[[threshold_label]]

  purrr::pmap_dfr(zone_meta, function(...){
    meta <- tibble::tibble(...)
    message(threshold_label, " / ", meta$zone)

    # -- plume properties ----------------------------------------------------
    df_area <- load_plume_ts(meta$zone, plume_dir = plume_dir, outlier_max = 20000)
    trend_area <- compute_monthly_trend(df_area$plume_area, df_area$date) |>
      dplyr::mutate(variable = "plume_area", category = "property")

    df_mass <- load_plume_ts(meta$zone, plume_dir = plume_dir, metric_col = mass_col, outlier_max = NULL)
    trend_mass <- compute_monthly_trend(df_mass$plume_area, df_mass$date) |>
      dplyr::mutate(variable = "SPM_mass", category = "property")

    # func/compute_plume_shape.py now covers both thresholds (2026-08-07), so
    # this file should always exist -- the guard is kept as defensive
    # robustness (skip compactness rather than failing the whole script) in
    # case that script is ever re-run for only one threshold in the future.
    shape_path <- paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv")
    if(file.exists(shape_path)){
      df_shape <- read_csv(shape_path, show_col_types = FALSE) |>
        dplyr::mutate(date = as.Date(date)) |>
        complete(date = seq(min(date), max(date), by = "day")) |>
        zoo::na.trim()
      trend_shape <- compute_monthly_trend(df_shape$compactness, df_shape$date) |>
        dplyr::mutate(variable = "compactness", category = "property")
    } else {
      message("compute_seasonal_trend: no PlumeShape.csv for ", threshold_label, "/", meta$zone,
             " -- skipping compactness for this threshold.")
      trend_shape <- NULL
    }

    df_coast <- compute_alongcoast_ts(meta$zone, meta, plume_dir)
    trend_coast <- compute_monthly_trend(df_coast$value, df_coast$date) |>
      dplyr::mutate(variable = "alongcoast_km", category = "property")

    # -- drivers ---------------------------------------------------------------
    trend_drivers <- purrr::map_dfr(drivers, function(driver_name){
      df_driver <- load_driver(driver_name, meta)
      compute_monthly_trend(df_driver$value, df_driver$date) |>
        dplyr::mutate(variable = driver_name, category = "driver")
    })

    dplyr::bind_rows(trend_area, trend_mass, trend_shape, trend_coast, trend_drivers) |>
      dplyr::mutate(threshold = threshold_label, zone = meta$zone, mouth_name = meta$mouth_name, .before = 1)
  })
})

readr::write_csv(monthly_trend_detail, "output/STATS/monthly_trend_summary.csv")

# Compact summary (manuscript Table 7): one row per zone x variable, dynamic
# threshold rows are what the main text shows; static rows are retained in
# the same file for the Supplementary comparison.
monthly_trend_compact <- monthly_trend_detail |>
  dplyr::filter(weight_choice == "ar") |>
  dplyr::summarise(
    n_months_fit = sum(!is.na(slope_p)),
    n_significant_months = sum(slope_p < 0.05, na.rm = TRUE),
    n_positive = sum(slope > 0, na.rm = TRUE),
    n_negative = sum(slope < 0, na.rm = TRUE),
    slope_min = if(all(is.na(slope))) NA_real_ else min(slope, na.rm = TRUE),
    slope_max = if(all(is.na(slope))) NA_real_ else max(slope, na.rm = TRUE),
    .by = c(threshold, zone, mouth_name, category, variable)
  ) |>
  dplyr::mutate(majority_direction = dplyr::case_when(
    n_positive > n_negative ~ "+",
    n_negative > n_positive ~ "-",
    TRUE ~ "mixed"
  )) |>
  dplyr::select(threshold, zone, mouth_name, category, variable, n_months_fit,
                n_significant_months, majority_direction, slope_min, slope_max)

readr::write_csv(monthly_trend_compact, "output/STATS/monthly_trend_compact_summary.csv")
print(dplyr::filter(monthly_trend_compact, threshold == "dynamic"), n = Inf)
