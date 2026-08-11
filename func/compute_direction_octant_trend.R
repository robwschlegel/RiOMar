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
# Run from repo root: Rscript func/compute_direction_octant_trend.R
source("func/multi.R")
# util.R::load_tide_gauge() needs tide.R::.load_tide_raw() -- see
# compute_seasonal_trend.R's identical comment (load_driver() relies on
# util.R already having tide.R's QC functions in scope even for non-tide
# drivers).
source("func/tide.R")

direction_drivers <- c("wind", "wave", "current")

octant_trend_annual_detail <- purrr::map_dfr(direction_drivers, function(driver_name){
  dir_col <- paste0(driver_name, "_dir")

  purrr::pmap_dfr(zone_meta, function(...){
    meta <- tibble::tibble(...)
    message("annual / ", driver_name, " / ", meta$zone)

    df_driver <- load_driver(driver_name, meta)
    compute_octant_trend(df_driver[[dir_col]], df_driver$date) |>
      dplyr::mutate(driver = driver_name, zone = meta$zone, mouth_name = meta$mouth_name, .before = 1)
  })
})

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
