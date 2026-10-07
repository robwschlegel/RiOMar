# One-off: compute the long-term linear trend (same AR(1)/HAC estimator as
# the area and mass trends, func/multi.R::fit_wls_hac_trend()) for two plume
# metrics that don't exist anywhere else in the pipeline:
#   - shape: daily compactness (4*pi*area/perimeter^2), from
#     func/analysis/compute_plume_shape.py's PlumeShape.csv (Crofton perimeter on the
#     panache PlumeMasks.nc daily masks).
#   - along-coast extension: daily along-coast position of the SPM-weighted
#     plume centroid, relative to the river mouth. The along-coast direction
#     is not assumed from geography -- it is estimated per zone as the first
#     principal component of the centroid's own long-term scatter (in local
#     km, relative to the river mouth), since a plume's centroid is
#     physically constrained to vary mostly along the coastline it is
#     attached to. Daily metric = signed projection onto that axis (km).
# Run from repo root: Rscript func/analysis/compute_shape_alongcoast_trend.R
source("func/multi.R")

plume_dir <- "output/panache/dynamic"

results <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)

  df_shape <- read_csv(paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv"), show_col_types = FALSE) |>
    dplyr::mutate(date = as.Date(date)) |>
    complete(date = seq(min(date), max(date), by = "day")) |>
    zoo::na.trim()
  # De-seasoned before trend fitting (deseason_doy(), func/multi.R), matching
  # the same Sutton et al. (2022) day-of-year adjustment driver_plume_trend()
  # applies to area/mass
  shape_adj <- deseason_doy(df_shape$compactness, df_shape$date)
  shape_trend <- fit_wls_hac_trend("ar", shape_adj, df_shape$date) |>
    dplyr::mutate(mouth_name = meta$mouth_name, metric = "compactness",
                  mean_adj = mean(shape_adj, na.rm = TRUE), sd_adj = sd(shape_adj, na.rm = TRUE))

  df_coast <- compute_alongcoast_ts(meta$zone, meta, plume_dir)
  coast_adj <- deseason_doy(df_coast$value, df_coast$date)
  coast_trend <- fit_wls_hac_trend("ar", coast_adj, df_coast$date) |>
    dplyr::mutate(mouth_name = meta$mouth_name, metric = "alongcoast_km",
                  mean_adj = mean(coast_adj, na.rm = TRUE), sd_adj = sd(coast_adj, na.rm = TRUE))

  dplyr::bind_rows(shape_trend, coast_trend)
})

results_summary <- results |>
  dplyr::mutate(slope_annualised = slope * 365.25,
                # HAC standard error and 95% CI of the annualised slope (Table 1)
                slope_se_annualised = slope_se * 365.25,
                ci95_lower = slope_annualised - 1.96 * slope_se_annualised,
                ci95_upper = slope_annualised + 1.96 * slope_se_annualised) |>
  dplyr::select(mouth_name, metric, n, mean_adj, sd_adj, slope, slope_annualised, slope_p,
                slope_se_annualised, ci95_lower, ci95_upper)

readr::write_csv(results_summary, "output/STATS/shape_alongcoast_trend_summary.csv")
print(results_summary, n = Inf)
