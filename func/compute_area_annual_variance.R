# One-off: compute manuscript Table 4's "Variance" row (previously an
# unfilled \tocheck{} placeholder, per Robert's exact specification -- mean
# plume area per year per zone, then the standard deviation of that annual
# time series). Uses the same daily plume-area series as
# func/compute_area_trend.R (combine_plume_driver("flow", meta), default
# metric_col = area_of_the_plume_mask_in_km2, so the annual means are
# consistent with Table 4's existing "Surface area mean (km^2)" row and its
# ceiling-screened outlier handling (util.R::load_plume_ts()).
#
# Run from repo root: Rscript func/compute_area_annual_variance.R
source("func/multi.R")

annual_variance <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  df <- combine_plume_driver("flow", meta)  # default metric_col = area_of_the_plume_mask_in_km2

  annual_means <- df |>
    dplyr::mutate(year = lubridate::year(date)) |>
    dplyr::summarise(annual_mean_area_km2 = mean(plume_area, na.rm = TRUE), .by = "year")

  tibble::tibble(zone = meta$zone, mouth_name = meta$mouth_name,
                 n_years = nrow(annual_means),
                 sd_annual_mean_area_km2 = stats::sd(annual_means$annual_mean_area_km2, na.rm = TRUE))
})

readr::write_csv(annual_variance, "output/STATS/area_annual_variance_summary.csv")
print(annual_variance, n = Inf)
