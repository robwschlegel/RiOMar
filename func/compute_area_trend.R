# One-off: compute the long-term linear trend in daily plume area, using the
# same Sutton et al. (2022)-style AR(1)/HAC estimator (func/multi.R::
# driver_plume_trend()/fit_wls_hac_trend()) already used for the SPM mass
# trend (func/compute_mass_spm_trend.R) and the abstract's plume-area trend
# claim. Fills manuscript.tex Table 5's still-empty "Surface area" row
# (mean km^2, trend km^2/yr) and provides the intercept/slope drawn as the
# trend line on manuscript Figure 4.
# Run from repo root: Rscript func/compute_area_trend.R
source("func/multi.R")

area_stats <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  df <- combine_plume_driver("flow", meta)  # default metric_col = area_of_the_plume_mask_in_km2
  trend <- driver_plume_trend(df, "flow", meta$mouth_name, save_plot = FALSE) |>
    dplyr::filter(variable == "plume", weight_choice == "ar", timestep == "daily")

  # mean/SD on the de-seasoned daily series (deseason_doy(), func/multi.R),
  # matching the same day-of-year adjustment used for the trend itself, and
  # the SD convention Robert asked Table 5 to use consistently across all
  # four plume variables (2026-08-05) -- not SD of annual means.
  area_adj <- deseason_doy(df$plume_area, df$date)

  tibble::tibble(zone = meta$zone, mouth_name = meta$mouth_name,
                 mean_area_km2 = mean(area_adj, na.rm = TRUE),
                 sd_area_km2 = sd(area_adj, na.rm = TRUE),
                 n = trend$n, intercept = trend$intercept, slope = trend$slope,
                 slope_annualised = trend$slope_annualised, slope_p = trend$slope_p)
})

readr::write_csv(area_stats, "output/STATS/area_trend_summary.csv")
print(area_stats, n = Inf)
