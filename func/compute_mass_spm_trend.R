# One-off: compute the long-term linear trend in SPM mass per plume, using
# the same Sutton et al. (2022)-style AR(1)/HAC estimator (func/multi.R::
# deseason_doy()) already used for the plume-area trend reported in
# the abstract. Run from repo root: Rscript func/compute_mass_spm_trend.R
source("func/multi.R")

mass_col <- "mass_SPM_in_the_plume_area_in_t"  # tonnes; sum_SPM [g/m3] * pixel_area [m2] * 1 m assumed depth, see plume_algorithm.py

mass_trends <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  df <- combine_plume_driver("flow", meta, metric_col = mass_col, outlier_max = NULL)

  # mean/SD on the de-seasoned daily series (deseason_doy(), func/multi.R),
  # same convention as compute_area_trend.R/compute_shape_alongcoast_trend.R
  mass_adj <- deseason_doy(df$plume_area, df$date)
  trend <- driver_plume_trend(df, "flow", meta$mouth_name, save_plot = FALSE) |>
    dplyr::mutate(zone = meta$zone,
                  mean_mass_t = mean(mass_adj, na.rm = TRUE), sd_mass_t = sd(mass_adj, na.rm = TRUE))
  trend
})

mass_trends_summary <- mass_trends |>
  dplyr::filter(variable == "plume", weight_choice == "ar") |>
  dplyr::mutate(slope_annualised_t_yr = slope_annualised,   # already tonnes
                slope_t = slope,                            # already tonnes, per native time unit
                intercept_t = intercept) |>                 # already tonnes, for the plume_area_timeseries figure's geom_abline trend line
  dplyr::select(zone, mouth_name, timestep, n, mean_mass_t, sd_mass_t, intercept_t, slope_t, slope_annualised_t_yr, slope_p)

readr::write_csv(mass_trends_summary, "output/STATS/mass_SPM_trend_summary.csv")
