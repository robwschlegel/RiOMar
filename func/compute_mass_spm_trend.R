# One-off: compute the long-term linear trend in SPM mass per plume, using
# the same Sutton et al. (2022)-style AR(1)/HAC estimator (func/multi.R::
# driver_plume_trend()) already used for the plume-area trend reported in
# the abstract. Run from repo root: Rscript func/compute_mass_spm_trend.R
source("func/multi.R")

mass_col <- "mass_SPM_in_the_plume_area_in_g_m"  # despite the name, this column is in kg (mean_SPM [g/m3] * area [km2] * 1000)

mass_trends <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  df <- combine_plume_driver("flow", meta, metric_col = mass_col, outlier_max = NULL)
  driver_plume_trend(df, "flow", meta$mouth_name, save_plot = FALSE)
})

mass_trends_summary <- mass_trends |>
  dplyr::filter(variable == "plume", weight_choice == "ar") |>
  dplyr::mutate(slope_annualised_t_yr = slope_annualised / 1000,   # kg/yr -> t/yr
                slope_t = slope / 1000) |>                        # kg -> t, per native time unit
  dplyr::select(mouth_name, timestep, n, slope_t, slope_annualised_t_yr, slope_p)

readr::write_csv(mass_trends_summary, "output/STATS/mass_SPM_trend_summary.csv")
