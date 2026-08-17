# One-off: compute, per zone, the best-lag correlation between each driver
# and plume area, and the driver's own AR(1)/HAC-weighted linear trend
# (func/multi.R::fit_wls_hac_trend(), the same estimator used for Table 5's
# plume-area/SPM-mass rows -- func/compute_area_trend.R), for manuscript
# Table 6's River discharge / Wind / Tide / Wave height rows.
#
# Best-lag search over 0-14 days (func/compute_driver_x11_correlation_table.R
# ::daily_flow_r_best_lag()/category_r_best_lag(), manuscript Figure 5) --
# applied here to all four drivers for consistency, plus a p-value
# (cor.test() at the identified best lag) alongside the lag itself.
# Run from repo root: Rscript func/compute_driver_correlation_trend.R
source("func/multi.R")
# util.R::load_tide_gauge() needs tide.R::.load_tide_raw() -- see func/figure.R's identical comment.
source("func/tide.R")

drivers <- c(flow = "River discharge", wind = "Wind", tide = "Tide", wave = "Wave height", current = "Current speed")
driver_units <- c(flow = "m^3 s^-1 yr^-1", wind = "m s^-1 yr^-1", tide = "m yr^-1", wave = "m yr^-1", current = "m s^-1 yr^-1")
MAX_LAG_DAILY <- 14

results <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  purrr::imap_dfr(drivers, function(driver_label, driver_name){
    df <- combine_plume_driver(driver_name, meta)

    cor_df <- driver_plume_correlation(df, max_lag_daily = MAX_LAG_DAILY) |> dplyr::filter(timestep == "daily")
    peak <- cor_df |> dplyr::slice_max(cor, n = 1)
    lagged_value <- dplyr::lag(df$value, peak$lag)
    peak_test <- cor.test(df$plume_area, lagged_value)

    trend <- fit_wls_hac_trend("ar", df$value, df$date)

    # mean/SD on the de-seasoned daily series (deseason_doy(), func/multi.R),
    # matching the same convention Table 5's generators already use
    # (func/compute_area_trend.R etc.), added 2026-08-11 per manuscript/TODO.md.
    value_adj <- deseason_doy(df$value, df$date)

    tibble::tibble(zone = meta$zone, driver = driver_name, driver_label = driver_label,
                   n = nrow(df), mean_value = mean(value_adj, na.rm = TRUE), sd_value = sd(value_adj, na.rm = TRUE),
                   r = peak$cor, lag_days = peak$lag, r_p = peak_test$p.value,
                   trend_annual = trend$slope * 365.25, trend_p = trend$slope_p,
                   unit = driver_units[[driver_name]])
  })
})

readr::write_csv(results, "output/STATS/driver_correlation_trend_summary.csv")
print(results, n = Inf)
