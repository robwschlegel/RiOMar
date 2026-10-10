# One-off: compute, per zone, the best-lag correlation between each driver
# and plume area, and the driver's own AR(1)/HAC-weighted linear trend
# (func/multi.R::fit_wls_hac_trend(), the same estimator used for the
# panache_stats_table's plume-area/SPM-mass rows -- func/analysis/compute_area_trend.R),
# for the driver_stats_table's River discharge / Wind / Tide / Wave height
# rows. See metadata/figure_table_registry.csv for current numbers.
#
# Best-lag search over 0-14 days (func/analysis/compute_driver_x11_correlation_table.R
# ::daily_flow_r_best_lag()/category_r_best_lag(), the
# daily_flow_lagged_correlation figure) -- applied here to all four drivers
# for consistency, plus a p-value (cor.test() at the identified best lag)
# alongside the lag itself.
# Run from repo root: Rscript func/analysis/compute_driver_correlation_trend.R
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

    # mean/SD on the de-seasoned daily series (deseason_doy(), func/multi.R),
    # matching the same convention the panache_stats_table's generators
    # already use (func/analysis/compute_area_trend.R etc.), added 2026-08-11 per metadata/TODO.md.
    value_adj <- deseason_doy(df$value, df$date)

    # Best-lag correlation between the de-seasoned driver and de-seasoned plume
    # area (2026-10-10: was the raw series, whose r mostly reflected the
    # seasonal cycle the two share; the X11 seasonal correlations already
    # cover that, so these measure day-to-day anomaly coupling instead)
    df_adj <- dplyr::mutate(df, value = value_adj, plume_area = deseason_doy(plume_area, date))
    cor_df <- driver_plume_correlation(df_adj, max_lag_daily = MAX_LAG_DAILY) |> dplyr::filter(timestep == "daily")
    peak <- cor_df |> dplyr::slice_max(cor, n = 1)
    lagged_value <- dplyr::lag(df_adj$value, peak$lag)
    peak_test <- cor.test(df_adj$plume_area, lagged_value)

    # Trend fit to the de-seasoned series, as for the plume trends and as
    # sec:linear_trends states (2026-10-07: was fit to the raw df$value, so
    # driver and plume trends were not on the same footing).
    trend <- fit_wls_hac_trend("ar", value_adj, df$date)

    tibble::tibble(zone = meta$zone, driver = driver_name, driver_label = driver_label,
                   n = nrow(df), mean_value = mean(value_adj, na.rm = TRUE), sd_value = sd(value_adj, na.rm = TRUE),
                   r = peak$cor, lag_days = peak$lag, r_p = peak_test$p.value,
                   trend_annual = trend$slope * 365.25, trend_p = trend$slope_p,
                   unit = driver_units[[driver_name]])
  })
})

readr::write_csv(results, "output/STATS/driver_correlation_trend_summary.csv")
print(results, n = Inf)
