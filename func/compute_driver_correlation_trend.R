# One-off: compute, per zone, the Pearson correlation between each driver
# and plume area, and the driver's own AR(1)/HAC-weighted linear trend
# (func/multi.R::fit_wls_hac_trend(), the same estimator used for Table 5's
# plume-area/SPM-mass rows -- func/compute_area_trend.R), for manuscript
# Table 6's River discharge / Wind / Tide / Wave height rows.
# Run from repo root: Rscript func/compute_driver_correlation_trend.R
source("func/multi.R")

drivers <- c(flow = "River discharge", wind = "Wind", tide = "Tide", wave = "Wave height")
driver_units <- c(flow = "m^3 s^-1 yr^-1", wind = "m s^-1 yr^-1", tide = "m yr^-1", wave = "m yr^-1")

results <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  purrr::imap_dfr(drivers, function(driver_label, driver_name){
    df <- combine_plume_driver(driver_name, meta)
    r <- cor(df$plume_area, df$value, use = "complete.obs")
    trend <- fit_wls_hac_trend("ar", df$value, df$date)
    tibble::tibble(zone = meta$zone, driver = driver_name, driver_label = driver_label,
                   n = nrow(df), r = r,
                   trend_annual = trend$slope * 365.25, trend_p = trend$slope_p,
                   unit = driver_units[[driver_name]])
  })
})

readr::write_csv(results, "output/STATS/driver_correlation_trend_summary.csv")
print(results, n = Inf)
