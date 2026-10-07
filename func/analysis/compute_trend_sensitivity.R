# Sensitivity checks on the linear trends (sec:linear_trends), added for the
# third round of co-author review (D. Doxaran, 2026-10):
#   1. raw vs de-seasoned: the same AR(1)/HAC trend (fit_wls_hac_trend()) fit
#      to the raw daily series and to the de-seasoned one (deseason_doy()),
#      for the four plume properties and the five magnitude drivers.
#   2. flow-adjusted: de-seasoned plume area / SPM mass regressed on
#      de-seasoned river flow at its best lag (driver_correlation_trend_summary.csv),
#      then the trend fit to the residuals, i.e. the part of the plume trend
#      that a linear dependence on flow does not account for.
#   3. seasonal timing: day of year at which the climatological seasonal
#      cycle of plume area, river flow and wave height first rises above its
#      annual mean in autumn, and the day of its maximum (X11 seasonal signal
#      for area and flow, daily climatology for wave height).
# Writes output/STATS/trend_sensitivity_summary.csv and
# output/STATS/seasonal_onset_timing.csv.
# Run from repo root: Rscript func/analysis/compute_trend_sensitivity.R
source("func/multi.R")

plume_dir <- "output/panache/dynamic"
YR <- 365.25
flow_lags <- read_csv("output/STATS/driver_correlation_trend_summary.csv", show_col_types = FALSE) |>
  dplyr::filter(driver == "flow") |>
  dplyr::distinct(zone, lag_days)

# Trend summary row for one series, fit both raw and de-seasoned
trend_pair <- function(value, date){
  value_adj <- deseason_doy(value, date)
  dplyr::bind_rows(
    fit_wls_hac_trend("ar", value, date) |> dplyr::mutate(series = "raw"),
    fit_wls_hac_trend("ar", value_adj, date) |> dplyr::mutate(series = "deseasoned")
  ) |>
    dplyr::mutate(mean_value = mean(value_adj, na.rm = TRUE))
}

sensitivity <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)

  # Plume area and SPM mass, paired with river flow
  df_area <- combine_plume_driver("flow", meta)
  df_mass <- combine_plume_driver("flow", meta, metric_col = "mass_SPM_in_the_plume_area_in_t")
  lag_d <- flow_lags$lag_days[flow_lags$zone == meta$zone]

  flow_adjusted <- function(df, metric){
    plume_adj <- deseason_doy(df$plume_area, df$date)
    flow_adj  <- dplyr::lag(deseason_doy(df$value, df$date), lag_d)
    fit <- lm(plume_adj ~ flow_adj, na.action = na.exclude)
    keep <- !is.na(residuals(fit))
    fit_wls_hac_trend("ar", residuals(fit)[keep], df$date[keep]) |>
      dplyr::mutate(series = "flow_adjusted", metric = metric,
                    mean_value = mean(plume_adj, na.rm = TRUE), flow_beta = coef(fit)[["flow_adj"]],
                    flow_lag_days = lag_d)
  }

  # Shape and along-coast centroid (same inputs as compute_shape_alongcoast_trend.R)
  df_shape <- read_csv(paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv"), show_col_types = FALSE) |>
    dplyr::mutate(date = as.Date(date)) |>
    complete(date = seq(min(date), max(date), by = "day")) |>
    zoo::na.trim()
  df_coast <- compute_alongcoast_ts(meta$zone, meta, plume_dir)

  plume_rows <- dplyr::bind_rows(
    trend_pair(df_area$plume_area, df_area$date) |> dplyr::mutate(metric = "plume_area"),
    trend_pair(df_mass$plume_area, df_mass$date) |> dplyr::mutate(metric = "spm_mass"),
    trend_pair(df_shape$compactness, df_shape$date) |> dplyr::mutate(metric = "compactness"),
    trend_pair(df_coast$value, df_coast$date) |> dplyr::mutate(metric = "alongcoast_km"),
    flow_adjusted(df_area, "plume_area"),
    flow_adjusted(df_mass, "spm_mass")
  )

  # Drivers from the cached daily matrix (func/driver_interactions.R) rather
  # than load_driver(), whose tide branch reruns the raw tide-gauge QC.
  # Raw and de-seasoned fits use the same rows, which is all this check needs.
  driver_mat <- read_csv(paste0("output/STATS/daily_driver_matrix_", meta$zone, ".csv"), show_col_types = FALSE)
  driver_cols <- c(flow = "flow", tide = "tide_range", wind = "wind_spd", current = "current", wave = "wave_height")
  driver_rows <- purrr::imap_dfr(driver_cols, function(col, d){
    df <- dplyr::filter(driver_mat, !is.na(.data[[col]]))
    trend_pair(df[[col]], df$date) |> dplyr::mutate(metric = d)
  })

  dplyr::bind_rows(plume_rows, driver_rows) |>
    dplyr::mutate(zone = meta$zone, mouth_name = meta$mouth_name, .before = 1)
})

sensitivity_summary <- sensitivity |>
  dplyr::mutate(slope_annualised = slope * YR,
                slope_se_annualised = slope_se * YR,
                pct_per_year = 100 * slope_annualised / abs(mean_value),
                pct_over_record = pct_per_year * (as.numeric(diff(range(as.Date(c("1998-01-01", "2025-12-31"))))) / YR)) |>
  dplyr::select(zone, mouth_name, metric, series, n, mean_value, slope_annualised, slope_se_annualised,
                slope_p, pct_per_year, pct_over_record, flow_beta, flow_lag_days)

readr::write_csv(sensitivity_summary, "output/STATS/trend_sensitivity_summary.csv")
print(sensitivity_summary, n = Inf, width = Inf)


# Seasonal timing ------------------------------------------------------------
# Climatological cycle by day of year, smoothed with a 31-day circular running
# mean so single storms don't set the onset date.
onset_timing <- function(doy, value, variable, zone){
  clim <- tibble::tibble(doy = doy, value = value) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .by = "doy") |>
    tidyr::complete(doy = 1:366) |>
    dplyr::arrange(doy) |>
    dplyr::mutate(value = zoo::na.approx(value, rule = 2))
  padded <- c(tail(clim$value, 15), clim$value, head(clim$value, 15))
  smooth <- zoo::rollmean(padded, 31, align = "center")
  above <- smooth > mean(smooth)
  # first upward crossing of the annual mean after 1 July (doy 182)
  up <- which(!above[-length(above)] & above[-1]) + 1
  tibble::tibble(zone = zone, variable = variable,
                 autumn_onset_doy = up[up >= 182][1], peak_doy = which.max(smooth))
}

x11_dir <- file.path("figures/ARTICLE", get_registry_row("x11_interannual_river_flow")$output_subdir, "DATA")
x11_plume <- read_csv(file.path(x11_dir, "ts_plume_data.csv"), show_col_types = FALSE)
x11_river <- read_csv(file.path(x11_dir, "ts_river_data.csv"), show_col_types = FALSE)

timing <- purrr::map_dfr(unique(zone_meta$zone), function(z){
  p <- dplyr::filter(x11_plume, Zone == z)
  r <- dplyr::filter(x11_river, Zone == z)
  w <- read_csv(paste0("output/STATS/driver_x11_inputs/", z, "_wave.csv"), show_col_types = FALSE)
  dplyr::bind_rows(
    onset_timing(yday(p$dates), p$Seasonal_signal, "plume_area_x11_seasonal", z),
    onset_timing(yday(r$dates), r$Seasonal_signal, "river_flow_x11_seasonal", z),
    onset_timing(yday(w$date), w$value, "wave_height_climatology", z)
  )
})

readr::write_csv(timing, "output/STATS/seasonal_onset_timing.csv")
print(timing, n = Inf)
