# func/compute_driver_x11_correlation_table.R
#
# NOTE (2026-08-11): the X11 trend-cycle correlation rows (wind/tide/wave/
# current vs plume area's weekly Census-X11 Interannual signal) moved to
# func/compute_driver_x11_correlation_table.py, which uses func/X11.py's
# weekly Census-Pezzulli decomposition instead of the monthly-only
# seasonal::seas()/X-13ARIMA-SEATS this file used to call directly (X-13
# rejects weekly-frequency input outright). This file now computes only the
# category correlations and the river-flow reference column below -- none of
# which are X11 decomposition -- and writes them to
# output/STATS/driver_x11_correlation_summary_category_rows.csv, which the
# Python script reads and merges with its own X11 rows into the final
# output/STATS/driver_x11_correlation_summary.csv. Called as an Rscript
# subprocess from compute_driver_x11_correlation_table.py; not meant to be
# run standalone for the final table anymore (its own CSV output is an
# intermediate file, not the finished one).
#
# Per Robert's original follow-up: the wave-height rows and the river-flow
# reference column use whichever 0-14 day lag maximises the correlation
# (the same lag search as manuscript Figure 5 / func/compute_wave_category_lagged_correlation.R),
# rather than a same-day (lag 0) correlation. Wind's own category rows are
# left at lag 0 -- only wave height and river flow were asked to use best-lag.
#
# Run from repo root: Rscript func/compute_driver_x11_correlation_table.R

source("func/multi.R")

CATEGORY_LEVELS <- c("calm (<3 m/s)", "onshore", "offshore")
MAX_LAG_DAILY <- 14

# Same wind-category derivation as compute_driver_x11_figures.R::plot_by_category().
add_wind_category <- function(df, meta){
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, cat_wind_spd = value, cat_wind_direction = direction)
  df <- df |> dplyr::left_join(df_wind, by = "date") |> tidyr::drop_na(cat_wind_spd, cat_wind_direction)
  df$wind_category <- dplyr::case_when(
    df$cat_wind_spd < 3            ~ "calm (<3 m/s)",
    df$cat_wind_direction == "off" ~ "offshore",
    TRUE                           ~ "onshore"
  )
  df
}

# Same-day (lag 0) category correlation -- used for wind's own category rows.
category_r_same_day <- function(driver_name, meta, category){
  df <- combine_plume_driver(driver_name, meta) |> add_wind_category(meta)
  df_cat <- dplyr::filter(df, wind_category == category)
  list(r = cor(df_cat$plume_area, df_cat$value, use = "complete.obs"), lag_days = 0L)
}

# Best-lag (0-14 days) category correlation -- used for wave's category rows.
# See func/compute_wave_category_lagged_correlation.R for why the lagged
# driver value must come from the full continuous series (not the
# category-filtered subset) for "N days" to mean real calendar days.
category_r_best_lag <- function(driver_name, meta, category, max_lag = MAX_LAG_DAILY){
  df <- combine_plume_driver(driver_name, meta) |> add_wind_category(meta)
  in_category <- df$wind_category == category
  cor_df <- purrr::map_dfr(0:max_lag, function(l){
    lagged_value <- dplyr::lag(df$value, l)
    tibble::tibble(lag = l, cor = cor(df$plume_area[in_category], lagged_value[in_category], use = "pairwise.complete.obs"))
  })
  peak <- cor_df |> dplyr::slice_max(cor, n = 1)
  list(r = peak$cor, lag_days = peak$lag)
}

# River flow reference column: the same whole-series best-lag search
# Supplementary Fig. daily_flow itself reports (func/figure.R::Figure_S_daily_flow(),
# formerly manuscript Figure 5 before the seasonal-analysis section
# repurposed that slot; driver_plume_correlation() + slice_max()), not
# restricted to any category.
daily_flow_r_best_lag <- function(meta, max_lag = MAX_LAG_DAILY){
  df <- combine_plume_driver("flow", meta)
  cor_df <- driver_plume_correlation(df, max_lag_daily = max_lag) |> dplyr::filter(timestep == "daily")
  peak <- cor_df |> dplyr::slice_max(cor, n = 1)
  list(r = peak$cor, lag_days = peak$lag)
}

message("Computing category + river-flow-reference rows...")
results <- purrr::map_dfr(zones, function(zone_name){
  meta <- get_zone_meta(zone_name = zone_name)
  flow <- daily_flow_r_best_lag(meta)

  wind_cat_rows <- purrr::map_dfr(CATEGORY_LEVELS, function(cat){
    res <- category_r_same_day("wind", meta, cat)
    tibble::tibble(zone = zone_name, driver = "wind", subset = paste0(cat, " (raw daily)"),
                   r = res$r, lag_days = res$lag_days,
                   river_flow_daily_r = flow$r, river_flow_lag_days = flow$lag_days)
  })

  wave_cat_rows <- purrr::map_dfr(CATEGORY_LEVELS, function(cat){
    res <- category_r_best_lag("wave", meta, cat)
    tibble::tibble(zone = zone_name, driver = "wave", subset = paste0(cat, " (raw daily, best lag)"),
                   r = res$r, lag_days = res$lag_days,
                   river_flow_daily_r = flow$r, river_flow_lag_days = flow$lag_days)
  })

  dplyr::bind_rows(wind_cat_rows, wave_cat_rows)
})

results <- results |>
  dplyr::mutate(zone = factor(zone, levels = ZONE_ORDER),
               r = round(r, 2), river_flow_daily_r = round(river_flow_daily_r, 2)) |>
  dplyr::arrange(zone, driver, subset)

readr::write_csv(results, "output/STATS/driver_x11_correlation_summary_category_rows.csv")
message("Wrote output/STATS/driver_x11_correlation_summary_category_rows.csv")
print(results, n = Inf)
