# func/compute_driver_x11_correlation_table.R
#
# Summarises every r value produced by func/compute_driver_x11_figures.R
# (the X11 monthly trend-cycle comparisons for wind/tide/wave/current, and
# the calm/onshore/offshore raw-daily category comparisons for wind/wave)
# into one table, alongside river discharge's own daily correlation with
# plume area as a reference column, so the two can be compared directly.
#
# Per Robert's follow-up: the wave-height rows and the river-flow reference
# column use whichever 0-14 day lag maximises the correlation (the same
# lag search as manuscript Figure 5 / func/compute_wave_category_lagged_correlation.R),
# rather than a same-day (lag 0) correlation. This only applies to daily-
# resolution numbers -- the X11 rows are monthly trend-cycle correlations,
# where a day-scale lag search doesn't have a coherent meaning, so those are
# left as same-period correlations for every driver including wave. Wind's
# own category rows are also left at lag 0 -- only wave height and river
# flow were asked to use best-lag.
#
# Deliberately recomputes rather than parsing the annotated PNGs: the r
# values are cheap to recompute (each X-13 call takes ~1-2s) and this keeps
# the table script independent of the figure script's plotting internals.
#
# Run from repo root: Rscript func/compute_driver_x11_correlation_table.R

source("func/multi.R")

other_drivers <- c("wind", "tide", "wave", "current")
CATEGORY_LEVELS <- c("calm (<3 m/s)", "onshore", "offshore")
MAX_LAG_DAILY <- 14

# Same X-13 spec as compute_driver_x11_figures.R (fixed airline model --
# automatic identification failed to converge for at least one series, see
# that script's header for details).
x11_trend_r <- function(driver_name, meta){
  df_daily <- combine_plume_driver(driver_name, meta)
  df_month <- driver_plume_timesteps(df_daily)$monthly
  start_year  <- lubridate::year(min(df_month$date))
  start_month <- lubridate::month(min(df_month$date))

  ts_plume  <- ts(zoo::na.approx(df_month$plume_area), frequency = 12, start = c(start_year, start_month))
  ts_driver <- ts(zoo::na.approx(df_month$value),      frequency = 12, start = c(start_year, start_month))

  plume_seas  <- seasonal::seas(ts_plume,  x11 = list(), arima.model = "(0 1 1)(0 1 1)",
                                regression.aictest = NULL, outlier = NULL)
  driver_seas <- seasonal::seas(ts_driver, x11 = list(), arima.model = "(0 1 1)(0 1 1)",
                                regression.aictest = NULL, outlier = NULL)

  list(r = cor(as.numeric(seasonal::trend(plume_seas)), as.numeric(seasonal::trend(driver_seas)), use = "complete.obs"),
       lag_days = NA_integer_)
}

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

message("Computing correlation table (this calls X-13 32 times, ~1-2 min)...")
results <- purrr::map_dfr(zones, function(zone_name){
  meta <- get_zone_meta(zone_name = zone_name)
  flow <- daily_flow_r_best_lag(meta)

  x11_rows <- purrr::map_dfr(other_drivers, function(d){
    res <- x11_trend_r(d, meta)
    tibble::tibble(zone = zone_name, driver = d, subset = "all (X11 monthly trend)",
                   r = res$r, lag_days = res$lag_days,
                   river_flow_daily_r = flow$r, river_flow_lag_days = flow$lag_days)
  })

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

  dplyr::bind_rows(x11_rows, wind_cat_rows, wave_cat_rows)
})

results <- results |>
  dplyr::mutate(zone = factor(zone, levels = ZONE_ORDER),
               r = round(r, 2), river_flow_daily_r = round(river_flow_daily_r, 2)) |>
  dplyr::arrange(zone, driver, subset)

readr::write_csv(results, "output/STATS/driver_x11_correlation_summary.csv")
message("Wrote output/STATS/driver_x11_correlation_summary.csv")
print(results, n = Inf)
