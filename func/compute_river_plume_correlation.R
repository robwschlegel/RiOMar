# One-off: per-river discharge-vs-own-plume-area lagged correlation, now
# that panache v5.0.0+ tracks individual river mouths (Results.csv's `river`
# column) rather than only a zone-level union. Same daily/monthly/annual
# lagged-Pearson approach as func/ROFI.R::rofi_plume_lagged_correlation()
# and func/compute_driver_correlation_trend.R, applied per river instead of
# per zone-aggregate driver. River <-> discharge-gauge mapping (including the
# Gironde = Garonne + Dordogne and Sevre = Sevre Niortaise + Lay basin
# groupings) lives in metadata/river_discharge_mapping.csv.
# Run from repo root: Rscript func/compute_river_plume_correlation.R
source("func/multi.R")

river_map <- readr::read_csv("metadata/river_discharge_mapping.csv", show_col_types = FALSE)

MAX_LAG_DAILY <- 30

# One river's lagged correlation (daily/monthly/annual) for one threshold.
river_plume_correlation <- function(zone, panache_river, discharge_slugs, plume_dir){
  river_slugs <- strsplit(discharge_slugs, "\\|")[[1]]

  df_plume <- load_plume_ts(zone, plume_dir, "area_of_the_plume_mask_in_km2", river = panache_river)
  df_flow <- load_river_flow_single(zone, river_slugs)

  df <- df_plume |>
    dplyr::select(date, plume_area) |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::rename(value = flow) |>
    zoo::na.trim()

  driver_plume_correlation(df, max_lag_daily = MAX_LAG_DAILY) |>
    dplyr::mutate(zone = zone, river = panache_river, .before = "lag")
}

cor_stats <- purrr::pmap_dfr(river_map, function(zone, panache_river, discharge_slugs){
  purrr::map_dfr(c("dynamic", "static"), function(threshold){
    river_plume_correlation(zone, panache_river, discharge_slugs,
                            plume_dir = paste0("output/panache/", threshold)) |>
      dplyr::mutate(threshold = threshold, .after = "river")
  })
})

readr::write_csv(cor_stats, "output/STATS/river_plume_lagged_correlation.csv")

# Compact summary (best lag + correlation per river x zone x threshold x
# timestep) for the manuscript supplement -- same distillation pattern as
# func/driver_interactions.R::summarise_monthly_driver_dominance().
summary_stats <- cor_stats |>
  dplyr::slice_max(cor, n = 1, by = c(zone, river, threshold, timestep)) |>
  dplyr::select(zone, river, threshold, timestep, lag, cor) |>
  dplyr::arrange(zone, river, threshold, timestep)

readr::write_csv(summary_stats, "output/STATS/river_plume_lagged_correlation_summary.csv")
print(summary_stats, n = Inf)
