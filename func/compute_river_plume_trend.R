# One-off: per-river discharge-vs-own-plume-area linear trend, now that
# panache v5.0.0+ tracks individual river mouths (Results.csv's `river`
# column) rather than only a zone-level union. Reuses multi.R::
# driver_plume_trend() unchanged -- same WLS-HAC trend-fitting and plotting
# logic as the zone-level trends_plume_flow_adj_<ZONE>.png figures
# (run_driver_suite("flow")), just pointed at each river's own plume mask
# and its own discharge gauge(s) instead of the zone's 'ALL' union mask and
# primary gauge. River <-> discharge-gauge mapping (including the Gironde =
# Garonne + Dordogne and Sevre = Sevre Niortaise + Lay basin groupings)
# lives in metadata/river_discharge_mapping.csv -- same mapping
# compute_river_plume_correlation.R uses.
# Run from repo root: Rscript func/compute_river_plume_trend.R
source("func/multi.R")

river_map <- readr::read_csv("metadata/river_discharge_mapping.csv", show_col_types = FALSE)

# One river's flow-vs-plume-area trend figure (dynamic threshold only --
# matches the zone-level trend figures, which are also dynamic-only).
river_plume_trend <- function(zone, panache_river, discharge_slugs){
  river_slugs <- strsplit(discharge_slugs, "\\|")[[1]]
  # Single underscore so the mouth_name/stats identifier reads as one clean
  # phrase, same convention as X11.py::
  # Apply_X11_method_on_time_series_per_river()'s river_id. Display label
  # (plot title + sanitised filename) uses the zone's pretty name instead.
  river_id <- paste0(zone, "_", stringr::str_replace_all(panache_river, " ", "_"))
  plot_label <- paste0(zone_title(zone), " ", panache_river)

  df_plume <- load_plume_ts(zone, plume_dir = "output/panache/dynamic",
                            metric_col = "area_of_the_plume_mask_in_km2", river = panache_river)
  df_flow <- load_river_flow_single(zone, river_slugs)

  df <- df_plume |>
    dplyr::select(date, plume_area) |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::rename(value = flow) |>
    zoo::na.trim()

  driver_plume_trend(df, "flow", river_id, plot_label = plot_label)
}

trend_stats <- purrr::pmap_dfr(river_map, function(zone, panache_river, discharge_slugs){
  river_plume_trend(zone, panache_river, discharge_slugs)
})

readr::write_csv(trend_stats, "output/STATS/river_plume_trend.csv")
