# One-off: compute the long-term linear trend (same AR(1)/HAC estimator as
# the area and mass trends, func/multi.R::fit_wls_hac_trend()) for two plume
# metrics that don't exist anywhere else in the pipeline:
#   - shape: daily compactness (4*pi*area/perimeter^2), from
#     func/compute_plume_shape.py's PlumeShape.csv (Crofton perimeter on the
#     panache PlumeMasks.nc daily masks).
#   - along-coast extension: daily along-coast position of the SPM-weighted
#     plume centroid, relative to the river mouth. The along-coast direction
#     is not assumed from geography -- it is estimated per zone as the first
#     principal component of the centroid's own long-term scatter (in local
#     km, relative to the river mouth), since a plume's centroid is
#     physically constrained to vary mostly along the coastline it is
#     attached to. Daily metric = signed projection onto that axis (km).
# Run from repo root: Rscript func/compute_shape_alongcoast_trend.R
source("func/multi.R")

plume_dir <- "output/panache/dynamic"

compute_alongcoast <- function(zone, meta){
  df <- read_csv(paste0(plume_dir, "/", zone, "/Results.csv"), show_col_types = FALSE) |>
    dplyr::mutate(date = as.Date(date)) |>
    dplyr::filter(!is.na(lon_weighted_centroid_of_the_plume_area), !is.na(lat_weighted_centroid_of_the_plume_area))

  lat0 <- meta$mouth_lat; lon0 <- meta$mouth_lon
  x_km <- (df$lon_weighted_centroid_of_the_plume_area - lon0) * 111.32 * cos(pi / 180 * lat0)
  y_km <- (df$lat_weighted_centroid_of_the_plume_area - lat0) * 111.32

  pc <- prcomp(cbind(x_km, y_km))
  axis1 <- pc$rotation[, 1]
  # Fix an arbitrary PCA sign convention: orient axis1 so its dominant
  # component (whichever of east-west/north-south carries more of the
  # loading) is positive, so "positive slope" has a stable, reportable
  # geographic meaning (east or north) instead of flipping arbitrarily.
  dominant_is_x <- abs(axis1[1]) >= abs(axis1[2])
  if((dominant_is_x && axis1[1] < 0) || (!dominant_is_x && axis1[2] < 0)) axis1 <- -axis1
  direction_label <- if(dominant_is_x) "east" else "north"
  message(zone, ": along-coast axis loadings x(east)=", round(axis1[1], 3),
          " y(north)=", round(axis1[2], 3), " -> positive = more ", direction_label)

  alongcoast_km <- as.numeric(cbind(x_km, y_km) %*% axis1)

  tibble::tibble(date = df$date, alongcoast_km = alongcoast_km) |>
    complete(date = seq(min(date), max(date), by = "day")) |>
    zoo::na.trim() |>
    dplyr::mutate(zone = zone, .before = "date")
}

results <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)

  df_shape <- read_csv(paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv"), show_col_types = FALSE) |>
    dplyr::mutate(date = as.Date(date)) |>
    complete(date = seq(min(date), max(date), by = "day")) |>
    zoo::na.trim()
  # De-seasoned before trend fitting (deseason_doy(), func/multi.R), matching
  # the same Sutton et al. (2022) day-of-year adjustment driver_plume_trend()
  # applies to area/mass -- this metric used to skip that step and fit the
  # trend on the raw daily series directly, an inconsistency with the other
  # three Table 5 variables found 2026-08-05 while adding mean/SD here.
  # Mean/SD are reported on this same de-seasoned series (not the raw one),
  # per Robert's clarification that these should reflect the treated data.
  shape_adj <- deseason_doy(df_shape$compactness, df_shape$date)
  shape_trend <- fit_wls_hac_trend("ar", shape_adj, df_shape$date) |>
    dplyr::mutate(mouth_name = meta$mouth_name, metric = "compactness",
                  mean_adj = mean(shape_adj, na.rm = TRUE), sd_adj = sd(shape_adj, na.rm = TRUE))

  df_coast <- compute_alongcoast(meta$zone, meta)
  coast_adj <- deseason_doy(df_coast$alongcoast_km, df_coast$date)
  coast_trend <- fit_wls_hac_trend("ar", coast_adj, df_coast$date) |>
    dplyr::mutate(mouth_name = meta$mouth_name, metric = "alongcoast_km",
                  mean_adj = mean(coast_adj, na.rm = TRUE), sd_adj = sd(coast_adj, na.rm = TRUE))

  dplyr::bind_rows(shape_trend, coast_trend)
})

results_summary <- results |>
  dplyr::mutate(slope_annualised = slope * 365.25) |>
  dplyr::select(mouth_name, metric, n, mean_adj, sd_adj, slope, slope_annualised, slope_p)

readr::write_csv(results_summary, "output/STATS/shape_alongcoast_trend_summary.csv")
print(results_summary, n = Inf)
