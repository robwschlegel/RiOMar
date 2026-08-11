# One-off: spatial-consistency diagnostic for wind/current/wave direction and
# magnitude within each zone's bbox. Addresses the manuscript/TODO.md Section
# 2.8 deferral -- func/multi.R::load_driver() collapses every driver to a
# single zone-day scalar (spatial mean) before any consumer sees it, so no
# consumer downstream of it can check whether pixels within a zone actually
# agree with each other on a given day. This script re-reads the same raw
# files via util.R::.nc_read_box() (which already returns the per-pixel field
# -- load_driver()'s daily loaders just discard it one line later) and
# computes, per zone-day: spatial min/median/mean/max/sd of magnitude, and the
# circular-statistics equivalents for direction (see circular_spatial_stats()
# below for why raw min/median/max don't work on a 0-360 degree bearing).
#
# Direction here is the "direction FROM" decimal-degree value (matching
# .speed_direction()'s existing wind/current convention, and wave's raw VMDR).
# The 8-octant categorical classification (multi.R::compass_octant()) is not
# computed here -- per-pixel direction is kept as raw degrees throughout so
# that classification can be bolted on later without reworking the loader.
#
# Not wired into code/4_time_series.py or code/5_figures.py -- run by hand:
#   Rscript func/compute_driver_spatial_variance.R
source("func/multi.R")

driver_names <- c("wind", "current", "wave")


# Per-pixel loading -----------------------------------------------------

# Load one driver's full per-pixel-per-day field for a zone, cropped to the
# zone's bbox, with magnitude/direction already derived per pixel. Mirrors
# load_driver()'s dispatch/file-discovery exactly, but stops one step short of
# its spatial-mean collapse so the pixel-level spread survives.
load_driver_spatial <- function(driver_name, meta){
  driver_name <- match.arg(driver_name, driver_names)
  zone_box <- dplyr::filter(zones_bbox, zone == meta$zone)
  lon_range <- c(zone_box$lon_min, zone_box$lon_max)
  lat_range <- c(zone_box$lat_min, zone_box$lat_max)

  if(driver_name == "wind"){
    wind_files <- dir(paste0("~/pCloudDrive/data/WIND/", meta$zone), pattern = "_daily_", full.names = TRUE)
    df <- purrr::map_dfr(wind_files, function(file_name){
      pix <- .nc_read_box(file_name, c("eastward_wind", "northward_wind"), lon_range, lat_range) |>
        dplyr::rename(u = eastward_wind, v = northward_wind)
      # Remove final day of data -- artefact from creating daily integrals
      # from hourly data (same fix as load_wind_sub()).
      final_date <- max(pix$date)
      dplyr::filter(pix, date != final_date)
    })
    vec <- .speed_direction(df$u, df$v, convention = "from")
    df$magnitude <- vec$speed
    df$direction <- vec$direction

  } else if(driver_name == "current"){
    dir_name <- path.expand(paste0("~/pCloudDrive/data/GLORYS/", meta$zone))
    df <- dplyr::bind_rows(
      .nc_read_box(file.path(dir_name, "glorys_199301_202412.nc"), c("uo", "vo"), lon_range, lat_range),
      .nc_read_box(file.path(dir_name, "glorys_uo_vo_202501_202512.nc"), c("uo", "vo"), lon_range, lat_range)
    ) |> dplyr::rename(u = uo, v = vo)
    vec <- .speed_direction(df$u, df$v, convention = "to")
    df$magnitude <- vec$speed
    df$direction <- vec$direction

  } else if(driver_name == "wave"){
    wave_files <- dir(paste0("~/pCloudDrive/data/WAVE/", meta$zone), pattern = "_daily_", full.names = TRUE)
    df <- purrr::map_dfr(wave_files, function(file_name){
      nc <- nc_open(file_name)
      has_wave_dir <- "VMDR" %in% names(nc$var)
      nc_close(nc)
      var_names <- if(has_wave_dir) c("VHM0", "VMDR") else "VHM0"
      pix <- .nc_read_box(file_name, var_names, lon_range, lat_range) |>
        dplyr::rename(magnitude = VHM0)
      pix$direction <- if(has_wave_dir) pix$VMDR else NA_real_
      # Remove final day of data -- same daily-integral artefact as load_wave().
      final_date <- max(pix$date)
      dplyr::filter(pix, date != final_date)
    })
  }

  df |> dplyr::select(date, lon, lat, magnitude, direction)
}


# Circular statistics -----------------------------------------------------

# Spatial statistics for a set of same-day direction bearings (decimal
# degrees, "from"/"to" convention as loaded above). Ordinary min/median/max
# break at the 0/360 wraparound (359 and 1 degrees are 2 degrees apart, not
# 358) so:
#   - mean is the circular mean (atan2 of averaged sin/cos), matching the
#     convention util.R::load_wave() already uses for its own zone-day mean.
#   - R is the resultant vector length (0-1): 1 = every pixel points the same
#     way, 0 = directions fully scattered. sd is the circular standard
#     deviation derived from R.
#   - min/median/max are computed on each pixel's signed deviation from the
#     circular mean (wrapped to +-180 degrees), not on the raw bearings --
#     e.g. max_dev = 12 means the most-disagreeing pixel that day pointed 12
#     degrees away from the zone's typical direction that day.
circular_spatial_stats <- function(direction_degrees){
  theta <- direction_degrees * pi / 180
  sin_bar <- mean(sin(theta), na.rm = TRUE)
  cos_bar <- mean(cos(theta), na.rm = TRUE)
  R <- sqrt(sin_bar^2 + cos_bar^2)
  mean_deg <- (atan2(sin_bar, cos_bar) * 180 / pi) %% 360
  sd_deg <- sqrt(-2 * log(R)) * 180 / pi

  dev <- ((direction_degrees - mean_deg + 180) %% 360) - 180  # signed, +-180
  tibble::tibble(
    direction_mean = mean_deg,
    direction_sd = sd_deg,
    direction_R = R,
    direction_min_dev = min(dev, na.rm = TRUE),
    direction_median_dev = median(dev, na.rm = TRUE),
    direction_max_dev = max(dev, na.rm = TRUE)
  )
}

# Reduce one driver's per-pixel-per-day field (load_driver_spatial() output)
# to one row per zone-day: ordinary spatial stats for magnitude, circular
# spatial stats for direction.
spatial_variance_stats <- function(df_pixel){
  df_pixel |>
    dplyr::summarise(
      magnitude_min = min(magnitude, na.rm = TRUE),
      magnitude_median = median(magnitude, na.rm = TRUE),
      magnitude_mean = mean(magnitude, na.rm = TRUE),
      magnitude_max = max(magnitude, na.rm = TRUE),
      magnitude_sd = sd(magnitude, na.rm = TRUE),
      circular_spatial_stats(direction),
      .by = "date"
    )
}


# Main loop ---------------------------------------------------------------

driver_spatial_variance_daily <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)
  purrr::map_dfr(driver_names, function(driver_name){
    message(meta$zone, " / ", driver_name)
    df_pixel <- load_driver_spatial(driver_name, meta)
    spatial_variance_stats(df_pixel) |>
      dplyr::mutate(driver_name = driver_name, zone = meta$zone, .before = 1)
  })
})

if (!dir.exists("output/STATS")) dir.create("output/STATS", recursive = TRUE)
readr::write_csv(driver_spatial_variance_daily, "output/STATS/driver_spatial_variance_daily.csv")


# Timestep aggregation ------------------------------------------------------

# Aggregate the daily spatial-stats table to monthly/annual means for
# plotting. Mirrors multi.R::driver_plume_timesteps()'s date-rounding
# convention. Magnitude stats and the deviation/sd/R direction columns
# aggregate with an ordinary mean; direction_mean is re-aggregated circularly
# (from the underlying sin/cos, not the degree values) so the monthly/annual
# rollup doesn't reintroduce the same wraparound error the daily stat avoids.
spatial_variance_timesteps <- function(df){
  circular_mean_only <- function(direction_degrees){
    theta <- direction_degrees * pi / 180
    (atan2(mean(sin(theta), na.rm = TRUE), mean(cos(theta), na.rm = TRUE)) * 180 / pi) %% 360
  }

  agg_cols <- c("magnitude_min", "magnitude_median", "magnitude_mean", "magnitude_max", "magnitude_sd",
               "direction_sd", "direction_R", "direction_min_dev", "direction_median_dev", "direction_max_dev")

  df_daily <- df |> dplyr::mutate(timestep = "daily")

  df_monthly <- df |>
    dplyr::mutate(date = lubridate::round_date(date, "month") + lubridate::days(14)) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(agg_cols), ~ mean(.x, na.rm = TRUE)),
                     direction_mean = circular_mean_only(direction_mean),
                     .by = c("driver_name", "zone", "date")) |>
    dplyr::mutate(timestep = "monthly")

  df_annual <- df |>
    dplyr::mutate(date = as.Date(paste0(lubridate::year(date), "-07-01"))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(agg_cols), ~ mean(.x, na.rm = TRUE)),
                     direction_mean = circular_mean_only(direction_mean),
                     .by = c("driver_name", "zone", "date")) |>
    dplyr::mutate(timestep = "annual")

  dplyr::bind_rows(df_daily, df_monthly, df_annual) |>
    dplyr::mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))
}

driver_spatial_variance_timesteps <- spatial_variance_timesteps(driver_spatial_variance_daily)


# Plotting --------------------------------------------------------------

# Min-max ribbon + mean/median line, for either the magnitude or direction
# family of columns, one row per timestep (daily/monthly/annual). Direction
# panels plot direction_mean itself, not the deviation columns -- the ribbon
# there is +-1 circular SD around the circular mean rather than min/max
# (min/max deviation is about spread, not a directly plottable bearing).
# ggplot_theme()'s font sizes are tuned for single-panel, article-scale
# figures; this grid packs 3 facets (daily/monthly/annual) x 2 columns
# (magnitude/direction) x 4 zones into one PNG, so it needs its own much
# smaller theme (same reasoning as figure.R::Figure_S_daily_flow's panel_theme).
.spatial_variance_panel_theme <- theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
       strip.background = element_rect(fill = "grey90"),
       strip.text = element_text(size = 7),
       axis.title = element_text(size = 8, colour = "black"),
       axis.text = element_text(size = 6, colour = "black"),
       axis.text.x = element_text(angle = 45, hjust = 1),
       panel.grid.minor = element_blank(),
       panel.border = element_rect(fill = NA, colour = "black"),
       plot.margin = margin(t = 4, r = 6, b = 2, l = 2))

.variance_panel <- function(df_ts, value_type){
  disp <- dplyr::filter(driver_display, driver_name == unique(df_ts$driver_name))

  if(value_type == "magnitude"){
    df_ts <- df_ts |> dplyr::rename(lo = magnitude_min, mid = magnitude_mean, med = magnitude_median, hi = magnitude_max)
    y_lab <- disp$driver_label
  } else {
    df_ts <- df_ts |>
      dplyr::mutate(lo = direction_mean - direction_sd, hi = direction_mean + direction_sd) |>
      dplyr::rename(mid = direction_mean)
    df_ts$med <- NA_real_
    y_lab <- "Direction (deg, +-1 circular SD)"
  }

  ggplot(df_ts, aes(x = date)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = disp$driver_colour, alpha = 0.25) +
    geom_line(aes(y = mid), colour = disp$driver_colour, linewidth = 0.5) +
    { if(value_type == "magnitude") geom_line(aes(y = med), colour = "black", linetype = "dashed", linewidth = 0.3) } +
    facet_wrap(~timestep, ncol = 1, scales = "free_x") +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(x = NULL, y = y_lab) +
    .spatial_variance_panel_theme
}

# One PNG per driver: magnitude spread (left) + direction spread (right),
# one row per zone, saved as a diagnostic (not manuscript ARTICLE/) figure --
# manuscript/TODO.md Section 2.8 has not decided whether this becomes real
# Methods/Results content or a Discussion-limitation caveat.
plot_spatial_variance <- function(driver_name, output_dir = "figures/driver_spatial_variance"){
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  panels <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    df_ts <- dplyr::filter(driver_spatial_variance_timesteps, driver_name == !!driver_name, zone == meta$zone)
    list(magnitude = .variance_panel(df_ts, "magnitude") + labs(title = zone_title(meta$zone)),
        direction = .variance_panel(df_ts, "direction") + labs(title = zone_title(meta$zone)))
  })

  plotlist <- panels |> purrr::map(function(p) list(p$magnitude, p$direction)) |> purrr::flatten()
  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(panels), align = "v",
                                 labels = panel_labels, font.label = list(size = 14, face = "bold"))

  save_plot_as_png(full_plot, paste0("spatial_variance_", driver_name), width = 14, height = 16, path = output_dir)
}

purrr::walk(driver_names, plot_spatial_variance)
