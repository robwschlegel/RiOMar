# One-off: spatial-consistency diagnostic for wind/current/wave direction and
# magnitude within each zone's bbox. Addresses the manuscript/TODO.md Section
# 2.8 deferral -- func/multi.R::load_driver() collapses every driver to a
# single zone-day scalar (spatial mean) before any consumer sees it, so no
# consumer downstream of it can check whether pixels within a zone actually
# agree with each other on a given day. This script re-reads the same raw
# files via util.R::.nc_read_box() (which already returns the per-pixel field
# -- load_driver()'s daily loaders just discard it one line later).
#
# Two complementary views:
#   - magnitude: spatial min/median/mean/max/sd per zone-day (ordinary
#     numeric stats -- see spatial_variance_stats()).
#   - direction: rather than a continuous circular-statistics spread (tried
#     first, then dropped per Robert's steer 2026-08-11 towards a more
#     directly interpretable view), each pixel is binned into one of
#     multi.R::compass_octant()'s 8 categories and counted per zone-day, then
#     rolled up to monthly/annual totals -- see octant_pixel_counts_daily().
#     A dominant, stable octant means the zone's direction agrees spatially;
#     counts spread evenly across all 8 means it doesn't.
#
# Direction here is the "direction FROM" decimal-degree value (matching
# .speed_direction()'s existing wind/current convention, and wave's raw
# VMDR) before binning into octants.
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


# Magnitude spatial statistics -----------------------------------------

# Reduce one driver's per-pixel-per-day field (load_driver_spatial() output)
# to one row per zone-day: ordinary spatial min/median/mean/max/sd of
# magnitude.
spatial_variance_stats <- function(df_pixel){
  df_pixel |>
    dplyr::summarise(
      magnitude_min = min(magnitude, na.rm = TRUE),
      magnitude_median = median(magnitude, na.rm = TRUE),
      magnitude_mean = mean(magnitude, na.rm = TRUE),
      magnitude_max = max(magnitude, na.rm = TRUE),
      magnitude_sd = sd(magnitude, na.rm = TRUE),
      .by = "date"
    )
}


# Direction octant pixel-counts -----------------------------------------

# Bin each pixel's direction into one of compass_octant()'s 8 categories
# (multi.R, already sourced) and count pixels per octant per zone-day, plus
# that count as a proportion of the day's total pixels. tidyr::complete()
# fills octants with zero pixels that day (rather than a missing row) so the
# monthly/annual aggregation below doesn't silently drop them.
#
# Proportion, not just the raw count, matters here: wind's per-pixel grid
# resolution jumps 4x at a file boundary in 2008 (28 pixels/zone-day before,
# 112 after -- a real characteristic of the source ERA5 files, not a bug),
# which otherwise swamps any octant's raw pixel_count with a step change
# unrelated to actual wind direction. Dividing by that day's total pixel
# count removes the confound; current and wave don't have this issue (their
# grids are constant), but proportion is computed uniformly for all three so
# the panels stay directly comparable.
octant_pixel_counts_daily <- function(df_pixel){
  octant_labels <- levels(compass_octant(0))
  counts <- df_pixel |>
    dplyr::mutate(octant = compass_octant(direction)) |>
    dplyr::filter(!is.na(octant)) |>
    dplyr::summarise(pixel_count = dplyr::n(), .by = c("date", "octant")) |>
    tidyr::complete(date = unique(df_pixel$date), octant = octant_labels, fill = list(pixel_count = 0))

  totals <- counts |> dplyr::summarise(total_pixels = sum(pixel_count), .by = "date")

  counts |>
    dplyr::left_join(totals, by = "date") |>
    dplyr::mutate(proportion = pixel_count / total_pixels) |>
    dplyr::select(date, octant, pixel_count, proportion)
}


# Main loop ---------------------------------------------------------------

# Single load_driver_spatial() pass per zone/driver, shared by the magnitude
# spatial-stats and the octant pixel-count tables below (avoids re-reading
# the raw NetCDFs twice -- the wave files alone are 500+ MB per zone).
zone_driver_results <- purrr::pmap(zone_meta, function(...){
  meta <- tibble::tibble(...)
  purrr::map(driver_names, function(driver_name){
    message(meta$zone, " / ", driver_name)
    df_pixel <- load_driver_spatial(driver_name, meta)
    variance <- spatial_variance_stats(df_pixel) |>
      dplyr::mutate(driver_name = driver_name, zone = meta$zone, .before = 1)
    octant_counts <- octant_pixel_counts_daily(df_pixel) |>
      dplyr::mutate(driver_name = driver_name, zone = meta$zone, .before = 1)
    list(variance = variance, octant_counts = octant_counts)
  })
}) |> purrr::flatten()

driver_spatial_variance_daily <- purrr::map_dfr(zone_driver_results, "variance")
driver_octant_pixel_counts_daily <- purrr::map_dfr(zone_driver_results, "octant_counts")

if (!dir.exists("output/STATS")) dir.create("output/STATS", recursive = TRUE)
readr::write_csv(driver_spatial_variance_daily, "output/STATS/driver_spatial_variance_daily.csv")
readr::write_csv(driver_octant_pixel_counts_daily, "output/STATS/driver_octant_pixel_counts_daily.csv")


# Timestep aggregation ------------------------------------------------------

# Aggregate the daily magnitude spatial-stats table to monthly/annual means
# for plotting. Mirrors multi.R::driver_plume_timesteps()'s date-rounding
# convention.
spatial_variance_timesteps <- function(df){
  agg_cols <- c("magnitude_min", "magnitude_median", "magnitude_mean", "magnitude_max", "magnitude_sd")

  df_daily <- df |> dplyr::mutate(timestep = "daily")

  df_monthly <- df |>
    dplyr::mutate(date = lubridate::round_date(date, "month") + lubridate::days(14)) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(agg_cols), ~ mean(.x, na.rm = TRUE)),
                     .by = c("driver_name", "zone", "date")) |>
    dplyr::mutate(timestep = "monthly")

  df_annual <- df |>
    dplyr::mutate(date = as.Date(paste0(lubridate::year(date), "-07-01"))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(agg_cols), ~ mean(.x, na.rm = TRUE)),
                     .by = c("driver_name", "zone", "date")) |>
    dplyr::mutate(timestep = "annual")

  dplyr::bind_rows(df_daily, df_monthly, df_annual) |>
    dplyr::mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))
}

driver_spatial_variance_timesteps <- spatial_variance_timesteps(driver_spatial_variance_daily)

# Aggregate the daily octant proportion table to monthly/annual means (mean
# of the daily proportion, not a ratio-of-sums -- keeps every day weighted
# equally regardless of that day's underlying pixel-grid resolution). Same
# date-rounding convention as above.
octant_pixel_counts_timesteps <- function(df){
  df_daily <- df |> dplyr::mutate(timestep = "daily")

  df_monthly <- df |>
    dplyr::mutate(date = lubridate::round_date(date, "month") + lubridate::days(14)) |>
    dplyr::summarise(proportion = mean(proportion, na.rm = TRUE),
                     .by = c("driver_name", "zone", "octant", "date")) |>
    dplyr::mutate(timestep = "monthly")

  df_annual <- df |>
    dplyr::mutate(date = as.Date(paste0(lubridate::year(date), "-07-01"))) |>
    dplyr::summarise(proportion = mean(proportion, na.rm = TRUE),
                     .by = c("driver_name", "zone", "octant", "date")) |>
    dplyr::mutate(timestep = "annual")

  dplyr::bind_rows(df_daily, df_monthly, df_annual) |>
    dplyr::mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))
}

driver_octant_pixel_counts_timesteps <- octant_pixel_counts_timesteps(driver_octant_pixel_counts_daily)


# Plotting --------------------------------------------------------------

# ggplot_theme()'s font sizes are tuned for single-panel, article-scale
# figures; these grids pack several facets into one PNG, so they need their
# own much smaller theme (same reasoning as figure.R::Figure_S_daily_flow's
# panel_theme).
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

# Min-max ribbon + mean/median line, one row per timestep (daily/monthly/
# annual).
.variance_panel <- function(df_ts, driver_name){
  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)
  df_ts <- df_ts |> dplyr::rename(lo = magnitude_min, mid = magnitude_mean, med = magnitude_median, hi = magnitude_max)

  ggplot(df_ts, aes(x = date)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = disp$driver_colour, alpha = 0.25) +
    geom_line(aes(y = mid), colour = disp$driver_colour, linewidth = 0.5) +
    geom_line(aes(y = med), colour = "black", linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~timestep, ncol = 1, scales = "free_x") +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(x = NULL, y = disp$driver_label) +
    .spatial_variance_panel_theme
}

# One PNG per driver: magnitude spread, one row per zone, saved as a
# diagnostic (not manuscript ARTICLE/) figure -- manuscript/TODO.md Section
# 2.8 has not decided whether this becomes real Methods/Results content or a
# Discussion-limitation caveat.
plot_spatial_variance <- function(driver_name, output_dir = "figures/driver_spatial_variance"){
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    df_ts <- dplyr::filter(driver_spatial_variance_timesteps, driver_name == !!driver_name, zone == meta$zone)
    .variance_panel(df_ts, driver_name) + labs(title = zone_title(meta$zone))
  })

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 1, nrow = length(plotlist),
                                 labels = panel_labels, font.label = list(size = 14, face = "bold"))

  save_plot_as_png(full_plot, paste0("spatial_variance_", driver_name), width = 10, height = 16, path = output_dir)
}

# Octant pixel-count barplot: bars = proportion of that zone-day's pixels in
# each octant, black line = simple OLS trend (geom_smooth(method = "lm"))
# fit purely for visual inspection, not a manuscript trend claim (that's
# multi.R::compute_octant_trend(), fit on the zone-day-averaged direction via
# func/compute_direction_octant_trend.R -- a different, day-count-per-year
# view, not this pixel-proportion-per-period one).
#
# y_max is shared across every panel in the calling plot_octant_pixel_counts()
# image (all 4 zones x 8 octants) -- these are investigation figures, not for
# publication, so a free per-facet y-axis is exactly what NOT to do: a rare
# octant auto-scaling to fill its own panel would visually hide that it's
# rare. One fixed y-axis makes a small proportion immediately look small.
.octant_pixel_count_panel <- function(df_ts, driver_name, zone_name, y_max){
  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)

  ggplot(df_ts, aes(x = date, y = proportion)) +
    geom_col(fill = disp$driver_colour) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
    facet_wrap(~octant, nrow = 2) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(x = NULL, y = "Proportion of pixels", title = zone_title(zone_name)) +
    .spatial_variance_panel_theme
}

# One PNG per driver x timestep (daily/monthly/annual): 8-octant proportion
# barplots, one row per zone, all sharing one y-axis across the whole image.
plot_octant_pixel_counts <- function(driver_name, timestep_label, output_dir = "figures/driver_spatial_variance"){
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  df_driver_ts <- dplyr::filter(driver_octant_pixel_counts_timesteps, driver_name == !!driver_name,
                                timestep == timestep_label)
  y_max <- max(df_driver_ts$proportion, na.rm = TRUE)

  plotlist <- purrr::pmap(zone_meta, function(...){
    meta <- tibble::tibble(...)
    df_ts <- dplyr::filter(df_driver_ts, zone == meta$zone)
    .octant_pixel_count_panel(df_ts, driver_name, meta$zone, y_max)
  })

  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 1, nrow = length(plotlist),
                                 labels = panel_labels, font.label = list(size = 14, face = "bold"))

  save_plot_as_png(full_plot, paste0("octant_pixel_counts_", timestep_label, "_", driver_name),
                   width = 14, height = 16, path = output_dir)
}

purrr::walk(driver_names, plot_spatial_variance)

purrr::walk(driver_names, function(dn){
  purrr::walk(c("daily", "monthly", "annual"), function(ts) plot_octant_pixel_counts(dn, ts))
})
