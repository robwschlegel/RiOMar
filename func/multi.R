# func/multi.R
# Loads all drivers of plume size, performs stats, plots results
#
# REFACTOR (2026-07) -------------------------------------------------------
# func/flow.R, func/wind.R, func/tide.R, func/ROFI.R, and func/surface.R have
# been consolidated into this file. Each of those scripts independently
# reimplemented the same "resolve zone from river mouth -> load driver ->
# load plume -> combine -> STL -> correlate (incl. lagged) -> plot" pipeline
# for a single driver (river flow, wind, tide, ROFI/coastal-current extent).
# That pattern is now written once, below, and parameterised by driver name.
# util.R already centralised the raw *loading* functions (load_river_flow,
# load_tide_gauge, load_wind_sub, load_ROFI, load_plume_ts) and a few
# generic helpers (lagged_correlation, make_pretty_title, stl_single,
# sec_axis_adjustement_factors) -- the individual driver scripts had mostly
# stopped using those and re-derived their own inline versions instead. This
# refactor routes everything back through those existing util.R functions.
#
# What changed, concretely:
#   - The zone/gauge lookup if/else block that was copy-pasted 3x in flow.R,
#     once each in wind.R/tide.R/surface.R, and again (slightly differently)
#     in multi_stl() below, is now the single `zone_meta` table + get_zone_meta().
#   - The inline plume-CSV-loading block in flow.R/wind.R/tide.R (which used
#     the raw `area_of_the_plume_mask_in_km2` column name and skipped the
#     day-gap `complete()` step) is replaced by the existing
#     util.R::load_plume_ts(), which already does this correctly.
#   - The on-/off-shore wind classification duplicated in wind.R, surface.R,
#     and multi_stl() below is now wind_add_direction().
#   - flow.R/wind.R/tide.R's near-identical "4-panel" comparison plot
#     (raw driver ts / plume ts / scatter / lag-correlation) is now
#     plot_driver_comparison().
#   - wind.R/tide.R's inline dual-y-axis STL plot (a hand-copy of
#     figure.R::plot_x11_river_and_plume(), hard-coded there to
#     the "river_flow" column name) is now plot_driver_plume_dual_axis(),
#     which is the generalisation flagged as a TODO in
#     manuscript/make_figures_tables.R for Figures 7-9.
#   - flow.R's flow_plume_trend_plus() (the Sutton et al. 2022-style
#     AR(1)/STL-weighted HAC trend estimator manuscript.tex Sec. 2.6.1
#     refers to) was river-flow-specific; it is now driver_plume_trend(),
#     usable for any driver.
#   - ROFI.R's comp_ROFI_plume() had the richest correlation structure
#     (daily/monthly/annual lag correlations) of the four scripts; that
#     structure is now the shared driver_plume_correlation() /
#     driver_plume_timesteps(), used for every driver, not just currents.
#   - surface.R's surface_plot_daily_maps() is carried over with light
#     cleanup (now uses zone_meta / make_pretty_title instead of its own
#     inline zone lookup); this is conceptually different from the
#     time-series comparisons above (per-pixel maps, not per-day series) so
#     it keeps its own section rather than being folded into the
#     driver_plume_* functions.
#   - func/flow.R, func/wind.R, func/tide.R, func/ROFI.R, func/surface.R now
#     just source() this file for backwards compatibility -- nothing that
#     called them directly should need to change, but new work should target
#     the functions here instead.
#   - Not touched / out of scope: func/VOG.R (a separate side-study, not a
#     plume-driver comparison), func/sentinel.R (Sentinel NetCDF loading
#     utility), func/analyse_spatiotemporal.R (one-off exploratory script for
#     a non-RiOMar bounding box), func/river_flow_prep.R, func/glorys_download.R,
#     func/ODATIS-MR_expert_subset.R (upstream data prep/download, not
#     driver-vs-plume analysis).
#
# NB: This refactor has not been executed against real data (no R available
# in the environment it was written in) -- the composition of already-tested
# util.R functions was checked function-by-function for signature and
# semantic equivalence with the code it replaces, but please run it and spot
# check a zone or two against the old figures/*.png outputs before trusting
# it blindly, and keep the old flow.R/wind.R/tide.R/ROFI.R/surface.R content
# in git history in case anything needs to be recovered.


# Analysis ideas ----------------------------------------------------------

# Create GIFs

# Basic time series comparisons to get seasonal and interannual comparisons
## Perform seasonal smoothing with heatwaveR
## also look into fixing the tidal range time series
## ultimately the point is to reduce the data in such a way that the time series can be related to one another
## what does an extreme event analysis reveal?
### how do X11 and STL differ?
## Look at seasonal Trends per month, not long term
## Get correlations of interannual and seasonal time series
## Dynamic Linear Models (DLM) is another option: library(dlm)

# Treat each pixel like its own time series and see what is happening with the forces when the pixel is triggered as a panache
## also how high SPM is while all this is happening
## show primary wind direction when pixel is triggered
## also relationship with SPM and tide range or category
## number of times pixel is flagged related to the size of the total panache when it is flagged
### Would need to relate wind with time lag to this as well
## could also tally the shape of the panache whenever pixel is flagged
# Other analyses
## Get mean offshore distance of centroid
## when creating GIFs, would be cool to have the centroid visualised as a 21 dot
### it could leave a trail of 1 dots behind it per day
### or fill colour could be left to show tide, wind, etc. on that day

# More ideas
## nmds of mean characteristics during plume events
## extreme event analysis
### need to establish a reasonable baseline, and go from there
### could be interesting to use the time varying X11 seasonality
### but then is it relevant to calculate event stats from the median signal?
  ### one could modify the analysis to always take from the base, being zero
### how would one establish the 90th percentile? Or rather, just always take everything over the seasonal signal
## EMD - see Vincent email
### percent contribution of each component to time series
### changes over time as well
## like a CTD cast, figure out a way to measure from when a plume peaks and then goes down below a certain threshold as a way of determining if it is an individual event
### and then get statistics from that
### also add up number of days with onshore wind, neap tide, etc.
### the spatial threshold could be based on the percentile of the total plume area over the full time series
### start searching by creating a contour plot of every 10th percentile
## also account for the lon/lat of the centroid and how that relates to the other drivers
# be able to say if the plume is increasing or decreasing in size so that the drivers on that day can be categorised under what the plume is doing

## Ideas for driver decomposition
# A SOM analysis of the spatial footprint during given drivers may work.
## Average the primary drivers over the timespan of the given large plume
## Then see if the SOM organises them by surface signature. I.e. size, centroid, shape, direction, duration
# Or PCA/ordination of days of panache given certain values for variables
# Use rle() to determine contiguous events temporarily
# Look at articles in pCloud folder for other analysis ideas and in situ data sources

## NB: For the *interaction* between drivers specifically (not just each
# driver individually vs. plume, which is what this file covers), see
# manuscript/driver_interactions_review.md and code/6_driver_interactions.R,
# which implement the GLM-with-interactions / GAM / regime-stratification
# road map agreed from that literature review.

## Temperature analysis
# Panache should be colder than coastal water
# It should be possible to a priori define the colder surface temperatures based on the panache pixels
# Use the ODATIS-MR SST for this
# Look at validation results to see if one product is better than other, also by zone


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ncdf4) # For reading NetCDF sources (load_plume_surface, load_wind_sub, load_wave, etc. -- util.R)
library(heatwaveR) # For seasonal smoothing analysis
library(seasonal) # For X11 analysis (currently not used)
library(RcppRoll) # For running means to get STL interannual signals closer to X11
library(patchwork)
library(sandwich) # For HAC covariance tests (driver_plume_trend)
library(lmtest) # For more detailed linear model tests (driver_plume_trend)
library(doParallel); doParallel::registerDoParallel(cores = 14)

# Common function
source("func/util.R")

# Zones, north to south (see func/util.R::ZONE_ORDER/order_zones())
zones <- ZONE_ORDER


# Zone / gauge metadata -----------------------------------------------------

# Canonical river mouth -> zone -> tide gauge lookup. Replaces the identical
# if/else block that was previously copy-pasted in flow_comp(), flow_trend(),
# flow_plume_trend_plus() (all three in the old flow.R), tide_calc() (tide.R),
# spatial_wind_calc() (wind.R), surface_plot() (surface.R), and multi_stl()
# further down this file.
zone_meta <- river_mouths |>
  dplyr::mutate(
    zone = dplyr::case_when(
      mouth_name == "Seine"       ~ "BAY_OF_SEINE",
      mouth_name == "Gironde"     ~ "BAY_OF_BISCAY",
      mouth_name == "Loire"       ~ "SOUTHERN_BRITTANY",
      mouth_name == "Grand Rhone" ~ "GULF_OF_LION",
      TRUE ~ NA_character_
    ),
    gauge = dplyr::case_when(
      zone == "BAY_OF_SEINE"      ~ "LE_HAVRE",
      zone == "BAY_OF_BISCAY"     ~ "PORT-BLOC",
      zone == "SOUTHERN_BRITTANY" ~ "SAINT-NAZAIRE",
      zone == "GULF_OF_LION"      ~ "MARSEILLE",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::arrange(match(zone, ZONE_ORDER))

# Look up zone metadata by either the river mouth name (as used in
# river_mouths/zone_meta) or the zone code (as used in zones_bbox / output/
# paths). Exactly one of mouth_name / zone_name should be supplied.
# get_zone_meta(mouth_name = "Seine")
# get_zone_meta(zone_name = "GULF_OF_LION")
get_zone_meta <- function(mouth_name = NULL, zone_name = NULL){
  if(!is.null(mouth_name)){
    out <- dplyr::filter(zone_meta, mouth_name == !!mouth_name)
  } else if(!is.null(zone_name)){
    out <- dplyr::filter(zone_meta, zone == !!zone_name)
  } else {
    stop("Supply either mouth_name or zone_name to get_zone_meta().")
  }
  if(nrow(out) != 1) stop("Zone/mouth not recognised in zone_meta.")
  return(out)
}

# Human-facing zone display name (e.g. "BAY_OF_SEINE" -> "Bay of Seine"),
# vectorised. The single shared source of this mapping across the R side of
# the codebase (previously duplicated ad hoc, e.g. Figure_S3_seasonal_boxplots()'s
# local zone_labels vector) -- use this instead of hand-writing zone titles or
# leaving zone codes in ALL_CAPS on any human-facing plot/table/caption. See
# func/util.py::zone_title() for the Python-side equivalent.
zone_display_names <- c(
  BAY_OF_SEINE      = "Bay of Seine",
  BAY_OF_BISCAY     = "Bay of Biscay",
  SOUTHERN_BRITTANY = "Southern Brittany",
  GULF_OF_LION      = "Gulf of Lion"
)
zone_title <- function(zone_name){
  out <- unname(zone_display_names[zone_name])
  if(anyNA(out) && !anyNA(zone_name)) stop("zone_title(): unrecognised zone code(s): ",
                                            paste(setdiff(zone_name, names(zone_display_names)), collapse = ", "))
  out
}


# Driver loading ------------------------------------------------------------

# Shared on-/off-shore wind classification. Previously duplicated verbatim
# (down to the same TODO comment) in wind.R::spatial_wind_calc(),
# surface.R::surface_plot(), and multi_stl() further down this file.
# NB: wind_spd/wind_dir are now computed upstream in util.R::load_wind_sub();
# this only adds the on-/off-shore label.
# TODO (carried over from the originals): think of a more sophisticated way
# to classify on-/off-shore than a hard-coded sign check per zone.
wind_add_direction <- function(df_wind, zone_name){
  if(zone_name %in% c("BAY_OF_BISCAY", "SOUTHERN_BRITTANY")){
    df_wind <- df_wind |> dplyr::mutate(direction = ifelse(u < 0, "off", "on"))
  } else if(zone_name == "BAY_OF_SEINE"){
    df_wind <- df_wind |> dplyr::mutate(direction = ifelse(v > 0, "off", "on"))
  } else if(zone_name == "GULF_OF_LION"){
    df_wind <- df_wind |> dplyr::mutate(direction = ifelse(v < 0, "off", "on"))
  } else {
    stop("Zone not recognised for wind direction classification.")
  }
  df_wind
}

# Load one driver's daily time series for a zone, in a common two-column
# (date, value) shape so downstream functions don't need to know which
# driver they're looking at. Generalises the per-driver loading blocks
# previously duplicated across flow.R/wind.R/tide.R/ROFI.R.
# load_driver("flow", get_zone_meta(zone_name = "GULF_OF_LION"))
load_driver <- function(driver_name, meta){
  driver_name <- match.arg(driver_name, c("flow", "tide", "wind", "current", "rofi", "wave"))

  if(driver_name == "flow"){
    df <- load_river_flow(paste0("data/RIVER_FLOW/", meta$zone)) |>
      dplyr::select(date, value = flow)

  } else if(driver_name == "tide"){
    df <- load_tide_gauge(paste0("data/TIDES/", meta$gauge)) |>
      dplyr::select(date, value = tide_range, tide_mean)

  } else if(driver_name == "wind"){
    lon_round <- plyr::round_any(meta$mouth_lon, 0.5)
    lat_round <- plyr::round_any(meta$mouth_lat, 0.5)
    lon_range <- c(lon_round - 0.5, lon_round + 0.5)
    lat_range <- c(lat_round - 0.5, lat_round + 0.5)
    wind_files <- dir(paste0("~/pCloudDrive/data/WIND/", meta$zone), pattern = "_daily_", full.names = TRUE)
    df_wind <- purrr::map_dfr(wind_files, load_wind_sub, lon_range, lat_range) |>
      wind_add_direction(meta$zone)
    df <- df_wind |> dplyr::select(date, value = wind_spd, wind_dir, direction, u, v)

  } else if(driver_name == "current"){
    # "current" = GLORYS surface current (eastward/northward velocity) near the
    # river mouth -- the "coastal ocean currents" driver discussed in
    # manuscript/driver_interactions_review.md. Unlike ROFI (see "rofi" below),
    # GLORYS covers all four zones, including the Gulf of Lion.
    lon_round <- plyr::round_any(meta$mouth_lon, 0.5)
    lat_round <- plyr::round_any(meta$mouth_lat, 0.5)
    lon_range <- c(lon_round - 0.5, lon_round + 0.5)
    lat_range <- c(lat_round - 0.5, lat_round + 0.5)
    df_current <- load_surface_current(meta$zone, lon_range, lat_range)
    df <- df_current |> dplyr::select(date, value = current_spd, current_dir, u, v)

  } else if(driver_name == "rofi"){
    # "rofi" = region-of-freshwater-influence extent, a model output distinct
    # from the ambient coastal current above. NB: no ROFI data for the Gulf of
    # Lion (see comment on the old ROFI.R::walk(zones[1:3], ...) call).
    rofi_files <- dir("data/ROFI", full.names = TRUE)
    df_rofi <- purrr::map_dfr(rofi_files, load_ROFI) |>
      dplyr::filter(zone == meta$zone)
    if(nrow(df_rofi) == 0) stop("No ROFI data available for zone ", meta$zone)
    df <- df_rofi |> dplyr::select(date, value = ROFI_surface)

  } else if(driver_name == "wave"){
    lon_round <- plyr::round_any(meta$mouth_lon, 0.5)
    lat_round <- plyr::round_any(meta$mouth_lat, 0.5)
    lon_range <- c(lon_round - 0.5, lon_round + 0.5)
    lat_range <- c(lat_round - 0.5, lat_round + 0.5)
    wave_files <- dir(paste0("~/pCloudDrive/data/WAVE/", meta$zone), pattern = "_daily_", full.names = TRUE)
    df_wave <- purrr::map_dfr(wave_files, load_wave, lon_range, lat_range)
    df <- df_wave |> dplyr::select(date, value = wave_height, wave_dir)
  }

  return(df)
}


# Combine plume + one driver --------------------------------------------

# Load plume + a single driver for one zone, join on date, and add STL
# interannual columns for both. This is the core object every comparison
# function below operates on.
# combine_plume_driver("flow", get_zone_meta(mouth_name = "Seine"))
# combine_plume_driver("flow", get_zone_meta(mouth_name = "Seine"), metric_col = "mass_SPM_in_the_plume_area_in_g_m", outlier_max = NULL)
combine_plume_driver <- function(driver_name, meta, metric_col = "area_of_the_plume_mask_in_km2", outlier_max = 20000){

  df_plume  <- load_plume_ts(meta$zone, metric_col = metric_col, outlier_max = outlier_max)  # util.R -- already handles gap-filling + outlier removal
  df_driver <- load_driver(driver_name, meta)

  df <- dplyr::left_join(df_plume, df_driver, by = "date") |>
    zoo::na.trim()

  start_date <- min(df$date)
  df$plume_stl  <- stl_single(df$plume_area, "inter", start_date)  # util.R
  df$driver_stl <- stl_single(df$value,      "inter", start_date)

  df <- df |> dplyr::mutate(driver_name = driver_name, zone = meta$zone, .before = "date")
  return(df)
}


# Multi-timestep correlation -------------------------------------------------

# Aggregate a combine_plume_driver() data.frame to daily/monthly/annual means.
# Generalises the timestep aggregation from the old ROFI.R::comp_ROFI_plume().
driver_plume_timesteps <- function(df){
  df_daily <- df |> dplyr::mutate(timestep = "daily")

  df_monthly <- df |>
    dplyr::mutate(date = round_date(date, "month") + days(14)) |>
    dplyr::summarise(plume_area = mean(plume_area, na.rm = TRUE),
                      value = mean(value, na.rm = TRUE), .by = "date") |>
    dplyr::mutate(timestep = "monthly")

  df_annual <- df |>
    dplyr::mutate(date = as.Date(paste0(year(date), "-07-01"))) |>
    dplyr::summarise(plume_area = mean(plume_area, na.rm = TRUE),
                      value = mean(value, na.rm = TRUE), .by = "date") |>
    dplyr::mutate(timestep = "annual")

  list(daily = df_daily, monthly = df_monthly, annual = df_annual)
}

# Daily/monthly/annual (lagged) correlation between plume area and one driver.
# Uses util.R::lagged_correlation() (x = plume_area gets lagged, y = driver
# value is the fixed reference) instead of the three near-identical inline
# tibble(lag = ..., cor = map_dbl(...)) blocks previously in
# flow.R/wind.R/tide.R, and the daily/monthly/annual structure previously
# unique to ROFI.R::comp_ROFI_plume().
driver_plume_correlation <- function(df, max_lag_daily = 30){
  ts_list <- driver_plume_timesteps(df)

  cor_daily <- lagged_correlation(x = ts_list$daily$plume_area, y = ts_list$daily$value, max_lag_daily) |>
    dplyr::mutate(timestep = "daily")
  cor_monthly <- lagged_correlation(x = ts_list$monthly$plume_area, y = ts_list$monthly$value, 12) |>
    dplyr::mutate(timestep = "monthly")
  n_annual_lag <- max(0, min(10, nrow(ts_list$annual) - 1))
  cor_annual <- lagged_correlation(x = ts_list$annual$plume_area, y = ts_list$annual$value, n_annual_lag) |>
    dplyr::mutate(timestep = "annual")

  dplyr::bind_rows(cor_daily, cor_monthly, cor_annual) |>
    dplyr::mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))
}


# Plotting --------------------------------------------------------------

# Driver display metadata (label + colour), used by the plotting functions
# below so callers only need to pass a driver_name.
driver_display <- tibble::tribble(
  ~driver_name, ~driver_label,             ~driver_colour,
  "flow",       "River flow (m^3 s-1)",    "blue",
  "tide",       "Tidal range (m)",         "darkgreen",
  "wind",       "Wind speed (m s-1)",      "purple",
  "current",    "Current speed (m s-1)",   "orchid",
  "rofi",       "ROFI extent (km^2)",      "goldenrod",
  "wave",       "Wave height (m)",         "steelblue"
)

# The 4-panel (a-d) comparison plot previously duplicated in
# flow.R::flow_comp(), wind.R::spatial_wind_calc(), and tide.R::tide_calc():
#   a) raw driver time series
#   b) raw plume-area time series
#   c) driver vs. plume scatter (+ linear fit)
#   d) lagged correlation (plume lagged behind driver, 0-30 days)
# mouth_name is used only for the plot title / output file name.
plot_driver_comparison <- function(df, driver_name, mouth_name){

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)
  cor_df <- driver_plume_correlation(df) |> dplyr::filter(timestep == "daily")

  driver_plot <- ggplot(df, aes(x = date, y = value)) +
    geom_line() +
    labs(y = disp$driver_label, x = NULL) +
    scale_x_date(expand = c(0, 0)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  panache_plot <- ggplot(df, aes(x = date, y = plume_area)) +
    geom_line() +
    labs(y = "plume area (km^2)", x = NULL) +
    scale_x_date(expand = c(0, 0)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  driver_plume_cor_plot <- ggplot(df, aes(x = value, y = plume_area)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 2) +
    labs(y = "plume area (km^2)", x = disp$driver_label) +
    theme(panel.border = element_rect(fill = NA, colour = "black"), legend.position = "bottom")

  driver_plume_cor_lag_plot <- ggplot(cor_df, aes(x = lag, y = cor)) +
    geom_point() +
    labs(x = paste("lag plume after", disp$driver_name, "(days)"), y = "correlation (r)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  plot_title <- grid::textGrob(paste0(mouth_name, " : ", driver_name, " vs plume size"),
                               gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(driver_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
  cor_plot <- ggpubr::ggarrange(driver_plume_cor_plot, driver_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
  full_plot_title <- ggpubr::ggarrange(plot_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")

  ggsave(filename = paste0("figures/driver_comparison/cor_plot_", driver_name, "_plume_", mouth_name, ".png"),
         plot = full_plot_title, width = 12, height = 6, dpi = 600)
  invisible(full_plot_title)
}

# Dual-y-axis STL plot (plume_stl on the left axis, driver_stl scaled onto
# the right axis). Generalises the inline duplicate of
# figure.R::plot_x11_river_and_plume() that wind.R and tide.R had
# each hand-copied and hard-coded to their own driver; also closes the TODO
# left in manuscript/make_figures_tables.R for Figures 7-9 (wind/wave/tide
# X11 panels), which asked for exactly this generalisation.
plot_driver_plume_dual_axis <- function(df, driver_name, zone_name){

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)

  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = df$driver_stl, var_ref = df$plume_stl)
  df <- df |> dplyr::mutate(driver_scaled = driver_stl * scaling_factor$diff + scaling_factor$adjust)
  unique_years <- df$date |> year() |> unique()

  pl <- ggplot(data = df) +
    geom_point(aes(x = date, y = plume_stl), color = "brown") +
    geom_path(aes(x = date, y = plume_stl), color = "brown") +
    geom_point(aes(x = date, y = driver_scaled), color = disp$driver_colour) +
    geom_path(aes(x = date, y = driver_scaled), color = disp$driver_colour) +
    scale_x_date(name = "",
                 breaks = paste(unique_years, "01-01", sep = "-") |> as.Date(),
                 labels = unique_years |> str_extract_all('[0-9][0-9]$') |> unlist()) +
    scale_y_continuous(name = "Plume area (km²)",
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = disp$driver_label)) +
    labs(title = zone_name) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y.left = element_text(color = "brown"),
          axis.ticks.y.left = element_line(color = "brown"),
          axis.line.y.left = element_line(color = "brown"),
          axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
          axis.text.y.right = element_text(color = disp$driver_colour),
          axis.ticks.y.right = element_line(color = disp$driver_colour),
          axis.line.y.right = element_line(color = disp$driver_colour),
          axis.title.y.right = element_text(color = disp$driver_colour, margin = unit(c(0, 0, 0, 7.5), "mm")),
          panel.border = element_rect(linetype = "solid", fill = NA))

  ggsave(filename = paste0("figures/driver_comparison/dual_axis_", driver_name, "_plume_", zone_name, ".png"),
         plot = pl, width = 12, height = 6, dpi = 300)
  invisible(pl)
}

# Bin a compass bearing (degrees, "from" convention) into one of 8 ordered
# compass octants. Used by func/driver_interactions.R::build_driver_matrix()
# to give wind_dir/wave_dir a categorical form (Robert: direction matters as
# much or more than magnitude for these two drivers, so they should enter the
# GLM/GAM/RF as a factor, not a raw degree value) -- kept separate from
# plot_driver_rose()'s own inline sector binning below since that one uses a
# finer, plot-resolution-tuned n_sectors (16) for bar drawing, not a fixed
# categorical driver definition.
compass_octant <- function(degrees){
  labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  factor(labels[(round(degrees / 45) %% 8) + 1], levels = labels)
}

# Direction/magnitude rose for wind or wave, coloured by the flow-controlled
# plume-area response -- manuscript Figure 7. Replaces plotting wind/wave as
# a raw magnitude time series (the original Figs 7-9 concept): Robert's
# point was that wind and wave direction matter as much or more than
# magnitude, which a time series of speed/height alone can't show. Bar
# length/frequency = how often the driver comes from that compass sector;
# fill colour = mean plume-area residual (after removing flow's effect via
# a plain lm(plume_area ~ flow), same logic as
# rhone_wind_wave_effect()'s area_resid) for days in that sector --
# so the rose shows both the driver's climatology and whether the plume
# responds differently when it comes from a given direction, in one plot.
# plot_driver_rose("wind", get_zone_meta(mouth_name = "Grand Rhone"))
# Residual of plume_area after removing flow's linear effect (df must have
# both columns) -- the "what's left once discharge is accounted for" signal
# shared by plot_driver_rose() and plot_category_scatter() below.
flow_controlled_residual <- function(df){
  residuals(lm(plume_area ~ flow, data = df))
}

# df_flow (columns date/plume_area/flow, as combine_plume_driver("flow", meta)
# |> select(date, plume_area, flow = value) produces) lets a caller that
# already has it -- e.g. Figure_7_driver_rose(), calling this twice per zone
# for wind and wave -- pass it in once rather than recomputing the same
# join/STL for every call; defaults to computing it here for standalone use.
plot_driver_rose <- function(driver_name, meta, n_sectors = 16, df_flow = NULL){
  driver_name <- match.arg(driver_name, c("wind", "wave"))
  dir_col <- paste0(driver_name, "_dir")
  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)

  if(is.null(df_flow)){
    df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)
  }
  df_driver <- load_driver(driver_name, meta)

  df <- df_flow |>
    dplyr::left_join(df_driver, by = "date") |>
    tidyr::drop_na(plume_area, flow, dplyr::all_of(dir_col))

  # Defensive fallback for a zone/driver missing direction entirely (see
  # load_wave()'s NB 2 -- this specifically affected the Bay of Seine's wave
  # direction until the 2026-07-31 WAVE re-download added VMDR there too;
  # kept here in case of a future data regression) -- drop_na() above would
  # empty df entirely in that case. A blank rose would be misleading (looks
  # like "no data collected" rather than "direction not available"), so say
  # so explicitly instead.
  if(nrow(df) == 0){
    pl <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = paste0(toupper(driver_name), " direction\nnot available\nfor this zone"),
               size = 5, colour = "grey40") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = paste0(zone_title(meta$zone), ": ", tolower(disp$driver_label))) +
      theme_void() + theme(plot.title = element_text(hjust = 0.5))

    ggsave(filename = paste0("figures/driver_comparison/rose_", driver_name, "_plume_", meta$mouth_name, ".png"),
          plot = pl, width = 7, height = 6, dpi = 300)
    return(invisible(pl))
  }

  df$area_resid <- flow_controlled_residual(df)

  sector_width <- 360 / n_sectors
  df$sector <- (round(df[[dir_col]] / sector_width) %% n_sectors) * sector_width

  df_summary <- df |>
    dplyr::summarise(n_days = dplyr::n(),
                     mean_area_resid = mean(area_resid, na.rm = TRUE), .by = "sector")

  compass_breaks <- seq(0, 360 - sector_width, by = max(sector_width, 45))
  compass_labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")[seq_along(compass_breaks)]

  pl <- ggplot(df_summary, aes(x = sector, y = n_days, fill = mean_area_resid)) +
    geom_col(width = sector_width * 0.9, colour = "grey30", linewidth = 0.2) +
    coord_polar(start = -(sector_width / 2) * pi / 180) +
    scale_x_continuous(breaks = compass_breaks, labels = compass_labels, limits = c(0, 360)) +
    scale_fill_gradient2(low = "steelblue", mid = "grey90", high = "firebrick", midpoint = 0,
                        name = "Plume-area\nresidual (km²)") +
    labs(x = NULL, y = "Days", title = paste0(zone_title(meta$zone), ": ", tolower(disp$driver_label))) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = "grey85"), plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("figures/driver_comparison/rose_", driver_name, "_plume_", meta$mouth_name, ".png"),
        plot = pl, width = 7, height = 6, dpi = 300)
  invisible(pl)
}

# Flow-controlled plume-area residual vs. wave height, coloured by on/off-
# shore wind category -- manuscript Figure 8. Generalises
# rhone_wind_wave_effect()'s pl_wave panel (Gulf of Lion only, with a
# bespoke Mistral/onshore/calm classification tuned to that zone's
# geography -- left as-is above, since it's a specific side-study, not
# replaced) to all four zones, using the generic on/off-shore
# classification wind_add_direction() already computes for every zone. Uses
# wave HEIGHT only (not direction), so this works for the Bay of Seine too,
# despite that zone's missing VMDR (see load_wave()'s NB 2).
# plot_category_scatter(get_zone_meta(mouth_name = "Grand Rhone"))
plot_category_scatter <- function(meta){
  df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, wind_spd = value, direction)
  df_wave <- load_driver("wave", meta) |> dplyr::select(date, wave_height = value)

  df <- df_flow |>
    dplyr::left_join(df_wind, by = "date") |>
    dplyr::left_join(df_wave, by = "date") |>
    tidyr::drop_na(plume_area, flow, wind_spd, direction, wave_height)
  df$area_resid <- flow_controlled_residual(df)

  df$wind_category <- dplyr::case_when(
    df$wind_spd < 3       ~ "calm (<3 m/s)",
    df$direction == "off" ~ "offshore",
    TRUE                  ~ "onshore"
  )
  df$wind_category <- factor(df$wind_category, levels = c("calm (<3 m/s)", "onshore", "offshore"))
  category_colours <- c("calm (<3 m/s)" = "grey50", "onshore" = "steelblue", "offshore" = "firebrick")

  ggplot(df, aes(x = wave_height, y = area_resid, colour = wind_category)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    scale_colour_manual(values = category_colours) +
    labs(x = "Wave height (m)", y = "Flow-controlled plume-area residual (km²)",
        colour = NULL, title = zone_title(meta$zone)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"), legend.position = "bottom")
}

# AR(1)-weighted / STL-weighted / unweighted linear trend + Newey-West (HAC)
# standard error correction, following the monthly time series adjustment
# methodology of Sutton et al. (2022) referenced in manuscript.tex Sec. 2.6.1.
# Extracted from driver_plume_trend() (below) so it can be reused directly on
# a bare (date, value) series -- e.g. plume shape or centroid drift -- without
# needing a paired driver series. driver_plume_trend() calls this internally;
# behaviour there is unchanged.
# fit_wls_hac_trend("ar", df$compactness, df$date)
ar_weights_func <- function(val_col, start_year, time_step){
  ts_obj <- ts(zoo::na.approx(val_col), frequency = time_step, start = c(start_year, 1))
  ar_model <- ar(ts_obj, order.max = 1)
  phi_est <- ar_model$ar
  sigma_est <- sqrt(phi_est)
  error_variance <- sigma_est^2 / (1 - phi_est^2)
  rep((1 / (error_variance^2)), length(val_col))
}
stl_weights_func <- function(val_col, start_year, time_step){
  ts_obj <- ts(zoo::na.approx(val_col), frequency = time_step, start = c(start_year, 1))
  stl_ts <- stl(ts_obj, s.window = "periodic")
  stl_var <- as.vector(stl_ts$time.series[, "remainder"])
  1 / (stl_var^2)
}
fit_wls_hac_trend <- function(weight_choice, val_col, date_col){
  start_year <- year(min(date_col))
  time_step <- if(length(val_col) < 1000) 12 else 365
  weights <- switch(weight_choice,
                    ar = ar_weights_func(val_col, start_year, time_step),
                    stl = stl_weights_func(val_col, start_year, time_step),
                    rep(1, length(val_col)))
  lm_model <- lm(val_col ~ date_col, weights = weights)
  lm_model_HAC <- coeftest(lm_model, vcov = vcovHAC(lm_model))
  tibble::tibble(n = length(val_col), time_step = time_step, start_year = start_year,
                 weight_choice = weight_choice, intercept = lm_model_HAC[1, 1],
                 slope = lm_model_HAC[2, 1], slope_se = lm_model_HAC[2, 2],
                 slope_t = lm_model_HAC[2, 3], slope_p = lm_model_HAC[2, 4])
}

# Generalises flow.R::flow_plume_trend_plus(), which implemented this for
# river flow only; used here for any driver.
# driver_plume_trend(combine_plume_driver("flow", get_zone_meta(mouth_name = "Grand Rhone")), "flow", "Grand Rhone")
# save_plot = FALSE skips the ggsave below -- needed when `df`'s "plume_area" column
# actually holds a different metric_col (see combine_plume_driver()), since the plot
# filename is keyed only on driver_name/mouth_name and would otherwise silently
# overwrite the existing plume-area trend figure for that driver/mouth.
driver_plume_trend <- function(df, driver_name, mouth_name, end_date = NULL, save_plot = TRUE){

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)

  df <- df |>
    dplyr::mutate(doy = yday(date), month = month(date), year = year(date))
  if(!is.null(end_date)) df <- dplyr::filter(df, date <= as.Date(end_date))

  # De-trend to get day-of-year and monthly climatological adjustments
  df_resid <- df |>
    dplyr::mutate(driver_resid = residuals(lm(value ~ date, na.action = na.exclude)),
                  plume_resid  = residuals(lm(plume_area ~ date, na.action = na.exclude)))

  df_doy_clim <- df_resid |>
    dplyr::summarise(driver_resid_doy = mean(driver_resid, na.rm = TRUE),
                      plume_resid_doy = mean(plume_resid, na.rm = TRUE), .by = "doy") |>
    dplyr::mutate(driver_resid_doy_clim = driver_resid_doy - mean(driver_resid_doy, na.rm = TRUE),
                  plume_resid_doy_clim  = plume_resid_doy - mean(plume_resid_doy, na.rm = TRUE))
  df_month_clim <- df_resid |>
    dplyr::summarise(driver_resid_monthly = mean(driver_resid, na.rm = TRUE),
                      plume_resid_monthly = mean(plume_resid, na.rm = TRUE), .by = "month") |>
    dplyr::mutate(driver_resid_monthly_clim = driver_resid_monthly - mean(driver_resid_monthly, na.rm = TRUE),
                  plume_resid_monthly_clim  = plume_resid_monthly - mean(plume_resid_monthly, na.rm = TRUE))

  df_daily <- df |>
    dplyr::left_join(df_doy_clim, by = "doy") |>
    dplyr::mutate(driver_doy_adj = value - driver_resid_doy_clim,
                  plume_doy_adj  = plume_area - plume_resid_doy_clim) |>
    dplyr::mutate(date_int = seq_len(dplyr::n()), .after = "date")

  df_monthly <- df |>
    dplyr::mutate(date = floor_date(date, "month")) |>
    dplyr::summarise(driver_monthly = mean(value, na.rm = TRUE),
                      plume_monthly = mean(plume_area, na.rm = TRUE), .by = c("year", "month", "date")) |>
    dplyr::left_join(df_month_clim, by = "month") |>
    dplyr::mutate(driver_monthly_adj = driver_monthly - driver_resid_monthly_clim,
                  plume_monthly_adj  = plume_monthly - plume_resid_monthly_clim) |>
    dplyr::mutate(date_int = seq_len(dplyr::n()), .after = "date")

  wls_driver_daily   <- plyr::ldply(c("ar", "stl", "none"), fit_wls_hac_trend, val_col = df_daily$driver_doy_adj, date_col = df_daily$date)
  wls_driver_monthly <- plyr::ldply(c("ar", "stl", "none"), fit_wls_hac_trend, val_col = df_monthly$driver_monthly_adj, date_col = df_monthly$date)
  wls_plume_daily    <- plyr::ldply(c("ar", "stl", "none"), fit_wls_hac_trend, val_col = df_daily$plume_doy_adj, date_col = df_daily$date)
  wls_plume_monthly  <- plyr::ldply(c("ar", "stl", "none"), fit_wls_hac_trend, val_col = df_monthly$plume_monthly_adj, date_col = df_monthly$date)

  stats <- dplyr::bind_rows(
    dplyr::mutate(wls_driver_daily,   variable = "driver", timestep = "daily"),
    dplyr::mutate(wls_driver_monthly, variable = "driver", timestep = "monthly"),
    dplyr::mutate(wls_plume_daily,    variable = "plume",  timestep = "daily"),
    dplyr::mutate(wls_plume_monthly,  variable = "plume",  timestep = "monthly")
  ) |>
    dplyr::mutate(driver_name = driver_name, mouth_name = mouth_name,
                  slope_annualised = dplyr::case_when(timestep == "daily" ~ slope * 365.25,
                                                      timestep == "monthly" ~ slope * 365.25,
                                                      TRUE ~ slope), .before = "n")

  if(save_plot){
  # Plot (daily = ar-weighted line; monthly = ar-weighted line), labelled with slope + p-value
  trend_labels_driver <- dplyr::filter(stats, variable == "driver", weight_choice == "ar")
  trend_labels_plume  <- dplyr::filter(stats, variable == "plume",  weight_choice == "ar")
  x_daily <- min(df_daily$date) + days(round(0.05 * as.numeric(diff(range(df_daily$date)))))
  x_monthly <- min(df_daily$date) + days(round(0.45 * as.numeric(diff(range(df_daily$date)))))
  y_plume  <- round(max(df_daily$plume_area, na.rm = TRUE) - stats::quantile(df_daily$plume_area, 0.2, na.rm = TRUE), -2)
  y_driver <- round(max(df_daily$value, na.rm = TRUE) - stats::quantile(df_daily$value, 0.05, na.rm = TRUE), -2)

  pl_plume <- ggplot(data = df_daily, aes(x = date, y = plume_area)) +
    geom_point(colour = "sienna", alpha = 0.1) +
    geom_point(aes(y = plume_doy_adj), colour = "darkblue", alpha = 0.1) +
    geom_point(data = df_monthly, aes(y = plume_monthly), colour = "sienna", alpha = 0.6, size = 3) +
    geom_point(data = df_monthly, aes(y = plume_monthly_adj), colour = "darkred", size = 3) +
    geom_abline(data = dplyr::filter(trend_labels_plume, timestep == "daily"),
                aes(intercept = intercept, slope = slope), linewidth = 2, colour = "darkblue") +
    geom_abline(data = dplyr::filter(trend_labels_plume, timestep == "monthly"),
                aes(intercept = intercept, slope = slope), linewidth = 2, colour = "darkred") +
    geom_label(data = dplyr::filter(trend_labels_plume, timestep == "monthly"), size = 5, hjust = 0, colour = "darkred",
               aes(x = x_daily, y = y_plume, label = paste0("Plume area slope = ", round(slope_annualised, 2), " km^2 yr-1\n",
                                                            "p-value = ", round(slope_p, 2)))) +
    geom_label(data = dplyr::filter(trend_labels_plume, timestep == "daily"), size = 5, hjust = 0, colour = "darkblue",
               aes(x = x_monthly, y = y_plume, label = paste0("Plume area slope = ", round(slope_annualised, 2), " km^2 yr-1\n",
                                                              "p-value = ", round(slope_p, 2)))) +
    labs(x = NULL, y = "Plume area [km^2]",
         title = paste0(mouth_name, " : plume area after statistical treatment (vs. ", driver_name, ")"),
         subtitle = "Red = adjusted monthly values; blue = adjusted daily values; brown = original data") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  pl_driver <- ggplot(data = df_daily, aes(x = date, y = value)) +
    geom_point(colour = "purple", alpha = 0.1) +
    geom_point(aes(y = driver_doy_adj), colour = "darkblue", alpha = 0.1) +
    geom_point(data = df_monthly, aes(y = driver_monthly), colour = "purple", alpha = 0.6, size = 3) +
    geom_point(data = df_monthly, aes(y = driver_monthly_adj), colour = "darkred", size = 3) +
    geom_abline(data = dplyr::filter(trend_labels_driver, timestep == "daily"),
                aes(intercept = intercept, slope = slope), linewidth = 2, colour = "darkblue") +
    geom_abline(data = dplyr::filter(trend_labels_driver, timestep == "monthly"),
                aes(intercept = intercept, slope = slope), linewidth = 2, colour = "darkred") +
    geom_label(data = dplyr::filter(trend_labels_driver, timestep == "monthly"), size = 5, hjust = 0, colour = "darkred",
               aes(x = x_daily, y = y_driver, label = paste0(disp$driver_label, " slope = ", round(slope_annualised, 2), " yr-1\n",
                                                             "p-value = ", round(slope_p, 2)))) +
    geom_label(data = dplyr::filter(trend_labels_driver, timestep == "daily"), size = 5, hjust = 0, colour = "darkblue",
               aes(x = x_monthly, y = y_driver, label = paste0(disp$driver_label, " slope = ", round(slope_annualised, 2), " yr-1\n",
                                                               "p-value = ", round(slope_p, 2)))) +
    labs(x = NULL, y = disp$driver_label,
         title = paste0(mouth_name, " : ", driver_name, " after statistical treatment"),
         subtitle = "Red = adjusted monthly values; blue = adjusted daily values; purple = original data") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  pl_combi <- ggpubr::ggarrange(pl_plume, pl_driver, ncol = 1, nrow = 2)
  ggsave(filename = paste0("figures/driver_comparison/trends_plume_", driver_name, "_adj_", mouth_name, ".png"),
         pl_combi, width = 12, height = 10)
  }

  return(stats)
}


# Runners ---------------------------------------------------------------

# Run the full comparison + trend suite for one driver across all four river
# mouths, returning the combined correlation and trend statistics (these feed
# manuscript.tex Table 6 "driver_stats" once finalised). This is the direct
# replacement for the old:
#   plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_comp)
#   plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_trend)
#   plyr::d_ply(.data = river_mouths, .variables = "row_name", .fun = flow_plume_trend_plus, .parallel = TRUE)
# (and the wind.R / tide.R / ROFI.R equivalents), now driver-agnostic.
# run_driver_suite("flow")
run_driver_suite <- function(driver_name){

  mouths <- if(driver_name == "rofi") dplyr::filter(zone_meta, zone != "GULF_OF_LION") else zone_meta  # no ROFI data for the Gulf of Lion

  results <- purrr::pmap(mouths, function(...){
    meta <- tibble::tibble(...)
    df <- combine_plume_driver(driver_name, meta)

    plot_driver_comparison(df, driver_name, meta$mouth_name)
    plot_driver_plume_dual_axis(df, driver_name, meta$zone)
    cor_stats <- driver_plume_correlation(df) |> dplyr::mutate(mouth_name = meta$mouth_name, zone = meta$zone, driver_name = driver_name)
    trend_stats <- driver_plume_trend(df, driver_name, meta$mouth_name)

    list(correlation = cor_stats, trend = trend_stats)
  })

  list(
    correlation = purrr::map_dfr(results, "correlation"),
    trend = purrr::map_dfr(results, "trend")
  )
}

# Run all drivers for all zones. NB: this recreates every figure previously
# produced by running flow.R + wind.R + tide.R + ROFI.R in sequence, plus the
# new GLORYS current driver; it is not run automatically on source() (see
# bottom of file).
run_all_driver_suites <- function(){
  driver_names <- c("flow", "tide", "wind", "current", "rofi", "wave")
  purrr::map(driver_names, run_driver_suite) |>
    purrr::set_names(driver_names)
}


# Surface / pixel-level multi-driver maps ------------------------------
# Facet daily plume maps by year x month for a zone.
surface_plot_daily_maps <- function(zone_name){

  df_plume <- plyr::ldply(zone_name, load_plume_surface, .parallel = FALSE) |>  # NB: do not run in parallel, it's done elsewhere
    dplyr::mutate(year = year(date), month = month(date), doy = yday(date), day = day(date))

  plot_daily <- ggplot(df_plume, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = day), alpha = 0.3) +
    scale_fill_viridis_c(option = "A", na.value = "transparent") +
    labs(x = NULL, y = NULL, fill = "Day of month", title = paste0("Daily plume maps per month and year for ", zone_name)) +
    facet_grid(year ~ month) + theme_bw() + theme(legend.position = "bottom")
  ggsave(filename = paste0("figures/driver_comparison/surface_daily_maps_", zone_name, ".png"), plot = plot_daily, height = 34, width = 36)
  invisible(plot_daily)
}

# Plot all surface daily maps
# NB: Tgis takes a while and is pretty heavy
# walk(zones, surface_plot_daily_maps)


# STL ---------------------------------------------------------------------

# Load all plume and driver data and perform stl. Refactored to route through
# combine_plume_driver()-style loading (still bespoke here because this
# function needs all four drivers joined at once, not one at a time).
# zone <- zones[4]
multi_stl <- function(zone){

  meta <- get_zone_meta(zone_name = zone)

  df_plume <- load_plume_ts(zone) |>
    dplyr::rename(SPM_threshold = tidyselect::last_col())  # standardise last column name so it combines across sites

  df_river_flow <- load_driver("flow", meta) |> dplyr::rename(flow = value)
  df_tide <- load_driver("tide", meta) |> dplyr::rename(tide_range = value)
  df_wind_full <- load_driver("wind", meta) |> dplyr::rename(wind_spd = value)

  df_all <- dplyr::left_join(df_plume, df_river_flow, by = "date") |>
    dplyr::left_join(df_tide, by = "date") |>
    dplyr::left_join(df_wind_full, by = "date") |>
    dplyr::mutate(plume_seas = stl_single(plume_area, out_col = "seas", start_date = min(df_plume$date)),
                  plume_inter = stl_single(plume_area, out_col = "inter", start_date = min(df_plume$date)),
                  plume_resid = stl_single(plume_area, out_col = "remain", start_date = min(df_plume$date)),
                  flow_seas = stl_single(flow, out_col = "seas", start_date = min(df_plume$date)),
                  flow_inter = stl_single(flow, out_col = "inter", start_date = min(df_plume$date)),
                  flow_resid = stl_single(flow, out_col = "remain", start_date = min(df_plume$date)),
                  tide_seas = stl_single(tide_range, out_col = "seas", start_date = min(df_plume$date)),
                  tide_inter = stl_single(tide_range, out_col = "inter", start_date = min(df_plume$date)),
                  tide_resid = stl_single(tide_range, out_col = "remain", start_date = min(df_plume$date)),
                  wind_seas = stl_single(wind_spd, out_col = "seas", start_date = min(df_plume$date)),
                  wind_inter = stl_single(wind_spd, out_col = "inter", start_date = min(df_plume$date)),
                  wind_resid = stl_single(wind_spd, out_col = "remain", start_date = min(df_plume$date))) |>
    dplyr::mutate(zone = zone, .before = "date")

  return(df_all)
}

# Compute all STL stats and save
if(!file.exists("output/STATS/stl_all.RData")){
  message("Computing STL stats for all zones...")
  stl_all <- plyr::ldply(zones, multi_stl, .parallel = TRUE)
  save(stl_all, file = "output/STATS/stl_all.RData")
}


# Multi-driver comparison -------------------------------------------------

# Plot the results
# df_stl <- stl_all
multi_plot <- function(df_stl){

  # Make pretty plot titles -- now via util.R::make_pretty_title() instead of
  # an inline case_when() duplicate of the same logic used elsewhere in this
  # file and in ROFI.R.
  df_pretty <- make_pretty_title(df_stl)

  # One year of data for seasonal plots
  df_mean <- df_pretty |>
    summarise(plume_mean = mean(plume_inter, na.rm = TRUE),
              flow_mean = mean(flow_inter, na.rm = TRUE),
              wind_mean = mean(wind_inter, na.rm = TRUE),
              tide_mean = mean(tide_inter, na.rm = TRUE), .by = c(zone, plot_title))
  df_seas <- df_pretty |>
    filter(year(date) == 1999) |>
    mutate(month = month(date, label = TRUE, abbr = TRUE),
           doy = yday(date)) |>
    dplyr::select(zone, plot_title, month, doy, plume_seas, flow_seas, tide_seas, wind_seas) |>
    distinct() |>
    left_join(df_mean, by = c("zone", "plot_title")) |>
    mutate(plume_seas = plume_seas + plume_mean,
           flow_seas = flow_seas + flow_mean,
           tide_seas = tide_seas + tide_mean,
           wind_seas = wind_seas + wind_mean)

  # Convenience wrappers for daily, seasonal, and interannual plot
  plot_daily <- function(df, y_col, line_colour, y_label, file_stub){
    unique_years <- df$date |> year() |> unique()
    pl_daily <- ggplot(data = df) +
      geom_path(aes_string(x = "date", y = y_col), color = line_colour) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_date(name = "", expand = c(0,0),
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
      scale_y_continuous(name = y_label) +
      labs( x = NULL) +
      ggplot_theme()
    ggsave(filename = paste0("figures/",file_stub,"_daily.png"), plot = pl_daily, width = 24, height = 24, dpi = 300)
  }
  plot_seas <- function(df, y_col, line_colour, y_label, file_stub){
    df_sub <- df[,c("plot_title", "month", y_col)]
    colnames(df_sub)[3] = "val"
    df_sub <- df_sub |>
      summarise(val_min = min(val, na.rm = TRUE),
                val_mean = mean(val, na.rm = TRUE),
                val_max = max(val, na.rm = TRUE), .by = c("plot_title", "month")) |>
      mutate(month_int = as.integer(month))
    pl_seas <- ggplot(data = df_sub, aes(x = month_int)) +
      geom_ribbon(aes(ymin = val_min, ymax = val_max), fill = line_colour, alpha = 0.3) +
      geom_path(aes(y = val_mean), color = line_colour, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_y_continuous(name = y_label) +
      scale_x_continuous(expand = c(0, 0), breaks = 1:12, labels = month.abb) +
      labs(x = NULL) +
      ggplot_theme()
    ggsave(filename = paste0("figures/",file_stub,"_seas.png"), plot = pl_seas, width = 24, height = 24, dpi = 300)
  }
  plot_inter <- function(df, y_col, line_colour, y_label, file_stub){
    unique_years <- df$date |> year() |> unique()
    colnames(df)[which(colnames(df) == y_col)] <- "value"
    df_sub <- df |>
      dplyr::select(plot_title, date, value) |>
      mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
      filter(date >= min(df$date)) |>
      group_by(plot_title, date) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "keep") |>
      group_by(plot_title) |>
      mutate(running_mean = roll_mean(value, n = 48, fill = NA, align = "center")) |>
      ungroup()
    pl_inter <- ggplot(data = df_sub) +
      geom_path(aes(x = date, y = running_mean), color = line_colour, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_date(name = "",
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
      scale_y_continuous(name = y_label) +
      ggplot_theme()
    ggsave(filename = paste0("figures/",file_stub,"_inter.png"), plot = pl_inter, width = 24, height = 24, dpi = 300)
  }

  # Daily time series
  plot_daily(df_pretty, "plume_area", "brown", "Plume area (km^2)", "driver_comparison/plume")

  # Seasonal time series
  plot_seas(df_seas, "plume_seas", "brown", "Plume area (km^2)", "driver_comparison/plume")
  plot_seas(df_seas, "flow_seas", "blue", "River flow (m^3 s-1)", "driver_comparison/flow")
  plot_seas(df_seas, "tide_seas", "darkgreen", "Tidal range (m)", "driver_comparison/tide")
  plot_seas(df_seas, "wind_seas", "purple", "Wind speed (m s-1)", "driver_comparison/wind")

  # Interannual time series
  plot_inter(df_pretty, "plume_inter", "brown", "Plume area (km^2)", "driver_comparison/plume")
  plot_inter(df_pretty, "flow_inter", "blue", "River flow (m^3 s-1)", "driver_comparison/flow")
  plot_inter(df_pretty, "tide_inter", "darkgreen", "Tidal range (m)", "driver_comparison/tide")
  plot_inter(df_pretty, "wind_inter", "purple", "Wind speed (m s-1)", "driver_comparison/wind")

  # Seasonal comparison plots
  comparison_plot_save(df_seas, "plume_seas", "flow_seas", "brown", "blue", "Plume area (km^2)", "River flow (m^3 s-1)", "driver_comparison/comparison_plume_flow_seas")
  comparison_plot_save(df_seas, "plume_seas", "wind_seas", "brown", "purple", "Plume area (km^2)", "Wind speed (m s-1)", "driver_comparison/comparison_plume_wind_seas")
  comparison_plot_save(df_seas, "plume_seas", "tide_seas", "brown", "darkgreen", "Plume area (km^2)", "Tidal range (m)", "driver_comparison/comparison_plume_tide_seas")

  # Interannual comparison plots
  comparison_plot_save(df_pretty, "plume_inter", "flow_inter", "brown", "blue", "Plume area (km^2)", "River flow (m^3 s-1)", "driver_comparison/comparison_plume_flow_inter")
  comparison_plot_save(df_pretty, "plume_inter", "tide_inter", "brown", "darkgreen", "Plume area (km^2)", "Tidal range (m)", "driver_comparison/comparison_plume_tide_inter")
  comparison_plot_save(df_pretty, "plume_inter", "wind_inter", "brown", "purple", "Plume area (km^2)", "Wind speed (m s-1)", "driver_comparison/comparison_plume_wind_inter")

  # Everything on one plot
  df_all_scaled <- df_pretty |>
    group_by(plot_title) |>
    mutate(plum_scaled = plume_inter/max(plume_inter, na.rm = TRUE),
           flow_scaled = flow_inter/max(flow_inter, na.rm = TRUE),
           tide_scaled = tide_inter/max(tide_inter, na.rm = TRUE),
           wind_scaled = wind_inter/max(wind_inter, na.rm = TRUE)) |>
    mutate(plum_scaled = plum_scaled/mean(plum_scaled, na.rm = TRUE),
           flow_scaled = flow_scaled/mean(flow_scaled, na.rm = TRUE),
           tide_scaled = tide_scaled/mean(tide_scaled, na.rm = TRUE),
           wind_scaled = wind_scaled/mean(wind_scaled, na.rm = TRUE)) |>
    dplyr::select(plot_title, date, plum_scaled:wind_scaled) |>
    pivot_longer(plum_scaled:wind_scaled) |>
    mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
    filter(date >= min(df_pretty$date)) |>
    group_by(plot_title, name, date) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "keep") |>
    group_by(plot_title, name) |>
    mutate(running_mean = roll_mean(value, n = 48, fill = NA, align = "center")) |>
    ungroup()

  all_plot <- ggplot(df_all_scaled, aes(x = date, y = running_mean)) +
    geom_path(aes(colour = name), linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1) +
    ggplot_theme()
  ggsave(filename = "figures/driver_comparison/all_plot.png", plot = all_plot, width = 20, height = 20, dpi = 300)
}

# Load STL calculated above
# load("output/STATS/stl_all.RData")

# Create plots
# multi_plot(stl_all)


# Missing data ------------------------------------------------------------

# Get missing dates of
if(!file.exists("output/STATS/missing_SPM.csv") | !file.exists("output/STATS/missing_chla.csv")){
  message("Computing missing SEXTANT files...")
  SPM_files_NA <- data.frame(file_name = dir("~/pCloudDrive/data/SEXTANT/SPM/", pattern = ".nc", recursive = TRUE)) |>
    mutate(base_name = basename(file_name)) |>
    separate(base_name, "-", extra = "drop") |>
    dplyr::rename(date = `-`) |>
    mutate(date = as.Date(date, format = "%Y%m%d")) |>
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |>
    filter(is.na(file_name))
  write_csv(SPM_files_NA, "output/STATS/missing_SPM.csv")
  chla_files_NA <- data.frame(file_name = dir("~/pCloudDrive/data/SEXTANT/CHLA/", pattern = ".nc", recursive = TRUE)) |>
    mutate(base_name = basename(file_name)) |>
    separate(base_name, "-", extra = "drop") |>
    dplyr::rename(date = `-`) |>
    mutate(date = as.Date(date, format = "%Y%m%d")) |>
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |>
    filter(is.na(file_name))
write_csv(chla_files_NA, "output/STATS/missing_chla.csv")
  
  # Filter down to missing days
  SPM_files_NA_count <- SPM_files_NA |>
    mutate(year = year(date),
          month = month(date, label = TRUE, abbr = TRUE)) |>
    summarise(miss_count_month_year = n(), .by = c("year", "month"))
  chla_files_NA_count <- chla_files_NA |>
    mutate(year = year(date),
          month = month(date, label = TRUE, abbr = TRUE)) |>
    summarise(miss_count_month_year = n(), .by = c("year", "month"))

  # Plot
  ggplot(SPM_files_NA_count, aes(x = month, y = miss_count_month_year)) +
    geom_col() +
    facet_wrap(~year) +
    labs(x = NULL, y = "count", title = "Monthly count of missing SPM SEXTANT files") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  ggsave("figures/validation/missng_SPM.png", width = 9, height = 9, dpi = 600)
  ggplot(chla_files_NA_count, aes(x = month, y = miss_count_month_year)) +
    geom_col() +
    facet_wrap(~year) +
    labs(x = NULL, y = "count", title = "Monthly count of missing chl a SEXTANT files") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  ggsave("figures/validation/missng_chla.png", width = 9, height = 9, dpi = 600)
}


# Run everything -----------------------------------------------------------

# NB: not run automatically on source() -- call explicitly
# run_all_driver_suites()
# purrr::walk(zones, surface_plot_daily_maps)


# Rhone only analyses ------------------------------------------------------

# Response to Claude's e-mail on the Grand Rhone plume-area
# trend. The e-mail raised four ideas; disposition of each below:
#   1. Concentration-detrending sensitivity test -- IMPLEMENTED,
#      rhone_detrend_test().
#   2. Has the balance of northern-tributary vs Cevennes ("cevenol")
#      floods shifted, explaining the concentration rise? -- PARTIALLY
#      implemented, rhone_flood_timing_shift(). RiOMar only has
#      combined Grand+Petit Rhone discharge at the mouth, not upstream
#      sub-basin gauges, so this cannot attribute the shift to a specific
#      tributary -- that needs the Observatoire des Sediments du Rhone
#      data / Olivier Radakovitch the e-mail points to. What's implemented
#      is a first-pass check of whether the *calendar timing* of high-flow
#      days has moved, using data already in the repo.
#   3. Does wind/wave forcing move the plume beyond what flow explains
#      (reframing "does the plume have its own seasonality" as "does it
#      respond to wind/wave directly"), including the specific Mistral vs
#      onshore/calm hypothesis -- IMPLEMENTED, rhone_wind_wave_effect().
#   4. Would a more stratified surface layer let the plume slide further
#      offshore? -- NOT implemented. The e-mail itself frames this as a
#      question for idealised numerical modelling of academic cases, not
#      something to test against the RiOMar time series, and RiOMar has no
#      surface density/stratification product to test it against even if
#      it were. Left as a discussion point for the reply, not a script.
#
# All three functions below are GULF_OF_LION / Grand Rhone only, as
# requested -- they are not generalised to the other three zones.

## 1. Concentration-detrending sensitivity test -------------------------------
# The plume-area trend could be, at least partly, an artefact of the SPM
# concentration near the Rhone mouth having risen over the record (see
# output/STATS/mass_SPM_trend_summary.csv and the earlier e-mail's claim)
# rather than the plume physically growing. The e-mail's proposed check:
# imagine the concentration had NO trend (a "trend-free" counterfactual,
# built "at least as a first approximation" by flattening the trend line),
# recompute the plume area under that counterfactual, and see whether the
# area trend survives.
#
# Re-running the full pixel-level panache detection under a hypothetical,
# de-trended SPM field is out of scope here -- it needs the raw daily SPM
# maps (not exposed to R) and a re-run of the Python flood-fill algorithm
# (func/plume.py::find_SPM_threshold()) for every day. An earlier version of
# this function used panache's own daily SPM_threshold (Results.csv) as a
# stand-in for concentration, but that is panache's classification cutoff
# for that day (func/plume.py:1254, `SPM_criterion = ds_reduced >
# SPM_threshold`) -- a stricter cutoff mechanically shrinks the classified
# area regardless of any real dilution physics, which confounds exactly the
# effect this test is trying to isolate. Instead, this uses real, independent
# in-situ SPM concentration measured at Arles on the Rhone itself
# (data/INSITU_data/OSR/ARLES_CMES-2.txt, hourly, mg/L, 2005-2023), from the
# Observatoire des Sediments du Rhone (OSR) -- the exact data source Claude's
# e-mail pointed to. Only quality flag "v" (validated/measured) rows are
# kept; flag "e" rows are the OSR's own gap-fill, estimated from a fixed
# discharge-to-concentration rating curve (see
# data/INSITU_data/OSR/Report_ARLES_CMES-2.txt), which would reintroduce a
# collinearity with flow that using real in-situ data is meant to avoid.
# Because the OSR record only covers 2005-2023, the joined record below is
# shorter than the full 1998-2025 satellite plume record.
# rhone_detrend_test()
rhone_detrend_test <- function(){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")

  # Daily-mean in-situ SPM concentration at Arles, validated readings only.
  df_conc <- read_delim("data/INSITU_data/OSR/ARLES_CMES-2.txt", delim = ";", skip = 3,
                        col_names = c("datetime", "value", "quality", "min", "max"),
                        col_types = "cdccc", show_col_types = FALSE) |>
    dplyr::filter(quality == "v") |>
    dplyr::mutate(date = as.Date(datetime, format = "%d/%m/%Y %H:%M:%S")) |>
    dplyr::summarise(conc = mean(value, na.rm = TRUE), .by = "date")

  # Standard plume-area + flow object used by every other driver in this
  # file. inner_join (not left_join) restricts the record to where the
  # satellite plume record and the OSR in-situ record overlap.
  df_flow <- combine_plume_driver("flow", meta)
  df <- df_flow |>
    dplyr::inner_join(df_conc, by = "date") |>
    zoo::na.trim()  # drops leading/trailing NA across plume_area, flow, and conc together

  # Interior gaps (e.g. days with no validated hourly readings) are linearly
  # interpolated so every remaining day has a value -- the same treatment
  # STL/trend fitting already gets elsewhere in this file (stl_single(),
  # ar_weights_func()).
  df$conc <- as.numeric(zoo::na.approx(df$conc, x = df$date, na.rm = FALSE))

  # a) Linear trend of the concentration proxy, and its "trend-free"
  #    counterfactual: the same series with the fitted slope subtracted
  #    back out (mean unchanged, so it stays on the original scale).
  conc_lm <- lm(conc ~ date, data = df)
  df$conc_fit <- predict(conc_lm)
  df$conc_detrended <- df$conc - (df$conc_fit - mean(df$conc_fit))

  # b) Empirical area ~ concentration sensitivity, controlling for flow so
  #    the concentration coefficient isn't just re-capturing flow's own,
  #    already-established effect on area (see driver_plume_trend("flow")).
  #    Caveat worth keeping in mind when reading beta_conc: flow and
  #    concentration are themselves correlated (bigger floods carry more
  #    sediment), so this is a partial, not fully independent, effect.
  area_conc_lm <- lm(plume_area ~ value + conc, data = df)
  beta_conc <- coef(area_conc_lm)["conc"]

  # c) Counterfactual area: subtract, day by day, the modelled contribution
  #    of that day's concentration sitting above (or below) its trend-free
  #    counterfactual.
  df$plume_area_conc_adj <- df$plume_area - beta_conc * (df$conc - df$conc_detrended)

  # d) Refit the trend on the raw area, the concentration-adjusted area,
  #    and the concentration proxy itself, all with the same AR-weighted /
  #    HAC estimator used throughout this file (fit_wls_hac_trend()), so
  #    the comparison is apples-to-apples with every other trend reported
  #    for this zone.
  trend_conc <- fit_wls_hac_trend("ar", df$conc, df$date) |> dplyr::mutate(series = "SPM concentration (Arles, OSR in-situ)")
  trend_raw  <- fit_wls_hac_trend("ar", df$plume_area, df$date) |> dplyr::mutate(series = "plume area (raw)")
  trend_adj  <- fit_wls_hac_trend("ar", df$plume_area_conc_adj, df$date) |> dplyr::mutate(series = "plume area (concentration-adjusted)")

  stats <- dplyr::bind_rows(trend_conc, trend_raw, trend_adj) |>
    dplyr::mutate(slope_yearly = slope * 365.25, beta_conc = beta_conc) |>
    dplyr::select(series, n, slope, slope_yearly, slope_se, slope_p, beta_conc)
  write_csv(stats, "output/STATS/rhone_conc_detrend.csv")

  # Plot: a) in-situ concentration with its fitted trend and trend-free
  # counterfactual; b) raw vs. concentration-adjusted plume area with their
  # respective AR-weighted trend lines, so the effect of the adjustment
  # (if any) is visible directly.
  # NB: conc is a high-variance daily series (floods vs. base flow), so
  # connecting every day with a line (geom_line) paints a near-solid band
  # that hides everything underneath it. Points + a clean geom_abline for
  # the fitted trend stays legible; the trend-free counterfactual is, at
  # this scale, visually indistinguishable from the raw series (the whole
  # record's drift is small next to the daily noise), so it is reported in
  # `stats` rather than over-plotted here.
  # SPM concentration is strongly right-skewed (flood spikes into the
  # thousands of mg/L vs. a base-flow median around 14 mg/L), so a linear
  # y-axis crushes almost every point flat against zero -- log10 is needed
  # to actually see the day-to-day/seasonal structure. The fitted trend
  # line itself is still the linear-scale fit used for the detrending math
  # above (conc_lm); log10() only changes how it is displayed here.
  # NB: geom_abline() does NOT work here -- it applies intercept/slope in
  # the already-log10-transformed coordinate space, not the original mg/L
  # space conc_lm was fit in, so the line would be drawn many orders of
  # magnitude off the visible panel. Instead, get the two endpoint values
  # by predicting conc_lm (linear space) at the start/end dates, then plot
  # those as ordinary data with geom_line() so they go through the same
  # log10 transform as the points.
  conc_trend_line <- tibble::tibble(date = range(df$date)) |>
    dplyr::mutate(conc_fit_line = predict(conc_lm, newdata = tibble::tibble(date = date)))

  pl_conc <- ggplot(df, aes(x = date)) +
    geom_point(aes(y = conc), colour = "grey50", alpha = 0.15, size = 0.5) +
    geom_line(data = conc_trend_line, aes(y = conc_fit_line), colour = "firebrick", linewidth = 1.2) +
    scale_y_log10() +
    labs(x = NULL, y = "SPM concentration, Arles (OSR in-situ, mg/L, log scale)",
         title = "Grand Rhone: in-situ SPM concentration at Arles",
         subtitle = "Red = fitted linear trend (see `stats` for the trend-free counterfactual comparison)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  fit_area_raw <- coef(lm(plume_area ~ date, data = df))
  fit_area_adj <- coef(lm(plume_area_conc_adj ~ date, data = df))

  pl_area <- ggplot(df, aes(x = date)) +
    geom_point(aes(y = plume_area), colour = "sienna", alpha = 0.15) +
    geom_point(aes(y = plume_area_conc_adj), colour = "darkblue", alpha = 0.15) +
    geom_abline(intercept = fit_area_raw[1], slope = fit_area_raw[2], colour = "sienna", linewidth = 1.2) +
    geom_abline(intercept = fit_area_adj[1], slope = fit_area_adj[2], colour = "darkblue", linewidth = 1.2) +
    labs(x = NULL, y = "Plume area (km^2)",
         title = "Grand Rhone: plume area, raw vs. concentration-adjusted",
         subtitle = "Brown = raw area (+ OLS trend); blue = area with the concentration-trend contribution removed (+ OLS trend)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  pl_combi <- ggpubr::ggarrange(pl_conc, pl_area, ncol = 1, nrow = 2)
  ggsave(filename = "figures/rhone_side_analyses/rhone_detrend_test.png", plot = pl_combi, width = 12, height = 10)

  return(list(stats = stats, data = df, plot = pl_combi))
}


## 2. Has the seasonal timing of Rhone floods shifted? ------------------------
# The e-mail's second idea asks *why* Rhone SPM concentration would be
# rising: has the balance between northern-tributary floods (winter/spring,
# snowmelt-and-rain driven) and Cevennes floods ("episodes cevenols",
# roughly Sep-Nov) shifted? A proper attribution needs the sub-basin gauge
# data the e-mail points to (Observatoire des Sediments du Rhone / Olivier
# Radakovitch) -- RiOMar only has the combined Grand+Petit Rhone discharge
# at the mouth (data/RIVER_FLOW/GULF_OF_LION), so this cannot say which
# tributary is responsible. What follows is a first-pass check of whether
# the *calendar timing* of high-flow days at the combined gauge has moved
# over the record -- worth reporting even though it can't yet attribute a
# cause.
# rhone_flood_timing_shift()
rhone_flood_timing_shift <- function(high_flow_quantile = 0.90){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")
  df_flow <- load_driver("flow", meta) |>
    dplyr::mutate(year = year(date), month = month(date), doy = yday(date))

  # Meteorological seasons (DJF/MAM/JJA/SON), extended from the Cevenol
  # (SON) focus to all four so winter/spring/summer can be compared on the
  # same footing. December is relabelled into the *following* year
  # (season_year = year + 1) so winter isn't split across two different
  # true winters -- e.g. Dec 2010 + Jan/Feb 2011 are one winter, both
  # tagged season_year 2011. For the other three seasons season_year is
  # just the calendar year. This has to happen *before* the high-flow
  # threshold below, because that threshold is also grouped by season_year
  # (see next comment) rather than raw calendar year.
  df_flow <- df_flow |>
    dplyr::mutate(season = dplyr::case_when(
                    month %in% 3:5  ~ "spring (MAM)",
                    month %in% 6:8  ~ "summer (JJA)",
                    month %in% 9:11 ~ "autumn (SON, Cevenol)",
                    TRUE            ~ "winter (DJF)"),
                  season_year = ifelse(month == 12, year + 1, year))
  season_levels <- c("winter (DJF)", "spring (MAM)", "summer (JJA)", "autumn (SON, Cevenol)")
  df_flow$season <- factor(df_flow$season, levels = season_levels)

  # "High flow" is defined relative to each season_year's own range (that
  # season_year's 90th percentile), not a single record-wide threshold -- a
  # fixed absolute threshold would mechanically flag more/fewer days in
  # years where mean flow itself has trended, which is exactly the
  # ambiguity this check is trying to avoid. Grouped by season_year, not
  # raw calendar year, so a single true winter always gets judged against
  # one consistent threshold -- grouping by calendar year instead would
  # let December and the following January/February of the same winter be
  # compared against two different years' flow distributions.
  df_flow <- df_flow |>
    dplyr::mutate(year_thresh = quantile(value, high_flow_quantile, na.rm = TRUE), .by = "season_year") |>
    dplyr::mutate(is_high_flow = value >= year_thresh)

  # Per (calendar) year: how many high-flow days, and their mean day-of-year
  # (NB: a circular statistic tied to the calendar year specifically,
  # because day-of-year itself resets at Jan 1 -- averaging doy within a
  # Dec-Feb season_year group would run straight across that reset and
  # produce a meaningless mid-year value, so this stays on raw "year").
  df_year <- df_flow |>
    dplyr::filter(is_high_flow) |>
    dplyr::summarise(n_high_flow = dplyr::n(), mean_doy_high_flow = mean(doy), .by = "year")

  # Per season-year: each season's share of that season-year's high-flow
  # days -- the four seasons' shares sum to ~1 within a season_year, so
  # this is the direct four-season generalisation of the single prop_autumn
  # check from the first pass. tidyr::complete() adds an explicit 0 (not a
  # missing row) for any season_year x season combination that had no
  # high-flow days at all that year -- "this season had none of the year's
  # extreme days" is real information (and changes the trend below
  # noticeably), not something to silently drop from the regression.
  df_season_year <- df_flow |>
    dplyr::filter(is_high_flow) |>
    dplyr::summarise(n_season = dplyr::n(), .by = c("season_year", "season")) |>
    tidyr::complete(season_year, season, fill = list(n_season = 0)) |>
    dplyr::mutate(n_total = sum(n_season), .by = "season_year") |>
    dplyr::mutate(prop_season = n_season / n_total)

  # Trend of each season's share over the years, and of the mean
  # day-of-year. One value per year/season-year here (not the
  # daily/monthly autocorrelated case fit_wls_hac_trend() is built for), so
  # a plain OLS trend is enough.
  season_trend <- function(season_name){
    d <- dplyr::filter(df_season_year, season == season_name)
    fit <- summary(lm(prop_season ~ season_year, data = d))$coefficients
    tibble::tibble(metric = paste0("proportion_of_high_flow_days_in_", season_name),
                   slope_per_year = fit["season_year", "Estimate"], p_value = fit["season_year", "Pr(>|t|)"])
  }
  stats_season <- purrr::map_dfr(season_levels, season_trend)

  fit_doy <- summary(lm(mean_doy_high_flow ~ year, data = df_year))$coefficients
  stats_doy <- tibble::tibble(metric = "mean_doy_of_high_flow_days",
                              slope_per_year = fit_doy["year", "Estimate"], p_value = fit_doy["year", "Pr(>|t|)"])

  stats <- dplyr::bind_rows(stats_doy, stats_season)
  write_csv(stats, "output/STATS/rhone_flood_seasonality.csv")

  pl <- ggplot(df_season_year, aes(x = season_year, y = prop_season)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
    facet_wrap(~season, ncol = 2) +
    labs(x = NULL, y = paste0("Proportion of top ", round((1 - high_flow_quantile) * 100), "% flow days falling in that season"),
         title = "Grand Rhone: has the seasonal timing of high-flow days shifted?",
         subtitle = "A rising trend means that season is becoming relatively more prominent for floods (Cevenol = autumn/SON panel)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  ggsave(filename = "figures/rhone_side_analyses/rhone_flood_timing_shift.png", plot = pl, width = 10, height = 8)

  return(list(stats = stats, data = df_season_year, plot = pl))
}


## 3. Does the plume respond to wind/wave forcing beyond flow? ---------------
# The e-mail's third idea reframes "does plume size have a seasonal cycle
# beyond river flow" as "does plume size respond directly to wind/wave
# forcing" -- the hypothesis being that winter's stronger wind/waves push
# the plume further offshore (larger detected area) for the same
# discharge, and that THIS, not a genuine seasonal cycle in the plume's own
# dynamics, is what a naive flow-vs-season comparison would pick up. It
# also raises a specific, testable version of that: Mistral (a strong,
# roughly N/NNW wind for the Gulf of Lion) should push the plume clear of
# the coast (larger detected area); weak or onshore/easterly wind should
# leave it hugging the coast, where it may be under-detected because of
# confounding coastal resuspension.
# rhone_wind_wave_effect()
rhone_wind_wave_effect <- function(){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")
  df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, wind_spd = value, wind_dir)
  df_wave <- load_driver("wave", meta) |> dplyr::select(date, wave_height = value)

  df <- df_flow |>
    dplyr::left_join(df_wind, by = "date") |>
    dplyr::left_join(df_wave, by = "date") |>
    tidyr::drop_na(plume_area, flow, wind_spd, wave_height)  # complete cases only, so the nested models below are fit to identical rows

  # a) Does wind/wave explain plume-area variance beyond flow alone? Nested
  #    model comparison (F-test via anova()) -- the same "does adding this
  #    term help" logic as code/6_driver_interactions.R::fit_interaction_glm(),
  #    scoped here to a plain additive comparison since the question is
  #    "does it matter at all", not "how do drivers interact".
  m_flow      <- lm(plume_area ~ flow, data = df)
  m_flow_wind <- lm(plume_area ~ flow + wind_spd, data = df)
  m_flow_wave <- lm(plume_area ~ flow + wave_height, data = df)

  anova_wind <- anova(m_flow, m_flow_wind)
  anova_wave <- anova(m_flow, m_flow_wave)

  model_comparison <- tibble::tibble(
    added_term = c("wind_spd", "wave_height"),
    f_statistic = c(anova_wind$F[2], anova_wave$F[2]),
    p_value = c(anova_wind$`Pr(>F)`[2], anova_wave$`Pr(>F)`[2]),
    r2_flow_only = summary(m_flow)$r.squared,
    r2_with_term = c(summary(m_flow_wind)$r.squared, summary(m_flow_wave)$r.squared)
  )
  write_csv(model_comparison, "output/STATS/rhone_wind_wave_models.csv")

  # b) Mistral vs onshore/easterly vs calm, holding flow's effect fixed by
  #    working with the *residual* plume area after removing the flow
  #    relationship (m_flow above) -- so what gets compared across wind
  #    categories is "area given today's flow", not raw area (which would
  #    just reflect whichever category happened to have wetter days).
  #    Sector definitions are an approximation, not a precise meteorological
  #    classification: Mistral ~ NNW-N (e.g. Guenard et al. 2006 use
  #    ~300-340 deg); the onshore/"Levant" sector ~E-SE is the rough
  #    opposite. wind_dir follows the meteorological "from" convention (see
  #    util.R::.speed_direction()), so these are compass bearings the wind
  #    is blowing FROM.
  df$area_resid <- residuals(m_flow)
  df$wind_category <- dplyr::case_when(
    df$wind_spd < 3                          ~ "calm (<3 m/s)",
    df$wind_dir >= 300 & df$wind_dir <= 350  ~ "Mistral (NNW-N)",
    df$wind_dir >= 90  & df$wind_dir <= 150  ~ "onshore/easterly",
    TRUE                                     ~ "other"
  )
  # Ordered calm -> onshore/easterly -> Mistral -> other (weakest/most
  # coast-hugging to strongest/most offshore-pushing, per the e-mail's
  # hypothesis), rather than the alphabetical default, so this order is
  # what every plot legend, boxplot axis, and summarise() output below uses.
  df$wind_category <- factor(df$wind_category,
                              levels = c("calm (<3 m/s)", "onshore/easterly", "Mistral (NNW-N)", "other"))
  wind_category_colours <- c("calm (<3 m/s)" = "grey50", "onshore/easterly" = "steelblue",
                              "Mistral (NNW-N)" = "firebrick", "other" = "goldenrod")

  category_summary <- df |>
    dplyr::summarise(n_days = dplyr::n(),
                      mean_area_resid = mean(area_resid, na.rm = TRUE),
                      sd_area_resid = sd(area_resid, na.rm = TRUE), .by = "wind_category") |>
    dplyr::arrange(wind_category)  # .by = doesn't preserve factor level order, arrange() does
  write_csv(category_summary, "output/STATS/rhone_wind_wave_categories.csv")
  category_anova <- summary(aov(area_resid ~ wind_category, data = df))

  pl_terms <- ggplot(df, aes(x = wind_spd, y = plume_area)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm", colour = "purple") +
    labs(x = "Wind speed (m/s)", y = "Plume area (km^2)",
         title = "Grand Rhone: plume area vs. wind speed",
         subtitle = "Raw relationship, not controlled for flow (see model_comparison)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  # c) Same flow-controlled residual as the category summary above, but as a
  #    scatter against wind speed (not a categorical boxplot) with a
  #    separate per-category trend line -- shows whether wind speed's
  #    relationship with (flow-controlled) plume area actually differs by
  #    wind category, rather than just comparing category means.
  pl_category <- ggplot(df, aes(x = wind_spd, y = area_resid, colour = wind_category)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    scale_colour_manual(values = wind_category_colours) +
    labs(x = "Wind speed (m/s)", y = "Plume area residual after removing the flow effect (km^2)",
         colour = NULL, title = "Grand Rhone: flow-controlled area vs. wind, by category") +
    theme(panel.border = element_rect(fill = NA, colour = "black"), legend.position = "none")

  # d) The e-mail's "or even wave, which may be simpler?" alternative to
  #    wind -- this plot did not exist before. Plume area vs. wave height,
  #    with wind category as colour so the wind effect is visible alongside
  #    the wave effect on the same panel, per the e-mail's suggestion.
  pl_wave <- ggplot(df, aes(x = wave_height, y = plume_area, colour = wind_category)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    scale_colour_manual(values = wind_category_colours) +
    labs(x = "Wave height (m)", y = "Plume area (km^2)", colour = NULL,
         title = "Grand Rhone: plume area vs. wave height, by wind category",
         subtitle = "Raw relationship, not controlled for flow") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  # NB: ggpubr's common.legend = TRUE (sharing one legend for pl_category
  # and pl_wave) renders a malformed black panel with this ggplot2 version
  # -- avoided here by just dropping pl_category's legend (identical
  # categories/colours to pl_wave's) and keeping pl_wave's own.
  pl_combi <- ggpubr::ggarrange(pl_terms, pl_category, pl_wave, ncol = 3, nrow = 1)
  ggsave(filename = "figures/rhone_side_analyses/rhone_wind_wave_effect.png", plot = pl_combi, width = 18, height = 6)

  return(list(model_comparison = model_comparison, wind_category_summary = category_summary,
              wind_category_anova = category_anova, data = df, plot = pl_combi))
}


# NB: not run automatically on source() -- call explicitly
# rhone_detrend_test()
# rhone_flood_timing_shift()
# rhone_wind_wave_effect()
