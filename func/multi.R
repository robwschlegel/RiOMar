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
#     plot_driver_plume_comparison().
#   - wind.R/tide.R's inline dual-y-axis STL plot (a hand-copy of
#     figure.R::make_the_X11_plot_of_river_and_plume(), hard-coded there to
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
#   - surface.R's two pixel-level, multi-driver spatial functions
#     (surface_plot, surface_plot_daily_maps) are carried over with light
#     cleanup (now use zone_meta / make_pretty_title instead of their own
#     inline zone lookups); these are conceptually different from the
#     time-series comparisons above (per-pixel maps, not per-day series) so
#     they keep their own section rather than being folded into the
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
# NB: surface_plot() below (merged in from func/surface.R) is a first pass at this idea.

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
library(tidync)
library(heatwaveR) # For seasonal smoothing analysis
library(seasonal) # For X11 analysis (currently not used)
library(RcppRoll) # For running means to get STL interannual signals closer to X11
library(patchwork)
library(sandwich) # For HAC covariance tests (driver_plume_trend)
library(lmtest) # For more detailed linear model tests (driver_plume_trend)
library(doParallel); doParallel::registerDoParallel(cores = 14)

# Common function
source("func/util.R")

# Zones
zones <- c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION")


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
  )

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


# Driver loading ------------------------------------------------------------

# Shared on-/off-shore wind classification. Previously duplicated verbatim
# (down to the same TODO comment) in wind.R::spatial_wind_calc(),
# surface.R::surface_plot(), and multi_stl() further down this file.
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
  df_wind |>
    dplyr::mutate(wind_spd = round(sqrt(u^2 + v^2), 2),
                  # NB: wind_dir is where the wind is coming FROM, not going to
                  wind_dir = round((270 - (atan2(v, u) * (180 / pi))) %% 360))
}

# Load one driver's daily time series for a zone, in a common two-column
# (date, value) shape so downstream functions don't need to know which
# driver they're looking at. Generalises the per-driver loading blocks
# previously duplicated across flow.R/wind.R/tide.R/ROFI.R.
# load_driver("flow", get_zone_meta(zone_name = "GULF_OF_LION"))
load_driver <- function(driver_name, meta){
  driver_name <- match.arg(driver_name, c("flow", "tide", "wind", "current"))

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
    # "current" = coastal-current / ROFI extent, i.e. the closest observational
    # proxy this project has for the "coastal ocean currents" driver discussed
    # in manuscript/driver_interactions_review.md. NB: no ROFI data for the
    # Gulf of Lion (see comment on the old ROFI.R::walk(zones[1:3], ...) call).
    rofi_files <- dir("data/ROFI", full.names = TRUE)
    df_rofi <- purrr::map_dfr(rofi_files, load_ROFI) |>
      dplyr::filter(zone == meta$zone)
    if(nrow(df_rofi) == 0) stop("No ROFI/current data available for zone ", meta$zone)
    df <- df_rofi |> dplyr::select(date, value = ROFI_surface)
  }

  return(df)
}


# Combine plume + one driver --------------------------------------------

# Load plume + a single driver for one zone, join on date, and add STL
# interannual columns for both. This is the core object every comparison
# function below operates on.
# combine_plume_driver("flow", get_zone_meta(mouth_name = "Seine"))
combine_plume_driver <- function(driver_name, meta){

  df_plume  <- load_plume_ts(meta$zone)                    # util.R -- already handles gap-filling + 20000km^2 outlier removal
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
  # TODO: Current should be the U and V vectors for the surface currents from the GLORYS model
  # ROFI is a different variable. This comes from a model output of the region of freshwater influence
  # Once GLORYS data, this should be corrected and ROFI should be it's own driver_name
  "current",    "ROFI / current extent (km^2)", "orchid"
)

# The 4-panel (a-d) comparison plot previously duplicated in
# flow.R::flow_comp(), wind.R::spatial_wind_calc(), and tide.R::tide_calc():
#   a) raw driver time series
#   b) raw plume-area time series
#   c) driver vs. plume scatter (+ linear fit)
#   d) lagged correlation (plume lagged behind driver, 0-30 days)
# mouth_name is used only for the plot title / output file name.
plot_driver_plume_comparison <- function(df, driver_name, mouth_name){

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

  ggsave(filename = paste0("figures/cor_plot_", driver_name, "_plume_", mouth_name, ".png"),
         plot = full_plot_title, width = 12, height = 6, dpi = 600)
  invisible(full_plot_title)
}

# Dual-y-axis STL plot (plume_stl on the left axis, driver_stl scaled onto
# the right axis). Generalises the inline duplicate of
# figure.R::make_the_X11_plot_of_river_and_plume() that wind.R and tide.R had
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

  ggsave(filename = paste0("figures/dual_axis_", driver_name, "_plume_", zone_name, ".png"),
         plot = pl, width = 12, height = 6, dpi = 300)
  invisible(pl)
}

# Weighted-least-squares trend estimator with AR(1)/STL-based weights and a
# Newey-West (HAC) standard error correction, following the monthly time
# series adjustment methodology of Sutton et al. (2022) referenced in
# manuscript.tex Sec. 2.6.1. Generalises flow.R::flow_plume_trend_plus(),
# which implemented this for river flow only; used here for any driver.
# driver_plume_trend(combine_plume_driver("flow", get_zone_meta(mouth_name = "Grand Rhone")), "flow", "Grand Rhone")
driver_plume_trend <- function(df, driver_name, mouth_name, end_date = NULL){

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

  # AR(1)-weighted / STL-weighted / unweighted linear models + HAC SEs
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
  lm_HAC_weights_func <- function(weight_choice, val_col, date_col){
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

  wls_driver_daily   <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, val_col = df_daily$driver_doy_adj, date_col = df_daily$date)
  wls_driver_monthly <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, val_col = df_monthly$driver_monthly_adj, date_col = df_monthly$date)
  wls_plume_daily    <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, val_col = df_daily$plume_doy_adj, date_col = df_daily$date)
  wls_plume_monthly  <- plyr::ldply(c("ar", "stl", "none"), lm_HAC_weights_func, val_col = df_monthly$plume_monthly_adj, date_col = df_monthly$date)

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

  # Plot (daily = ar-weighted line; monthly = ar-weighted line), labelled with slope + p-value
  trend_labels_driver <- dplyr::filter(stats, variable == "driver", weight_choice == "ar")
  trend_labels_plume  <- dplyr::filter(stats, variable == "plume",  weight_choice == "ar")
  x_daily <- min(df_daily$date) + days(round(0.05 * as.numeric(diff(range(df_daily$date)))))
  x_monthly <- min(df_daily$date) + days(round(0.45 * as.numeric(diff(range(df_daily$date)))))

  pl_plume <- ggplot(data = df_daily, aes(x = date, y = plume_area)) +
    geom_point(colour = "sienna", alpha = 0.1) +
    geom_point(aes(y = plume_doy_adj), colour = "darkblue", alpha = 0.1) +
    geom_point(data = df_monthly, aes(y = plume_monthly), colour = "sienna", alpha = 0.6, size = 3) +
    geom_point(data = df_monthly, aes(y = plume_monthly_adj), colour = "darkred", size = 3) +
    geom_abline(data = dplyr::filter(trend_labels_plume, timestep == "daily"),
                aes(intercept = intercept, slope = slope), linewidth = 2, colour = "darkblue") +
    geom_abline(data = dplyr::filter(trend_labels_plume, timestep == "monthly"),
                aes(intercept = intercept, slope = slope), linewidth = 2, colour = "darkred") +
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
    labs(x = NULL, y = disp$driver_label,
         title = paste0(mouth_name, " : ", driver_name, " after statistical treatment"),
         subtitle = "Red = adjusted monthly values; blue = adjusted daily values; purple = original data") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  pl_combi <- ggpubr::ggarrange(pl_plume, pl_driver, ncol = 1, nrow = 2)
  ggsave(filename = paste0("figures/trends_plume_", driver_name, "_adj_", mouth_name, ".png"),
         pl_combi, width = 12, height = 10)

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

  mouths <- if(driver_name == "current") dplyr::filter(zone_meta, zone != "GULF_OF_LION") else zone_meta  # no ROFI data for the Gulf of Lion

  results <- purrr::pmap(mouths, function(...){
    meta <- tibble::tibble(...)
    df <- combine_plume_driver(driver_name, meta)

    plot_driver_plume_comparison(df, driver_name, meta$mouth_name)
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

# Run all four drivers for all zones. NB: this recreates every figure
# previously produced by running flow.R + wind.R + tide.R + ROFI.R in
# sequence; it is not run automatically on source() (see bottom of file).
run_all_driver_suites <- function(){
  purrr::map(c("flow", "tide", "wind", "current"), run_driver_suite) |>
    purrr::set_names(c("flow", "tide", "wind", "current"))
}


# Surface / pixel-level multi-driver maps ------------------------------
# Carried over from func/surface.R with light cleanup (zone_meta /
# make_pretty_title instead of an inline zone lookup); conceptually distinct
# from the time-series work above since these operate per-pixel rather than
# per-day, so they keep their own section.

# Load all daily plume pixels for a zone, join on-date driver values, and
# summarise per pixel (count, driver ranges, on-/off-shore wind proportion).
# zone_name <- "GULF_OF_LION"
surface_plot <- function(zone_name){

  meta <- get_zone_meta(zone_name = zone_name)

  # Zone-specific plotting bounding boxes (wider than the driver-search box
  # used elsewhere, purely for map extent) -- unchanged from surface.R
  bbox_wide <- switch(zone_name,
    "BAY_OF_SEINE"      = list(lon_W = -1.0, lon_E = 0.5, lat_N = 0.5, lat_S = -0.2),
    "BAY_OF_BISCAY"     = list(lon_W = -3.0, lon_E = 0.1, lat_N = 1.0, lat_S = -0.5),
    "SOUTHERN_BRITTANY" = list(lon_W = -3.0, lon_E = 0.2, lat_N = 0.5, lat_S = -1.0),
    "GULF_OF_LION"      = list(lon_W = -2.0, lon_E = 3.0, lat_N = 0.0, lat_S = -2.5),
    stop("Zone not recognised.")
  )
  lon_round <- plyr::round_any(meta$mouth_lon, 0.5); lat_round <- plyr::round_any(meta$mouth_lat, 0.5)
  lon_range_wide <- c(lon_round + bbox_wide$lon_W, lon_round + bbox_wide$lon_E)
  lat_range_wide <- c(lat_round + bbox_wide$lat_S, lat_round + bbox_wide$lat_N)

  # Load all daily plume maps for the zone
  plume_dir <- paste0("output/panache/", zone_name)
  plume_files <- dir(plume_dir, pattern = ".csv", recursive = TRUE, full.names = TRUE)
  df_plume <- plyr::ldply(plume_files, load_plume_surface, .parallel = TRUE)  # util.R

  # Load drivers via the shared loader
  df_flow <- load_driver("flow", meta) |> dplyr::rename(flow = value)
  df_tide <- load_driver("tide", meta) |> dplyr::rename(tide_range = value)
  df_wind <- load_driver("wind", meta) |> dplyr::rename(wind_spd = value)

  df_full <- df_plume |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::left_join(df_tide, by = "date") |>
    dplyr::left_join(df_wind, by = "date")

  suppressWarnings(
    df_pixel <- dplyr::summarise(df_full,
                                 count = dplyr::n(),
                                 flow_min = min(flow, na.rm = TRUE), flow_mean = mean(flow, na.rm = TRUE), flow_max = max(flow, na.rm = TRUE),
                                 tide_range_min = min(tide_range, na.rm = TRUE), tide_range_mean = mean(tide_range, na.rm = TRUE), tide_range_max = max(tide_range, na.rm = TRUE),
                                 wind_spd_min = min(wind_spd, na.rm = TRUE), wind_spd_mean = mean(wind_spd, na.rm = TRUE), wind_spd_max = max(wind_spd, na.rm = TRUE),
                                 count_on = sum(direction == "on", na.rm = TRUE), count_off = sum(direction == "off", na.rm = TRUE),
                                 .by = c("lon", "lat")) |>
    dplyr::mutate(prop_n = count / length(unique(df_full$date)),
                  prop_on = count_on / count, prop_off = count_off / count,
                  dplyr::across(dplyr::everything(), ~ ifelse(is.finite(.), ., NA)))
  )

  plot_count <- ggplot(df_pixel, aes(x = lon, y = lat)) +
    annotation_borders(regions = "France", fill = "grey70") +
    geom_tile(aes(fill = log10(count))) + scale_fill_viridis_c() +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Plume count (n)") + theme_bw() + theme(legend.position = "bottom")

  plot_count_prop <- ggplot(df_pixel, aes(x = lon, y = lat)) +
    annotation_borders(regions = "France", fill = "grey70") +
    geom_tile(aes(fill = prop_n)) + scale_fill_viridis_c() +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Plume count proportion (n/all_days)") + theme_bw() + theme(legend.position = "bottom")

  plot_flow <- df_pixel |>
    dplyr::select(lon, lat, flow_min, flow_mean, flow_max) |>
    tidyr::pivot_longer(cols = c(flow_min, flow_mean, flow_max), names_to = "var") |>
    dplyr::mutate(var = factor(var, levels = c("flow_min", "flow_mean", "flow_max"))) |>
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = value)) + scale_fill_viridis_c(option = "A") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "River flow range") + facet_wrap(~var) + theme_bw() + theme(legend.position = "bottom")

  plot_wind_spd <- df_pixel |>
    dplyr::select(lon, lat, wind_spd_min, wind_spd_mean, wind_spd_max) |>
    tidyr::pivot_longer(cols = c(wind_spd_min, wind_spd_mean, wind_spd_max), names_to = "var") |>
    dplyr::mutate(var = factor(var, levels = c("wind_spd_min", "wind_spd_mean", "wind_spd_max"))) |>
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = value)) + scale_fill_viridis_c(option = "C") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Wind speed range") + facet_wrap(~var) + theme_bw() + theme(legend.position = "bottom")

  plot_wind_dir <- df_pixel |>
    dplyr::select(lon, lat, prop_on, prop_off) |>
    tidyr::pivot_longer(cols = c(prop_on, prop_off), names_to = "var") |>
    dplyr::mutate(var = factor(var, levels = c("prop_on", "prop_off"))) |>
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = value)) +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Proportion of on- or off-shore winds") + facet_wrap(~var) + theme_bw() + theme(legend.position = "bottom")

  plot_tide_range <- df_pixel |>
    dplyr::select(lon, lat, tide_range_min, tide_range_mean, tide_range_max) |>
    tidyr::pivot_longer(cols = c(tide_range_min, tide_range_mean, tide_range_max), names_to = "var") |>
    dplyr::mutate(var = factor(var, levels = c("tide_range_min", "tide_range_mean", "tide_range_max"))) |>
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = value)) + scale_fill_viridis_c(option = "E") +
    coord_quickmap(xlim = lon_range_wide, ylim = lat_range_wide) +
    labs(x = NULL, y = NULL, title = "Tidal range range") + facet_wrap(~var) + theme_bw() + theme(legend.position = "bottom")

  plot_multi <- (plot_count | plot_count_prop) / (plot_flow | plot_tide_range) / (plot_wind_spd | plot_wind_dir)
  ggsave(filename = paste0("figures/surface_stats_", zone_name, ".png"), plot = plot_multi, height = 14, width = 26)
  invisible(df_pixel)
}

# Facet daily plume maps by year x month for a zone.
surface_plot_daily_maps <- function(zone_name){

  df_plume <- plyr::ldply(zone_name, load_plume_surface, .parallel = FALSE) |>  # NB: do not run in parallel, it's done elsewhere
    dplyr::mutate(year = year(date), month = month(date), doy = yday(date), day = day(date))

  plot_daily <- ggplot(df_plume, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = day), alpha = 0.3) +
    scale_fill_viridis_c(option = "A", na.value = "transparent") +
    labs(x = NULL, y = NULL, fill = "Day of month", title = paste0("Daily plume maps per month and year for ", zone_name)) +
    facet_grid(year ~ month) + theme_bw() + theme(legend.position = "bottom")
  ggsave(filename = paste0("figures/surface_daily_maps_", zone_name, ".png"), plot = plot_daily, height = 34, width = 36)
  invisible(plot_daily)
}


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
  plot_daily(df_pretty, "plume_area", "brown", "Plume area (km^2)", "plume")

  # Seasonal time series
  plot_seas(df_seas, "plume_seas", "brown", "Plume area (km^2)", "plume")
  plot_seas(df_seas, "flow_seas", "blue", "River flow (m^3 s-1)", "flow")
  plot_seas(df_seas, "tide_seas", "darkgreen", "Tidal range (m)", "tide")
  plot_seas(df_seas, "wind_seas", "purple", "Wind speed (m s-1)", "wind")

  # Interannual time series
  plot_inter(df_pretty, "plume_inter", "brown", "Plume area (km^2)", "plume")
  plot_inter(df_pretty, "flow_inter", "blue", "River flow (m^3 s-1)", "flow")
  plot_inter(df_pretty, "tide_inter", "darkgreen", "Tidal range (m)", "tide")
  plot_inter(df_pretty, "wind_inter", "purple", "Wind speed (m s-1)", "wind")

  # Seasonal comparison plots
  comparison_plot_save(df_seas, "plume_seas", "flow_seas", "brown", "blue", "Plume area (km^2)", "River flow (m^3 s-1)", "comparison_plume_flow_seas")
  comparison_plot_save(df_seas, "plume_seas", "wind_seas", "brown", "purple", "Plume area (km^2)", "Wind speed (m s-1)", "comparison_plume_wind_seas")
  comparison_plot_save(df_seas, "plume_seas", "tide_seas", "brown", "darkgreen", "Plume area (km^2)", "Tidal range (m)", "comparison_plume_tide_seas")

  # Interannual comparison plots
  comparison_plot_save(df_pretty, "plume_inter", "flow_inter", "brown", "blue", "Plume area (km^2)", "River flow (m^3 s-1)", "comparison_plume_flow_inter")
  comparison_plot_save(df_pretty, "plume_inter", "tide_inter", "brown", "darkgreen", "Plume area (km^2)", "Tidal range (m)", "comparison_plume_tide_inter")
  comparison_plot_save(df_pretty, "plume_inter", "wind_inter", "brown", "purple", "Plume area (km^2)", "Wind speed (m s-1)", "comparison_plume_wind_inter")

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
  ggsave(filename = "figures/all_plot.png", plot = all_plot, width = 20, height = 20, dpi = 300)
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
  ggsave("figures/missng_SPM.png", width = 9, height = 9, dpi = 600)
  ggplot(chla_files_NA_count, aes(x = month, y = miss_count_month_year)) +
    geom_col() +
    facet_wrap(~year) +
    labs(x = NULL, y = "count", title = "Monthly count of missing chl a SEXTANT files") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  ggsave("figures/missng_chla.png", width = 9, height = 9, dpi = 600)
}


# Run everything -----------------------------------------------------------

# NB: not run automatically on source() -- call explicitly
# run_all_driver_suites()
# purrr::walk(zones, surface_plot)
# purrr::walk(zones, surface_plot_daily_maps)
  