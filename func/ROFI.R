# func/ROFI.R
# The analyses used to compare ROFI to plume size

source("func/multi.R")

# Flow / plume / ROFI succession analysis, implementing the points raised by
# Maud about why the SPM turbid plume and the ROFI differ in
# size and in the timing of their annual maximum:
#   1. She measures the ROFI only from the estuary mouth seaward, so ~600 km^2
#      of the Gironde estuary that our panache plume area *does* include is
#      never counted in her Bay of Biscay ROFI -- part of why it looks
#      smaller than the plume there. load_plume_no_estuary() below recomputes
#      the plume area with the same estuary excluded, using the lat/lon boxes
#      from the masking code she sent (second email), so the two are
#      comparable on the same footing.
#   2. Her ROFI threshold (33 psu) is high enough that the ROFI keeps
#      expanding well after the SPM plume has started sedimenting out, so she
#      expects the SPM plume maximum 1-2 months before the ROFI maximum, and
#      the succession flow peak -> plume peak -> ROFI peak to be visible if
#      flow is added to the same time series. plot_flow_plume_rofi_succession()
#      and rofi_plume_peak_lag() below do exactly that.
#   3. Her river sediment-concentration-seasonality hypothesis (more sediment
#      load in early winter after the dry season, letting the plume travel
#      further before dropping below the SPM threshold) is not tested because 
#      this project has no in-river SPM concentration time series, only discharge.

# Estuary bounding boxes, transcribed as-is from the mask_estuary Python code
# Maud sent. Applied here to the panache lon/lat grid so the exclusion matches.
in_estuary <- function(lon, lat, zone_name){
  if(zone_name == "BAY_OF_BISCAY"){            # Gironde
    lat < 45.75 & lon > -1.06
  } else if(zone_name == "BAY_OF_SEINE"){      # Seine
    (lat < 49.51 & lon > 0.21) | (lat < 49.34 & lon > 0.04)
  } else if(zone_name == "SOUTHERN_BRITTANY"){ # Loire
    (lat > 47.12 & lon > -2.11) | (lat > 47.275 & lon > -2.26)
  } else {
    stop("No estuary mask defined for zone ", zone_name, " (Maud's ROFI only covers Seine/Gironde/Loire).")
  }
}

# Area of one PlumeMasks.nc pixel in km^2, same formula as
# func/plume.py::return_stats_dictionnary() and func/compute_plume_shape.py
# (degrees -> km via 111.32 km/deg, longitude scaled by cos(latitude)).
pixel_area_km2 <- function(zone_name, plume_dir = "output/panache/dynamic"){
  nc_dat <- nc_open(paste0(plume_dir, "/", zone_name, "/PlumeMasks.nc"))
  lon <- ncvar_get(nc_dat, "lon")
  lat <- ncvar_get(nc_dat, "lat")
  nc_close(nc_dat)
  dlon <- mean(diff(lon))
  dlat <- mean(diff(lat))
  (dlon * 111.32 * cos(pi / 180 * mean(lat))) * (dlat * 111.32)
}

# Daily plume area with Maud's estuary boxes excluded, built by re-counting
# the per-pixel plume hits already in PlumeMasks.nc (via util.R's
# load_plume_surface()) rather than reloading the raw SEXTANT SPM
# concentration.
# NB: load_plume_surface() only returns rows for pixels flagged as plume, so
# a day where every plume pixel happens to fall inside the estuary box would
# otherwise silently disappear instead of counting as area = 0. Results.csv
# (loaded here only for its date column) is used as the reference set of
# "day actually processed by panache" to catch that -- it already records
# zero-plume days explicitly (area = 0), never by omission.
# This recomputes area independently of load_plume_ts(), so util.R's
# plume_area_ceiling guard (the per-zone hard cap on implausible
# dynamic-threshold spillover days) is reapplied here by hand -- otherwise a
# spillover day would be filtered out of the standard plume area but not out
# of this one.
load_plume_no_estuary <- function(zone_name, plume_dir = "output/panache/dynamic"){
  pixel_area_km2 <- pixel_area_km2(zone_name, plume_dir)

  df_ref <- read_csv(paste0(plume_dir, "/", zone_name, "/Results.csv"), show_col_types = FALSE) |>
    dplyr::transmute(date = as.Date(date))

  df_ref |>
    dplyr::left_join(
      load_plume_surface(zone_name, plume_dir) |>
        dplyr::filter(!in_estuary(lon, lat, zone_name)) |>
        dplyr::summarise(plume_area = dplyr::n() * pixel_area_km2, .by = "date"),
      by = "date") |>
    dplyr::mutate(plume_area = tidyr::replace_na(plume_area, 0),
                  plume_area = ifelse(plume_area > plume_area_ceiling[[zone_name]], NA, plume_area)) |>
    tidyr::complete(date = seq(min(date), max(date), by = "day")) |>  # calendar gaps stay NA, same convention as load_plume_ts()
    dplyr::mutate(zone = zone_name, .before = "date")
}

# Join flow + plume + ROFI on date for one of the three zones with ROFI data.
# exclude_estuary = TRUE uses load_plume_no_estuary() so the plume area is
# comparable to Maud's estuary-excluded ROFI; FALSE uses the standard
# load_plume_ts() plume area (the one used everywhere else in the pipeline).
combine_flow_plume_rofi <- function(meta, exclude_estuary = TRUE){
  df_plume <- if(exclude_estuary) load_plume_no_estuary(meta$zone) else load_plume_ts(meta$zone)
  df_flow  <- load_driver("flow", meta) |> dplyr::select(date, flow = value)
  df_rofi  <- load_driver("rofi", meta) |> dplyr::select(date, rofi = value)

  df <- df_plume |>
    dplyr::select(date, plume_area) |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::left_join(df_rofi, by = "date") |>
    zoo::na.trim()

  start_date <- min(df$date)
  df$flow_stl  <- stl_single(df$flow,       "inter", start_date)
  df$plume_stl <- stl_single(df$plume_area, "inter", start_date)
  df$rofi_stl  <- stl_single(df$rofi,       "inter", start_date)

  df |> dplyr::mutate(zone = meta$zone, exclude_estuary = exclude_estuary, .before = "date")
}

# Single-zone dual-axis panel: plume area and ROFI extent share the left axis
# (same units, km^2, so no rescaling between them), river flow gets the right
# axis (rescaled onto the left axis's range, same convention as
# plot_driver_plume_dual_axis()). Raw daily values only so the 
# flow -> plume -> ROFI succession Maud described is visible exactly
# as observed, not as the interannual signal.
plot_flow_plume_rofi_panel <- function(zone_name, show_axis_titles = TRUE){
  meta <- get_zone_meta(zone_name = zone_name)
  df <- combine_flow_plume_rofi(meta, exclude_estuary = TRUE)

  scaling_factor <- sec_axis_adjustement_factors(var_to_scale = df$flow, var_ref = c(df$plume_area, df$rofi))
  df <- df |> dplyr::mutate(flow_scaled = flow * scaling_factor$diff + scaling_factor$adjust)

  ggplot(df, aes(x = date)) +
    geom_line(aes(y = plume_area, colour = "Plume area"), alpha = 0.8) +
    geom_line(aes(y = rofi, colour = "ROFI extent"), alpha = 0.8) +
    geom_line(aes(y = flow_scaled, colour = "River flow"), alpha = 0.8) +
    scale_colour_manual(name = NULL, values = c("Plume area" = "sienna", "ROFI extent" = "goldenrod", "River flow" = "blue")) +
    scale_y_continuous(name = if(show_axis_titles) "Plume area / ROFI extent (km^2)" else NULL,
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff,
                                           name = if(show_axis_titles) "River flow (m^3 s-1)" else NULL)) +
    labs(x = NULL, title = zone_title(zone_name)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          axis.title.y.right = element_text(color = "blue"),
          axis.text.y.right = element_text(color = "blue"),
          axis.ticks.y.right = element_line(color = "blue"),
          axis.line.y.right = element_line(color = "blue"))
}

# Stack the three zones' dual-axis panels into a single figure (one file, not
# three), sharing one legend, so the flow -> plume -> ROFI succession can be
# compared across all three rivers at a glance.
plot_flow_plume_rofi_succession <- function(){
  rofi_zones <- order_zones(c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY"))
  panels <- purrr::map(rofi_zones, plot_flow_plume_rofi_panel, show_axis_titles = FALSE)

  # Per-panel axis titles suppressed in favour of one shared left/right label
  # on the assembled composite (Figure 4/6 convention), via
  # ggpubr::annotate_figure() -- ggarrange()'s common.legend achieves the
  # same shared-legend effect patchwork::plot_layout(guides = "collect") did.
  full_plot <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = panels, ncol = 1, nrow = length(panels),
                      common.legend = TRUE, legend = "bottom", align = "v"),
    left = ggpubr::text_grob("Plume area / ROFI extent (km²)", rot = 90, size = 14),
    right = ggpubr::text_grob("River flow (m³ s⁻¹)", rot = -90, size = 14, color = "blue"))

  # manuscript Figure S8 (flow -> plume -> ROFI succession)
  if (!dir.exists("figures/ARTICLE/FIGURE_S4")) dir.create("figures/ARTICLE/FIGURE_S4", recursive = TRUE)
  ggsave(filename = "figures/ARTICLE/FIGURE_S4/flow_plume_rofi_succession.png", plot = full_plot, width = 12, height = 10, dpi = 300)
  invisible(full_plot)
}

# Lagged correlation between two of {plume_area, flow, rofi} at daily/
# monthly/annual timesteps, for one zone.
# lagged_correlation(x, y, lag)'s "lag" is: the shift k (days) at which x
# from k days earlier best correlates with y today -- e.g. x_col =
# "plume_area", y_col = "rofi" asks at what k plume_area(t-k) best matches
# rofi(t), i.e. x_col leads y_col. max_lag_daily defaults to 90 (vs. the
# pipeline's usual 30) to give Maud's hypothesized 1-2 month (~30-60 day)
# lag headroom on both sides; checked out to 120 days for BAY_OF_SEINE
# plume-vs-ROFI specifically, whose correlation is weak (r < 0.28
# throughout) and multi-modal rather than one clean peak -- so its reported
# peak lag is less trustworthy than the other pairs/zones, which do show a
# single clear maximum well inside this window.
rofi_plume_lagged_correlation <- function(zone_name, x_col, y_col, max_lag_daily = 90){
  meta <- get_zone_meta(zone_name = zone_name)
  df_wide <- combine_flow_plume_rofi(meta, exclude_estuary = TRUE)
  df <- tibble::tibble(date = df_wide$date, value = df_wide[[x_col]], plume_area = df_wide[[y_col]])

  driver_plume_correlation(df, max_lag_daily) |>
    dplyr::mutate(zone = zone_name, x = x_col, y = y_col, .before = "lag")
}

# TODO: Run this for each calendar month in order to show how the lagged correlations
# change by month (i.e. season). Use different coloured lines to denote the months.
# Visualise the daily-timestep lagged correlations across zones and variable
# pairs: correlation-vs-lag curve per panel, with the lag of maximum r marked in red
plot_rofi_plume_lagged_correlation <- function(cor_stats){
  df_daily <- dplyr::filter(cor_stats, timestep == "daily") |>
    dplyr::mutate(zone = factor(zone_title(zone), levels = zone_title(order_zones(unique(zone)))))
  df_peak <- df_daily |> dplyr::slice_max(cor, n = 1, by = c(zone, pair_label))

  pl <- ggplot(df_daily, aes(x = lag, y = cor)) +
    geom_line(colour = "steelblue") +
    geom_point(data = df_peak, colour = "firebrick", size = 2) +
    geom_text(data = df_peak, aes(label = paste0(lag, "d")), vjust = -1, size = 3, colour = "firebrick") +
    facet_grid(zone ~ pair_label) +
    labs(x = "Lag (days)", y = "Correlation (r)",
         title = "Lagged correlation between river flow, plume area and ROFI extent",
         subtitle = "Daily resolution; estuary excluded from plume area; red = lag of maximum r") +
    theme_bw()

  # manuscript Figure S9
  if (!dir.exists("figures/ARTICLE/FIGURE_S4")) dir.create("figures/ARTICLE/FIGURE_S4", recursive = TRUE)
  ggsave(filename = "figures/ARTICLE/FIGURE_S4/rofi_plume_lagged_correlation.png", plot = pl, width = 12, height = 9, dpi = 300)
  invisible(pl)
}

# Run the succession plot + lagged-correlation stats for all three zones with
# ROFI data (no Gulf of Lion ROFI data), estuary excluded from the plume area
# throughout (the only version we're keeping -- see combine_flow_plume_rofi()),
# and save the combined correlation table + figures.
run_rofi_plume_succession <- function(){
  rofi_zones <- order_zones(c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY"))
  pairs <- tibble::tribble(
    ~x_col,       ~y_col,  ~pair_label,
    "plume_area", "flow",  "Plume vs flow",
    "plume_area", "rofi",  "Plume vs ROFI",
    "rofi",       "flow",  "ROFI vs flow"
  )

  cor_stats <- purrr::map_dfr(rofi_zones, function(zone_name){
    purrr::pmap_dfr(pairs, function(x_col, y_col, pair_label){
      rofi_plume_lagged_correlation(zone_name, x_col, y_col) |>
        dplyr::mutate(pair_label = pair_label)
    })
  })

  write_csv(cor_stats, "output/STATS/rofi_plume_lagged_correlation.csv")
  plot_flow_plume_rofi_succession()
  plot_rofi_plume_lagged_correlation(cor_stats)
  invisible(cor_stats)
}
