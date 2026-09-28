# Plotting --------------------------------------------------------------

# Driver display metadata (label + colour), used by the plotting functions
# below so callers only need to pass a driver_name.
driver_display <- tibble::tribble(
  ~driver_name, ~driver_label,             ~driver_colour,
  "flow",       "River flow (m³ s⁻¹)",     "blue",
  "tide",       "Tidal range (m)",         "darkgreen",
  "wind",       "Wind speed (m s⁻¹)",      "purple",
  "current",    "Current speed (m s⁻¹)",   "orchid",
  "wave",       "Wave height (m)",         "steelblue"
)

# The 4-panel (a-d) comparison plot
#   a) raw driver time series
#   b) raw plume-area time series
#   c) driver vs. plume scatter (+ linear fit)
#   d) lagged correlation (plume lagged behind driver, 0-30 days)
# zone_name is used only for the plot title / output file name.
plot_driver_comparison <- function(df, driver_name, zone_name){

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)
  cor_df <- driver_plume_correlation(df) |> dplyr::filter(timestep == "daily")

  driver_plot <- ggplot(df, aes(x = date, y = value)) +
    geom_line() +
    labs(y = disp$driver_label, x = NULL) +
    scale_x_date(expand = c(0, 0)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  panache_plot <- ggplot(df, aes(x = date, y = plume_area)) +
    geom_line() +
    labs(y = "plume area (km²)", x = NULL) +
    scale_x_date(expand = c(0, 0)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  driver_plume_cor_plot <- ggplot(df, aes(x = value, y = plume_area)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 2) +
    labs(y = "plume area (km²)", x = disp$driver_label) +
    theme(panel.border = element_rect(fill = NA, colour = "black"), legend.position = "bottom")

  driver_plume_cor_lag_plot <- ggplot(cor_df, aes(x = lag, y = cor)) +
    geom_point() +
    labs(x = paste("lag plume after", disp$driver_name, "(days)"), y = "correlation (r)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  plot_title <- grid::textGrob(paste0(zone_title(zone_name), " : ", driver_name, " vs plume size"),
                               gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(driver_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
  cor_plot <- ggpubr::ggarrange(driver_plume_cor_plot, driver_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
  full_plot_title <- ggpubr::ggarrange(plot_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")

  ggsave(filename = paste0("figures/driver_comparison/cor_plot_", driver_name, "_plume_", zone_name, ".png"),
         plot = full_plot_title, width = 12, height = 6, dpi = 600)
  invisible(full_plot_title)
}

# Bin a compass bearing (degrees, "from" convention) into one of 8 ordered
# compass octants. 
# Kept separate from plot_driver_rose()'s own inline sector binning below 
# since that one uses a finer, plot-resolution-tuned n_sectors (16) for bar drawing, 
# not a fixed categorical driver definition.
compass_octant <- function(degrees){
  labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  factor(labels[(round(degrees / 45) %% 8) + 1], levels = labels)
}

# Extract residual for driver accounting for river flow
flow_controlled_residual <- function(df){
  residuals(lm(plume_area ~ flow, data = df))
}

# Data prep shared by plot_driver_rose() below and by
# plot_driver_rose_diagram()'s shared-colour-scale computation (figure.R):
# join driver direction/magnitude onto the flow-controlled plume series, bin
# into compass sectors, and winsorize the per-sector mean residual. Returns
# a zero-row tibble (same columns) when direction data isn't available for
# this zone/driver, so callers can detect that case uniformly.
compute_driver_rose_summary <- function(driver_name, meta, n_sectors = 8, df_flow = NULL){
  driver_name <- match.arg(driver_name, c("wind", "wave", "current"))
  dir_col <- paste0(driver_name, "_dir")

  if(is.null(df_flow)){
    df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)
  }
  df_driver <- load_driver(driver_name, meta)

  df <- df_flow |>
    dplyr::left_join(df_driver, by = "date") |>
    tidyr::drop_na(plume_area, flow, dplyr::all_of(dir_col))

  # Defensive fallback for a zone/driver missing direction entirely.
  if(nrow(df) == 0){
    return(tibble::tibble(sector = numeric(0), n_days = integer(0),
                          mean_area_resid = numeric(0), pct_days = numeric(0),
                          mean_area_resid_plot = numeric(0)))
  }

  df$area_resid <- flow_controlled_residual(df)

  sector_width <- 360 / n_sectors
  df$sector <- (round(df[[dir_col]] / sector_width) %% n_sectors) * sector_width

  df_summary <- df |>
    dplyr::summarise(n_days = dplyr::n(),
                     mean_area_resid = mean(area_resid, na.rm = TRUE), .by = "sector") |>
    dplyr::mutate(pct_days = 100 * n_days / sum(n_days))

  # Sectors with few days give a noisy mean residual that can be an extreme
  # outlier purely from small-sample luck, which stretches the fill colour
  # scale and washes out the contrast between the well-sampled, more
  # meaningful sectors. Winsorize: clamp any under-sampled sector's fill
  # value to the min/max range spanned by the well-sampled sectors, rather
  # than letting it set the scale's extremes. This only affects the fill
  # colour, never bar height/radius (n_days/pct_days are untouched).
  #
  # "Well-sampled" is each sector's expected/uniform share of days
  # (100/n_sectors, e.g. 12.5% for 8 sectors), not a fixed 5% -- a fixed 5%
  # cutoff let a sector just above it (e.g. 7% of days) both define the
  # clamp range *and* be exempt from clamping, even when far less sampled
  # than the actually-dominant sectors. Using the same threshold both
  # ways is symmetric: every sector below it is subject to clamping,
  # every sector at or above it both defines the range and is left unclamped.
  uniform_share <- 100 / n_sectors
  well_sampled_range <- df_summary |>
    dplyr::filter(pct_days >= uniform_share) |>
    dplyr::pull(mean_area_resid) |>
    range(na.rm = TRUE)

  # A single well-sampled sector (found 2026-08-10: Southern Brittany's wave
  # direction is 80% from due west, leaving only one sector at/above
  # uniform_share) collapses well_sampled_range to one repeated value, and
  # every under-sampled sector then gets clamped TO that value -- every
  # sector's plotted fill becomes identical, rendering as a single-colour
  # legend instead of a gradient. Skip clamping in that case too (not just
  # the non-finite case) and fall back to each sector's own raw value.
  well_sampled_range_valid <- all(is.finite(well_sampled_range)) && diff(well_sampled_range) > 0

  df_summary |>
    dplyr::mutate(mean_area_resid_plot = if (well_sampled_range_valid) {
      dplyr::case_when(
        pct_days < uniform_share & mean_area_resid > well_sampled_range[2] ~ well_sampled_range[2],
        pct_days < uniform_share & mean_area_resid < well_sampled_range[1] ~ well_sampled_range[1],
        TRUE ~ mean_area_resid
      )
    } else {
      mean_area_resid
    })
}

# Direction/magnitude rose for wind or wave or currents,
# coloured by the flow-controlled plume-area response.
# fill_limits lets a caller (plot_driver_rose_diagram()) impose a shared
# colour scale across several calls (e.g. one zone's wind/wave/current
# roses); NULL keeps this call's own auto-scaled range. show_legend = FALSE
# suppresses this panel's own legend, for callers that draw one shared
# legend externally instead.
plot_driver_rose <- function(driver_name, meta, n_sectors = 8, df_flow = NULL, fill_limits = NULL, show_legend = TRUE){
  driver_name <- match.arg(driver_name, c("wind", "wave", "current"))

  df_summary <- compute_driver_rose_summary(driver_name, meta, n_sectors, df_flow)

  # A blank rose would be misleading (looks like "no data collected" rather
  # than "direction not available"), so say so explicitly instead.
  if(nrow(df_summary) == 0){
    pl <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = paste0(toupper(driver_name), " direction\nnot available\nfor this zone"),
               size = 7, colour = "grey40") +
      xlim(0, 1) + ylim(0, 1) +
      theme_void()

    ggsave(filename = paste0("figures/driver_comparison/rose_", driver_name, "_plume_", meta$mouth_name, ".png"),
          plot = pl, width = 7, height = 6, dpi = 300)
    return(invisible(pl))
  }

  sector_width <- 360 / n_sectors
  compass_breaks <- seq(0, 360 - sector_width, by = max(sector_width, 45))
  compass_labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")[seq_along(compass_breaks)]

  pl <- ggplot(df_summary, aes(x = sector, y = pct_days, fill = mean_area_resid_plot)) +
    geom_col(width = sector_width * 0.9, colour = "grey30", linewidth = 0.2) +
    # geom_col() already centres each bar on its own x value (sector's
    # midpoint, e.g. 0 degrees for N), so no extra rotation is needed to
    # align bars with the compass labels below -- the half-sector `start`
    # offset previously here was based on the opposite (wrong) assumption
    # that `sector` was a bin's left edge, and rotated the whole rose one
    # half-sector counter-clockwise, putting N to the left of 12 o'clock
    # instead of at it (found 2026-08-10, visible in every panel of the
    # driver_rose_diagram figure).
    coord_polar(start = 0) +
    scale_x_continuous(breaks = compass_breaks, labels = compass_labels, limits = c(0, 360)) +
    scale_fill_gradient2(low = "steelblue", mid = "grey90", high = "firebrick", midpoint = 0,
                        name = "Plume-area\nresidual (km²)", limits = fill_limits) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = "grey85"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(1.1, "cm"),
          legend.position = if (show_legend) "right" else "none")

  ggsave(filename = paste0("figures/driver_comparison/rose_", driver_name, "_plume_", meta$mouth_name, ".png"),
        plot = pl, width = 7, height = 6, dpi = 300)
  invisible(pl)
}

# Flow-controlled plume-area residual vs. wave height, coloured by on/off-
# shore wind category
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
    df$wind_spd < 3       ~ "calm (<3 m s⁻¹)",
    df$direction == "off" ~ "offshore",
    TRUE                  ~ "onshore"
  )
  df$wind_category <- factor(df$wind_category, levels = c("calm (<3 m s⁻¹)", "onshore", "offshore"))
  category_colours <- c("calm (<3 m s⁻¹)" = "grey50", "onshore" = "steelblue", "offshore" = "firebrick")

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
# methodology of Sutton et al. (2022).
# Extracted from driver_plume_trend() (below) so it can be reused directly on
# a bare (date, value) series -- e.g. plume shape or centroid drift -- without
# needing a paired driver series. driver_plume_trend() calls this internally;
# behaviour there is unchanged.
# fit_wls_hac_trend("ar", df$compactness, df$date)
ar_weights_func <- function(val_col, start_year, time_step){
  ts_obj <- ts(zoo::na.approx(val_col), frequency = time_step, start = c(start_year, 1))
  ar_model <- ar(ts_obj, order.max = 1)
  phi_est <- ar_model$ar
  # order.max = 1 lets ar() select order 0 (no AR component) by AIC when the
  # series shows no significant AR(1) structure, leaving phi_est empty --
  # not observed on the daily series this was originally always called with
  # (autocorrelation there is essentially always significant), but real for
  # short annual series (e.g. compute_octant_trend()). Fall back to
  # unweighted in that case rather than erroring.
  # sqrt(phi_est) below assumes a positive AR(1) coefficient -- true for
  # every persistent daily geophysical series this was originally used on,
  # but a short annual series (e.g. compute_octant_trend()) can legitimately
  # show negative/oscillating autocorrelation, which would otherwise produce
  # NaN weights and crash lm(). Same unweighted fallback as the order-0 case
  # above.
  if(length(phi_est) == 0 || phi_est <= 0) return(rep(1, length(val_col)))
  sigma_est <- sqrt(phi_est)
  error_variance <- sigma_est^2 / (1 - phi_est^2)
  weights <- rep((1 / (error_variance^2)), length(val_col))
  if(!all(is.finite(weights))) return(rep(1, length(val_col)))
  weights
}
# Day-of-year climatological de-seasoning, extracted as a standalone 
# single-series helper so it can be applied to any daily plume-property 
# series (e.g. compactness, along-coast centroid position) that doesn't 
# have a paired driver series to deseason alongside it. 
# driver_plume_trend() above does the same day-of-year adjustment
# inline for the paired driver/plume case, but couldn't be reused directly
# for a single series. Returns the de-seasoned series, same length/order as input.
deseason_doy <- function(value, date){
  doy <- yday(date)
  doy <- ifelse(!leap_year(date) & doy >= 60, doy + 1L, doy)
  resid <- residuals(lm(value ~ date, na.action = na.exclude))
  doy_clim <- tibble::tibble(doy = doy, resid = resid) |>
    dplyr::summarise(resid_doy = mean(resid, na.rm = TRUE), .by = "doy") |>
    dplyr::mutate(resid_doy_clim = resid_doy - mean(resid_doy, na.rm = TRUE))
  tibble::tibble(doy = doy, value = value) |>
    dplyr::left_join(doy_clim, by = "doy") |>
    dplyr::mutate(value_adj = value - resid_doy_clim) |>
    dplyr::pull(value_adj)
}

fit_wls_hac_trend <- function(weight_choice, val_col, date_col, time_step = NULL){
  start_year <- year(min(date_col))
  if(is.null(time_step)) time_step <- if(length(val_col) < 1000) 12 else 365
  weights <- switch(weight_choice,
                    ar = ar_weights_func(val_col, start_year, time_step),
                    rep(1, length(val_col)))
  lm_model <- lm(val_col ~ date_col, weights = weights)
  lm_model_HAC <- coeftest(lm_model, vcov = vcovHAC(lm_model))
  tibble::tibble(n = length(val_col), time_step = time_step, start_year = start_year,
                 weight_choice = weight_choice, intercept = lm_model_HAC[1, 1],
                 slope = lm_model_HAC[2, 1], slope_se = lm_model_HAC[2, 2],
                 slope_t = lm_model_HAC[2, 3], slope_p = lm_model_HAC[2, 4])
}

# Per-calendar-month linear trend (sec:seasonal_methods), reusing the same
# de-seasoning (deseason_doy()) and AR(1)/HAC estimator (fit_wls_hac_trend())
# as the annual trend analysis (sec:linear_trends), refit separately within
# each of the 12 calendar months rather than once across the full year.
compute_monthly_trend <- function(value, date, min_n = 30){
  value_adj <- deseason_doy(value, date)
  month_vec <- month(date)
  purrr::map_dfr(1:12, function(m){
    idx <- which(month_vec == m & !is.na(value_adj))
    if(length(idx) < min_n){
      return(tibble::tibble(month = m, n = length(idx), time_step = NA_real_, start_year = NA_real_,
                            weight_choice = "ar", intercept = NA_real_, slope = NA_real_,
                            slope_se = NA_real_, slope_t = NA_real_, slope_p = NA_real_))
    }
    fit_wls_hac_trend("ar", value_adj[idx], date[idx]) |>
      dplyr::mutate(month = m, .before = 1)
  })
}

# Annual compass-octant occupancy trend (sec:linear_trends). Direction is
# circular (0-360 degrees), so a raw-angle OLS trend has no principled
# meaning across the 0/360 wrap and would mask bimodal regime shifts -- this
# instead bins each day into compass_octant()'s 8 categories (the same
# categorical treatment already used for the GAM/GLM/RF direction
# predictors, driver_interactions.R::build_driver_matrix()) and trends each
# octant's annual proportion-of-days using the same fit_wls_hac_trend()
# engine as every other trend in this pipeline, with an explicit time_step =
# 1 (the auto-detected 12/365 heuristic assumes a daily-resolution series,
# wrong for this annual one-point-per-year series). Octants below min_years
# of data or min_occurrence average share are returned as NA rows rather
# than fit, mirroring compute_monthly_trend()'s min_n guard -- a rarely
# occurring octant (e.g. Southern Brittany's wave direction is ~80% west,
# multi.R:444-446) would otherwise produce a spurious-looking "significant"
# trend from a handful of days.
compute_octant_trend <- function(degrees, date, min_years = 15, min_occurrence = 0.01){
  octant <- compass_octant(degrees)
  labels <- levels(octant)
  daily <- tibble::tibble(year = year(date), octant = octant) |>
    dplyr::filter(!is.na(octant))

  annual_totals <- daily |> dplyr::summarise(n_total = dplyr::n(), .by = "year")

  annual_counts <- daily |>
    dplyr::summarise(n_days = dplyr::n(), .by = c("year", "octant")) |>
    tidyr::complete(year = annual_totals$year, octant = labels, fill = list(n_days = 0)) |>
    dplyr::left_join(annual_totals, by = "year") |>
    dplyr::mutate(proportion = n_days / n_total)

  purrr::map_dfr(labels, function(oct){
    df_oct <- dplyr::filter(annual_counts, octant == oct) |> dplyr::arrange(year)
    n_years <- nrow(df_oct)
    mean_occurrence <- mean(df_oct$proportion, na.rm = TRUE)

    if(n_years < min_years || is.na(mean_occurrence) || mean_occurrence < min_occurrence){
      return(tibble::tibble(octant = oct, n_years = n_years, mean_occurrence = mean_occurrence,
                            time_step = NA_real_, start_year = NA_real_, weight_choice = "ar",
                            intercept = NA_real_, slope = NA_real_, slope_se = NA_real_,
                            slope_t = NA_real_, slope_p = NA_real_))
    }

    fit_wls_hac_trend("ar", df_oct$proportion, as.Date(paste0(df_oct$year, "-07-01")), time_step = 1) |>
      dplyr::mutate(octant = oct, n_years = n_years, mean_occurrence = mean_occurrence, .before = 1)
  })
}

# Monthly companion to compute_octant_trend() (sec:seasonal_methods), for
# the calendar-month-grouped seasonal analysis convention used throughout
# this pipeline (balanced year-count per month, unlike a calendar-boundary-
# crossing season). Unlike compute_monthly_trend(), which de-seasons a
# continuous daily series once and then fits OLS directly on that series'
# per-month day-level subset, direction is categorical -- there is no
# continuous magnitude to de-season -- so this aggregates first (each
# calendar month's annual proportion-of-days-in-octant, one point per year)
# and trends that aggregated series, reusing the same fit_wls_hac_trend()
# engine. Same min_years/min_occurrence guards as compute_octant_trend().
compute_monthly_octant_trend <- function(degrees, date, min_years = 10, min_occurrence = 0.01){
  octant <- compass_octant(degrees)
  labels <- levels(octant)
  daily <- tibble::tibble(year = year(date), month = month(date), octant = octant) |>
    dplyr::filter(!is.na(octant))

  monthly_totals <- daily |> dplyr::summarise(n_total = dplyr::n(), .by = c("year", "month"))

  monthly_counts <- daily |>
    dplyr::summarise(n_days = dplyr::n(), .by = c("year", "month", "octant")) |>
    tidyr::complete(tidyr::nesting(year, month), octant = labels, fill = list(n_days = 0)) |>
    dplyr::left_join(monthly_totals, by = c("year", "month")) |>
    dplyr::mutate(proportion = n_days / n_total)

  purrr::map_dfr(1:12, function(m){
    purrr::map_dfr(labels, function(oct){
      df_oct <- dplyr::filter(monthly_counts, month == m, octant == oct) |> dplyr::arrange(year)
      n_years <- nrow(df_oct)
      mean_occurrence <- mean(df_oct$proportion, na.rm = TRUE)

      if(n_years < min_years || is.na(mean_occurrence) || mean_occurrence < min_occurrence){
        return(tibble::tibble(month = m, octant = oct, n_years = n_years, mean_occurrence = mean_occurrence,
                              time_step = NA_real_, start_year = NA_real_, weight_choice = "ar",
                              intercept = NA_real_, slope = NA_real_, slope_se = NA_real_,
                              slope_t = NA_real_, slope_p = NA_real_))
      }

      fit_wls_hac_trend("ar", df_oct$proportion,
                        as.Date(paste0(df_oct$year, "-", sprintf("%02d", m), "-15")), time_step = 1) |>
        dplyr::mutate(month = m, octant = oct, n_years = n_years, mean_occurrence = mean_occurrence, .before = 1)
    })
  })
}

# Along-coast projection of the SPM-weighted plume centroid
# The along-coast direction is estimated per zone as
# the first principal component of the centroid's own long-term scatter (in
# local km, relative to the river mouth), then each day's centroid is
# projected onto that axis. Shared by func/analysis/compute_seasonal_trend.R,
# func/figure.R::plot_seasonal_boxplot_heatmap(), and
# func/analysis/compute_shape_alongcoast_trend.R (which used to carry its own
# near-identical local compute_alongcoast() before it was deduplicated onto
# this shared version -- verified to reproduce identical output first,
# since that script feeds the published panache_stats_table).
# compute_alongcoast_ts("GULF_OF_LION", get_zone_meta(zone_name = "GULF_OF_LION"), "output/panache/dynamic")
compute_alongcoast_ts <- function(zone, meta, plume_dir){
  df <- read_csv(paste0(plume_dir, "/", zone, "/Results.csv"), show_col_types = FALSE) |>
    dplyr::mutate(date = as.Date(date)) |>
    dplyr::filter(.data$river == "ALL") |>
    dplyr::filter(!is.na(lon_weighted_centroid_of_the_plume_area), !is.na(lat_weighted_centroid_of_the_plume_area))

  lat0 <- meta$mouth_lat; lon0 <- meta$mouth_lon
  x_km <- (df$lon_weighted_centroid_of_the_plume_area - lon0) * 111.32 * cos(pi / 180 * lat0)
  y_km <- (df$lat_weighted_centroid_of_the_plume_area - lat0) * 111.32

  pc <- prcomp(cbind(x_km, y_km))
  axis1 <- pc$rotation[, 1]
  # Fix an arbitrary PCA sign convention: orient axis1 so its dominant
  # component (whichever of east-west/north-south carries more of the
  # loading) is positive, so "positive slope"/"positive value" has a stable,
  # reportable geographic meaning (east or north) instead of flipping
  # arbitrarily between calls.
  dominant_is_x <- abs(axis1[1]) >= abs(axis1[2])
  if((dominant_is_x && axis1[1] < 0) || (!dominant_is_x && axis1[2] < 0)) axis1 <- -axis1

  alongcoast_km <- as.numeric(cbind(x_km, y_km) %*% axis1)

  tibble::tibble(date = df$date, value = alongcoast_km) |>
    complete(date = seq(min(date), max(date), by = "day")) |>
    zoo::na.trim()
}

# Needed when `df`'s "plume_area" column
# actually holds a different metric_col (see combine_plume_driver()), since the plot
# filename is keyed only on driver_name/mouth_name and would otherwise silently
# overwrite the existing plume-area trend figure for that driver/mouth.
driver_plume_trend <- function(df, driver_name, mouth_name, end_date = NULL, save_plot = TRUE, plot_label = mouth_name){

  disp <- dplyr::filter(driver_display, driver_name == !!driver_name)

  # doy normalised to a common 1-366 scale across leap/non-leap years
  df <- df |>
    dplyr::mutate(doy = yday(date), doy = ifelse(!leap_year(date) & doy >= 60, doy + 1L, doy),
                  month = month(date), year = year(date))
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

  wls_driver_daily   <- plyr::ldply(c("ar", "none"), fit_wls_hac_trend, val_col = df_daily$driver_doy_adj, date_col = df_daily$date)
  wls_driver_monthly <- plyr::ldply(c("ar", "none"), fit_wls_hac_trend, val_col = df_monthly$driver_monthly_adj, date_col = df_monthly$date)
  wls_plume_daily    <- plyr::ldply(c("ar", "none"), fit_wls_hac_trend, val_col = df_daily$plume_doy_adj, date_col = df_daily$date)
  wls_plume_monthly  <- plyr::ldply(c("ar", "none"), fit_wls_hac_trend, val_col = df_monthly$plume_monthly_adj, date_col = df_monthly$date)

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
                aes(x = x_daily, y = y_plume, label = paste0("Plume area slope = ", round(slope_annualised, 2), " km² yr⁻¹\n",
                                                              "p-value = ", round(slope_p, 2)))) +
      geom_label(data = dplyr::filter(trend_labels_plume, timestep == "daily"), size = 5, hjust = 0, colour = "darkblue",
                aes(x = x_monthly, y = y_plume, label = paste0("Plume area slope = ", round(slope_annualised, 2), " km² yr⁻¹\n",
                                                                "p-value = ", round(slope_p, 2)))) +
      labs(x = NULL, y = "Plume area [km²]",
          title = paste0(plot_label, " : plume area after statistical treatment (vs. ", driver_name, ")"),
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
                aes(x = x_daily, y = y_driver, label = paste0(disp$driver_label, " slope = ", round(slope_annualised, 2), " yr⁻¹\n",
                                                              "p-value = ", round(slope_p, 2)))) +
      geom_label(data = dplyr::filter(trend_labels_driver, timestep == "daily"), size = 5, hjust = 0, colour = "darkblue",
                aes(x = x_monthly, y = y_driver, label = paste0(disp$driver_label, " slope = ", round(slope_annualised, 2), " yr⁻¹\n",
                                                                "p-value = ", round(slope_p, 2)))) +
      labs(x = NULL, y = disp$driver_label,
          title = paste0(plot_label, " : ", driver_name, " after statistical treatment"),
          subtitle = "Red = adjusted monthly values; blue = adjusted daily values; purple = original data") +
      theme(panel.border = element_rect(fill = NA, colour = "black"))

    pl_combi <- ggpubr::ggarrange(pl_plume, pl_driver, ncol = 1, nrow = 2)
    plot_label_file <- stringr::str_replace_all(plot_label, " ", "_")
    ggsave(filename = paste0("figures/driver_comparison/trends_plume_", driver_name, "_adj_", plot_label_file, ".png"),
          pl_combi, width = 12, height = 10)
  }

  return(stats)
}


