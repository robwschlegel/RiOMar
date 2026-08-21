# One-off exploratory figure: zone x month heatmap of each plume property's
# and driver's per-calendar-month linear trend, expressed as a percent-of-
# mean change per year rather than raw units -- lets properties/drivers on
# very different scales (e.g. plume area in km^2 vs. tidal range in m) sit on
# one shared colour scale. Red = increasing, blue = decreasing.
#
# Trend: output/STATS/monthly_trend_summary.csv (func/compute_seasonal_trend.R),
# dynamic threshold, AR(1)/HAC-weighted per-calendar-month slope
# (compute_monthly_trend(), func/multi.R) -- the same slope behind Table 7 and
# the Supplementary month-by-month table. Percent-per-year here is
# 100 * (slope * 365.25) / mean, where mean is that zone x variable's overall
# mean value pooled across all months and years (dynamic threshold), not a
# per-month mean -- so all 12 monthly cells for a given zone x variable share
# one baseline and stay comparable to each other in absolute, not just
# relative, terms.
#
# The mean is recomputed here from the same source loaders
# (load_plume_ts()/load_driver()/PlumeShape.csv/compute_alongcoast_ts())
# compute_seasonal_trend.R used to produce the slopes, rather than read from
# Figure_5's monthly_boxplot_data.csv -- that file stores SPM mass in tonnes
# (Figure_5_seasonal_analysis's plume_area/1000 conversion) while
# monthly_trend_summary.csv's SPM_mass slope is fit on raw kg, and it may not
# include every variable (e.g. compactness) if it predates a PlumeShape.csv
# run. Recomputing the mean directly from the same raw series as the slope
# guarantees matching units and full variable coverage.
#
# Run from repo root: Rscript func/generate_monthly_trend_pct_heatmap.R
source("func/multi.R")
# util.R::load_tide_gauge() needs tide.R::.load_tide_raw() -- see func/figure.R's identical comment.
source("func/tide.R")

YR <- 365.25
mass_col <- "mass_SPM_in_the_plume_area_in_g_m"  # kg, see compute_mass_spm_trend.R
drivers <- c("flow", "wind", "tide", "wave", "current")

variable_display <- c(
  plume_area    = "Plume area (km²)",
  SPM_mass      = "SPM mass (t)",
  compactness   = "Compactness",
  alongcoast_km = "Along-coast drift (km)",
  flow          = "River discharge (m³ s⁻¹)",
  wind          = "Wind speed (m s⁻¹)",
  tide          = "Tidal range (m)",
  wave          = "Wave height (m)",
  current       = "Current speed (m s⁻¹)"
)

plume_dir <- "output/panache/dynamic"

overall_mean <- purrr::pmap_dfr(zone_meta, function(...){
  meta <- tibble::tibble(...)

  df_area <- load_plume_ts(meta$zone, plume_dir = plume_dir, outlier_max = 20000) |>
    dplyr::transmute(variable = "plume_area", value = plume_area)

  df_mass <- load_plume_ts(meta$zone, plume_dir = plume_dir, metric_col = mass_col, outlier_max = NULL) |>
    dplyr::transmute(variable = "SPM_mass", value = plume_area)  # kg, matches monthly_trend_summary.csv's unconverted slope

  shape_path <- paste0(plume_dir, "/", meta$zone, "/PlumeShape.csv")
  df_shape <- if(file.exists(shape_path)){
    readr::read_csv(shape_path, show_col_types = FALSE) |>
      dplyr::transmute(variable = "compactness", value = compactness)
  } else NULL

  df_coast <- compute_alongcoast_ts(meta$zone, meta, plume_dir) |>
    dplyr::transmute(variable = "alongcoast_km", value = value)

  df_drivers <- purrr::map_dfr(drivers, function(driver_name){
    load_driver(driver_name, meta) |>
      dplyr::transmute(variable = driver_name, value = value)
  })

  dplyr::bind_rows(df_area, df_mass, df_shape, df_coast, df_drivers) |>
    dplyr::mutate(zone = meta$zone, .before = 1)
}) |>
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .by = c(zone, variable))

# "Southern Brittany" abbreviated to "S. Brittany" for this axis label only
# (not zone_title() itself, which other figures/tables still use
# unabbreviated) to save left-margin white space in this 3x3 panel grid.
zone_labels <- zone_title(rev(zones))
zone_labels[zone_labels == "Southern Brittany"] <- "S. Brittany"

trend <- readr::read_csv("output/STATS/monthly_trend_summary.csv", show_col_types = FALSE) |>
  dplyr::filter(threshold == "dynamic", weight_choice == "ar") |>
  dplyr::left_join(overall_mean, by = c("zone", "variable")) |>
  dplyr::mutate(pct_per_year = 100 * (slope * YR) / mean_value,
                # geom_tile's discrete y-axis places the first factor level at the
                # bottom, so levels are reversed from ZONE_ORDER here to read
                # north (Bay of Seine) at top -> south (Gulf of Lion) at bottom,
                # top-to-bottom -- matching the project's north-to-south panel
                # convention (see ZONE_ORDER, func/util.R) used by every
                # facet_wrap multi-zone figure elsewhere.
                zone = factor(zone, levels = rev(zones), labels = zone_labels),
                month = factor(month, levels = 1:12, labels = month.abb),
                variable = factor(variable, levels = names(variable_display), labels = unname(variable_display)))

p_heatmap <- ggplot(trend, aes(x = month, y = zone, fill = pct_per_year)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = dplyr::filter(trend, slope_p < 0.05), fill = NA, colour = "black",
           linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~variable, ncol = 3) +
  scale_fill_gradient2(name = "% change\nper year", low = "steelblue4", mid = "white", high = "firebrick3", midpoint = 0) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
       axis.text.y = element_text(size = 10), panel.grid = element_blank())

if(!dir.exists("figures/driver_comparison")) dir.create("figures/driver_comparison", recursive = TRUE)
ggsave(filename = "figures/driver_comparison/monthly_trend_pct_heatmap.png", plot = p_heatmap, width = 12, height = 10, dpi = 300)
message("Wrote figures/driver_comparison/monthly_trend_pct_heatmap.png")
