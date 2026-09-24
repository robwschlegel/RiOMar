# func/generate_x11_driver_correlation_heatmap.R
#
# Manuscript slot "x11_driver_correlation_heatmap" -- see
# manuscript/figure_table_registry.csv. Replaces the originally-planned
# ~10-page per-driver seasonal/interannual X11 time-series comparison
# Supplementary section (manuscript/TODO.md, superseded 2026-09) with a
# single compact heatmap of Pearson r between plume area's X11 Seasonal (and,
# separately, Interannual) component and each of 5 drivers' own component,
# per zone.
#
# Follows the same split already used for river flow's own X11 comparison
# (func/X11.R::make_the_plot(), func/figure.R::plot_x11_river_and_plume()):
# Python computes and persists the X11 decomposition, R reads the CSVs,
# computes the correlation via cor(), and does all the plotting.
#
# wind/tide/wave/current: reads output/panache/dynamic/<Zone>/X11_ANALYSIS/
# <driver>/<driver>_WEEKLY.csv, written by
# func/compute_x11_driver_signals.py -- run that first (or here
# automatically if missing).
# flow: reads the same weekly X11 output CSVs the main-text Fig 6/7
# r-value annotations already use, via X11.R::get_X11_data() -- not
# recomputed, so the heatmap's flow cells match the main-text prose exactly.
#
# Run from repo root: Rscript func/generate_x11_driver_correlation_heatmap.R

source("func/X11.R")  # loads tidyverse etc. and sources multi.R (which sources util.R) itself; also provides get_X11_data()

OTHER_DRIVERS <- c("wind", "tide", "wave", "current")
X11_DIR <- "output/panache/dynamic"

driver_csv_path <- function(zone_name, driver_name){
  file.path(X11_DIR, zone_name, "X11_ANALYSIS", driver_name, paste0(driver_name, "_WEEKLY.csv"))
}

if (!all(file.exists(unlist(lapply(zones, function(z) driver_csv_path(z, OTHER_DRIVERS)))))) {
  system("python func/compute_x11_driver_signals.py")
}

# One row per zone x driver x component (seasonal/interannual); a missing
# CSV (X11 decomposition failed its quality cutoff for that zone/driver, see
# compute_x11_driver_signals.py) yields NA rather than an error, rendered as
# a distinct grey cell below rather than silently dropped.
heat_data <- purrr::map_dfr(zones, function(zone_name){

  # River flow: same source as the main-text Fig 6/7 r-value annotations.
  X11_data <- get_X11_data(where_are_saved_X11_results = X11_DIR, Zone = zone_name,
                           Data_source = "SEXTANT", sensor_name = "merged",
                           atmospheric_correction = "Standard", Temporal_resolution = "WEEKLY")
  flow_row <- tibble::tibble(
    zone = zone_name, driver = "flow",
    component = c("seasonal", "interannual"),
    r = c(cor(X11_data$Seasonal_signal_plume_area, X11_data$Seasonal_signal_river_flow, use = "complete.obs"),
         cor(X11_data$Interannual_signal_plume_area, X11_data$Interannual_signal_river_flow, use = "complete.obs")))

  other_rows <- purrr::map_dfr(OTHER_DRIVERS, function(driver_name){
    plume_path <- file.path(X11_DIR, zone_name, "X11_ANALYSIS", "area_of_the_plume_mask_in_km2",
                            "SEXTANT_merged_Standard_WEEKLY.csv")
    csv_path <- driver_csv_path(zone_name, driver_name)
    if (!file.exists(plume_path) || !file.exists(csv_path)) {
      return(tibble::tibble(zone = zone_name, driver = driver_name,
                            component = c("seasonal", "interannual"), r = NA_real_))
    }
    # plume_path's CSV (written by the frozen X11.py::apply_X11_method_and_save_results())
    # uses "dates" (plural); csv_path's CSV (func/compute_x11_driver_signals.py,
    # this session's own new script) uses "date" (singular) -- align here
    # rather than change either writer.
    plume_ts <- readr::read_csv(plume_path, show_col_types = FALSE) |> dplyr::rename(date = dates)
    driver_ts <- readr::read_csv(csv_path, show_col_types = FALSE)
    merged <- dplyr::inner_join(plume_ts, driver_ts, by = "date", suffix = c("_plume", "_driver"))
    tibble::tibble(
      zone = zone_name, driver = driver_name,
      component = c("seasonal", "interannual"),
      r = c(cor(merged$Seasonal_signal_plume, merged$Seasonal_signal_driver, use = "complete.obs"),
           cor(merged$Interannual_signal_plume, merged$Interannual_signal_driver, use = "complete.obs")))
  })

  dplyr::bind_rows(flow_row, other_rows)
})

readr::write_csv(heat_data, "output/STATS/driver_x11_seasonal_interannual_r.csv")

driver_display <- c(flow = "River flow", wind = "Wind speed", tide = "Tidal range",
                    wave = "Wave height", current = "Current speed")
zone_labels <- zone_title(rev(zones))
zone_labels[zone_labels == "Southern Brittany"] <- "S. Brittany"

DRIVER_ORDER <- c("flow", "wave", "wind", "current", "tide")  # Robert's call, 2026-09-24

heat_data <- heat_data |>
  dplyr::mutate(
    zone = factor(zone, levels = rev(zones), labels = zone_labels),
    driver = factor(driver, levels = DRIVER_ORDER, labels = unname(driver_display[DRIVER_ORDER])),
    component = factor(component, levels = c("seasonal", "interannual"),
                       labels = c("Seasonal component", "Interannual component")))

# Symmetric-about-zero diverging scale, same as the driver-rose figure
# (func/multi.R's scale_fill_gradient2(low="steelblue", ..., midpoint=0)) --
# the natural fit for a correlation, unlike Figure 4's ratio-centred
# purple/orange scale (midpoint=1).
p_heatmap <- ggplot(heat_data, aes(x = driver, y = zone, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_label(aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))),
            size = 3.5, fill = "white", linewidth = 0, label.padding = unit(0.12, "lines")) +
  facet_wrap(~component, ncol = 2) +
  scale_fill_gradient2(name = "Pearson r", low = "steelblue", mid = "grey90", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey80") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
       axis.text.y = element_text(size = 10), panel.grid = element_blank())

output_subdir <- get_registry_row("x11_driver_correlation_heatmap")$output_subdir
figure_dir <- file.path("figures", "ARTICLE", output_subdir)
save_plot_as_png(p_heatmap, registry_basename(output_subdir), width = 9, height = 4.2, path = figure_dir)
message("Wrote ", registry_filename(output_subdir))
