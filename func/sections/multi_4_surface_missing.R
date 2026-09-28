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


# Missing data ------------------------------------------------------------

# Get missing dates of
if(!file.exists("output/STATS/missing_SPM.csv") | !file.exists("output/STATS/missing_chla.csv")){
  message("Computing missing SEXTANT files...")
  SPM_files_NA <- data.frame(file_name = dir(riomar_data_path("SEXTANT", "SPM"), pattern = ".nc", recursive = TRUE)) |>
    mutate(base_name = basename(file_name)) |>
    separate(base_name, "-", extra = "drop") |>
    dplyr::rename(date = `-`) |>
    mutate(date = as.Date(date, format = "%Y%m%d")) |>
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |>
    filter(is.na(file_name))
  write_csv(SPM_files_NA, "output/STATS/missing_SPM.csv")
  chla_files_NA <- data.frame(file_name = dir(riomar_data_path("SEXTANT", "CHLA"), pattern = ".nc", recursive = TRUE)) |>
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
  ggsave("figures/validation/missing_SPM.png", width = 9, height = 9, dpi = 600)
  ggplot(chla_files_NA_count, aes(x = month, y = miss_count_month_year)) +
    geom_col() +
    facet_wrap(~year) +
    labs(x = NULL, y = "count", title = "Monthly count of missing chl a SEXTANT files") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  ggsave("figures/validation/missing_chla.png", width = 9, height = 9, dpi = 600)
}


# Run everything -----------------------------------------------------------

# NB: not run automatically on source() -- call explicitly
# purrr::walk(zones, surface_plot_daily_maps)


