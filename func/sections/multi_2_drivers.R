# Driver loading ------------------------------------------------------------

# Shared on-/off-shore wind classification
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

# Study period every data source is meant to be bounded to (Table~\ref{tab:metadata}).
# CMEMS/HydroPortail downloads request exactly this range, so their files are
# already bounded on disk; SHOM's tide web form has no such request parameter
# and returns each gauge's entire digitised archive instead (e.g. Marseille
# back to 1849, all four gauges' 2026 files already present) -- load_driver()
# below enforces the cap explicitly so a raw archive's actual span can never
# silently widen a driver's trend-fitting window.
STUDY_DATE_RANGE <- as.Date(c("1998-01-01", "2025-12-31"))

# Load one driver's daily time series for a zone, in a common two-column
# (date, value) shape so downstream functions don't need to know which
# driver they're looking at.
load_driver <- function(driver_name, meta){
  driver_name <- match.arg(driver_name, c("flow", "tide", "wind", "current", "wave"))

  if(driver_name == "flow"){
    df <- load_river_flow(paste0("data/RIVER_FLOW/", meta$zone)) |>
      dplyr::select(date, value = flow)

  } else if(driver_name == "tide"){
    df <- load_tide_gauge(paste0("data/TIDES/", meta$gauge)) |>
      dplyr::select(date, value = tide_range, tide_mean)

  } else if(driver_name == "wind"){
    zone_box <- dplyr::filter(zones_bbox, zone == meta$zone)
    lon_range <- c(zone_box$lon_min, zone_box$lon_max)
    lat_range <- c(zone_box$lat_min, zone_box$lat_max)
    wind_files <- dir(riomar_data_path("WIND", meta$zone), pattern = "_daily_", full.names = TRUE)
    df_wind <- purrr::map_dfr(wind_files, load_wind_sub, lon_range, lat_range) |>
      wind_add_direction(meta$zone)
    df <- df_wind |> dplyr::select(date, value = wind_spd, wind_dir, direction, u, v)

  } else if(driver_name == "current"){
    zone_box <- dplyr::filter(zones_bbox, zone == meta$zone)
    lon_range <- c(zone_box$lon_min, zone_box$lon_max)
    lat_range <- c(zone_box$lat_min, zone_box$lat_max)
    df_current <- load_surface_current(meta$zone, lon_range, lat_range)
    df <- df_current |> dplyr::select(date, value = current_spd, current_dir, u, v)

  } else if(driver_name == "wave"){
    zone_box <- dplyr::filter(zones_bbox, zone == meta$zone)
    lon_range <- c(zone_box$lon_min, zone_box$lon_max)
    lat_range <- c(zone_box$lat_min, zone_box$lat_max)
    wave_files <- dir(riomar_data_path("WAVE", meta$zone), pattern = "_daily_", full.names = TRUE)
    df_wave <- purrr::map_dfr(wave_files, load_wave, lon_range, lat_range)
    df <- df_wave |> dplyr::select(date, value = wave_height, wave_dir)
  }

  dplyr::filter(df, date >= STUDY_DATE_RANGE[1], date <= STUDY_DATE_RANGE[2])
}


# Combine plume + one driver --------------------------------------------

# Load plume + a single driver for one zone, join on date. This is the core
# object every comparison function below operates on.
# combine_plume_driver("flow", get_zone_meta(mouth_name = "Seine"))
# combine_plume_driver("flow", get_zone_meta(mouth_name = "Seine"), metric_col = "mass_SPM_in_the_plume_area_in_t", outlier_max = NULL)
combine_plume_driver <- function(driver_name, meta, metric_col = "area_of_the_plume_mask_in_km2", outlier_max = NULL,
                                 plume_dir = "output/panache/dynamic"){

  df_plume  <- load_plume_ts(meta$zone, plume_dir = plume_dir, metric_col = metric_col,
                             outlier_max = outlier_max)  # util.R -- already handles gap-filling + outlier removal
  df_driver <- load_driver(driver_name, meta)

  df <- dplyr::left_join(df_plume, df_driver, by = "date") |>
    zoo::na.trim()

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
# Uses util.R::lagged_correlation() (x = driver value gets lagged, y = plume
# area is the fixed reference, so a positive "lag" means the driver leads the plume
driver_plume_correlation <- function(df, max_lag_daily = 30){
  ts_list <- driver_plume_timesteps(df)

  cor_daily <- lagged_correlation(x = ts_list$daily$value, y = ts_list$daily$plume_area, max_lag_daily) |>
    dplyr::mutate(timestep = "daily")
  cor_monthly <- lagged_correlation(x = ts_list$monthly$value, y = ts_list$monthly$plume_area, 12) |>
    dplyr::mutate(timestep = "monthly")
  n_annual_lag <- max(0, min(10, nrow(ts_list$annual) - 1))
  cor_annual <- lagged_correlation(x = ts_list$annual$value, y = ts_list$annual$plume_area, n_annual_lag) |>
    dplyr::mutate(timestep = "annual")

  dplyr::bind_rows(cor_daily, cor_monthly, cor_annual) |>
    dplyr::mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))
}


