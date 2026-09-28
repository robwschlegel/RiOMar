# Loading -----------------------------------------------------------------

# Zone-specific hard ceiling on classified plume area (km^2), applied
# unconditionally in load_plume_ts() below regardless of which metric_col is
# requested.
plume_area_ceiling <- c(BAY_OF_BISCAY = 12000, BAY_OF_SEINE = 2500,
                        GULF_OF_LION = 10000, SOUTHERN_BRITTANY = 6000)

# Load time series of plume values. `metric_col` selects which Results.csv
# column becomes `plume_area` downstream (every driver_plume_* function in
# multi.R is written against that column name regardless of what it holds,
# so this is the one place a different plume metric -- e.g. mass_SPM_in_the_
# plume_area_in_t -- needs to be wired in). `outlier_max` screens the
# `plume_area` column only, dropping values above the threshold as NA; it
# defaults to NULL (no screening). Pass an explicit value (e.g. 20000, a
# physically implausible km^2 area for these zones) at call sites that need
# the guard -- do not rely on an implicit default, since a metric other than
# area may not share the same physically sensible ceiling.
load_plume_ts <- function(zone, plume_dir = "output/panache/dynamic",
                          metric_col = "area_of_the_plume_mask_in_km2", outlier_max = NULL,
                          river = "ALL"){
  file_name <- paste0(plume_dir, "/", zone, "/Results.csv")
  suppressMessages({
    df_plume <- read_csv(file_name) |>
      dplyr::mutate(date = as.Date(date)) |>
      dplyr::filter(.data$river == !!river) |>
      dplyr::select(-river) |>
      dplyr::select(date:confidence_index_in_perc) |>
      complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |>
      dplyr::mutate(plume_area_ceiling_exceeded = area_of_the_plume_mask_in_km2 > plume_area_ceiling[[zone]],
                    dplyr::across(-c(date, plume_area_ceiling_exceeded),
                                  ~ ifelse(plume_area_ceiling_exceeded, NA, .))) |>
      dplyr::select(-plume_area_ceiling_exceeded) |>
      dplyr::rename(plume_area = !!rlang::sym(metric_col))
    if(!is.null(outlier_max)) df_plume <- dplyr::mutate(df_plume, plume_area = ifelse(plume_area > outlier_max, NA, plume_area))
    df_plume <- df_plume |>
      zoo::na.trim() |>
      mutate(zone = zone, .before = "date")
  })
  return(df_plume)
}

# Get the (lon, lat) footprint of the daily plume mask for a zone
load_plume_surface <- function(zone, plume_dir = "output/panache/dynamic"){
  file_name <- paste0(plume_dir, "/", zone, "/PlumeMasks.nc")
  nc_dat <- nc_open(file_name)
  lon  <- ncvar_get(nc_dat, "lon")
  lat  <- ncvar_get(nc_dat, "lat")
  time <- .nc_time_to_date(nc_dat, "time")
  # plume_mask dims are [lon, lat, river, time]; keep only the combined 'ALL'
  # river layer (panache v5.0.0+ adds one mask layer per individual river
  # plus this union layer -- see plume_algorithm.py::_stack_river_masks())
  river_names <- ncvar_get(nc_dat, "river")
  idx_all <- which(river_names == "ALL")
  mask <- ncvar_get(nc_dat, "plume_mask",
                    start = c(1, 1, idx_all, 1),
                    count = c(-1, -1, 1, -1))  # [lon, lat, time]
  nc_close(nc_dat)

  hit <- which(mask == 1, arr.ind = TRUE)
  df_plume_surface <- tibble::tibble(zone = zone,
                                     date = time[hit[, 3]],
                                     lon  = lon[hit[, 1]],
                                     lat  = lat[hit[, 2]])
  message(paste0("Loaded ", nrow(df_plume_surface), " plume surface points for zone ", zone))
  return(df_plume_surface)
}

# Load river flow data
# dir_name <- "data/RIVER_FLOW/BAY_OF_SEINE"
load_river_flow <- function(dir_name){

  # Get file list (non-recursive: only files directly in dir_name, not old/ or HydroPortail/)
  files_to_load <- list.files(path = dir_name, pattern = "\\.(csv)$", full.names = TRUE)

  # All current files have a two-column header: date,debit
  data_list <- lapply(files_to_load, function(file){
    df <- read.csv(file, header = TRUE)
    df$date <- as.Date(df$date)
    df[, c("flow", "date")] <- list(df$debit, df$date)
    df[, c("flow", "date")]
  })

  names(data_list) <- basename(files_to_load)

  # Combine all dataframes and aggregate by date
  dplyr::bind_rows(data_list) |>
    dplyr::summarise(flow = sum(flow, na.rm = TRUE),
                     n_rivers = dplyr::n(), .by = "date") |>
    dplyr::filter(n_rivers == length(data_list))
}

# Load one river's flow, or a small named group of rivers sharing a mouth/
# estuary (e.g. Gironde = Garonne + Dordogne, Sevre = Sevre Niortaise + Lay)
# summed together -- same per-file read/sum logic as load_river_flow(),
# just restricted to specific files instead of every CSV in the zone
# directory. Used to match panache's now-individual river-mouth plume output
# (see metadata/river_discharge_mapping.csv) against its own discharge.
# load_river_flow_single("GULF_OF_LION", "grand_rhone")
# load_river_flow_single("BAY_OF_BISCAY", c("garonne", "dordogne"))
load_river_flow_single <- function(zone, river_slugs){
  dir_name <- file.path("data/RIVER_FLOW", zone)
  files_to_load <- file.path(dir_name, paste0(river_slugs, ".csv"))

  data_list <- lapply(files_to_load, function(file){
    df <- read.csv(file, header = TRUE)
    df$date <- as.Date(df$date)
    df[, c("flow", "date")] <- list(df$debit, df$date)
    df[, c("flow", "date")]
  })

  dplyr::bind_rows(data_list) |>
    dplyr::summarise(flow = sum(flow, na.rm = TRUE),
                     n_rivers = dplyr::n(), .by = "date") |>
    dplyr::filter(n_rivers == length(files_to_load)) |>
    dplyr::select(-n_rivers)
}

# Load tide gauge data
load_tide_gauge <- function(dir_name){

  station <- basename(dir_name)
  df_tide <- .load_tide_raw(dir_name)

  # Flag calendar days whose sub-daily curve does not look tidal (see
  # qc_tide_days() / func/tide.R) before tide_mean/tide_range are computed
  df_flags <- qc_tide_days(df_tide, station)

  df_tide_daily <- df_tide |>
    mutate(date = as.Date(t)) |>
    dplyr::summarise(tide_mean = round(mean(tide, na.rm = TRUE), 2),
                     tide_range = max(tide, na.rm = TRUE)-min(tide, na.rm = TRUE), .by = "date") |>
    dplyr::left_join(df_flags, by = "date") |>
    dplyr::mutate(tide_bad = dplyr::coalesce(tide_bad, TRUE),
                  tide_mean = ifelse(tide_bad, NA, tide_mean),
                  tide_range = ifelse(tide_bad, NA, tide_range)) |>
    dplyr::select(date, tide_mean, tide_range, tide_qc_reason = reason)
  return(df_tide_daily)
}

# Compute speed and compass bearing from eastward (u) and northward (v) vector
# components. convention = "from" reports the direction the vector originates
# from (meteorological convention, e.g. wind); convention = "to" reports the
# direction the vector is heading towards (oceanographic convention, e.g.
# currents). Used by load_wind_sub() and load_surface_current().
.speed_direction <- function(u, v, convention = c("from", "to")){
  convention <- match.arg(convention)
  speed <- sqrt(u^2 + v^2)
  bearing_to <- (90 - atan2(v, u) * (180 / pi)) %% 360
  direction <- if(convention == "from") (bearing_to + 180) %% 360 else bearing_to
  list(speed = speed, direction = direction)
}

# Decode a NetCDF file's numeric time variable into Date, following the CF
# "<seconds|hours|days> since <origin>" units convention used by every
# NetCDF source read below.
.nc_time_to_date <- function(nc, time_var = "time"){
  units_str <- ncatt_get(nc, time_var, "units")$value
  m <- regmatches(units_str, regexec("^(seconds|hours|days) since (.+)$", units_str))[[1]]
  if(length(m) != 3) stop("Unrecognised time units string: ", units_str)
  multiplier <- switch(m[2], seconds = 1, hours = 3600, days = 86400)
  raw <- as.vector(ncvar_get(nc, time_var))
  as.Date(as.POSIXct(raw * multiplier, origin = m[3], tz = "UTC"))
}

# Read one or more variables from a NetCDF file, cropped to a lon/lat box, as
# a long (lon, lat, date, <var>...) tibble
.nc_read_box <- function(file_name, var_names, lon_range, lat_range,
                         lon_var = "longitude", lat_var = "latitude", time_var = "time"){
  nc <- nc_open(file_name)
  on.exit(nc_close(nc))

  lon <- ncvar_get(nc, lon_var)
  lat <- ncvar_get(nc, lat_var)
  lon_idx <- which(lon >= lon_range[1] & lon <= lon_range[2])
  lat_idx <- which(lat >= lat_range[1] & lat <= lat_range[2])
  date <- .nc_time_to_date(nc, time_var)

  grid <- expand.grid(lon = lon[lon_idx], lat = lat[lat_idx], date = date)

  for(v in var_names){
    n_dims <- length(nc$var[[v]]$dim)
    if(n_dims < 3) stop("Unexpected variable dimensionality (", n_dims, "D) for ", v, " in ", file_name)
    n_middle <- n_dims - 3  # any dims between lat and time (e.g. GLORYS's single-level depth)
    start <- c(min(lon_idx), min(lat_idx), rep(1, n_middle), 1)
    count <- c(length(lon_idx), length(lat_idx), rep(1, n_middle), -1)
    arr <- ncvar_get(nc, v, start = start, count = count)
    grid[[v]] <- as.vector(arr)
  }

  tibble::as_tibble(grid)
}

# Load wind data
load_wind_sub <- function(file_name, lon_range, lat_range){
  wind_df <- .nc_read_box(file_name, c("eastward_wind", "northward_wind"), lon_range, lat_range) |>
    dplyr::rename(u = eastward_wind, v = northward_wind) |>
    dplyr::select(date, lon, lat, u, v) |>
    dplyr::summarise(u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE), .by = "date")

  # Remove final day of data
  ## it is an artefact from creating daily integrals from hourly data
  final_date <- max(wind_df$date)
  wind_df <- filter(wind_df, date != final_date)

  # Zone-average wind speed and direction (direction = where the wind is coming FROM)
  wind_vec <- .speed_direction(wind_df$u, wind_df$v, convention = "from")
  wind_df$wind_spd <- round(wind_vec$speed, 2)
  wind_df$wind_dir <- round(wind_vec$direction)
  return(wind_df)
}

# Load wave data (significant wave height + mean direction)
load_wave <- function(file_name, lon_range, lat_range){
  nc <- nc_open(file_name)
  has_wave_dir <- "VMDR" %in% names(nc$var)
  nc_close(nc)

  var_names <- if(has_wave_dir) c("VHM0", "VMDR") else "VHM0"
  wave_df_raw <- .nc_read_box(file_name, var_names, lon_range, lat_range) |>
    dplyr::rename(wave_height = VHM0)
  if(has_wave_dir) wave_df_raw <- dplyr::rename(wave_df_raw, wave_dir = VMDR) else wave_df_raw$wave_dir <- NA_real_

  wave_df <- wave_df_raw |>
    dplyr::select(date, lon, lat, wave_height, wave_dir) |>
    dplyr::summarise(wave_height = mean(wave_height, na.rm = TRUE),
                     wave_dir_x = mean(sin(wave_dir * pi / 180), na.rm = TRUE),
                     wave_dir_y = mean(cos(wave_dir * pi / 180), na.rm = TRUE), .by = "date") |>
    dplyr::mutate(wave_height = round(wave_height, 2),
                  wave_dir = if(has_wave_dir) round(atan2(wave_dir_x, wave_dir_y) * (180 / pi)) %% 360 else NA_real_) |>
    dplyr::select(date, wave_height, wave_dir)

  # Remove final day of data
  ## it is an artefact from creating daily integrals from hourly data
  final_date <- max(wave_df$date)
  wave_df <- filter(wave_df, date != final_date)
  return(wave_df)
}

# Load one GLORYS file's surface currents (uo/vo), reduced to the zone/box daily mean.
.load_current_sub <- function(file_name, lon_range, lat_range){
  .nc_read_box(file_name, c("uo", "vo"), lon_range, lat_range) |>
    dplyr::rename(u = uo, v = vo) |>
    dplyr::select(date, lon, lat, u, v) |>
    dplyr::summarise(u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE), .by = "date")
}

# Load GLORYS surface current (eastward/northward velocity) data for a zone
load_surface_current <- function(zone_name, lon_range, lat_range){
  dir_name <- riomar_data_path("GLORYS", zone_name)
  current_df <- dplyr::bind_rows(
    .load_current_sub(file.path(dir_name, "glorys_199301_202412.nc"), lon_range, lat_range),
    .load_current_sub(file.path(dir_name, "glorys_uo_vo_202501_202512.nc"), lon_range, lat_range)
  )

  # Zone-average current speed and direction (direction = where the current is flowing TO)
  current_vec <- .speed_direction(current_df$u, current_df$v, convention = "to")
  current_df$current_spd <- round(current_vec$speed, 2)
  current_df$current_dir <- round(current_vec$direction)
  return(current_df)
}

