# func/util.R
# The storage point for many functions re-used by other scripts


# Meta-data ---------------------------------------------------------------

# Canonical zone order, north to south, for consistent panel ordering across
# every multi-zone figure/table (top to bottom, or left to right): Bay of
# Seine, Southern Brittany, Bay of Biscay, Gulf of Lion. See
# func/util.py::ZONE_ORDER/order_zones() for the Python equivalent. Use
# order_zones() to sort any vector/data frame of zone codes -- e.g. replacing
# a plain sort()/unique() that would otherwise come out alphabetical --
# rather than hand-ordering zones at each call site.
ZONE_ORDER <- c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION")

order_zones <- function(zone_vector){
  zone_vector[order(match(zone_vector, ZONE_ORDER))]
}

# In the future this will be taken from define_parameters() in func/util.py
river_mouths <- data.frame(row_name = 1:4,
                           mouth_name = c("Seine", "Gironde", "Loire", "Grand Rhone"),
                           mouth_lon = c(0.145, -1.05, -2.10, 4.83),
                           mouth_lat = c(49.43, 45.59, 47.29, 43.41))

# The zones en large
zones_list <- c("GULF_OF_LION", "BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY")

# Zone bounding boxes
zones_bbox <- data.frame(zone = zones_list,
                         lon_min = c(3.50, -1.50, -4.00, -5.00),
                         lon_max = c(6.00, 0.50, -0.50, -1.50),
                         lat_min  = c(42.25, 49.25, 44.50, 46.5),
                         lat_max = c(44.00, 50.25, 46.50, 48.00)) |>
  mutate(zone_pretty = factor(zone,
                              levels = ZONE_ORDER,
                              labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")), .after = "zone") |>
  dplyr::arrange(match(zone, ZONE_ORDER))

# Create FRANCE bounding box with same structure as zones_bbox
france_bbox <- data.frame(zone = "FRANCE",
                         lon_min = c(-7.8),
                         lon_max = c(10.3),
                         lat_min  = c(41.2),
                         lat_max = c(51.5)) 


# Tide gauge sub-daily QC ---------------------------------------------------
# Flags calendar days whose sub-daily curve does not look tidal, so that
# load_tide_gauge() (below) can null out tide_mean/tide_range for those days
# before they reach func/multi.R's driver comparisons. See func/tide.R for
# the diagnostics (example good/bad days, per-station counts) behind this and
# for the reasoning behind the threshold choices below.
#
# Harmonic model: M2 + S2 + K1 + O1 (dominant semi-diurnal + diurnal tidal
# constituents). Confirmed empirically for all six gauges via the tidal form
# number F = (K1+O1)/(M2+S2) (func/tide.R): all six are semidiurnal or
# mixed-mainly-semidiurnal, so the same 4-constituent model is used
# everywhere -- only the pass/fail thresholds differ by gauge.
tide_periods_hr <- c(M2 = 12.4206012, S2 = 12.0, K1 = 23.93447213, O1 = 25.81933871)

# R^2 floor for the local harmonic fit (see .tide_day_qc()), below which a
# day's tide_range/tide_mean is considered unreliable. Two tiers, not one:
# the Atlantic/Channel gauges (Le Havre, Port-Bloc, Saint-Nazaire) have a
# strong, clean semidiurnal signal (M2 amplitude 1.4-2.5 m) -- even their
# worst 1st-percentile day still fits the harmonic model with R^2 > 0.92, so
# 0.85 only catches genuine failures. The Mediterranean gauges (Fos-sur-Mer,
# Marseille, Port-de-Bouc) have a tiny tidal signal (M2 amplitude ~0.06 m)
# that is easily swamped by non-tidal sea-level noise (storm surge, seiches,
# and at Fos-sur-Mer, harbour-scale wind chop) -- plenty of genuinely good
# days sit at R^2 0.7-0.9, so 0.30 is used instead: below that, the day's
# water level is no longer tide-dominated, whatever the physical cause, so
# it is not a reliable tidal-range value regardless.
tide_r2_threshold <- c("LE_HAVRE" = 0.85, "PORT-BLOC" = 0.85, "SAINT-NAZAIRE" = 0.85,
                       "FOS-SUR-MER" = 0.30, "MARSEILLE" = 0.30, "PORT_DE_BOUC" = 0.30)

# Count local extrema in a sub-daily curve, merging turning points whose
# prominence (jump to the neighbouring turning point) is below `prom` --
# removes sensor-noise-driven wiggles that are not real tidal extrema.
.count_extrema <- function(vals, prom){
  n <- length(vals)
  if(n < 3) return(0)
  signs <- sign(diff(vals))
  signs[signs == 0] <- NA
  signs <- zoo::na.locf(signs, na.rm = FALSE)
  keep <- unique(c(1, which(diff(signs) != 0) + 1, n))
  while(length(keep) > 2){
    vdiffs <- abs(diff(vals[keep]))
    if(all(vdiffs >= prom, na.rm = TRUE)) break
    rm_i <- which.min(vdiffs) + 1
    if(rm_i <= 1 || rm_i >= length(keep)) break
    keep <- keep[-rm_i]
  }
  max(length(keep) - 2, 0)
}

# QC for one calendar day. t_all/y_all are the full (t, tide) vectors for one
# station (hourly, already cleaned -- see .load_tide_raw()); spike_thresh is
# that station's own historical 99.9th-percentile |rate of change| (m/hr).
# Checks, in order (first failure wins):
#   1. sparse/gappy sampling: < 20 of ~24 hourly obs, or a single gap > 4h --
#      the day's max/min cannot be trusted
#   2. spike/step: an hour-to-hour jump beyond this gauge's own historical
#      rate of change (calibrated per-station: what counts as an impossible
#      jump scales directly with each gauge's tidal range)
#   3. extrema mismatch: the number of local extrema in the raw day differs
#      by >= 3 from the number predicted by the local harmonic fit (wrong
#      shape -- e.g. a day that is essentially flat when the tide should
#      have turned twice, or wildly over-oscillating)
#   4. low R^2: the local harmonic fit explains too little of the day's own
#      variance (residuals too large relative to the local tidal signal --
#      see tide_r2_threshold for the per-station floor)
.tide_day_qc <- function(day, t_all, y_all, station, spike_thresh){
  # NB: no explicit tz -- t_all (from .load_tide_raw()) is parsed with
  # as.POSIXct()'s system-default tz too, so day boundaries here must match
  # it or every check below silently slices the wrong hours into "today"
  day_start <- as.POSIXct(paste(day, "00:00:00"))
  day_end   <- day_start + 24*3600

  idx_day <- which(t_all >= day_start & t_all < day_end)
  n_day <- length(idx_day)
  if(n_day < 5) return(data.frame(date = day, tide_bad = TRUE, reason = "sparse"))

  t_day <- t_all[idx_day]; y_day <- y_all[idx_day]
  dt_day <- as.numeric(diff(t_day), units = "hours")
  max_gap_hr <- if(length(dt_day) > 0) max(dt_day) else Inf
  if(n_day < 20 || max_gap_hr > 4) return(data.frame(date = day, tide_bad = TRUE, reason = "sparse"))

  # Spike/step: a single-point reversal (jumps beyond this gauge's own
  # historical rate, then springs most of the way back within the next
  # step) inconsistent with tidal physics. Deliberately NOT just "a fast
  # rate of change" -- a genuine spring-tide flood/ebb can be just as fast,
  # but moves monotonically in one direction; only requiring the rate
  # threshold flagged several textbook-clean spring tides at the
  # high-amplitude Atlantic/Channel gauges (see func/tide.R diagnostics),
  # whereas a real glitch shows up as a jump immediately undone.
  ok <- dt_day >= 0.5 & dt_day <= 1.5
  rate_signed <- ifelse(ok, diff(y_day)/dt_day, NA)
  n_r <- length(rate_signed)
  if(n_r >= 2){
    is_reversal <- abs(rate_signed[-n_r]) > spike_thresh & abs(rate_signed[-1]) > spike_thresh &
      sign(rate_signed[-n_r]) != sign(rate_signed[-1])
    if(any(is_reversal, na.rm = TRUE)) return(data.frame(date = day, tide_bad = TRUE, reason = "spike"))
  }

  # Local harmonic fit on a +-24h buffer window (2-3 full tidal cycles either
  # side for a stable fit), evaluated against just this day's own points
  win_start <- day_start - 24*3600; win_end <- day_end + 24*3600
  idx_win <- which(t_all >= win_start & t_all < win_end)
  if(length(idx_win) < 12) return(data.frame(date = day, tide_bad = TRUE, reason = "sparse"))
  tw <- t_all[idx_win]; yw <- y_all[idx_win]
  hrs <- as.numeric(difftime(tw, day_start, units = "hours"))
  X <- cbind(1, cos(2*pi*hrs/tide_periods_hr[1]), sin(2*pi*hrs/tide_periods_hr[1]),
                cos(2*pi*hrs/tide_periods_hr[2]), sin(2*pi*hrs/tide_periods_hr[2]),
                cos(2*pi*hrs/tide_periods_hr[3]), sin(2*pi*hrs/tide_periods_hr[3]),
                cos(2*pi*hrs/tide_periods_hr[4]), sin(2*pi*hrs/tide_periods_hr[4]))
  coefs <- .lm.fit(X, yw)$coefficients
  yhat_w <- as.vector(X %*% coefs)
  idx_d_in_w <- which(tw >= day_start & tw < day_end)
  y_d <- yw[idx_d_in_w]; yhat_d <- yhat_w[idx_d_in_w]
  ss_tot <- sum((y_d - mean(y_d))^2)
  r2 <- if(ss_tot > 0) 1 - sum((y_d - yhat_d)^2)/ss_tot else 1

  # Expected extrema from the fitted curve at 6-min resolution vs. observed
  # extrema in the raw data (0.03 m prominence filter, ~3x the gauges'
  # 0.01 m reading resolution, to ignore sensor jitter)
  hrs_fine <- seq(0, 24, by = 1/10)
  Xf <- cbind(1, cos(2*pi*hrs_fine/tide_periods_hr[1]), sin(2*pi*hrs_fine/tide_periods_hr[1]),
                 cos(2*pi*hrs_fine/tide_periods_hr[2]), sin(2*pi*hrs_fine/tide_periods_hr[2]),
                 cos(2*pi*hrs_fine/tide_periods_hr[3]), sin(2*pi*hrs_fine/tide_periods_hr[3]),
                 cos(2*pi*hrs_fine/tide_periods_hr[4]), sin(2*pi*hrs_fine/tide_periods_hr[4]))
  extrema_exp <- sum(diff(sign(diff(as.vector(Xf %*% coefs)))) != 0)
  extrema_obs <- .count_extrema(y_day, prom = 0.03)
  if(abs(extrema_obs - extrema_exp) >= 3) return(data.frame(date = day, tide_bad = TRUE, reason = "extrema_mismatch"))

  if(r2 < tide_r2_threshold[[station]]) return(data.frame(date = day, tide_bad = TRUE, reason = "low_r2"))

  return(data.frame(date = day, tide_bad = FALSE, reason = NA_character_))
}

# Per-day tidal QC flags for one station's cleaned sub-daily (t, tide)
# series (see .load_tide_raw()). Returns one row per calendar date with
# tide_bad (logical) and reason (NA when good).
# qc_tide_days(.load_tide_raw("data/TIDES/MARSEILLE"), "MARSEILLE")
qc_tide_days <- function(df_tide, station){
  if(!station %in% names(tide_r2_threshold)){
    stop("Unrecognised tide gauge '", station, "' -- add an R^2 threshold to tide_r2_threshold in func/util.R.")
  }

  # This gauge's own 99.9th-percentile *daily* extreme rate of change --
  # i.e. the 99.9th percentile of each day's own worst hour-to-hour jump
  # (transitions within the same calendar day only, matching the within-day
  # slicing .tide_day_qc() uses for its own spike check).
  daily_max_rate <- df_tide |>
    mutate(date = as.Date(t)) |>
    dplyr::group_by(date) |>
    dplyr::summarise(max_rate = {
      dt <- as.numeric(diff(t), units = "hours")
      ok <- dt >= 0.5 & dt <= 1.5
      rate <- ifelse(ok, abs(diff(tide)/dt), NA)
      if(all(is.na(rate))) NA_real_ else max(rate, na.rm = TRUE)
    }, .groups = "drop")
  spike_thresh <- quantile(daily_max_rate$max_rate, 0.999, na.rm = TRUE)

  days <- unique(as.Date(df_tide$t))
  plyr::ldply(days, .tide_day_qc, t_all = df_tide$t, y_all = df_tide$tide,
              station = station, spike_thresh = spike_thresh, .parallel = TRUE)
}

# Read + clean one tide gauge's raw sub-daily record (all years available,
# source == 4 "hourly validated" readings only). Shared by load_tide_gauge()
# below and the QC diagnostics in func/tide.R.
# NB: Saint-Nazaire also carries source == 6 ("Pleines et basses mers")
# rows -- sparse high/low-water-only readings interleaved a few minutes off
# the hourly grid. These are dropped so the sub-daily curve QC'd above sits
# on a clean, evenly-sampled hourly grid.
.load_tide_raw <- function(dir_name){
  tide_files <- dir(dir_name, pattern = ".txt", full.names = TRUE)
  suppressMessages(
    df_tide <- map_dfr(tide_files, read_delim, col_names = c("t", "tide", "source"), skip = 14, delim = ";", col_select = c("t", "tide", "source"))
  )
  df_tide |>
    dplyr::filter(source == 4) |>
    mutate(t = as.POSIXct(t, format = "%d/%m/%Y %H:%M:%S")) |>
    dplyr::distinct(t, .keep_all = TRUE) |>
    dplyr::arrange(t) |>
    dplyr::select(t, tide)
}


# Pixels ------------------------------------------------------------------

# Simple wrapper to extract start and end times of an ODATIS-MR file
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_CDOM-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20230424__FRANCE_03_MOD_T-FNU-NS_DAY_00.nc"
# file_name <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/02/01/19980201-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
get_start_end_time <- function(file_name){
  
  # Get global info
  df_info <- ncdump::NetCDF(file_name)$attribute[1]$global
  
  # Extract dates accordingly
  if(grepl("SEXTANT", file_name)){
    df_time <- df_info |> 
      mutate(start_time = ymd_hms(paste(gsub("UTC", "", start_date), gsub("UTC", "", start_time)), tz = "GMT"),
             end_time = ymd_hms(paste(gsub("UTC", "", stop_date), gsub("UTC", "", stop_time)), tz = "GMT")) |>
      dplyr::select(start_time, end_time) 
  } else{
    df_time <- df_info |> 
      mutate(start_time = as.POSIXct(gsub("T|Z", " ", start_time), format = "%Y%m%d %H%M%S", tz = "GMT"),
             end_time = as.POSIXct(gsub("T|Z", " ", end_time), format = "%Y%m%d %H%M%S", tz = "GMT")) |> 
      dplyr::select(start_time, end_time)
  }

  # Exit
  return(df_time)
}

# Simple wrapper to extract and save full satellite coord grids
get_sat_grid <- function(file_name){
  nc_data <- nc_open(file_name)
  nc_lon <- as.vector(ncvar_get(nc_data, "lon"))
  nc_lat <- as.vector(ncvar_get(nc_data, "lat"))
  nc_close(nc_data)
  coords_sat <- expand.grid(lon = nc_lon, lat = nc_lat, KEEP.OUT.ATTRS = FALSE)
  return(coords_sat)
}

# Create indexes of which pixels match the 1 km range grid around the in situ sites
# target_site = zone_site_df[10,]; n_pixels = 9; dist_range = 1; sat_rast = rast_sat
# rm(site_sp, site_buffer, cropped_rast, pixel_coords)
get_pixels <- function(target_site, sat_grid, sat_rast, n_pixels, dist_range = 1){
  
  # Filter the sat_grid by the distance from the target site
  # NB: dist_tange/100*2 is a quick and dirty way to reduce the number of pixels to calculate the distance for
  pixel_coords <- sat_grid |> 
    filter(lon >= target_site$lon - dist_range/100*2, lon <= target_site$lon + dist_range/100*2,
           lat >= target_site$lat - dist_range/100*2, lat <= target_site$lat + dist_range/100*2)
                             
  # Calculate distance of pixels from the target in km
  pixel_coords$dist <- round(distHaversine(pixel_coords, target_site[,c("lon", "lat")])/1000, 2)
  
  # Select the nearest 49 (7x7) or 9 (3x3) pixels depending on product
  pixel_coords <- pixel_coords |> arrange(dist) |> slice_head(n = n_pixels)
  
  # Get the pixel IDs and exit
  pixel_coords$cell_numbers <- as.integer(raster::cellFromXY(sat_rast, xy = pixel_coords))

  # Clean up and exit
  pixel_coords <- pixel_coords |> 
    mutate(zone = target_site$zone, source = target_site$source, site = target_site$site) |> 
    dplyr::select(zone, source, site, lon, lat, dist, cell_numbers)
  return(pixel_coords)
}

# Wrapper for get_pixels() to write the pixels to .csv per sat product
# zone_site_df <- filter(zone_sites, zone == "BAY_OF_BISCAY")
# sat_name <- "MODIS"; var_name_full <- "SPM-G-NS_mean"
# sat_name <- "SEXTANT"; var_name_full <- "SPM"
# nc_file_base <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc"
write_pixel_coords <- function(sat_name, var_name_full){
  
  # Load the in situ site metadata
  zone_sites <- read_csv("metadata/in_situ_site_list.csv", show_col_types = FALSE) |> 
    filter(!is.na(zone))

  # Determine file and variable names
  # Logic gate handles the different dates of the fil used for pixels
  if(sat_name == "SEXTANT"){
    if(var_name_full == "SPM"){
      var_name_stub <- "SPIM"
      sat_var <- "analysed_spim"
    } else if(var_name_full == "CHLA"){
      var_name_stub <- "CHL"
      sat_var <- "analysed_chl_a"
    }
    nc_file_base <- paste0("~/pCloudDrive/data/SEXTANT/",var_name_full,
      "/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-",var_name_stub,"-ATL-v01-fv01-OI.nc")
  } else if(sat_name == "MODIS"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/MODIS/L3m_20200101__FRANCE_03_MOD_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "MERIS"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/MERIS/L3m_20090101__FRANCE_03_MER_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "OLCI-A"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/OLCI-A/L3m_20200101__FRANCE_03_OLA_",var_name_full,"_DAY_00.nc")
  } else if(sat_name == "OLCI-B"){
    nc_file_base <- paste0("ODATIS-MR/FRANCE/OLCI-B/L3m_20200101__FRANCE_03_OLB_",var_name_full,"_DAY_00.nc")
  }
  if(sat_name != "SEXTANT"){
      sat_var <- paste0(var_name_full,"_mean")
  }
  
  # Get number of pixels to extract based on product type
  n_pixels <- ifelse(sat_name == "SEXTANT", 9, 49)

  # Get the grid and raster bases
  grid_sat <- get_sat_grid(nc_file_base)
  rast_sat <- raster::raster(nc_file_base, varname = sat_var)
  
  # Extract all pixels
  plan(multisession, workers = parallel::detectCores() - 4)
  zone_pixels <- future_map_dfr(1:nrow(zone_sites),
                                function(i) get_pixels(target_site = zone_sites[i,], 
                                                       sat_grid = grid_sat, sat_rast = rast_sat, 
                                                       n_pixels = n_pixels, dist_range = 1),
                                .options = furrr_options(seed = TRUE))
  plan(sequential)
  
  # Save and exit
  write_csv(zone_pixels, paste0("metadata/zone_pixels_",sat_name,"_",var_name_full,".csv"))
}

# Once the pixels have been determined, use this to extract the data
# file_name <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily/L3m_20020704__FRANCE_03_MOD_CHL-OC5-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/BAY_OF_SEINE/daily/L3m_20020619__FRANCE_03_MER_SPM-G-PO_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/SOUTHERN_BRITTANY/daily/L3m_20120908__FRANCE_03_MOD_CDOM-NS_DAY_00.nc"
# file_name <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/SOUTHERN_BRITTANY/daily/L3m_20040825__FRANCE_03_MOD_NRRS555-NS_DAY_00.nc"
# ncdump::NetCDF(file_name)
# df <- zone_pixels
extract_pixels <- function(file_name, df){
  
  # Determine variable names from file pathway
  if(grepl("SEXTANT", file_name)){
    if(grepl("SPM", file_name)){
      var_base_name <- "SPM"
      var_nc_name <- "analysed_spim"
    } else if(grepl("CHL", file_name)){
      var_base_name <- "CHLA"
      var_nc_name <- "analysed_chl_a"
    } else {
      stop("File structure not recognised")
    }
  } else {
    var_base_name <- str_split(basename(file_name), "_")[[1]][7]
    # var_col_name <- str_split(var_base_name, "-")[[1]][1]
    var_nc_name <- paste0(var_base_name,"_mean")
  }
  
  # Get date values
  if(grepl("SEXTANT", file_name)){
    nc_date <- as.Date(str_split(basename(file_name), "-")[[1]][1], format = "%Y%m%d")
  } else {
    nc_date <- as.Date(str_split(basename(file_name), "_")[[1]][2], format = "%Y%m%d")
  }
  
  # Legacy code kept here for testing purposes
  # Get raster
  # sat_rast <- raster::raster(file_name, varname = var_nc_name)
  
  # Get values
  # sat_vals <- raster::extract(sat_rast, df$cell_numbers)
  
  # Add to data.frame
  # df_res <- df |> 
  #   mutate(date = nc_date,
  #          variable = var_col_name,
  #          value = sat_vals) |> 
  #   dplyr::select(-cell_numbers)
  
  # Extract data and exit
  # NB: A couple of NetCDF files are mysteriously corrupt
  df_res <- tryCatch({
    df |>
      mutate(date = nc_date,
             variable = var_base_name,
             value = raster::extract(raster::raster(file_name, varname = var_nc_name), cell_numbers)) |>
      dplyr::select(-cell_numbers)
  }, error = function(e) {
    message(file_name, " could not be read : ",e$message)
    df |>
      mutate(date = nc_date,
             variable = var_base_name,
             value = NA) |>
      dplyr::select(-cell_numbers)
  })
  return(df_res)
  
  # Load the full nc file to test the raster extraction method
  # df_nc <- tidync::tidync(file_name) |> tidync::hyper_tibble() |>
  #   mutate(lon = plyr::round_any(as.numeric(lon), 0.005),
  #          lat = round(as.numeric(lat), 2)) |>
  #   dplyr::select(lon, lat, analysed_spim)
  
  # Merge and test similarity
  # df_test <- df_res |>
  #   mutate(lon = plyr::round_any(lon, 0.005),
  #          lat = round(lat, 2)) |>
  #   left_join(df_nc, by = c("lon", "lat")) |>
  #   mutate(diff = mean(value-analysed_spim, na.rm = TRUE))
  
  # Test in situ stations with all missing data
  # Eyrac, pk 30, pk 52, Anse du Piquet, Baie d'Yves (a), Cotard, Ile d'Aix, Truscat, Antoine, Luc-sur-Mer
  # is_site <- "Luc-sur-Mer"
  # df_is_test <- df_res |> filter(site == is_site)
  # df_is_test_mean <- summarise(df_is_test, lon = mean(lon), lat = mean(lat), .by = c("zone", "source", "site"))
  # df_nc_test <- df_nc |> filter(lon >= min(df_is_test$lon)-1, lon <= max(df_is_test$lon)+1,
  #                               lat >= min(df_is_test$lat)-1, lat <= max(df_is_test$lat)+1)
  # ggplot() +
  #   annotation_borders() +
  #   geom_raster(data = df_nc_test, aes( x = lon, y = lat, fill = analysed_spim)) +
  #   geom_raster(data = df_is_test, aes( x = lon, y = lat), fill = "red") +
  #   geom_point(data = df_is_test_mean, aes(x = lon, y = lat), colour = "darkred") +
  #   coord_quickmap(xlim = range(df_nc_test$lon), ylim = range(df_nc_test$lat)) +
  #   scale_fill_viridis_c()
  
  # Clean up
  # rm(file_name, df, var_col_name, var_nc_name, df_res, nc_date, df_nc, df_test, df_is_test, df_is_test_mean, df_nc_test, is_site)
}

# Wrapper function to be called per variable
# sat_name <- "SEXTANT"; var_name <- "SPM"
process_pixels <- function(sat_name, var_name){
  
  # Load zone pixels
  file_stub <- paste0(sat_name,"_",var_name)
  
  # Get file pathways
  if(sat_name == "SEXTANT"){
    files_path <- file.path("~/pCloudDrive/data/SEXTANT",var_name)
    var_name_file <- ifelse(var_name == "SPM", "SPIM", ifelse(var_name == "CHLA", "CHL"))
  } else {
    files_path <- file.path("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR", sat_name, zone_name)
    var_name_file <- var_name
  }
  files_var <- dir(files_path, recursive = TRUE, full.names = TRUE, pattern = var_name_file)

  # Filter out .png files from SEXTANT folders
  if(sat_name == "SEXTANT"){
    files_var <- files_var[!grepl(".png", files_var)]
  }

  # Extract data for all files and save
  if(length(files_var) > 0){

    file_name_var <- paste0("output/MATCH_UP_DATA/FRANCE/zone_data_",file_stub,".csv")

    if(!file.exists(file_name_var)){

      # Extract data from all pixels
      message("Started ",var_name," extraction at : ", Sys.time())
      zone_pixels <- read_csv(paste0("metadata/zone_pixels_",file_stub,".csv"), show_col_types = FALSE)

      plan(multisession, workers = parallel::detectCores() - 2)
      zone_data_var <- future_map_dfr(files_var, extract_pixels, df = zone_pixels,
                                      .options = furrr_options(seed = TRUE))
      
      # Save results
      message("Saving ",var_name," extraction at : ", Sys.time())
      data.table::fwrite(zone_data_var, file_name_var)
      plan(sequential)
      gc()
      gc()
    }

    # Create median value time series
    file_name_median_all <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,"_all.csv")
    file_name_median_small <- paste0("output/MATCH_UP_DATA/FRANCE/zone_median_",file_stub,"_small.csv")
    
    # Create medians etc. from all pixels
    if(!file.exists(file_name_median_all)){
      # Load data from .csv
      if(!exists("zone_data_var")){
        message("Loading ",var_name," extraction at : ", Sys.time())
        zone_data_var <- data.table::fread(file_name_var)
      }
      message("Started median all calculations at : ", Sys.time())
      zone_median_all <- zone_data_var |> 
        filter(value > 0) |>
        summarise(median = median(value, na.rm = TRUE), 
                  mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE),
                  n = n(), 
                  .by = c("zone", "source", "site", "date", "variable"))
      data.table::fwrite(zone_median_all, file_name_median_all)
      rm(zone_median_all); gc()
    }

    # Create medians etc. from 'small' pixels
    if(!file.exists(file_name_median_small)){
      # Load data from .csv
      if(!exists("zone_data_var")){
        message("Loading ",var_name," extraction at : ", Sys.time())
        zone_data_var <- data.table::fread(file_name_var)
      }
      message("Started median small calculations at : ", Sys.time())
      if(sat_name == "SEXTANT"){
        slice_n <- 1
      } else{
        slice_n <- 9
      }
      zone_median_small <- zone_data_var |> 
        filter(value > 0) |> # NB: This removes NA pixels before selecting the nearest pixel, which may be incorrect
        group_by(zone, source, site, date, variable) |> 
        arrange(dist) |> 
        slice_head(n = slice_n) |> 
        ungroup() |> 
        summarise(median = median(value, na.rm = TRUE), 
                  mean = mean(value, na.rm = TRUE),
                  sd = sd(value, na.rm = TRUE),
                  n = n(), 
                  .by = c("zone", "source", "site", "date", "variable"))
      data.table::fwrite(zone_median_small, file_name_median_small)
      rm(zone_median_small); gc()
    }

    # Cleanup and exit
    rm(zone_data_var, file_name_var, file_name_median_all, file_name_median_small, file_stub, zone_pixels)
    gc()
    gc()
  }
}

# Function that extracts all of the data identified in the product/zone lon/lat metadata files
# sat_name <- "SEXTANT"
# sat_name <- "MODIS"; zone_name <- "SOUTHERN_BRITTANY"
extract_pixels_all <- function(sat_name, zone_name = NULL){#, overwrite = FALSE){
  
  message("Started run on ", sat_name, " : ", Sys.time())
  
  # Get file pathways
  if(sat_name == "SEXTANT"){
    process_pixels(sat_name, "SPM")
    process_pixels(sat_name, "CHLA")
  } else {

    ## CHL1
    process_pixels(sat_name, "CHL1-AC"); gc(); gc()
    process_pixels(sat_name, "CHL1-NS"); gc(); gc()
    process_pixels(sat_name, "CHL1-PO"); gc(); gc()

    ## CHL-GONS - No NS products
    process_pixels(sat_name, "CHL-GONS-AC"); gc(); gc()
    process_pixels(sat_name, "CHL-GONS-PO"); gc(); gc()
    
    ## CHL-OC5
    process_pixels(sat_name, "CHL-OC5-AC"); gc(); gc()
    process_pixels(sat_name, "CHL-OC5-PO"); gc(); gc()
    process_pixels(sat_name, "CHL-OC5-NS"); gc(); gc()

    ## SPM-G
    process_pixels(sat_name, "SPM-G-AC"); gc(); gc()
    process_pixels(sat_name, "SPM-G-PO"); gc(); gc()
    process_pixels(sat_name, "SPM-G-NS"); gc(); gc()
    
    ## SPM-R
    process_pixels(sat_name, "SPM-R-AC"); gc(); gc()
    process_pixels(sat_name, "SPM-R-PO"); gc(); gc()
    process_pixels(sat_name, "SPM-R-NS"); gc(); gc()

    ## TUR
    process_pixels(sat_name, "T-FNU-AC"); gc(); gc()
    process_pixels(sat_name, "T-FNU-PO"); gc(); gc()
    process_pixels(sat_name, "T-FNU-NS"); gc(); gc()
    
    ## SST - Only NS and PO
    process_pixels(sat_name, "SST-NS"); gc(); gc()
    process_pixels(sat_name, "SST-PO"); gc(); gc()

    ## SST-NIGHT - Only for NS
    process_pixels(sat_name, "SST-NIGHT-NS"); gc(); gc()
  }
  
  # Clean and exit
  message("Finished ", sat_name," at : ", Sys.time())
}


# Loading -----------------------------------------------------------------

# Zone-specific hard ceiling on classified plume area (km^2), applied
# unconditionally in load_plume_ts() below regardless of which metric_col is
# requested. Replaces a single flat 20000 km^2 guard with one ceiling per
# zone, chosen by hand (2026-07) from func/surface.R's suspicious-plume-day
# scan against the daily plume-shape grids -- each zone's plausible maximum
# extent differs by an order of magnitude (e.g. the Gulf of Lion's shelf vs.
# the Seine's), so a single global cutoff was either too loose for the Seine
# or too tight for the Gulf of Lion. BAY_OF_SEINE's ceiling in particular is
# deliberately tight (2500 km^2): this zone is the most turbid, so its
# dynamic threshold is the most prone to spilling into open-ocean pixels,
# and it's better to lose several real high-area days than keep the
# spillover ones.
plume_area_ceiling <- c(BAY_OF_BISCAY = 12000, BAY_OF_SEINE = 2500,
                        GULF_OF_LION = 10000, SOUTHERN_BRITTANY = 6000)

# Load time series of plume values. `metric_col` selects which Results.csv
# column becomes `plume_area` downstream (every driver_plume_* function in
# multi.R is written against that column name regardless of what it holds,
# so this is the one place a different plume metric -- e.g. mass_SPM_in_the_
# plume_area_in_g_m -- needs to be wired in). `outlier_max` is only a
# sensible guard for the area column (20000 km^2 is physically implausible
# for these zones); pass NULL to skip it for other metrics. The
# plume_area_ceiling cap above is applied to area_of_the_plume_mask_in_km2
# first and unconditionally, before metric_col is even selected -- so it
# takes effect no matter which metric is ultimately requested, while every
# other Results.csv column (mean_SPM, mass, centroids, confidence, n_pixel)
# is left exactly as read, even on a day whose area gets nulled here.
load_plume_ts <- function(zone, plume_dir = "output/panache/dynamic",
                          metric_col = "area_of_the_plume_mask_in_km2", outlier_max = 20000){
  file_name <- paste0(plume_dir, "/", zone, "/Results.csv")
  suppressMessages({
    df_plume <- read_csv(file_name) |>
      dplyr::mutate(date = as.Date(date)) |>
      dplyr::select(date:confidence_index_in_perc) |>
      complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |>
      dplyr::mutate(area_of_the_plume_mask_in_km2 = ifelse(area_of_the_plume_mask_in_km2 > plume_area_ceiling[[zone]],
                                                            NA, area_of_the_plume_mask_in_km2)) |>
      dplyr::rename(plume_area = !!rlang::sym(metric_col))
    if(!is.null(outlier_max)) df_plume <- dplyr::mutate(df_plume, plume_area = ifelse(plume_area > outlier_max, NA, plume_area))
    df_plume <- df_plume |>
      zoo::na.trim() |>
      mutate(zone = zone, .before = "date")
  })
  return(df_plume)
}

# Get the (lon, lat) footprint of the daily plume mask for a zone. Panache's
# 2026-07-17 rewrite replaced the old per-day CSV output (one file per date,
# lon/lat columns) with a single PlumeMasks.nc cube per zone/threshold
# (plume_mask(time, lat, lon), 0/1) -- read that instead, keeping only pixels
# flagged as plume. plume_dir mirrors load_plume_ts()'s static/dynamic switch.
load_plume_surface <- function(zone, plume_dir = "output/panache/dynamic"){
  file_name <- paste0(plume_dir, "/", zone, "/PlumeMasks.nc")
  nc_dat <- nc_open(file_name)
  lon  <- ncvar_get(nc_dat, "lon")
  lat  <- ncvar_get(nc_dat, "lat")
  time <- .nc_time_to_date(nc_dat, "time")
  mask <- ncvar_get(nc_dat, "plume_mask")  # [lon, lat, time]
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
# NetCDF source read below. Replaces tidync's automatic CF time decoding
# (see .nc_read_box() below for why tidync was dropped).
.nc_time_to_date <- function(nc, time_var = "time"){
  units_str <- ncatt_get(nc, time_var, "units")$value
  m <- regmatches(units_str, regexec("^(seconds|hours|days) since (.+)$", units_str))[[1]]
  if(length(m) != 3) stop("Unrecognised time units string: ", units_str)
  multiplier <- switch(m[2], seconds = 1, hours = 3600, days = 86400)
  raw <- ncvar_get(nc, time_var)
  as.Date(as.POSIXct(raw * multiplier, origin = m[3], tz = "UTC"))
}

# Read one or more variables from a NetCDF file, cropped to a lon/lat box, as
# a long (lon, lat, date, <var>...) tibble -- replaces
# tidync::tidync() |> hyper_filter() |> hyper_tibble(). Switched from tidync
# to ncdf4 throughout this repo (2026-07-31) because tidync's RNetCDF
# dependency fails to load under rpy2's embedded R in this environment
# ("RNetCDF.so: undefined symbol: nc_reclaim_data"), even though it loads
# fine under a plain Rscript session; ncdf4 has no such conflict.
# ncdf4::ncvar_get() returns arrays ordered (lon, lat, [depth,] time) for
# every product read here. start/count reads only the lon/lat box directly
# off disk -- a contiguous hyperslab, since lon_idx/lat_idx come from a
# monotonic coordinate range -- instead of pulling the full spatial/time
# extent into memory first and subsetting in R; tidync's hyper_filter() used
# to push this same filter down to the read. Any dimension between lat and
# time (e.g. GLORYS current files' single-level depth) is read as a bare
# size-1 slice, which ncvar_get()'s default collapse_degen = TRUE then drops
# automatically, leaving the same (lon, lat, time) shape regardless of
# whether the file has that extra dimension.
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

# Load wave data (significant wave height + mean direction).
# NB: unlike wind/current, VMDR (mean wave direction) is a raw angle, not
# derived from vector components -- the pixel-box mean below is a plain
# circular (unit-vector) mean over the lon/lat box, which avoids the 0/360
# wrap-around problem a naive arithmetic mean would have *within this
# function*. The same problem still exists one step upstream: the "_daily_"
# source files are produced by cdo daymean directly on hourly VMDR, which
# is an arithmetic (not circular) mean, so any day where wave direction
# crosses the 0/360 wrap can already carry a biased wave_dir before it
# reaches this function.
# NB 2 (found 2026-07-31, pre-dates the tidync->ncdf4 migration; RESOLVED
# 2026-07-31 by the WAVE data re-download): the Bay of Seine's wave file
# (IBI_daily_199801_202601.nc) used to have no VMDR variable at all (VHM0 +
# VTM02 only, unlike the other three zones' VHM0 + VMDR), an upstream
# download gap. The re-download confirmed in that day's driver-interactions
# rerun added VMDR for the Bay of Seine too, so wave_dir is no longer
# expected to be all-NA for that zone -- the has_wave_dir fallback below
# (and plot_driver_rose()'s all-NA handling) is kept as defensive code in
# case a future re-download regresses, not because the gap is current.
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

# Load GLORYS surface current (eastward/northward velocity) data for a zone.
# NB: current data through 2024-12-31 lives in the combined historical file
# ("glorys_199301_202412.nc"); 2025 onward is a separate uo/vo-only file
# ("glorys_uo_vo_202501_202512.nc") -- both live under
# ~/pCloudDrive/data/GLORYS/<zone_name>/ and are combined here.
load_surface_current <- function(zone_name, lon_range, lat_range){
  dir_name <- path.expand(paste0("~/pCloudDrive/data/GLORYS/", zone_name))
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

# Load ROFI surface NetcDF
load_ROFI <- function(file_name){
  # Get zone from file name
  if(grepl("Seine", file_name)){
    zone <- "BAY_OF_SEINE"
  } else if (grepl("Gironde", file_name)){
    zone <- "BAY_OF_BISCAY"
  } else if (grepl("Loire", file_name)){
    zone <- "SOUTHERN_BRITTANY"
  } else if (grepl("Rhone", file_name)){
    zone <- "GULF_OF_LION"
  } else {
    stop("Zone not recognised from ROFI file name")
  }
  nc <- nc_open(file_name)
  date <- .nc_time_to_date(nc, "time_counter")
  rofi_surface <- ncvar_get(nc, "ROFI_surface")
  nc_close(nc)

  df_ROFI <- tibble::tibble(zone = zone, date = date, ROFI_surface = rofi_surface) |>
    dplyr::summarise(ROFI_surface = mean(ROFI_surface, na.rm = TRUE), .by = c("zone", "date"))
  return(df_ROFI)
}

# Load X11 results
load_X11 <- function(zone_name, plume = TRUE){
  if(plume){
    file_name <- paste0("output/FIXED_THRESHOLD/GULF_OF_LION/X11_ANALYSIS/area_of_the_plume_mask_in_km2/SEXTANT_merged_Standard_WEEKLY.csv")
  } else {
    file_name <- paste0("output/FIXED_THRESHOLD/GULF_OF_LION/X11_ANALYSIS/river_flow/River_flow___WEEKLY.csv")
  }
  suppressMessages(
    df_X11 <- read_csv(file_name) |> 
      dplyr::rename(date = dates) |> 
      complete(date = seq(min(date), max(date), by = "day")) |> 
      mutate(plume_area = ifelse(plume_area > 20000, NA, plume_area)) |> 
      zoo::na.trim() |> 
      mutate(zone = zone, .before = "date")
  )
  return(df_X11)
}


# Statistics --------------------------------------------------------------

# Check for leap year
is_leap_year <- function(year) {
  (year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)
}

# Function to adjust DOY for non-leap years
adjust_doy <- function(year, doy) {
  if (!is_leap_year(year) && doy > 59) {
    doy + 1
  } else {
    doy
  }
}

# Lagged correlations
# NB: The x vector is lagged, meaning it is effectively pushed forward in time
lagged_correlation <- function(x, y, max_lag) {
  df_cor <- tibble(
    lag = 0:max_lag,
    cor = map_dbl(0:max_lag, ~ cor(y, lag(x, .), use = "pairwise.complete.obs"))
  )
  return(df_cor)
}

# Calculate STL
stl_single <- function(x_col, out_col, start_date, ts_freq = 365.25){
  
  # Create ts object and calculate stl
  ts_x <- ts(zoo::na.approx(x_col), frequency = ts_freq, start = c(year(start_date), quarter(start_date)))
  stl_x <- stl(ts_x, s.window = "periodic", t.window = 11, inner = 5, outer = 0)
  #, robust = TRUE) # rather allow outliers to influence results
  
  # Add NA to end if necessary
  if(length(stl_x$time.series[,2]) != length(x_col)){

    # Create a matrix of NA values
    na_matrix <- matrix(NA, nrow = length(x_col)-length(stl_x$time.series[,2]), ncol = ncol(stl_x$time.series))
    
    # Append the NA matrix to the original matrix
    stl_x$time.series <- rbind(stl_x$time.series, na_matrix)
  }
  
  # Decide which column to take
  if(out_col == "seas"){
    return(as.vector(stl_x$time.series[,1]))
  } else if(out_col == "inter"){
    return(as.vector(stl_x$time.series[,2]))
  } else if(out_col == "remain"){
    return(as.vector(stl_x$time.series[,3]))
  } else {
    stop("'out_col' not recognised")
  }
}

# Load panache time series based on river mouth name
plume_clim_calc <- function(zone){
  file_name <- paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")
  suppressMessages(
    df_plume <- read_csv(file_name) |> 
      dplyr::select(date:confidence_index_in_perc) |>
      complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
      dplyr::rename(plume_area = area_of_the_plume_mask_in_km2) |> 
      mutate(plume_area = ifelse(plume_area > 20000, NA, plume_area),) |> 
      zoo::na.trim() |> 
      mutate(zone = zone, .before = "date")
  )
  df_clim <- heatwaveR::ts2clm(df_plume, x = date, y = plume_area, climatologyPeriod = c(min(df_plume$date), max(df_plume$date))) |> 
    dplyr::select(zone, date, doy, plume_area, seas, thresh) |> 
    dplyr::rename(plume_seas = seas, plume_thresh = thresh)
}

# The full suite of stats to calculate
compute_stats <- function(x_vec, y_vec){
  
  if(!is.numeric(x_vec)) stop("x_vec is not numeric")
  if(!is.numeric(y_vec)) stop("y_vec is not numeric")
  
  if(length(x_vec) < 3){
    return(data.frame(row.names = NULL,
                      n = length(x_vec),
                      Slope = NA,
                      Slope_log = NA,
                      RMSE = NA,
                      MSA = NA,
                      MAPE = NA,
                      Bias = NA,
                      Error = NA))
  }
  
  # Calculate RMSE (Root Mean Square Error)
  rmse <- sqrt(mean((y_vec - x_vec)^2, na.rm = TRUE))
  
  # Calculate MAPE (Mean Absolute Percentage Error)
  mape <- mean(abs((y_vec - x_vec) / x_vec), na.rm = TRUE) * 100
  
  # Calculate MSA (Mean Squared Adjustment)
  msa <- mean(abs(y_vec - x_vec), na.rm = TRUE)
  
  # Calculate linear slope
  lin_fit <- lm(y_vec ~ x_vec)
  slope <- coef(lin_fit)[2]
  
  # Calculate log-log linear slope
  log_lin_fit <- lm(log10(y_vec) ~ log10(x_vec))
  log_slope <- coef(log_lin_fit)[2]
  
  # Calculate Bias
  log_ratio <- log10(y_vec / x_vec)
  log_ratio_median <- median(log_ratio, na.rm = TRUE)
  bias_perc <- 100 * (sign(log_ratio_median) * (10^abs(log_ratio_median) - 1))
  
  # Calculate error
  log_ratio_median_abs <- median(abs(log_ratio), na.rm = TRUE)
  error_perc <- 100 * (10^log_ratio_median_abs - 1)
  
  # Combine int data.frame and exit
  return(data.frame(row.names = NULL,
                    n = length(x_vec),
                    Slope = round(slope, 2),
                    Slope_log = round(log_slope, 2),
                    RMSE = round(rmse, 6),
                    MSA = round(msa, 6),
                    MAPE = round(mape, 2),
                    Bias = round(bias_perc, 2),
                    Error = round(error_perc, 2)))
}


# Plotting ----------------------------------------------------------------

# Create pretty plot labels from zone values
make_pretty_title <- function(df){
  df <- df |> 
    mutate(plot_title = case_when(zone == "BAY_OF_SEINE" ~ "Bay of Seine",
                                  zone == "SOUTHERN_BRITTANY" ~ "Southern Brittany",
                                  zone == "BAY_OF_BISCAY" ~ "Bay of Biscay",
                                  zone == "GULF_OF_LION" ~ "Gulf of Lion"), .after = "zone") |> 
    mutate(plot_title = factor(plot_title, levels = c("Bay of Seine", "Southern Brittany", "Bay of Biscay", "Gulf of Lion")))
  return(df)
}

# Scale one value to another for tidier double-y-axis plots
sec_axis_adjustement_factors <- function(var_to_scale, var_ref){
  
  index_to_keep <- which(is.finite(var_ref))
  var_ref <- var_ref[index_to_keep]
  
  index_to_keep <- which(is.finite(var_to_scale))
  var_to_scale <- var_to_scale[index_to_keep]
  
  max_var_to_scale <- max(var_to_scale, na.rm = T) 
  min_var_to_scale <- min(var_to_scale, na.rm = T) 
  max_var_ref <- max(var_ref, na.rm = T) 
  min_var_ref <- min(var_ref, na.rm = T) 
  
  diff_to_scale <- max_var_to_scale - min_var_to_scale
  diff_to_scale <- ifelse(diff_to_scale == 0, 1 , diff_to_scale)
  diff_ref <- max_var_ref - min_var_ref
  diff <- diff_ref / diff_to_scale
  
  adjust <- (max_var_ref - max_var_to_scale*diff) 
  
  return(data.frame(diff = diff, adjust = adjust, operation = "scaled var = (var_to_scale * diff) + adjust",
                    trans_axis_operation = "var_to_scale = {scaled_var - adjust} / diff)"))
}

# Consistent theme for project. Gridlines-off variant, matching every
# current manuscript figure (this copy used to also show gridlines plus a
# plot.subtitle style that nothing else used; func/figure.R's and
# func/X11.R's copies already agreed on gridlines-off, so this copy was
# brought in line with those rather than the reverse -- see the refactor
# that deduplicated the three near-identical copies of this function).
ggplot_theme <-   function(){
  theme(text = element_text(size = 35, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 55),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(angle = 0),
        axis.ticks.length = unit(.25, "cm"))
}

# Convenience wrapper for saving png files
save_plot_as_png <- function(plot, name = c(), width = 14, height = 8.27, path, res = 150){
  
  graphics.off()
  
  if (name %>% length() == 1) {
    if (dir.exists(path) == FALSE) {dir.create(path, recursive = TRUE)}
    path <- file.path(path, paste(name, ".png", sep = ""))
  } else {
    path <- paste(path, ".png", sep = "")
  }
  
  if (grepl(pattern = ".png.png", path)) {path <- path %>% gsub(pattern = ".png.png", replacement = ".png", x = .)}
  
  png(path, width = width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
  
}

# Comparison plots
# df <- df_pretty; var_1 <- "plume_inter"; var_2 <- "flow_inter"
# df <- df_seas; var_1 <- "plume_seas"; var_2 <- "flow_seas"
# colour_1 <- "brown"; colour_2 <- "blue"; label_1 <- "Plume area (km^2)"; label_2 <- "River flow (m^3 s-1)"; file_stub <- "comparison_plume_flow_inter"
comparison_plot <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2){
  
  if(grepl("seas", var_1)){
    df_sub <- df[,c("plot_title", "month", var_1, var_2)]
    colnames(df_sub) <- c("plot_title", "month", "var_1", "var_2")
  } else {
    df_sub <- df[,c("plot_title", "date", var_1, var_2)]
    colnames(df_sub) <- c("plot_title", "date", "var_1", "var_2")
  }
  
  # Plot base
  if(grepl("seas", var_1)){
    
    # Scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_sub$var_2, df_sub$var_1)
    df_scaling <- summarise(df_sub, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_sub, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    
    # Get range for ribbon plot
    df_scale_sub <- df_scale |> 
      summarise(var_1_min = min(var_1, na.rm = TRUE),
                var_1_mean = mean(var_1, na.rm = TRUE),
                var_1_max = max(var_1, na.rm = TRUE),
                var_2_min = min(var_2_scaled, na.rm = TRUE),
                var_2_mean = mean(var_2_scaled, na.rm = TRUE),
                var_2_max = max(var_2_scaled, na.rm = TRUE), .by = c("plot_title", "month")) |> 
      mutate(month_int = as.integer(month))
    
    # Plot them
    pl_base <- ggplot(data = df_scale_sub, aes(x = month_int)) + 
      # Var 1
      geom_ribbon(aes(ymin = var_1_min, ymax = var_1_max), fill = colour_1, alpha = 0.3) +
      geom_path(aes(y = var_1_mean), color = colour_1, linewidth = 2) +
      # Var 2
      geom_ribbon(aes(ymin = var_2_min, ymax = var_2_max), fill = colour_2, alpha = 0.3) +
      geom_path(aes(y = var_2_mean), color = colour_2, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_continuous(expand = c(0, 0), breaks = 1:12, labels = month.abb)
  } else {
    
    # Perform rolling mean
    df_roll_mean <- df_sub |>
      mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
      # mutate(date = round_date(date, unit = "months")) |>
      filter(date >= min(df$date)) |>
      group_by(plot_title, date) |>
      summarise(var_1 = mean(var_1, na.rm = TRUE),
                var_2 = mean(var_2, na.rm = TRUE), .groups = "keep") |>
      group_by(plot_title) |>
      mutate(var_1 = roll_mean(var_1, n = 48, fill = NA, align = "center"),
             var_2 = roll_mean(var_2, n = 48, fill = NA, align = "center")) |>
      ungroup()
    
    # Then get the scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_roll_mean$var_2, df_roll_mean$var_1)
    df_scaling <- summarise(df_roll_mean, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_roll_mean, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    unique_years <- df_scale$date |> year() |> unique()
    
    pl_base <- ggplot(data = df_scale) +
      # Var 1 data
      geom_point(aes(x = date, y = var_1), color = colour_1) +
      geom_path(aes(x = date, y = var_1), color = colour_1) +
      # geom_smooth(aes(x = date, y = var_1), method = "lm", se = FALSE, color = colour_1) +
      # Var 2 data
      geom_point(aes(x = date, y = var_2_scaled), color = colour_2) +
      geom_path(aes(x = date, y = var_2_scaled), color = colour_2) +
      # geom_smooth(aes(x = date, y = var_2_scaled), method = "lm", se = FALSE, color = colour_2) +
      # Facet
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      # X-axis labels
      scale_x_date(name = "", expand = c(0, 0),
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) 
  }
  
  # Finish up the comparison plot
  pl_comp <- pl_base +
    # Y-axis labels
    scale_y_continuous(name = label_1,
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
                                           name = label_2)) +
    labs( x = NULL) +
    # Extra bits
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y.left = element_text(color = colour_1),
          axis.ticks.y.left = element_line(color = colour_1),
          axis.line.y.left = element_line(color = colour_1),
          axis.title.y.left = element_text(color = colour_1, margin = unit(c(0, 7.5, 0, 0), "mm")),
          axis.text.y.right = element_text(color = colour_2),
          axis.ticks.y.right = element_line(color = colour_2),
          axis.line.y.right = element_line(color = colour_2),
          axis.title.y.right = element_text(color = colour_2, margin = unit(c(0, 0, 0, 7.5), "mm")),
          panel.border = element_rect(linetype = "solid", fill = NA))
  return(pl_comp)
}

# Convenience wrapper to run and save comparison plots
# NB: this is hard coded to work with four plots
comparison_plot_save <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2, file_stub){
  comp_list <- plyr::dlply(df, c("zone"), comparison_plot, var_1 = var_1, var_2 = var_2, 
                           colour_1 = colour_1, colour_2 = colour_2, label_1 = label_1, label_2 = label_2)
  comp_fig <- comp_list[[2]] + comp_list[[4]] + comp_list[[1]] + comp_list[[3]] + patchwork::plot_layout(ncol = 1, axes = "collect")
  ggsave(filename = paste0("figures/",file_stub,".png"), plot = comp_fig, width = 24, height = 24, dpi = 300)
}

# It does what it says on the tin
var_labels <- function(var_name){
  
  if(grepl("CHL|RRS", var_name)){
    var_colour <- "green4"
    if(grepl("CHL", var_name)){
      axis_limits <- c(10^-2, 10^2)
      unit <- expression(mg~m^-3)
    } else {
      unit <- expression(sr-1)
      axis_limits <- c(10^-5, 10^-1)
    }
  } else if(grepl("TEMP|SST", var_name)) {
    unit <- expression('"°C"')
    axis_limits <- c(3, 30)
    var_colour <- "orange4"
    # NB: Must be run after SST because the T for turbidity is an issue
  } else if(grepl("SPM|TUR|T|POC|CDOM", var_name)){
    axis_limits <- c(10^-2, 10^3)
    var_colour<- "brown4"
    if(grepl("TUR|T", var_name)){
      unit <- expression(NTU)
    } else if(var_name == "CDOM"){
      unit <- expression(m-1)
    } else if(var_name == "POC"){
      unit <- expression(mg~m^-3)
    } else {
      unit <- expression(g~m^-3)
    }
  } else {
    stop(paste("Could not find ", var_name))
  }
  
  # List and exit
  to_return <- list("unit" = unit,
                    "var_colour" = var_colour, 
                    "axis_limits" = axis_limits)
  return(to_return)
}

# Convenience colour wrapper
colours_of_stations <- function(){
  
  # colour_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
  #                  "Point C" = mako(n = 1,begin = 0.85,end = 0.85), 
  #                  "Luc-sur-Mer" = mako(n = 1,begin = 0.90,end = 0.9), 
  #                  "Smile" = mako(n = 1,begin = 0.95,end = 0.95),
  #                  
  #                  "Bizeux" = viridis(n = 1,begin = 1,end = 1), # Bretagne
  #                  "Le Buron" = viridis(n = 1,begin = 0.925,end = 0.925), 
  #                  "Cézembre" = viridis(n = 1,begin = 0.85,end = 0.85), 
  #                  "Estacade" = viridis(n = 1,begin = 0.775,end = 0.775), 
  #                  "Astan" = viridis(n = 1,begin = 0.7,end = 0.7), 
  #                  "Portzic" = viridis(n = 1,begin = 0.625,end = 0.625), 
  #                  
  #                  "Antioche" = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
  #                  "pk 86" = plasma(n = 1,begin = 0.1,end = 0.1), 
  #                  "pk 52" = plasma(n = 1,begin = 0.15,end = 0.15),
  #                  "pk 30" = plasma(n = 1,begin = 0.2,end = 0.2),
  #                  "Comprian" = plasma(n = 1,begin = 0.25,end = 0.25), 
  #                  "Eyrac" = plasma(n = 1,begin = 0.3,end = 0.3), 
  #                  "Bouee 13" = plasma(n = 1,begin = 0.35,end = 0.35), 
  #                  
  #                  "Sola" = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
  #                  "Sete" = rocket(n = 1,begin = 0.85,end = 0.85), 
  #                  "Frioul" = rocket(n = 1,begin = 0.90,end = 0.90),
  #                  "Point B" = rocket(n = 1,begin = 0.95,end = 0.95)
  # ) 
  
  # colour_values = c('Manche orientale - Mer du Nord' = mako(n = 1,begin = 0.8,end = 0.8), # Manche
  #                  'Baie de Seine' = mako(n = 1,begin = 0.95,end = 0.95),
  # 
  #                  'Manche occidentale' = viridis(n = 1,begin = 1,end = 1), # Bretagne
  #                  'Bretagne Sud' = viridis(n = 1,begin = 0.625,end = 0.625),
  # 
  #                  'Pays de la Loire - Pertuis' = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
  #                  'Sud Golfe de Gascogne' = plasma(n = 1,begin = 0.35,end = 0.35),
  # 
  #                  'Golfe du Lion' = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
  #                  'Mer ligurienne - Corse' = rocket(n = 1,begin = 0.95,end = 0.95)
  # )
  
  # colour_values = c('BAY OF SEINE' = viridis::mako(n = 1, begin = 0.8,end = 0.8),
  #                  'SOUTHERN BRITTANY' = viridis::viridis(n = 1, begin = 0.9,end = 0.9),
  #                  'GULF OF BISCAY' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
  #                  'GULF OF LION' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  colour_values = c('Bay of Seine' = viridis::mako(n = 1, begin = 0.8,end = 0.8),
                    'S. Brittany' = viridis::viridis(n = 1, begin = 0.9,end = 0.9),
                    'Bay of Biscay' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
                    'Gulf of Lion' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  
  return(colour_values)
}

# Plot linear trends and stats for matched data
# var_combi = zone_all_in_situ_monthly$variable_combi[1]; df = zone_all_in_situ_monthly; df_stats = zone_all_monthly_lm
validation_lm_plots <- function(var_combi, sat_name, median_base, df, df_stats){
  
  # Filter and pivot data
  df_var_sub <- df |> 
    filter(variable_combi == var_combi) |> 
    pivot_longer(cols = value_in_situ:value_satellite) |> 
    mutate(name = gsub("value|_", "", name),
           name = gsub("insitu", "in situ", name)) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  
    # Get y-axis label and units
  if(df_var_sub$variable[1] == "TEMP"){
    y_lab <- "Temperature (°C)"
    y_lim <- c(4, 32)
    var_sat <- "SST"
  } else if(df_var_sub$variable[1] == "CHLA"){
    y_lab <- "Chlorophyll-a (mg m-3)"
    y_lim <- c(0, 20)
    var_sat <- "CHL"
  } else if(df_var_sub$variable[1] == "TUR"){
    y_lab <- "Turbidity (NTU)"
    y_lim <- c(0, 20)
  } else if(df_var_sub$variable[1] == "SPM"){
    y_lab <- "SPM (g m-3)"
    y_lim <- c(0, 20)
  } else {
    # Not worried about other variables for the moment
    y_lab <- NA
    y_lim <- c(0, 20)
  }
  
  # Filter stats labels
  df_stats_var_sub <- df_stats |> 
    filter(variable_combi == var_combi) |> 
    mutate(y = max(y_lim)-2) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  
  # Create title and subtitle
  var_name_sat <- df_var_sub$variable_sat[1]
  var_name_is <- df_var_sub$variable[1]
  cor_name <- df_var_sub$correction[1]
  proc_name <- df_var_sub$processing[1]
  grid_name <- df_var_sub$grid_size[1]
  plot_title = paste(sat_name, var_name_sat, "vs.", "in situ", var_name_is)
  plot_subtitle = paste("Correction:", cor_name, "| Processing:", proc_name, "| Grid size:", grid_name)
  
  # Create plot
  plot_zone_var_TS <- ggplot(data = df_var_sub, aes(x = date, y = value)) +
    geom_point(aes(y = value, colour = name, shape = source), alpha = 0.4) +
    geom_line(aes(y = value, colour = name, linetype = source), alpha = 0.4) +
    geom_smooth(aes(colour = name, linetype = source), method = "lm", se = FALSE) +
    # In situ values
    geom_label(data = filter(df_stats_var_sub, source == "REPHY"), 
               aes(x = as.Date("2005-01-01"), y = y, label = paste0(source," : ", round(slope_is, 2)," / year"), vjust = 1.0)) +
    geom_label(data = filter(df_stats_var_sub, source == "SOMLIT"), 
               aes(x = as.Date("2005-01-01"), y = y, label = paste0(source," : ", round(slope_is, 2)," / year"), vjust = -0.5)) +
    # Satellite values
    geom_label(data = filter(df_stats_var_sub, source == "REPHY"), colour = "red",
               aes(x = as.Date("2020-01-01"), y = y, label = paste0(source," : ", round(slope_sat, 2)," / year"), vjust = 1.0)) +
    geom_label(data = filter(df_stats_var_sub, source == "SOMLIT"),  colour = "red",
               aes(x = as.Date("2020-01-01"), y = y, label = paste0(source," : ", round(slope_sat, 2)," / year"), vjust = -0.5)) +
    facet_wrap(~zone_pretty,nrow = 2) +
    scale_colour_manual(values = c("black", "red")) +
    # scale_y_continuous(limits = c(min(df_var_sub$value, na.rm = TRUE), max(df_var_sub$value, na.rm = TRUE)*1.1)) +
    coord_cartesian(ylim = y_lim) +
    labs(title = plot_title, subtitle = plot_subtitle,
         y = y_lab, x = NULL, colour = "Type", shape = "Source", linetype = "Source") +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 25, face = "bold"), 
          strip.text = element_text(size = 20),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 20),
          axis.title = element_text(size = 23),
          axis.text = element_text(size = 20),
          legend.position = "bottom")
  # plot_zone_var_TS
  ggsave(paste0("figures/validation/ts/ts_",sat_name,"_",var_combi,"_",median_base,".png"), plot_zone_var_TS, width = 14, height = 8)
  return()
}

# The figure code wrapper
# var_name = "SPM"; match_up_df = zone_in_situ_SEXTANT; match_up_stats = stats_SEXTANT; sat_name = "SEXTANT"
# sat_name = sat_name; match_up_df = zone_all_in_situ; match_up_stats = zone_all_in_situ_stats; var_combi = unique(zone_all_in_situ_stats$variable_combi)[1]
validation_plots <- function(var_combi, sat_name, median_base, match_up_df, match_up_stats){

  # Subset datasets for chosen variable
  match_up_df_var <- match_up_df |> 
    filter(variable_combi == var_combi) |>
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  match_up_stats_var <- match_up_stats |> 
    filter(variable_combi == var_combi,
           zone == "GLOBAL",
           source == "ALL",
           site == "ALL",
           season == "ALL")
  
  # Get the variable names etc. for plotting
  var_name_is <- match_up_df_var$variable[1]
  var_name_sat <- match_up_df_var$variable_sat[1]
  cor_name <- match_up_df_var$correction[1]
  proc_name <- match_up_df_var$processing[1]
  grid_name <- match_up_df_var$grid_size[1]
  
  # Get axis labels
  plot_meta_is <- var_labels(var_name_is)
  plot_meta_sat <- var_labels(var_name_sat)

  # Get 1:1 line limits
  identity_line <- data.frame(x = plot_meta_is$axis_limits, 
                              y = plot_meta_sat$axis_limits)

  # Get colour values
  colour_values <- colours_of_stations()
  
  # Get stats to plot
  if(var_name_sat %in% c("SST", "SST-NIGHT")){
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope
    # R2_value <- match_up_stats_var$r2 # Not currently calculated
  } else {
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope_log
    # R2_value <- match_up_stats_var$r2_log # Not currently calculated
  }
  
  # Create title and subtitle
  plot_title = paste(sat_name, var_name_sat, "vs.", "in situ", var_name_is)
  plot_subtitle = paste("Correction:", cor_name, "| Processing:", proc_name, "| Grid size:", grid_name)
  
  # Create the plot
  scatterplot <- ggplot(data = match_up_df_var, aes(x = value_in_situ, y = value_satellite)) + 
    
    geom_point(aes(colour = zone_pretty, shape = source), size = 6, show.legend = TRUE) + 
    # TODO: Add 2:1 and 1:2 dashed lines
    geom_line(data = identity_line, aes(x = x, y = y), linetype = "dashed", show.legend = FALSE) +
    
    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('In~situ~measurements~(', plot_meta_is$unit, ')'))) +
    
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Satellite~estimates~(', plot_meta_sat$unit, ')'))) + 
    
    coord_equal(xlim = plot_meta_is$axis_limits, 
                ylim = plot_meta_sat$axis_limits) +
    
    annotate(geom = 'text', x = plot_meta_is$axis_limits[1], y = plot_meta_sat$axis_limits[2], 
             hjust = 0, vjust = 1, color = "black", size = 12.5,
             label = paste('Error = ', round(ifelse(Error_value |> is.numeric(), Error_value, NA), 1), "%\n",
                           'Bias = ', round(ifelse(Bias_value |> is.numeric(), Bias_value, NA), 1), "%\n",
                           # TODO: Change this to show Slope_log or Slope depending on the variable tested (e.g. SST or not)
                           # 'R²_log = ', round(R2_value, 2),"\n", # Not currently calculated
                           'Slope = ', round(ifelse(Slope_value |> is.numeric(), Slope_value, NA), 2),"\n",
                           'n = ', nrow(match_up_df_var), sep = "")) +
            #  label = paste("Slope = ", round(ifelse(Slope_value |> is.numeric(), Slope_value, NA), 2),"\n",
            #                # 'R² = ', round(statistics_values$r2_log, 2),"\n",
            #                "n = ", nrow(match_up_df_var), sep = "")) +
    
    labs(title = plot_title, subtitle = plot_subtitle) +
    
    scale_color_manual(name = "zone", values = colour_values, drop = FALSE) +
    
    # scale_linetype_manual(values = c("Identity line" = "dashed",
    #                                  "Linear regression" = "solid"), name = "") +
    
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 10), order = 1),
           shape = guide_legend(ncol = 1, override.aes = list(size = 10), order = 2),
           # linetype = guide_legend(override.aes = list(color = c("black"), shape = c(NA), linetype = c("dashed")), ncol = 2, order = 2)
    ) +
    
    ggplot_theme() + 
    
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 30),
          legend.margin = margin(5, 10, 5, 5),
          plot.subtitle = element_text(size = 20, hjust = 0.5, color = "black", face = "bold.italic"),
          plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # scatterplot
  
  if(var_name_sat == "SST"){ 
    scatterplot <- scatterplot +  
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      scale_x_continuous(name = parse(text = paste0('In~situ~measurements~(', plot_meta_is$unit, ')'))) +
      scale_y_continuous(name = parse(text = paste0('Satellite~estimates~(', plot_meta_sat$unit, ')')))
  } else {
    scatterplot <- scatterplot + 
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      annotation_logticks()
  }
  # scatterplot
  
  # Add histograms to x and y axes
  scatterplot_with_side_hist <- ggMarginal(scatterplot, type = "histogram", groupFill = TRUE, alpha = 1)
  # scatterplot_with_side_hist
  ggsave(paste0("figures/validation/scatterplot/",sat_name,"_",var_combi,".png"), 
         plot = scatterplot_with_side_hist, height = 16, width = 16, bg = "white")
  
  # Barplot of frequency per year
  bar_plot_freq_per_year <- match_up_df_var |> 
    mutate(Year = year(date)) |> 
    dplyr::count(Year) |> 
    ggplot(aes(x = Year)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_is$var_colour, linewidth = 1.5) + 
    scale_x_continuous(breaks = seq(1998, 2025, by = 3), name = "", labels = function(x) substring(x, 3, 4)) +
    scale_y_continuous(expand = c(0,0), name = "n per year") +
    labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_year
  ggsave(paste0("figures/validation/barplot/annual_",sat_name,"_",var_combi,".png"), 
         plot = bar_plot_freq_per_year, height = 10, width = 14, bg = "white")
  
  # Barplot of monthly counts
  bar_plot_freq_per_month <- match_up_df_var |> 
    mutate(Month = month(date)) |> 
    dplyr::count(Month) |> 
    ggplot(aes(x = Month)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_is$var_colour, linewidth = 2) + 
    scale_x_continuous(breaks = 1:12, name = "", labels = function(x) month.abb[x]) +
    scale_y_continuous(expand = c(0,0), name = "n per month") +
    coord_cartesian(xlim = c(1,12)) +
    labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_month
  ggsave(paste0("figures/validation/barplot/monthly_",sat_name,"_",var_combi,".png"), 
         plot = bar_plot_freq_per_month, height = 10, width = 14, bg = "white")
  return()
}

# Run all validation stats and produce the plots
# sat_name = "SEXTANT"; median_base = "all"
# sat_name = "MODIS"; median_base = "all"
# sat_name = "OLCI-B"; median_base = "all"
validate_sensor <- function(sat_name, median_base){
  
  # Get the pixel cutoff based on product
  if(median_base == "small"){
    if(sat_name == "SEXTANT"){
      pixel_n_cut <- 1
    } else {
      pixel_n_cut <- 3
    }
  } else {
    if(sat_name == "SEXTANT"){
      pixel_n_cut <- 3
    } else {
      pixel_n_cut <- 16
    }
  }
  
  # Load prepped data
  sat_files <- dir("output/MATCH_UP_DATA/FRANCE", full.names = TRUE, 
                   pattern = paste0("zone_median_",sat_name))
  sat_files <- sat_files[grepl(paste0("_",median_base), sat_files)]
  zone_median <- map_dfr(sat_files, data.table::fread) |> mutate(date = as.Date(date)) |> 
    # Remove all rows that are below the pixel cutoff and CV cutoff of 20%
    mutate(sd = case_when(n <= 2 ~ 0, TRUE ~ sd), # Necessary for CV for 1 pixel count matchups for SEXTANT 'small'
           cv = sd / median) |> 
    filter(n >= pixel_n_cut, cv <= 0.20) #|> 
    # Complete all dates
    # complete(date = seq(min(date), max(date), by = "day"), fill = list(median = NA), 
    #          nesting(zone, source, site, variable))    
  
  # Make variable name conversions as necessary
  if(sat_name == "SEXTANT"){
    # Create a TUR data.frame of SPM data and re-add to the main dataset
    zone_median_tur <- zone_median |> 
      filter(variable == "SPM") |> 
      mutate(variable = "TUR")
    zone_median <- bind_rows(zone_median, zone_median_tur) |> 
      distinct() |>
      mutate(variable_is = variable) |> 
      mutate(variable = case_when(variable == "CHLA" ~ "CHL",
                                  variable == "TUR" ~ "SPM",
                                  TRUE ~ variable)) |> 
      mutate(correction = "Standard",
             processing = "OC5",
             grid_size = case_when(median_base == "small" ~ "1x1",
                                   median_base == "all" ~ "3x3"))
    rm(zone_median_tur); gc()
  } else {
    zone_median <- zone_median |> 
      mutate(variable_is = case_when(grepl("CHL", variable) ~ "CHLA",
                                     grepl("SST", variable) ~ "TEMP",
                                     grepl("SPM", variable) ~ "SPM",
                                     grepl("T-FNU", variable) ~ "TUR",
                                     TRUE ~ variable)) |> 
      mutate(correction = case_when(grepl("-AC", variable) ~ "AC",
                                     grepl("-PO", variable) ~ "PO",
                                     grepl("-NS", variable) ~ "NS")) |> 
      mutate(processing = case_when(grepl("-GONS-", variable) ~ "GONS",
                                     grepl("-OC5-", variable) ~ "OC5",
                                     grepl("-G-", variable) ~ "G",
                                     grepl("-R-", variable) ~ "R",
                                     grepl("-FNU-", variable) ~ "FNU")) |> 
      mutate(grid_size = case_when(median_base == "small" ~ "3x3",
                                  median_base == "all" ~ "7x7")) |> 
      # Clean up variable names for further use
      mutate(variable = gsub("-AC|-PO|-NS", "", variable)) |> 
      mutate(variable = gsub("-GONS|-OC5|-G|-R|-FNU", "", variable))
  }
  
  # Load in situ data and complete the date column
  zone_data_in_situ <- read_csv("data/INSITU_data/zone_data_in_situ.csv", show_col_types = FALSE) |> 
    dplyr::select(-lon, -lat) #|> 
    # complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA), 
    #          nesting(zone, zone_pretty, source, site, time_class, variable))
  
  # Combine extracted sat data with in situ
  zone_all_in_situ_base <- zone_data_in_situ |> 
    left_join(zone_median, by = c("zone", "source", "site", "date", "variable" = "variable_is"),
              relationship = "many-to-many") |> # For day and night time classes for SST
    dplyr::rename(value_in_situ = value, value_satellite = median, variable_sat = variable.y) |> 
    filter(!is.na(variable_sat)) |> 
    filter(!(variable_sat == "SST" & time_class == "night")) |>
    filter(!(variable_sat == "SST-NIGHT" & time_class == "day")) |>
    mutate(season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter", month(date) %in% 3:5  ~ "Spring",
      month(date) %in% 6:8  ~ "Summer", month(date) %in% 9:11 ~ "Autumn"), .after = "date") |> 
    dplyr::select(zone, dplyr::everything(), -time_class) |> 
    mutate(variable_combi = paste(variable, variable_sat, correction, processing, grid_size, sep = "_"))

  # Create monthly average TS for lm analysis
  zone_all_in_situ_monthly <- zone_all_in_situ_base |> 
    mutate(date = floor_date(date, unit = "month")) |> 
    summarise(value_in_situ = mean(value_in_situ, na.rm = TRUE),
              value_satellite = mean(value_satellite, na.rm = TRUE),
              .by = c("variable_combi", "correction", "processing", "grid_size", "zone", "zone_pretty", "source", "season", "date", "variable", "variable_sat"))
  
  # Calculate linear model stats to look at change over time
  zone_in_situ_monthly_lm <- zone_all_in_situ_monthly |> 
    filter(value_in_situ > 0) |> 
    group_by(variable_combi, correction, processing, grid_size, zone, zone_pretty, source, variable, variable_sat) |> 
    do(broom::tidy(lm(value_in_situ ~ date, data = .))) |> 
    filter(term == "date") |> 
    dplyr::rename(slope_is = estimate, p_is = p.value) |> 
    dplyr::select(variable_combi:variable_sat, slope_is, p_is)
  zone_sat_monthly_lm <- zone_all_in_situ_monthly |> 
    filter(value_satellite > 0) |> 
    group_by(variable_combi, correction, processing, grid_size, zone, zone_pretty, source, variable, variable_sat) |> 
    do(broom::tidy(lm(value_satellite ~ date, data = .))) |> 
    filter(term == "date") |> 
    dplyr::rename(slope_sat = estimate, p_sat = p.value) |> 
    dplyr::select(variable_combi:variable_sat, slope_sat, p_sat)
  
  # Combine results
  zone_all_monthly_lm <- left_join(zone_in_situ_monthly_lm, zone_sat_monthly_lm,
                                   by = join_by(variable_combi, correction, processing, grid_size, zone, zone_pretty, 
                                                source, variable, variable_sat)) |> 
    # Convert to values / year
    mutate(slope_is = slope_is*365.25, slope_sat = slope_sat*365.25)
  
  # Save statistics
  write_csv(zone_all_monthly_lm, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/",sat_name,"_lm_stats_",median_base,".csv"))

  # Plot the linear TS matchups
  plan(multicore, workers = 10)
  future_walk(unique(zone_all_in_situ_monthly$variable_combi), validation_lm_plots, 
              sat_name = sat_name, median_base = median_base,
              df = zone_all_in_situ_monthly, df_stats = zone_all_monthly_lm)
  plan(sequential)
  
  # Remove any missing values for stats matchups
  zone_all_in_situ <- zone_all_in_situ_base |> 
    filter(value_in_situ > 0, value_satellite > 0)
  
  # Create big grid of sites to look for obvious outliers
  # ggplot(data = zone_in_situ, aes(x = value_in_situ, y = value_satellite)) +
  #   geom_point(aes(colour = source, shape = variable)) +
  #   facet_wrap(~site, scales = "free")
  
  # Stats for all groups together by variable
  zone_in_situ_stats_01 <- zone_all_in_situ |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all groups together by variable and season
  zone_in_situ_stats_02 <- zone_all_in_situ |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones by variable
  zone_in_situ_stats_03 <- zone_all_in_situ |> mutate(source = "ALL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones by variable and season
  zone_in_situ_stats_04 <- zone_all_in_situ |> mutate(source = "ALL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sources by variable
  zone_in_situ_stats_05 <- zone_all_in_situ |> mutate(zone = "GLOBAL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sources by variable and season
  zone_in_situ_stats_06 <- zone_all_in_situ |> mutate(zone = "GLOBAL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones and sources by variable
  zone_in_situ_stats_07 <- zone_all_in_situ |> mutate(site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones and sources by variable and season
  zone_in_situ_stats_08 <- zone_all_in_situ |> mutate(site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sites by variable
  zone_in_situ_stats_09 <- zone_all_in_situ |> mutate(site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sites by variable and season
  zone_in_situ_stats_10 <- zone_all_in_situ |> mutate(site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for each site by variable
  zone_in_situ_stats_11 <- zone_all_in_situ |> mutate(season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for each site by variable and season
  zone_in_situ_stats_12 <- zone_all_in_situ |>
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Bind all together
  zone_all_in_situ_stats <- bind_rows(zone_in_situ_stats_01, zone_in_situ_stats_02, zone_in_situ_stats_03,
                                      zone_in_situ_stats_04, zone_in_situ_stats_05, zone_in_situ_stats_06,
                                      zone_in_situ_stats_07, zone_in_situ_stats_08, zone_in_situ_stats_09,
                                      zone_in_situ_stats_10, zone_in_situ_stats_11, zone_in_situ_stats_12) |> 
    mutate(sensor = sat_name, .before = "correction")
  
  # Save results
  write_csv(zone_all_in_situ_stats, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/",sat_name,"_stats_",median_base,".csv"))
  rm(zone_in_situ_stats_01, zone_in_situ_stats_02, zone_in_situ_stats_03,
     zone_in_situ_stats_04, zone_in_situ_stats_05, zone_in_situ_stats_06,
     zone_in_situ_stats_07, zone_in_situ_stats_08, zone_in_situ_stats_09,
     zone_in_situ_stats_10, zone_in_situ_stats_11, zone_in_situ_stats_12); gc()
  
  # Run all of the plots per variable pairing
  zone_all_in_situ_stats$variable_combi <- paste0(zone_all_in_situ_stats$variable,"_",
                                                 zone_all_in_situ_stats$variable_sat, "_",
                                                 zone_all_in_situ_stats$correction, "_",
                                                 zone_all_in_situ_stats$processing, "_",
                                                 zone_all_in_situ_stats$grid_size)
  plan(multicore, workers = 10)
  future_walk(unique(zone_all_in_situ_stats$variable_combi), validation_plots,
                     sat_name = sat_name, median_base = median_base,
                     match_up_df = zone_all_in_situ, match_up_stats = zone_all_in_situ_stats)
  plan(sequential)
}


# Tables ---------------------------------------------------------------

# Create pretty tables
# file_path = "output/MATCH_UP_DATA/FRANCE/STATISTICS/"; sat_name = "SEXTANT"
# TODO: This will need to be modified to work with other satellite products
validation_tables <- function(file_path, sat_name) {
  
  # Load all output stats
  # TODO: Clean up zone names and order by latitude
  stat_files <- map_dfr(dir(file_path, pattern = paste0(sat_name,"_stats_all"), full.names = TRUE), read_csv) |> 
    filter(zone != "GLOBAL", source == "ALL", site == "ALL", season == "ALL") |> 
    dplyr::select(zone, Bias, Error, n, variable) |> 
    mutate_at(vars(-zone, -variable), ~ round(., 1)) |> 
    mutate(var_label = paste("When compared with in situ", ifelse(variable == 'TUR', 'TURB', variable))) |>
    dplyr::select(-variable)
  # dplyr::rename(`Bias (%)` = Bias, `Error (%)` = Error)
  
  desired_colnames <- names(stat_files) |> 
    str_replace_all("Bias", "Bias (%)") |> 
    str_replace_all("Error", "Error (%)") |> 
    str_remove_all("var_label")
  names(desired_colnames) <- names(stat_files)
  
  # Create the SPM/TUR table
  # TODO: Wrap this into a function to be called per variable
  table_SPM_TUR <- stat_files |>
    filter(!grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite SPM', subtitle = 'Compared with SPM and Turbidity in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    # Pretty tabs
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3), # A bit more padding for summaries
                row_group.padding = px(4)) |> # And even more for our groups
    # More styling
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Create the CHLA table
  table_CHLA <- stat_files |>
    filter(grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite CHL', subtitle = 'Compared with CHLA in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3),
                row_group.padding = px(4)) |> 
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Save and exit
  ## SPM/TUR
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.html"), inline_css = TRUE)
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.png"), expand = 10)
  ## CHLA
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.html"), inline_css = TRUE)
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.png"), expand = 10)
}

