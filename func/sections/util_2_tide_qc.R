# Tide gauge sub-daily QC ---------------------------------------------------
# Flags calendar days whose sub-daily curve does not look tidal, so that
# load_tide_gauge() (below) can null out tide_mean/tide_range for those days
# before they reach func/multi.R's driver comparisons. See func/tide.R for
# the diagnostics (example good/bad days, per-station counts) behind this and
# for the reasoning behind the threshold choices below.
#
# Harmonic model: M2 + S2 + K1 + O1 (dominant semi-diurnal + diurnal tidal
# constituents). Confirmed empirically for all four gauges via the tidal form
# number F = (K1+O1)/(M2+S2) (func/tide.R): all four are semidiurnal or
# mixed-mainly-semidiurnal, so the same 4-constituent model is used
# everywhere -- only the pass/fail thresholds differ by gauge.
tide_periods_hr <- c(M2 = 12.4206012, S2 = 12.0, K1 = 23.93447213, O1 = 25.81933871)

# R^2 floor for the local harmonic fit (see .tide_day_qc()), below which a
# day's tide_range/tide_mean is considered unreliable. Two tiers, not one:
# the Atlantic/Channel gauges (Le Havre, Port-Bloc, Saint-Nazaire) have a
# strong, clean semidiurnal signal (M2 amplitude 1.4-2.5 m) -- even their
# worst 1st-percentile day still fits the harmonic model with R^2 > 0.92, so
# 0.85 only catches genuine failures. The Mediterranean gauge (Marseille) has
# a tiny tidal signal (M2 amplitude ~0.06 m) that is easily swamped by
# non-tidal sea-level noise (storm surge, seiches) -- plenty of genuinely
# good days sit at R^2 0.7-0.9, so 0.30 is used instead: below that, the
# day's water level is no longer tide-dominated, whatever the physical
# cause, so it is not a reliable tidal-range value regardless.
tide_r2_threshold <- c("LE_HAVRE" = 0.85, "PORT-BLOC" = 0.85, "SAINT-NAZAIRE" = 0.85,
                       "MARSEILLE" = 0.30)

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
  future::plan(future::multisession, workers = parallel::detectCores() - 4)
  out <- furrr::future_map_dfr(days, .tide_day_qc, t_all = df_tide$t, y_all = df_tide$tide,
                               station = station, spike_thresh = spike_thresh,
                               .options = furrr::furrr_options(seed = TRUE))
  future::plan(future::sequential)
  out
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
    # tz = "UTC" matches the raw file's own declared "Fuseau horaire : UTC"
    # header. Without it, as.POSIXct() uses the system's local timezone,
    # which silently produces NA for the ~2 hours/year that don't exist
    # locally across a spring-forward DST transition (harmless on a UTC- or
    # non-DST-observing machine, but real on e.g. Europe/Paris) -- these NAs
    # then propagate into qc_tide_days()'s unique(as.Date(t)) as a literal
    # NA "day", which .tide_day_qc() can't format into a parseable string.
    mutate(t = as.POSIXct(t, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")) |>
    dplyr::distinct(t, .keep_all = TRUE) |>
    dplyr::arrange(t) |>
    dplyr::select(t, tide)
}


