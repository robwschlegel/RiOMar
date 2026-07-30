# func/tide.R
# Tide-gauge sub-daily QC diagnostics
#
# DEPRECATED (2026-07) -- the driver-vs-plume comparison workflow this file
# used to hold (tide_calc(), etc.) was merged into func/multi.R as part of
# consolidating flow.R/wind.R/tide.R/ROFI.R/surface.R -- see the
# "REFACTOR (2026-07)" comment block at the top of func/multi.R. That
# workflow now lives at:
#   tide_calc(mouth_row) -> combine_plume_driver("tide", meta) + plot_driver_plume_comparison(df, "tide", meta$mouth_name) + plot_driver_plume_dual_axis(df, "tide", meta$zone)
# The pre-refactor version of this file is preserved in git history
# (`git log -- func/tide.R`).
#
# This file is now repurposed (2026-07) for a second, unrelated job: the
# actual per-day tidal QC *logic* (harmonic fit, spike/extrema/sparse
# checks, per-station thresholds) lives in qc_tide_days() / .tide_day_qc()
# in func/util.R, because that is what load_tide_gauge() -- and therefore
# every downstream driver comparison in multi.R -- consumes. This file is
# the diagnostic companion: it re-derives the tidal regime of each gauge
# (form number), plots concrete good/bad day examples per station, and
# tabulates how many days were flagged and why. Nothing here feeds back into
# the pipeline; it exists purely to document/justify the thresholds in
# util.R. Re-run it whenever those thresholds change.

source("func/multi.R")

tide_stations <- c("FOS-SUR-MER", "LE_HAVRE", "MARSEILLE", "PORT_DE_BOUC", "PORT-BLOC", "SAINT-NAZAIRE")


# Tidal regime per gauge (form number) --------------------------------------

# Form number F = (K1+O1)/(M2+S2) amplitude ratio, from a 4-constituent
# harmonic fit over each gauge's last 3 years of hourly data:
#   F < 0.25              semidiurnal
#   0.25 <= F < 1.5       mixed, mainly semidiurnal
#   1.5 <= F < 3.0        mixed, mainly diurnal
#   F >= 3.0              diurnal
# This is what tide_periods_hr / tide_r2_threshold in func/util.R are
# calibrated against: all six gauges turn out semidiurnal or
# mixed-mainly-semidiurnal, so the same M2+S2+K1+O1 model is used
# everywhere, and only the per-gauge R^2 floor differs (Atlantic/Channel
# vs. Mediterranean, see below).
tide_form_number <- function(station){
  df <- .load_tide_raw(paste0("data/TIDES/", station)) |>
    dplyr::filter(t >= max(t) - 3*365*24*3600)
  hrs <- as.numeric(difftime(df$t, df$t[1], units = "hours"))
  X <- cbind(1, cos(2*pi*hrs/tide_periods_hr[1]), sin(2*pi*hrs/tide_periods_hr[1]),
                cos(2*pi*hrs/tide_periods_hr[2]), sin(2*pi*hrs/tide_periods_hr[2]),
                cos(2*pi*hrs/tide_periods_hr[3]), sin(2*pi*hrs/tide_periods_hr[3]),
                cos(2*pi*hrs/tide_periods_hr[4]), sin(2*pi*hrs/tide_periods_hr[4]))
  co <- .lm.fit(X, df$tide)$coefficients
  amp <- c(M2 = sqrt(co[2]^2+co[3]^2), S2 = sqrt(co[4]^2+co[5]^2),
           K1 = sqrt(co[6]^2+co[7]^2), O1 = sqrt(co[8]^2+co[9]^2))
  data.frame(station = station, t(amp), F = round((amp["K1"]+amp["O1"])/(amp["M2"]+amp["S2"]), 3))
}

# purrr::map_dfr(tide_stations, tide_form_number)


# Per-day QC summary, all gauges ---------------------------------------------

# Re-derives qc_tide_days() output for all six gauges and tabulates how many
# days were flagged, and why. Saved to output/STATS/ alongside the other
# pipeline QC summaries (missing_SPM.csv, missing_chla.csv in multi.R).
tide_qc_all <- function(){
  purrr::map_dfr(tide_stations, function(st){
    df <- load_tide_gauge(paste0("data/TIDES/", st))
    df$station <- st
    df
  })
}

tide_qc_summary_table <- function(df_qc){
  df_qc |>
    dplyr::group_by(station) |>
    dplyr::summarise(n_days = dplyr::n(),
                     n_flagged = sum(is.na(tide_range)),
                     pct_flagged = round(100*n_flagged/n_days, 2),
                     n_sparse = sum(tide_qc_reason == "sparse", na.rm = TRUE),
                     n_spike = sum(tide_qc_reason == "spike", na.rm = TRUE),
                     n_extrema_mismatch = sum(tide_qc_reason == "extrema_mismatch", na.rm = TRUE),
                     n_low_r2 = sum(tide_qc_reason == "low_r2", na.rm = TRUE),
                     .groups = "drop")
}


# Good/bad day plots ----------------------------------------------------

# Plots one calendar day's raw hourly curve (+-1 day of context) against the
# same local harmonic fit qc_tide_days() judged it by, so the reason a day
# was flagged (or not) is visible rather than asserted.
# plot_tide_qc_day(.load_tide_raw("data/TIDES/MARSEILLE"), as.Date("2015-01-02"), "MARSEILLE", NA)
plot_tide_qc_day <- function(df_tide, day, station, reason){

  # NB: no explicit tz -- must match the (system-default) tz df_tide$t was
  # parsed with in .load_tide_raw(), same reasoning as .tide_day_qc() in util.R
  day_start <- as.POSIXct(paste(day, "00:00:00"))
  day_end   <- day_start + 24*3600
  win_start <- day_start - 24*3600
  win_end   <- day_end + 24*3600

  df_win <- dplyr::filter(df_tide, t >= win_start, t < win_end)
  hrs <- as.numeric(difftime(df_win$t, day_start, units = "hours"))
  X <- cbind(1, cos(2*pi*hrs/tide_periods_hr[1]), sin(2*pi*hrs/tide_periods_hr[1]),
                cos(2*pi*hrs/tide_periods_hr[2]), sin(2*pi*hrs/tide_periods_hr[2]),
                cos(2*pi*hrs/tide_periods_hr[3]), sin(2*pi*hrs/tide_periods_hr[3]),
                cos(2*pi*hrs/tide_periods_hr[4]), sin(2*pi*hrs/tide_periods_hr[4]))
  coefs <- .lm.fit(X, df_win$tide)$coefficients
  df_win$fit <- as.vector(X %*% coefs)

  status <- if(is.na(reason)) "GOOD" else paste0("FLAGGED: ", reason)
  plot_title <- paste0(station, " -- ", format(day), " (", status, ")")

  ggplot(df_win, aes(x = t)) +
    annotate("rect", xmin = day_start, xmax = day_end, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
    geom_line(aes(y = fit), colour = "firebrick", linetype = "dashed", linewidth = 0.6) +
    geom_line(aes(y = tide), colour = "black") +
    geom_point(aes(y = tide), size = 1) +
    labs(title = plot_title, x = NULL, y = "tide height (m)",
        subtitle = "black = observed hourly record | red dashed = local M2+S2+K1+O1 harmonic fit | grey band = the flagged calendar day") +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          plot.subtitle = element_text(size = 9))
}

# One good + one bad example day per station, saved to figures/qc/tide_qc/.
# The bad-day reason is chosen per station from whichever check actually
# fired for it (see tide_qc_summary_table()), so each plot demonstrates a
# real, representative failure mode rather than a cherry-picked edge case.
tide_qc_example_days <- data.frame(
  station = tide_stations,
  good_day = as.Date(c("2015-01-02", "2015-01-02", "2015-01-02", "2023-04-29", "2015-01-02", "2015-01-02")),
  bad_day  = as.Date(c("2012-05-17", "2023-03-31", "1849-11-07", "2023-04-28", "2000-11-06", "1992-12-02")),
  stringsAsFactors = FALSE
)

save_tide_qc_examples <- function(){
  dir.create("figures/qc/tide_qc", showWarnings = FALSE, recursive = TRUE)
  for(i in seq_len(nrow(tide_qc_example_days))){
    st <- tide_qc_example_days$station[i]
    raw <- .load_tide_raw(paste0("data/TIDES/", st))
    flags <- qc_tide_days(raw, st)

    good_day <- tide_qc_example_days$good_day[i]
    bad_day  <- tide_qc_example_days$bad_day[i]
    good_reason <- flags$reason[flags$date == good_day]
    bad_reason  <- flags$reason[flags$date == bad_day]

    p_good <- plot_tide_qc_day(raw, good_day, st, good_reason)
    p_bad  <- plot_tide_qc_day(raw, bad_day, st, bad_reason)

    ggsave(paste0("figures/qc/tide_qc/", st, "_good_", good_day, ".png"), p_good, width = 8, height = 4.5, dpi = 200)
    ggsave(paste0("figures/qc/tide_qc/", st, "_bad_", bad_day, ".png"), p_bad, width = 8, height = 4.5, dpi = 200)
  }
}

# Run the diagnostics:
# df_qc <- tide_qc_all()
# tide_qc_summary_table(df_qc)
# write_csv(tide_qc_summary_table(df_qc), "output/STATS/tide_qc_summary.csv")
# save_tide_qc_examples()
