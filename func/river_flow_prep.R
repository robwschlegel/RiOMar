# func/river_flow_prep.R

# This script combines and cleans river discharge data downloaded from the
# Hub'Eau hydrometrie API (see download_hubeau_river_flow() in func/dl.py,
# wired into code/0_download_data.py). Station-to-zone mapping lives in
# metadata/HydroPortail_station_list.csv. For each river, a "target" gauge
# (the one whose scale/location we want the output series to represent,
# usually the one closest to the river mouth) is extended using a "donor"
# gauge with better temporal coverage, via combine_gauges() below. Every
# output CSV carries a `flag` column ("measured" / "reconstructed" /
# "measured_low_flow_uncertain") so reconstructed values are never presented
# as directly observed.
#
# Exploratory/diagnostic plots for these fits live in
# func/river_flow_diagnostics.R, not here.


# Setup -------------------------------------------------------------------

library(tidyverse)

zones_list <- c("GULF_OF_LION", "BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY")


# Record-extension methods -------------------------------------------------

# MOVE.1 (Hirsch 1982, WRR 18(4):1081-1088): a reduced-major-axis (RMA) fit
# through log-transformed values, not OLS. OLS minimises vertical residuals
# only, which attenuates the slope and shrinks the extended series' variance
# toward its mean (flood peaks under-predicted, low flows over-predicted);
# RMA instead preserves the variance ratio between the two gauges, which is
# what a record-extension transfer needs. Good default whenever the two
# gauges' relationship is a reasonably stable power law across the flow
# range they share.
#
# pred_x: the target series (may be gappy) used to calibrate the fit.
# pred_y: the donor series (long/complete); predictions are returned for
#         every date pred_y has, not just the overlap with pred_x.
predict_move1 <- function(pred_x, pred_y) {

  df <- data.frame(pred_x = pred_x, pred_y = pred_y) |>
    filter(pred_x > 1, pred_y > 1) # Log transform falls over under 1

  slope <- sign(cor(log(df$pred_x), log(df$pred_y))) * sd(log(df$pred_y)) / sd(log(df$pred_x))
  intercept <- mean(log(df$pred_y)) - slope * mean(log(df$pred_x))

  # Invert log_y = intercept + slope*log_x for log_x, over the full pred_y range
  exp((log(pred_y) - intercept) / slope)
}

# QPPQ (Fennessey 1994): quantile / flow-duration-curve transfer. Every
# donor value is mapped to its percentile within the donor's own long-term
# flow-duration curve, then that percentile is inverted through the
# target's flow-duration curve (estimated from whatever target data
# exists) to predict the target-equivalent value. Unlike MOVE.1 this does
# not assume one global log-linear relationship holds across the full flow
# range, which makes it a better fit where the gauge relationship is
# regime-dependent rather than a stable power law -- e.g. a low-flow
# rating-curve breakdown (Vilaine/Rieux, officially flagged unreliable
# below 20 m3/s) or tidal/lock backwater effects that only distort part of
# the range (Charente, Sevre Niortaise).
#
# pred_x, pred_y: same convention as predict_move1.
predict_qppq <- function(pred_x, pred_y) {
  donor_ecdf <- ecdf(pred_y)
  target_vals <- pred_x[!is.na(pred_x)]
  p <- donor_ecdf(pred_y)
  as.numeric(quantile(target_vals, probs = p, na.rm = TRUE, type = 7, names = FALSE))
}


# Data loading --------------------------------------------------------------

# Load a single HydroPortail/Hub'Eau CSV export
load_HP <- function(file_name) {
  df <- read_csv(file_name, show_col_types = FALSE) |>
    mutate(date = as.Date(`Date (TU)`)) |>
    filter(`Valeur (en m³/s)` >= 0) |>
    summarise(debit = mean(`Valeur (en m³/s)`, na.rm = TRUE), .by = date)
  return(df)
}

# Load every file in dir_HP matching `pattern` (a station code, or several
# alternated with "|" to stitch successive station codes at the same site
# into one series, e.g. after equipment upgrades) and average by date.
load_gauge <- function(pattern, dir_HP) {
  map_dfr(dir(dir_HP, pattern = pattern, full.names = TRUE), load_HP) |>
    summarise(debit = mean(debit, na.rm = TRUE), .by = date)
}

# Extend a target gauge's record using a donor gauge with better temporal
# coverage. Real target values are always preferred; the donor-derived
# prediction (via `method`) only fills dates the target lacks. Every row is
# flagged "measured" or "reconstructed" so the two are never conflated
# downstream. If flag_below is set, real target values under that threshold
# are flagged "measured_low_flow_uncertain" instead (for gauges with an
# official low-flow rating-curve caveat), rather than being dropped -- the
# aim is to minimise gaps in the 1998-2025 series, not to hide uncertain
# but real readings.
combine_gauges <- function(target_pattern, donor_pattern, dir_HP,
                            method = "move1", flag_below = -Inf) {

  target <- load_gauge(target_pattern, dir_HP)
  donor <- load_gauge(donor_pattern, dir_HP)

  predict_fun <- switch(method,
    move1 = predict_move1,
    qppq = predict_qppq,
    stop("Unknown method: ", method)
  )

  full_join(target, donor, by = "date", suffix = c("_target", "_donor")) |>
    mutate(
      debit_pred = predict_fun(pred_x = debit_target, pred_y = debit_donor),
      debit = ifelse(!is.na(debit_target), debit_target, debit_pred),
      flag = case_when(
        !is.na(debit_target) & debit_target < flag_below ~ "measured_low_flow_uncertain",
        !is.na(debit_target) ~ "measured",
        TRUE ~ "reconstructed"
      )
    ) |>
    dplyr::select(date, debit, flag)
}


# Per-river configuration ----------------------------------------------------

# One row per output river. type "direct" loads target_pattern as-is (one
# station, or several alternated with "|" for a simple no-correction
# stitch/average of equally-trustworthy gauges, e.g. Garonne's two nearby
# stations). type "extend" runs combine_gauges(target_pattern, donor_pattern).
#
# Method choices: MOVE.1 is the default. QPPQ is used specifically where the
# gauge relationship is known to be regime-dependent rather than a stable
# power law -- see the note on predict_qppq() above -- namely the Vilaine
# (Rieux's low-flow rating-curve caveat), Charente and Sevre Niortaise (both
# tidal/lock-influenced at the relevant gauges, per station referential
# comments and regional hydrology documentation).
#
# The Rhone's Grand/Petit split does not need a bespoke partition model:
# fitting each branch (Arles / Fourques) directly against the Tarascon
# gauge via MOVE.1 on the full 1998-2022 triplet overlap already captures
# the flow-dependence (fitted log-log slopes of ~0.99 for Grand Rhone and
# ~1.12 for Petit Rhone, vs. a flat ratio that swings from ~0.90 in floods
# to ~0.98 in baseflow across the same data) with ~5-6% median error --
# the same two-gauge machinery used everywhere else in this table.
river_config <- tibble::tribble(
  ~zone,                ~river,        ~type,     ~target_pattern,          ~donor_pattern,   ~method,  ~flag_below, ~filename,
  "GULF_OF_LION",       "grand_rhone", "extend",  "V730000302",             "V720001002",     "move1",  -Inf,        "grand_rhone.csv",
  "GULF_OF_LION",       "petit_rhone", "extend",  "V730000202",             "V720001002",     "move1",  -Inf,        "petit_rhone.csv",
  "BAY_OF_SEINE",       "seine",       "extend",  "H320000101|H320000104", "H322011003",      "move1",  -Inf,        "seine.csv",
  "BAY_OF_SEINE",       "orne",        "direct",  "I362101001",             NA,               NA,       -Inf,        "orne.csv",
  "SOUTHERN_BRITTANY",  "loire",       "direct",  "L800001020",             NA,               NA,       -Inf,        "loire.csv",
  "SOUTHERN_BRITTANY",  "vilaine",     "extend",  "J930061101",             "J770061002|J770061003", "qppq", 20,     "vilaine.csv",
  "BAY_OF_BISCAY",      "lay",         "extend",  "N3511610",               "N330161010",      "move1", -Inf,        "lay.csv",
  "BAY_OF_BISCAY",      "garonne",     "direct",  "O909001001|O900001002",  NA,               NA,       -Inf,        "garonne.csv",
  "BAY_OF_BISCAY",      "dordogne",    "direct",  "P555001001",             NA,               NA,       -Inf,        "dordogne.csv",
  "BAY_OF_BISCAY",      "charente",    "extend",  "R520001001",             "R222001001",      "qppq",  -Inf,        "charente.csv"
)


# Preparation driver ----------------------------------------------------------

prep_HP <- function(zone) {

  dir_main <- file.path("data/RIVER_FLOW", zone)
  dir_HP <- file.path("data/RIVER_FLOW", zone, "HydroPortail")

  zone_rivers <- river_config[river_config$zone == zone, ]

  pwalk(zone_rivers, function(river, type, target_pattern, donor_pattern, method, flag_below, filename, ...) {
    df <- if (type == "direct") {
      load_gauge(target_pattern, dir_HP) |> mutate(flag = "measured")
    } else {
      combine_gauges(target_pattern, donor_pattern, dir_HP, method = method, flag_below = flag_below)
    }
    df <- df |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |>
      arrange(date) |>
      dplyr::select(date, debit, flag)
    write_csv(df, file.path(dir_main, filename))
  })

  # Sevre Niortaise: three gauges (Marans/mouth-most, Niort, Azay), all
  # tide/lock-influenced (per EPTB Charente-region hydrology documentation
  # on the wider basin) -- each of the two inland gauges is independently
  # transferred onto Marans' scale via QPPQ, then blended (real Marans value
  # preferred; otherwise the mean of whichever donor-derived predictions are
  # available that day), rather than the previous flat multipliers chained
  # through a straight rbind-and-average.
  if (zone == "BAY_OF_BISCAY") {
    sevre_niort <- combine_gauges("N611061001", "N430062201", dir_HP, method = "qppq")
    sevre_azay <- combine_gauges("N611061001", "N401061001", dir_HP, method = "qppq")

    sevre <- full_join(
      dplyr::select(sevre_niort, date, debit_niort = debit, flag_niort = flag),
      dplyr::select(sevre_azay, date, debit_azay = debit, flag_azay = flag),
      by = "date"
    ) |>
      mutate(
        debit = case_when(
          flag_niort == "measured" ~ debit_niort,
          !is.na(debit_niort) & !is.na(debit_azay) ~ (debit_niort + debit_azay) / 2,
          !is.na(debit_niort) ~ debit_niort,
          !is.na(debit_azay) ~ debit_azay,
          TRUE ~ NA_real_
        ),
        flag = ifelse(flag_niort == "measured", "measured", "reconstructed")
      ) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |>
      arrange(date) |>
      dplyr::select(date, debit, flag)

    write_csv(sevre, file.path(dir_main, "sevre_niortaise.csv"))
  }
}

# Run it
walk(zones_list, prep_HP)
