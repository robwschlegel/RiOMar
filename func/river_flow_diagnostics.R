# func/river_flow_diagnostics.R

# Diagnostic plots for the gauge-extension fits in func/river_flow_prep.R.
# Not part of the production pipeline -- run interactively to inspect how
# well a given river's target/donor correction performs before trusting its
# output CSV. Source river_flow_prep.R first (defines load_gauge,
# predict_move1, predict_qppq, river_config).

source("func/river_flow_prep.R")

# Time series + log-log scatter for one "extend"-type river in river_config.
# Shows the target (real), donor (rescaled onto the target's fitted scale),
# and the final combined/flagged series together.
plot_gauge_extension <- function(river_name) {

  cfg <- river_config[river_config$river == river_name, ]
  if (nrow(cfg) == 0) stop("Unknown river: ", river_name)
  if (cfg$type != "extend") stop(river_name, " is not an 'extend'-type river")

  dir_HP <- file.path("data/RIVER_FLOW", cfg$zone, "HydroPortail")
  target <- load_gauge(cfg$target_pattern, dir_HP)
  donor <- load_gauge(cfg$donor_pattern, dir_HP)
  combined <- combine_gauges(cfg$target_pattern, cfg$donor_pattern, dir_HP,
                              method = cfg$method, flag_below = cfg$flag_below)

  joined <- full_join(target, donor, by = "date", suffix = c("_target", "_donor"))
  predict_fun <- if (cfg$method == "move1") predict_move1 else predict_qppq
  joined$debit_pred <- predict_fun(pred_x = joined$debit_target, pred_y = joined$debit_donor)

  p_scatter <- ggplot(joined, aes(x = debit_donor, y = debit_target)) +
    geom_point(alpha = 0.3) +
    geom_line(aes(y = debit_pred), colour = "red") +
    scale_x_log10() + scale_y_log10() +
    labs(title = paste(river_name, "-", cfg$method, "fit"),
         x = paste("donor:", cfg$donor_pattern), y = paste("target:", cfg$target_pattern))

  p_series <- combined |>
    ggplot(aes(x = date, y = debit, colour = flag)) +
    geom_line() +
    labs(title = paste(river_name, "- combined series"), y = "River flow (m³ s⁻¹)")

  list(scatter = p_scatter, series = p_series, joined = joined, combined = combined)
}

# Residual summary (measured vs. reconstructed agreement in the overlap) for
# every "extend"-type river, as a quick before-you-trust-it check.
summarise_gauge_extensions <- function() {
  extend_rivers <- river_config$river[river_config$type == "extend"]
  map_dfr(extend_rivers, function(r) {
    cfg <- river_config[river_config$river == r, ]
    dir_HP <- file.path("data/RIVER_FLOW", cfg$zone, "HydroPortail")
    combined <- combine_gauges(cfg$target_pattern, cfg$donor_pattern, dir_HP,
                                method = cfg$method, flag_below = cfg$flag_below)
    tibble::tibble(
      river = r,
      method = cfg$method,
      n_total = nrow(combined),
      n_measured = sum(combined$flag == "measured", na.rm = TRUE),
      n_reconstructed = sum(combined$flag == "reconstructed", na.rm = TRUE),
      n_low_flow_flagged = sum(combined$flag == "measured_low_flow_uncertain", na.rm = TRUE),
      n_gap = sum(is.na(combined$debit))
    )
  })
}

# Sevre Niortaise: the two donor-based QPPQ transfers (Niort, Azay) plotted
# against Marans (the mouth-most/target gauge), plus how they're blended.
plot_sevre_niortaise <- function() {
  dir_HP <- file.path("data/RIVER_FLOW", "BAY_OF_BISCAY", "HydroPortail")
  sevre_niort <- combine_gauges("N611061001", "N430062201", dir_HP, method = "qppq")
  sevre_azay <- combine_gauges("N611061001", "N401061001", dir_HP, method = "qppq")

  bind_rows(
    mutate(sevre_niort, donor = "Niort"),
    mutate(sevre_azay, donor = "Azay")
  ) |>
    ggplot(aes(x = date, y = debit, colour = flag)) +
    geom_line() +
    facet_wrap(~donor, ncol = 1) +
    labs(title = "Sevre Niortaise - independent donor transfers onto Marans' scale")
}

# How far outside its MOVE.1 calibration range each donor value used to
# reconstruct post-overlap target values actually falls. predict_move1()
# fits its log-log slope/intercept only on dates where target and donor
# overlap, then applies that fit to every donor value going forward with no
# in-range check -- this quantifies how much of that forward application is
# pure extrapolation vs. still inside the calibrated donor range, and by how
# far in log10-donor-units. river_name defaults to both Rhone branches,
# whose donor (Tarascon) record continues past the target (Arles/Fourques)
# gauges' 2022 cutoff (see river_config's comment on the Rhone split).
summarise_move1_extrapolation <- function(river_names = c("grand_rhone", "petit_rhone")) {
  map_dfr(river_names, function(r) {
    cfg <- river_config[river_config$river == r, ]
    if (cfg$method != "move1") stop(r, " is not a move1-method river")
    dir_HP <- file.path("data/RIVER_FLOW", cfg$zone, "HydroPortail")
    target <- load_gauge(cfg$target_pattern, dir_HP)
    donor <- load_gauge(cfg$donor_pattern, dir_HP)

    joined <- full_join(target, donor, by = "date", suffix = c("_target", "_donor"))
    calib <- joined |> filter(!is.na(debit_target), !is.na(debit_donor),
                              debit_target > 1, debit_donor > 1)
    log_donor_calib_range <- range(log(calib$debit_donor))

    reconstructed <- joined |> filter(is.na(debit_target), !is.na(debit_donor), debit_donor > 1) |>
      mutate(log_donor = log(debit_donor),
             out_of_range = log_donor < log_donor_calib_range[1] | log_donor > log_donor_calib_range[2],
             log10_units_beyond_range = pmax(0, log_donor_calib_range[1] - log_donor, log_donor - log_donor_calib_range[2]) / log(10))

    tibble::tibble(
      river = r,
      calibration_period = paste(format(range(calib$date), "%Y"), collapse = " to "),
      calib_donor_range_m3s = paste(round(exp(log_donor_calib_range), 1), collapse = " to "),
      n_reconstructed_post_overlap = nrow(reconstructed),
      n_out_of_range = sum(reconstructed$out_of_range),
      pct_out_of_range = round(100 * mean(reconstructed$out_of_range), 1),
      max_log10_units_beyond_range = round(max(reconstructed$log10_units_beyond_range), 2),
      max_donor_value_m3s = round(max(reconstructed$debit_donor), 1)
    )
  })
}

# Example usage:
# plot_gauge_extension("grand_rhone")$series
# plot_gauge_extension("grand_rhone")$scatter
# summarise_gauge_extensions()
# plot_sevre_niortaise()
# summarise_move1_extrapolation()
