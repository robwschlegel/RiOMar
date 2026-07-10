#!/usr/bin/env Rscript
# code/6_driver_interactions.R
#
# Multi-driver interaction analysis for river plume size, implementing the
# road map agreed in manuscript/driver_interactions_review.md Section 3
# (steps 1-6). Feeds manuscript.tex Section 3.3 "Driver analysis"
# (subsections 3.3.1-3.3.5: River discharge / Wind stress / Wave height /
# Tidal range / Multi-driver) and Figure 10 (the GAM visualisation).
#
# Road map step -> what this script does -> where the output goes:
#   1. Baseline additive GLM              -> fit_baseline_glm()          -> output/STATS/driver_glm_comparison.csv
#   2. + pairwise interaction terms, LRT/AIC -> fit_interaction_glm() / compare_glms() -> output/STATS/driver_glm_comparison.csv
#   3. GAM with te() tensor smooths       -> fit_gam() / plot_gam_figure()  -> output/STATS/driver_gam_summary.csv, manuscript/figures/Figure_10.png
#   4. Regime stratification (Rossby proxy, 6 m/s wind threshold, tidal bin) -> refit_by_regime() -> output/STATS/driver_regime_glm.csv
#   5. Per-metric models (area, centroid lon/lat, not just area)         -> fit_metric_models() over multiple response variables -> output/STATS/driver_metric_models_<response>.csv
#   6. Exploratory random forest + iml H-statistic                       -> fit_rf_diagnostic() -> output/STATS/driver_rf_importance.csv, output/STATS/driver_rf_interaction_hstat.csv
#
# Step 7 (energy-budget / instability-mode decomposition) is explicitly
# out of scope per the review -- it needs a dynamical model or
# high-frequency current/density data RiOMar doesn't currently have, and is
# flagged there as future work for the Discussion/Conclusion, not this
# script.
#
# Known simplifications / gaps, so nothing here is mistaken for more
# rigorous than it is:
#   - "Rossby number" in add_regime_labels() is NOT the true discharge
#     Rossby number from Fong & Geyer (2002) (Ro = U / (f * L), needing
#     mouth width and outflow velocity); it is a per-zone standardised
#     discharge (z-score) used purely as a cheap, ordered stand-in for
#     "is today a high- or low-discharge day relative to this river's own
#     range". Replace with the real Ro if/when mouth geometry (width,
#     depth) is added to river_mouths / zone_meta in func/util.R.
#   - Wave height is one of the drivers named throughout
#     driver_interactions_review.md and manuscript.tex (Figure 8, Table 6
#     row "Wave height") but is not currently loaded anywhere in func/
#     (no load_wave_height()-equivalent exists). build_driver_matrix()
#     below is written so a wave_height column just needs to be left_join()-ed
#     in once that loader exists -- every GLM/GAM/RF formula here is built
#     from whatever driver columns are present in the data, not a hard-coded
#     list, so wave height will be picked up automatically.
#   - "Centroid" per-metric models use the unweighted mean pixel lon/lat of
#     the daily plume mask as a simple proxy for plume position (Ralston et
#     al. 2024's alongshore/cross-shore split is not reproduced here, since
#     that needs a coastline-following coordinate rotation per zone this
#     project doesn't yet have). SPM mass (also named in Ralston et al. as a
#     separate response worth modelling) is not modelled here because it is
#     not in the daily plume time series RiOMar currently outputs (only a
#     confidence_index_in_perc column) -- see manuscript.tex Table 5's
#     existing \tocheck{} placeholder for the same gap.
#   - No R interpreter was available when this script was written, so it has
#     been checked for syntactic and semantic consistency against func/
#     util.R and func/multi.R (which it depends on) but has NOT been
#     executed. Run it end to end and sanity-check a zone or two against
#     the plots before trusting the numbers in the manuscript.


# Setup ---------------------------------------------------------------------
# Like every other code/*.R and code/*.py stage script in this pipeline
# (see CLAUDE.md), this must be run from the repo root, e.g.:
#   Rscript code/6_driver_interactions.R

if(!dir.exists("func")) stop("Run this script from the repo root (e.g. `Rscript code/6_driver_interactions.R`), not from inside code/.")

library(tidyverse)
library(mgcv)       # GAM with tensor-product smooths (step 3)
library(broom)       # tidy() model summaries (steps 1,2,4,5)
library(patchwork)
library(gratia)
library(ranger)
library(iml)

# Run multi-driver analyses and load project common functions
source("func/multi.R")

if (!dir.exists("output/STATS")) dir.create("output/STATS", recursive = TRUE)

# Step 0: build the daily multi-driver data frame per zone ------------------
# One row per day, one column per driver, joined on date. Every downstream
# function below derives its formula from whichever driver columns are
# actually present, so adding a new driver (e.g. wave height) here is enough
# to have it flow through every step automatically.

build_driver_matrix <- function(zone_name){

  meta <- get_zone_meta(zone_name = zone_name)

  df_plume <- load_plume_ts(zone_name)
  df_flow  <- load_driver("flow", meta) |> dplyr::select(date, flow = value)
  df_tide  <- load_driver("tide", meta) |> dplyr::select(date, tide_range = value)
  df_wind  <- load_driver("wind", meta) |> dplyr::select(date, wind_spd = value, wind_dir, direction)

  df <- df_plume |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::left_join(df_tide, by = "date") |>
    dplyr::left_join(df_wind, by = "date")

  # No ROFI/current data for the Gulf of Lion (see func/multi.R::run_driver_suite())
  if(zone_name != "GULF_OF_LION"){
    df_current <- tryCatch(load_driver("current", meta) |> dplyr::select(date, current = value),
                           error = function(e) NULL)
    if(!is.null(df_current)) df <- dplyr::left_join(df, df_current, by = "date")
  }

  df |>
    dplyr::mutate(zone = zone_name, .before = "date") |>
    zoo::na.trim()
}

driver_matrices <- purrr::map(zones, build_driver_matrix) |> purrr::set_names(zones)

# Names of the driver columns available for each zone (varies: no "current"
# for GULF_OF_LION; will pick up "wave_height" automatically once available).
available_drivers <- function(df){
  intersect(c("flow", "wind_spd", "tide_range", "current", "wave_height"), names(df))
}


# Step 1: baseline additive GLM ----------------------------------------------
# manuscript.tex Sec. 2.7 "as currently written" -- the reference point every
# later step is compared against.

fit_baseline_glm <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  form <- stats::as.formula(paste(response, "~", paste(drivers, collapse = " + ")))
  stats::glm(form, data = df, family = gaussian())
}


# Step 2: + pairwise interaction terms ---------------------------------------
# plume_area ~ flow + wind_spd + tide_range (+ current, + wave_height) +
#              all pairwise interactions among available drivers.
# Compared against the Step 1 model with a likelihood-ratio test and AIC.

fit_interaction_glm <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  pair_terms <- utils::combn(drivers, 2, FUN = function(x) paste(x, collapse = ":"))
  form <- stats::as.formula(paste(response, "~", paste(c(drivers, pair_terms), collapse = " + ")))
  stats::glm(form, data = df, family = gaussian())
}

compare_glms <- function(zone_name, response = "plume_area"){
  df <- driver_matrices[[zone_name]]
  m0 <- fit_baseline_glm(df, response)
  m1 <- fit_interaction_glm(df, response)
  lrt <- stats::anova(m0, m1, test = "Chisq")
  tibble::tibble(zone = zone_name, response = response,
                 aic_additive = stats::AIC(m0), aic_interaction = stats::AIC(m1),
                 deviance_additive = m0$deviance, deviance_interaction = m1$deviance,
                 lrt_chisq = lrt$Deviance[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chi)`[2])
}

glm_comparison_stats <- purrr::map_dfr(zones, compare_glms)
readr::write_csv(glm_comparison_stats, "output/STATS/driver_glm_comparison.csv")
# -> Sec 3.3.5 "Multi-driver": "do drivers interact" answered per zone via lrt_p.


# Step 3: GAM with tensor-product smooths ------------------------------------
# One te(driver_i, driver_j) term per pairwise combination of available
# drivers. If a particular pair was non-significant in Step 2's interaction
# GLM, it is still included here (per driver_interactions_review.md: "or for
# all pairs if you'd rather not pre-select") -- restrict pair_terms below to
# only the significant pairs from glm_comparison_stats if a lighter-weight
# per-pair LRT is preferred instead.

fit_gam <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  pair_terms <- utils::combn(drivers, 2, simplify = FALSE)
  te_terms <- purrr::map_chr(pair_terms, ~ paste0("te(", .x[1], ", ", .x[2], ")"))
  form <- stats::as.formula(paste(response, "~", paste(te_terms, collapse = " + ")))
  mgcv::gam(form, data = df, method = "REML")
}

gam_models <- purrr::map(zones, ~ fit_gam(driver_matrices[[.x]])) |> purrr::set_names(zones)

gam_summary <- purrr::imap_dfr(gam_models, function(m, zone_name){
  s <- summary(m)
  tibble::tibble(zone = zone_name, term = rownames(s$s.table),
                 edf = s$s.table[, "edf"], ref.df = s$s.table[, "Ref.df"],
                 F = s$s.table[, "F"], p = s$s.table[, "p-value"],
                 r_sq_adj = s$r.sq, deviance_explained = s$dev.expl)
})
readr::write_csv(gam_summary, "output/STATS/driver_gam_summary.csv")

# Figure 10: composite tensor-smooth visualisation across all four zones.
# This replaces the make_figures_tables.R::generate_figure_10_gam() TODO
# stub; see the note added there.
plot_gam_figure <- function(){
  plots <- purrr::imap(gam_models, function(m, zone_name){
    gratia::draw(m) + patchwork::plot_annotation(title = zone_name)
  })
  combined <- patchwork::wrap_plots(plots, ncol = 2)
  ggplot2::ggsave("figures/Figure_10.png", combined, width = 16, height = 12, dpi = 300)
  invisible(combined)
}
plot_gam_figure()


# Step 4: regime stratification ----------------------------------------------
# Discharge "Rossby proxy" (see gap note at top of file), 6 m/s wind
# threshold (Fofonova et al. 2015), and a per-zone median-split tidal-range
# bin (recalibrating Verney et al. 2024's ~6 m absolute threshold, since
# tidal range differs enormously between the Seine and the Rhone -- see
# driver_interactions_review.md Sec. 3, step 4).

add_regime_labels <- function(df){
  df |>
    dplyr::mutate(
      flow_z = as.numeric(scale(flow)),
      discharge_regime = ifelse(flow_z >= 0, "high (>=zone median)", "low (<zone median)"),
      wind_regime = ifelse(wind_spd >= 6, "high (>=6 m/s)", "low (<6 m/s)"),
      tide_bin = ifelse(tide_range >= stats::median(tide_range, na.rm = TRUE),
                        "spring (>=zone median)", "neap (<zone median)")
    )
}

refit_by_regime <- function(zone_name, regime_col, response = "plume_area", min_n = 30){
  df <- add_regime_labels(driver_matrices[[zone_name]])
  drivers <- available_drivers(df)
  form <- stats::as.formula(paste(response, "~", paste(drivers, collapse = " + ")))

  purrr::map_dfr(split(df, df[[regime_col]]), function(sub){
    if(nrow(sub) < min_n) return(NULL)
    m <- stats::glm(form, data = sub, family = gaussian())
    broom::tidy(m) |> dplyr::mutate(n = nrow(sub), regime_value = unique(sub[[regime_col]]))
  }) |>
    dplyr::mutate(zone = zone_name, regime = regime_col, response = response, .before = 1)
}

regime_stats <- purrr::map_dfr(zones, function(z){
  dplyr::bind_rows(
    refit_by_regime(z, "discharge_regime"),
    refit_by_regime(z, "wind_regime"),
    refit_by_regime(z, "tide_bin")
  )
})
readr::write_csv(regime_stats, "output/STATS/driver_regime_glm.csv")
# -> "check whether driver coefficients genuinely differ by regime" (step 4);
#    a self-contained result for Sec 3.3.5 even before the combined GAM is final.


# Step 5: per-metric models (not just plume area) ----------------------------
# Ralston et al. (2024): area, alongshore position, and buoyancy/SPM signature
# respond to different driver combinations, so model each response separately
# rather than collapsing everything into one "plume area" GLM/GAM. Centroid
# lon/lat is the position metric available without a coastline-following
# coordinate transform (see gap note at top of file). mean_SPM_in_the_plume_area
# and mass_SPM_in_the_plume_area_in_g_m come directly from the Results.csv time
# series loaded in build_driver_matrix() and are therefore already in driver_matrices.

compute_daily_centroid <- function(zone_name){
  plume_dir <- paste0("output/panache/", zone_name)
  plume_files <- dir(plume_dir, pattern = ".csv", recursive = TRUE, full.names = TRUE)
  if(length(plume_files) == 0){
    warning("No pixel-level plume files found for ", zone_name, " -- skipping centroid metrics.")
    return(NULL)
  }
  df_plume <- plyr::ldply(plume_files, load_plume_surface, .parallel = TRUE)
  df_plume |>
    dplyr::summarise(centroid_lon = mean(lon, na.rm = TRUE),
                      centroid_lat = mean(lat, na.rm = TRUE),
                      n_pixels = dplyr::n(), .by = "date") |>
    dplyr::mutate(zone = zone_name, .before = "date")
}

# NB: this reloads every daily pixel-level CSV per zone (the same cost as
# func/multi.R::surface_plot()) -- expect this to be the slow part of the
# script. Cache to disk if re-running often.
centroid_matrices <- purrr::map(zones, function(z){
  cent <- compute_daily_centroid(z)
  if(is.null(cent)) return(NULL)
  dplyr::left_join(driver_matrices[[z]], cent |> dplyr::select(-zone), by = "date")
}) |> purrr::set_names(zones)
centroid_matrices <- purrr::compact(centroid_matrices)

fit_metric_models <- function(df, response){
  list(
    glm_additive    = fit_baseline_glm(df, response),
    glm_interaction = fit_interaction_glm(df, response),
    gam             = fit_gam(df, response)
  )
}

metric_responses <- c("plume_area", "mean_SPM_in_the_plume_area", "mass_SPM_in_the_plume_area_in_g_m",
                      "centroid_lon", "centroid_lat")
# plume_area and SPM metrics come from driver_matrices; centroid metrics need centroid_matrices.
ts_responses <- c("plume_area", "mean_SPM_in_the_plume_area", "mass_SPM_in_the_plume_area_in_g_m")
metric_model_stats <- purrr::map(metric_responses, function(resp){
  mats <- if(resp %in% ts_responses) driver_matrices else centroid_matrices
  purrr::imap_dfr(mats, function(df, zone_name){
    if(!(resp %in% names(df))) return(NULL)
    df_resp <- tidyr::drop_na(df, dplyr::all_of(c(resp, available_drivers(df))))
    if(nrow(df_resp) < 30) return(NULL)
    models <- fit_metric_models(df_resp, resp)
    tibble::tibble(zone = zone_name, response = resp,
                   aic_additive = stats::AIC(models$glm_additive),
                   aic_interaction = stats::AIC(models$glm_interaction),
                   gam_r_sq_adj = summary(models$gam)$r.sq,
                   gam_deviance_explained = summary(models$gam)$dev.expl)
  })
}) |> purrr::set_names(metric_responses)

purrr::iwalk(metric_model_stats, function(stats_df, resp){
  readr::write_csv(stats_df, paste0("output/STATS/driver_metric_models_", resp, ".csv"))
})
# -> Sec 3.3.5 "Multi-driver": compare aic_additive/aic_interaction/gam fit
#    quality across plume_area, mean_SPM_in_the_plume_area,
#    mass_SPM_in_the_plume_area_in_g_m, centroid_lon, and centroid_lat.


# Step 6: exploratory random forest + H-statistic ----------------------------
# Purely diagnostic (per driver_interactions_review.md: "treat its output as
# a hypothesis generator rather than a result to publish directly") -- run
# before finalising which pairs get a te() term in Step 3, or after, as a
# sanity check that Step 3 didn't miss an interaction.

fit_rf_diagnostic <- function(zone_name, response = "plume_area"){
  df <- driver_matrices[[zone_name]]
  drivers <- available_drivers(df)
  df_complete <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))

  rf <- ranger::ranger(stats::as.formula(paste(response, "~", paste(drivers, collapse = " + "))),
                       data = df_complete[, c(response, drivers)],
                       importance = "permutation", num.trees = 500)

  predictor <- iml::Predictor$new(rf, data = df_complete[, drivers], y = df_complete[[response]])
  interaction_hstat <- iml::Interaction$new(predictor)$results

  list(model = rf, importance = rf$variable.importance, interaction = interaction_hstat)
}

rf_results <- purrr::map(zones, fit_rf_diagnostic) |> purrr::set_names(zones) |> purrr::compact()

rf_importance <- purrr::imap_dfr(rf_results, function(res, zone_name){
  tibble::tibble(zone = zone_name, driver = names(res$importance), importance = res$importance)
})
readr::write_csv(rf_importance, "output/STATS/driver_rf_importance.csv")

rf_interaction <- purrr::imap_dfr(rf_results, function(res, zone_name){
  dplyr::mutate(res$interaction, zone = zone_name, .before = 1)
})
if(nrow(rf_interaction) > 0) readr::write_csv(rf_interaction, "output/STATS/driver_rf_interaction_hstat.csv")
# -> hypothesis-generation only: cross-check any driver pair with a high
#    H-statistic here against whether it already has a significant te() term
#    in gam_summary (Step 3) or a significant interaction in
#    glm_comparison_stats (Step 2); if not, that's a candidate to add.


message("code/6_driver_interactions.R complete. Outputs written to output/STATS/driver_*.csv and figures/Figure_10.png.")

