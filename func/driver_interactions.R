# func/driver_interactions.R
#
# Multi-driver interaction analysis for river plume size, implementing the
# road map agreed in manuscript/driver_interactions_review.md Section 3
# (steps 1-6). Feeds manuscript.tex Section 3.3 "Driver analysis"
# (subsections 3.3.1-3.3.5: River discharge / Wind stress / Wave height /
# Tidal range / Multi-driver) and Figure 10 (the GAM visualisation).
#
# Called from code/4_time_series.py via rpy2, the same way func/X11.R,
# func/figure.R, func/validate.R, and func/plume.R are called from their
# respective func/*.py modules -- source() this file, then call
# run_driver_interactions_analysis(). This keeps figure/stat generation
# centralised behind the numbered code/ pipeline (per CLAUDE.md) rather than
# split across a separate, unnumbered driver-interactions script.
#
# NB: ROFI extent is NOT one of the drivers modelled here. ROFI is a
# downstream, model-derived consequence of the same river discharge that
# drives the satellite-detected plume (func/ROFI.R's succession analysis
# shows the plume changing before the ROFI does, both following discharge),
# not an independent physical forcing on plume size -- it is used elsewhere
# in this project purely as an external accuracy check on panache's plume-area
# estimates (Supplementary Fig. S4 in manuscript.tex), never as a model input.
#
# Road map step -> what this script does -> where the output goes:
#   1. Baseline additive GLM              -> fit_baseline_glm()          -> output/STATS/driver_glm_comparison.csv
#   2. + pairwise interaction terms, LRT/AIC -> fit_interaction_glm() / compare_glms() -> output/STATS/driver_glm_comparison.csv
#   3. GAM with te() tensor smooths       -> fit_gam() / plot_gam_figure()  -> output/STATS/driver_gam_summary.csv, figures/ARTICLE/FIGURE_S7/Figure_S7.png (main-text Fig. 10, moved to Supplementary S7 2026-07-31)
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
#     row "Wave height"). It is now loaded via util.R::load_wave() (VHM0/VMDR
#     from the CMEMS IBI/MED wave reanalysis) and joined into
#     build_driver_matrix() as wave_height/wave_dir -- every GLM/GAM/RF
#     formula here is built from whatever driver columns are present in the
#     data, not a hard-coded list, so wave height is picked up automatically.
#     NB (2026-07): the underlying WAVE .nc files were still being downloaded
#     when this loader was wired in, so this path has not yet been run
#     end-to-end against real data -- re-check once the download completes.
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
# Like every other func/*.R module called via rpy2 (see CLAUDE.md), this
# assumes the R session's working directory is the repo root -- true when
# sourced from code/4_time_series.py, which is itself run as
# `python code/4_time_series.py` from the repo root.

if(!dir.exists("func")) stop("func/driver_interactions.R must be sourced with the repo root as the working directory.")

library(tidyverse)
library(mgcv)       # GAM with tensor-product smooths (step 3)
library(broom)       # tidy() model summaries (steps 1,2,4,5)
library(patchwork)
library(gratia)
library(ranger)
library(iml)

# iml::Interaction uses future.apply internally; the R6 predictor object can
# exceed the default 500 MiB per-worker globals limit.
options(future.globals.maxSize = 2 * 1024^3)   # 2 GiB

# Run multi-driver analyses and load project common functions
source("func/multi.R")


# Step 0: build the daily multi-driver data frame per zone ------------------
# One row per day, one column per driver, joined on date. Every downstream
# function below derives its formula from whichever driver columns are
# actually present, so adding a new driver (e.g. wave height) here is enough
# to have it flow through every step automatically.
#
# plume_dir: path to the panache output root (e.g. "output/panache/dynamic").
#            Passed to util.R::load_plume_ts() so the correct threshold run
#            is used.

build_driver_matrix <- function(zone_name, plume_dir = "output/panache/dynamic"){

  meta <- get_zone_meta(zone_name = zone_name)

  df_plume   <- load_plume_ts(zone_name, plume_dir = plume_dir)
  df_flow    <- load_driver("flow", meta) |> dplyr::select(date, flow = value)
  df_tide    <- load_driver("tide", meta) |> dplyr::select(date, tide_range = value)
  df_wind    <- load_driver("wind", meta) |> dplyr::select(date, wind_spd = value, wind_dir, direction)
  df_current <- load_driver("current", meta) |> dplyr::select(date, current = value, current_dir)
  df_wave    <- load_driver("wave", meta) |> dplyr::select(date, wave_height = value, wave_dir)

  df_plume |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::left_join(df_tide, by = "date") |>
    dplyr::left_join(df_wind, by = "date") |>
    dplyr::left_join(df_current, by = "date") |>
    dplyr::left_join(df_wave, by = "date") |>
    dplyr::mutate(zone = zone_name, .before = "date",
                  # Categorical (8-octant) form of direction, since direction
                  # matters as much or more than magnitude for wind/wave/
                  # current -- see multi.R::compass_octant(). Bay of Seine's
                  # wave_dir is all-NA (load_wave()'s NB 2), so wave_dir_cat
                  # is all-NA there too; available_drivers() below drops it
                  # for that zone automatically. current_dir_cat added
                  # 2026-08-05 (per Robert): current_dir was already loaded
                  # above but never turned into a category or added to
                  # available_drivers()'s candidates, so current direction
                  # fed none of the GLM/GAM/RF models -- only current speed
                  # did.
                  wind_dir_cat = compass_octant(wind_dir),
                  wave_dir_cat = compass_octant(wave_dir),
                  current_dir_cat = compass_octant(current_dir)) |>
    zoo::na.trim()
}

# Names of the driver columns available for each zone -- any candidate column
# that is entirely NA for this zone (e.g. wave_dir_cat for the Bay of Seine)
# is dropped rather than passed into a model formula, generalising the same
# per-zone gap-handling already done in util.R::load_wave() and
# multi.R::plot_driver_rose().
available_drivers <- function(df){
  candidates <- c("flow", "wind_spd", "tide_range", "current", "wave_height", "wind_dir_cat", "wave_dir_cat", "current_dir_cat")
  drivers <- intersect(candidates, names(df))
  Filter(function(d) any(!is.na(df[[d]])), drivers)
}


# Step 1: baseline additive GLM ----------------------------------------------

fit_baseline_glm <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  form <- stats::as.formula(paste(response, "~", paste(drivers, collapse = " + ")))
  stats::glm(form, data = df, family = gaussian())
}


# Step 2: + pairwise interaction terms ---------------------------------------

fit_interaction_glm <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  pair_terms <- utils::combn(drivers, 2, FUN = function(x) paste(x, collapse = ":"))
  form <- stats::as.formula(paste(response, "~", paste(c(drivers, pair_terms), collapse = " + ")))
  stats::glm(form, data = df, family = gaussian())
}

compare_glms <- function(zone_name, driver_matrices, response = "plume_area"){
  df <- driver_matrices[[zone_name]]
  drivers <- available_drivers(df)
  df_valid <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))
  if(nrow(df_valid) < 30 || stats::var(df_valid[[response]]) < 1e-6){
    message("compare_glms: skipping ", zone_name, " (insufficient variation in ", response, ")")
    return(NULL)
  }
  m0 <- fit_baseline_glm(df_valid, response)
  m1 <- fit_interaction_glm(df_valid, response)
  lrt <- stats::anova(m0, m1, test = "Chisq")
  tibble::tibble(zone = zone_name, response = response,
                 aic_additive = stats::AIC(m0), aic_interaction = stats::AIC(m1),
                 deviance_additive = m0$deviance, deviance_interaction = m1$deviance,
                 lrt_chisq = lrt$Deviance[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chi)`[2])
}


# Step 3: GAM with tensor-product smooths ------------------------------------

fit_gam <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  df_valid <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))
  if(nrow(df_valid) < 30 || stats::var(df_valid[[response]]) < 1e-6) return(NULL)

  # te() tensor-product smooths need numeric arguments -- categorical drivers
  # (wind_dir_cat/wave_dir_cat) enter as flat parametric main-effect terms
  # instead, alongside the te() smooths over every numeric pair.
  is_categorical <- purrr::map_lgl(drivers, ~ !is.numeric(df_valid[[.x]]))
  numeric_drivers <- drivers[!is_categorical]
  categorical_drivers <- drivers[is_categorical]

  pair_terms <- utils::combn(numeric_drivers, 2, simplify = FALSE)
  te_terms <- purrr::map_chr(pair_terms, ~ paste0("te(", .x[1], ", ", .x[2], ")"))
  form <- stats::as.formula(paste(response, "~", paste(c(te_terms, categorical_drivers), collapse = " + ")))
  mgcv::gam(form, data = df_valid, method = "REML")
}

# Moved 2026-07-31 from main-text Figure 10 to Supplementary Figure S7, per
# Robert's request (didn't contribute enough to the main narrative) and
# redesigned at the same time to fix the original layout: nesting one
# already-multi-panel gratia::draw() patchwork object inside a second
# patchwork::wrap_plots() call produced garbled, overlapping titles and
# illegibly small axis text -- the same nested-patchwork bug already fixed
# once for the X11 composite figures (see
# figure.R::save_x11_component_composite()) and fixed here the same way:
# render each zone's panel grid to its own PNG (5-column layout instead of
# 10-across, gratia::draw()'s redundant "Basis: Tensor product" caption
# suppressed, and a real title), then stack the four zone images with
# magick::image_append() rather than nesting patchwork objects.
plot_gam_figure <- function(gam_models, fig_path){
  if (!dir.exists(dirname(fig_path))) dir.create(dirname(fig_path), recursive = TRUE)
  tmp_dir <- tempfile("gam_zone_panels_")
  dir.create(tmp_dir)

  zone_files <- purrr::imap_chr(gam_models, function(m, zone_name){
    p <- gratia::draw(m, ncol = 5, caption = FALSE) +
      patchwork::plot_annotation(title = zone_title(zone_name),
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5)))
    out_file <- file.path(tmp_dir, paste0(zone_name, ".png"))
    ggplot2::ggsave(out_file, p, width = 20, height = 9, dpi = 200)
    out_file
  })

  composite <- magick::image_append(magick::image_read(zone_files), stack = TRUE)
  magick::image_write(composite, fig_path)
  unlink(tmp_dir, recursive = TRUE)
  invisible(composite)
}

# Partial-dependence curve for one driver from a fitted GAM: vary that
# driver over its observed range while holding every other driver at its
# median, predict plume_area from the model, and return the curve with a
# +/-2SE band (predict.gam(se.fit = TRUE)). Works with fit_gam() as-is
# (pairwise te() tensor smooths only, no univariate s() terms) since partial
# dependence is a post-hoc prediction technique -- it doesn't need a
# particular smooth-term structure, just a model to predict from. Used by
# func/figure.R::Figure_9_gam_partial() (manuscript Figure 9).
# gam_partial_effect(fit_gam(driver_matrices[["GULF_OF_LION"]]), "wind_spd", driver_matrices[["GULF_OF_LION"]])
gam_partial_effect <- function(gam_model, driver_name, df, n_points = 50){
  drivers <- available_drivers(df)
  df_valid <- tidyr::drop_na(df, dplyr::all_of(c("plume_area", drivers)))

  # Hold every driver but driver_name at a "typical" value: median for
  # numeric drivers, most-frequent level (as a plain string) for categorical
  # ones (median() errors on a factor). categorical drivers arrive here as
  # character, not factor (read straight from daily_driver_matrix_*.csv), so
  # levels(col) is NULL; factor(x, levels = NULL) creates a broken 0-level
  # factor -- not "use the default levels" the way omitting the argument
  # would -- which silently makes every predict.gam() call below return NA
  # (confirmed 2026-07-31: this was producing entirely blank Figure 9 panels
  # for the two zones/drivers with a categorical predictor in the model,
  # i.e. every zone, since wind_dir_cat/wave_dir_cat are always present).
  # Passing the modal value as a plain string instead lets predict.gam()
  # re-level it against the fitted model's own stored factor levels.
  newdata <- df_valid[1, drivers, drop = FALSE]
  for(d in drivers){
    col <- df_valid[[d]]
    newdata[[d]] <- if(is.numeric(col)) stats::median(col, na.rm = TRUE) else
      names(which.max(table(col)))
  }
  newdata <- newdata[rep(1, n_points), , drop = FALSE]
  newdata[[driver_name]] <- seq(min(df_valid[[driver_name]], na.rm = TRUE),
                                max(df_valid[[driver_name]], na.rm = TRUE), length.out = n_points)

  pred <- stats::predict(gam_model, newdata = newdata, se.fit = TRUE)
  tibble::tibble(driver_name = driver_name, x = newdata[[driver_name]], fit = pred$fit, se = pred$se.fit)
}


# Step 4: regime stratification ----------------------------------------------

add_regime_labels <- function(df){
  df |>
    dplyr::mutate(
      flow_z = as.numeric(scale(flow)),
      discharge_regime = ifelse(flow_z >= 0, "high (>=zone median)", "low (<zone median)"),
      wind_regime = ifelse(wind_spd >= 6, "high (>=6 m/s)", "low (<6 m/s)"),
      tide_bin = ifelse(tide_range >= stats::median(tide_range, na.rm = TRUE),
                        "spring (>=zone median)", "neap (<zone median)"),
      # NB: unlike wind_regime's literature-grounded 6 m/s threshold (Fofonova
      # et al. 2015), no equivalent fixed threshold for current speed or wave
      # height was found in the reviewed literature -- Poppeschi et al. (2024)
      # use a dynamic (climatological P90) threshold on wave orbital velocity,
      # not a fixed value on Hs, so these two regimes use the same zone-median
      # split as tide_bin rather than an invented absolute cutoff.
      current_regime = ifelse(current >= stats::median(current, na.rm = TRUE),
                              "high (>=zone median)", "low (<zone median)"),
      wave_regime = ifelse(wave_height >= stats::median(wave_height, na.rm = TRUE),
                           "high (>=zone median)", "low (<zone median)")
    )
}

refit_by_regime <- function(zone_name, regime_col, driver_matrices, response = "plume_area", min_n = 30){
  df <- add_regime_labels(driver_matrices[[zone_name]])
  drivers <- available_drivers(df)
  df_valid <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))
  if(nrow(df_valid) < min_n || stats::var(df_valid[[response]]) < 1e-6) return(NULL)
  df <- df_valid
  form <- stats::as.formula(paste(response, "~", paste(drivers, collapse = " + ")))

  purrr::map_dfr(split(df, df[[regime_col]]), function(sub){
    if(nrow(sub) < min_n) return(NULL)
    m <- stats::glm(form, data = sub, family = gaussian())
    broom::tidy(m) |> dplyr::mutate(n = nrow(sub), regime_value = unique(sub[[regime_col]]))
  }) |>
    dplyr::mutate(zone = zone_name, regime = regime_col, response = response, .before = 1)
}


# Step 5: per-metric models --------------------------------------------------

fit_metric_models <- function(df, response){
  list(
    glm_additive    = fit_baseline_glm(df, response),
    glm_interaction = fit_interaction_glm(df, response),
    gam             = fit_gam(df, response)
  )
}

metric_responses <- c("plume_area", "mean_SPM_in_the_plume_area", "mass_SPM_in_the_plume_area_in_g_m",
                      "lon_weighted_centroid_of_the_plume_area", "lat_weighted_centroid_of_the_plume_area")


# Step 6: exploratory random forest + H-statistic ----------------------------

fit_rf_diagnostic <- function(zone_name, driver_matrices, response = "plume_area"){
  df <- driver_matrices[[zone_name]]
  drivers <- available_drivers(df)
  df_complete <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))
  if(nrow(df_complete) < 30 || stats::var(df_complete[[response]]) < 1e-6){
    message("fit_rf_diagnostic: skipping ", zone_name, " (insufficient variation in ", response, ")")
    return(NULL)
  }

  rf <- ranger::ranger(stats::as.formula(paste(response, "~", paste(drivers, collapse = " + "))),
                       data = df_complete[, c(response, drivers)],
                       importance = "permutation", num.trees = 500)

  predictor <- iml::Predictor$new(rf, data = df_complete[, drivers], y = df_complete[[response]])
  interaction_hstat <- iml::Interaction$new(predictor)$results

  list(model = rf, importance = rf$variable.importance, interaction = interaction_hstat)
}


# Runner: execute all steps for one set of plume results ---------------------

run_full_analysis <- function(plume_dir, stats_dir, fig_path){

  if(!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

  message("[", Sys.time(), "] Step 0/6: building daily driver matrices for ", length(zones), " zones (", plume_dir, ")...")

  # Step 0: build driver matrices
  driver_matrices <- purrr::map(zones, ~ build_driver_matrix(.x, plume_dir = plume_dir)) |>
    purrr::set_names(zones)

  # Save daily combined tables
  purrr::iwalk(driver_matrices, function(df, zone){
    readr::write_csv(df, file.path(stats_dir, paste0("daily_driver_matrix_", zone, ".csv")))
  })

  message("[", Sys.time(), "] Step 1-2/6: fitting baseline + interaction GLMs and comparing via LRT/AIC...")

  # Step 2: GLM comparison
  glm_comparison_stats <- purrr::map(zones, compare_glms, driver_matrices = driver_matrices) |>
    purrr::compact() |> dplyr::bind_rows()
  readr::write_csv(glm_comparison_stats, file.path(stats_dir, "driver_glm_comparison.csv"))

  message("[", Sys.time(), "] Step 3/6: fitting GAMs with tensor-product smooths and drawing Figure 10...")

  # Step 3: GAM
  gam_models <- purrr::map(zones, ~ fit_gam(driver_matrices[[.x]])) |>
    purrr::set_names(zones) |> purrr::compact()

  gam_summary <- purrr::imap_dfr(gam_models, function(m, zone_name){
    s <- summary(m)
    tibble::tibble(zone = zone_name, term = rownames(s$s.table),
                   edf = s$s.table[, "edf"], ref.df = s$s.table[, "Ref.df"],
                   F = s$s.table[, "F"], p = s$s.table[, "p-value"],
                   r_sq_adj = s$r.sq, deviance_explained = s$dev.expl)
  })
  readr::write_csv(gam_summary, file.path(stats_dir, "driver_gam_summary.csv"))
  if(length(gam_models) > 0) plot_gam_figure(gam_models, fig_path)

  message("[", Sys.time(), "] Step 4/6: refitting GLMs by discharge/wind/tide regime...")

  # Step 4: regime stratification
  regime_stats <- purrr::map_dfr(zones, function(z){
    dplyr::bind_rows(
      refit_by_regime(z, "discharge_regime", driver_matrices),
      refit_by_regime(z, "wind_regime", driver_matrices),
      refit_by_regime(z, "tide_bin", driver_matrices),
      refit_by_regime(z, "current_regime", driver_matrices),
      refit_by_regime(z, "wave_regime", driver_matrices)
    )
  })
  readr::write_csv(regime_stats, file.path(stats_dir, "driver_regime_glm.csv"))

  message("[", Sys.time(), "] Step 5/6: fitting per-metric models for ", length(metric_responses), " response variables...")

  # Step 5: per-metric models
  metric_model_stats <- purrr::map(metric_responses, function(resp){
    purrr::imap_dfr(driver_matrices, function(df, zone_name){
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
    readr::write_csv(stats_df, file.path(stats_dir, paste0("driver_metric_models_", resp, ".csv")))
  })

  message("[", Sys.time(), "] Step 6/6: fitting random forest + H-statistic interaction diagnostics...")

  # Step 6: random forest
  rf_results <- purrr::map(zones, fit_rf_diagnostic, driver_matrices = driver_matrices) |>
    purrr::set_names(zones) |> purrr::compact()

  rf_importance <- purrr::imap_dfr(rf_results, function(res, zone_name){
    tibble::tibble(zone = zone_name, driver = names(res$importance), importance = res$importance)
  })
  readr::write_csv(rf_importance, file.path(stats_dir, "driver_rf_importance.csv"))

  rf_interaction <- purrr::imap_dfr(rf_results, function(res, zone_name){
    dplyr::mutate(res$interaction, zone = zone_name, .before = 1)
  })
  if(nrow(rf_interaction) > 0){
    readr::write_csv(rf_interaction, file.path(stats_dir, "driver_rf_interaction_hstat.csv"))
  }

  message("run_full_analysis() complete. Outputs written to ", stats_dir, " and ", fig_path,
          " at ", Sys.time())
}


# Entry point: called from code/4_time_series.py via rpy2 -------------------
# Runs both the dynamic-threshold (main results) and static-threshold
# (supplementary) driver analyses, mirroring the two run_full_analysis()
# calls previously made directly at the bottom of code/6_driver_interactions.R
# (now deleted -- this file replaces it, called in-sequence from Stage 4
# rather than as a separate, unnumbered Stage 6 script).

run_driver_interactions_analysis <- function(){

  # 2026-07-31: Figure 10 moved from main-text to Supplementary Figure S7
  # (Robert: didn't contribute enough to the main narrative), redesigned in
  # the same pass -- see plot_gam_figure()'s note.
  message("== Driver interactions: dynamic threshold (main results) ==")
  run_full_analysis(
    plume_dir = "output/panache/dynamic",
    stats_dir = "output/STATS",
    fig_path  = "figures/ARTICLE/FIGURE_S7/Figure_S7.png"
  )

  message("== Driver interactions: static threshold (supplementary) ==")
  run_full_analysis(
    plume_dir = "output/panache/static",
    stats_dir = "output/STATS/static",
    fig_path  = "figures/ARTICLE/FIGURE_S7/Figure_S7_static.png"
  )

  message("func/driver_interactions.R::run_driver_interactions_analysis() complete.")
  invisible(TRUE)
}
