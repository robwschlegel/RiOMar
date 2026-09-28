# func/driver_interactions.R
# Multi-driver interaction analysis for river plume size.

# Called from code/4_time_series.py via rpy2, the same way func/X11.R,
# func/figure.R, func/validate.R, and func/plume.R are called from their
# respective func/*.py modules. 
# source() this file, then call run_driver_interactions_analysis(). 
# This keeps figure/stat generation centralised behind the numbered 
# code/ pipeline (per CLAUDE.md) rather than split across a separate, 
# unnumbered driver-interactions script.

# What this script does
#   1. Baseline additive GLM
#   2. + pairwise interaction terms, LRT/AIC
#   3. Per-metric models: additive GLM, interaction GLM, and GAM with te()
#      tensor smooths, for each of 5 response metrics (area, mean/mass SPM,
#      centroid lon/lat)
#   4. Exploratory random forest importance
#
# Removed 2026-09-28 (none had a manuscript consumer left): a zone-level GAM
# fit for driver_gam_summary.csv + its figures (plot_gam_figure()), the
# discharge/wind/tide/current/wave regime-stratified GLMs
# (add_regime_labels()/refit_by_regime(), driver_regime_glm.csv), the random
# forest H-statistic diagnostic, and the per-month repeat of the whole
# sequence (run_monthly_driver_interactions_analysis() and
# summarise_monthly_driver_dominance()). The CSVs they last wrote are still
# on disk and cited as frozen results in manuscript.tex; the code is in git
# history (last present in commit 533e42c) and summarised in
# manuscript/reviewer_responses.md.

# Known simplifications / gaps, so nothing here is mistaken for more
# rigorous than it is:
#   - "Centroid" per-metric models use the SPM-weighted mean pixel lon/lat
#     of the daily plume mask (lon/lat_weighted_centroid_of_the_plume_area),
#     not an unweighted mean. Ralston et al. (2024)'s alongshore/cross-shore
#     split is still not reproduced here, since that needs a coastline-
#     following coordinate rotation per zone this project doesn't yet have.


# Setup ---------------------------------------------------------------------
# Like every other func/*.R module called via rpy2 (see CLAUDE.md), this
# assumes the R session's working directory is the repo root.

if(!dir.exists("func")) stop("func/driver_interactions.R must be sourced with the repo root as the working directory.")

library(tidyverse)
library(mgcv)       # fit_gam(): used by the per-metric models step and by
                    # func/figure.R's Figure S8 (gam_partial_effects)
library(ranger)

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

  df <- df_plume |>
    dplyr::left_join(df_flow, by = "date") |>
    dplyr::left_join(df_tide, by = "date") |>
    dplyr::left_join(df_wind, by = "date") |>
    dplyr::left_join(df_current, by = "date") |>
    dplyr::left_join(df_wave, by = "date")

  df |>
    dplyr::mutate(zone = zone_name, .before = "date",
                  # Categorical (8-octant) form of direction
                  wind_dir_cat = compass_octant(wind_dir),
                  wave_dir_cat = compass_octant(wave_dir),
                  current_dir_cat = compass_octant(current_dir)) |>
    zoo::na.trim()
}

# Names of the driver columns available for each zone.
# Any candidate column that is entirely NA for this zone
# is dropped rather than passed into a model formula.
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
  message("  [", Sys.time(), "] ", zone_name, ": fitting baseline + interaction GLMs...")
  m0 <- fit_baseline_glm(df_valid, response)
  m1 <- fit_interaction_glm(df_valid, response)
  lrt <- stats::anova(m0, m1, test = "Chisq")
  tibble::tibble(zone = zone_name, response = response,
                 aic_additive = stats::AIC(m0), aic_interaction = stats::AIC(m1),
                 deviance_additive = m0$deviance, deviance_interaction = m1$deviance,
                 lrt_chisq = lrt$Deviance[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chi)`[2])
}


# GAM with tensor-product smooths (used by step 3 and figure.R) ------------------------------------

fit_gam <- function(df, response = "plume_area"){
  drivers <- available_drivers(df)
  df_valid <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))
  if(nrow(df_valid) < 30 || stats::var(df_valid[[response]]) < 1e-6) return(NULL)

  # te() tensor-product smooths need numeric arguments.
  # Categorical drivers(wind_dir_cat/wave_dir_cat) enter as flat parametric 
  # main-effect terms instead, alongside the te() smooths over every numeric pair.
  is_categorical <- purrr::map_lgl(drivers, ~ !is.numeric(df_valid[[.x]]))
  numeric_drivers <- drivers[!is_categorical]
  categorical_drivers <- drivers[is_categorical]

  pair_terms <- utils::combn(numeric_drivers, 2, simplify = FALSE)
  te_terms <- purrr::map_chr(pair_terms, ~ paste0("te(", .x[1], ", ", .x[2], ")"))
  form <- stats::as.formula(paste(response, "~", paste(c(te_terms, categorical_drivers), collapse = " + ")))
  mgcv::gam(form, data = df_valid, method = "REML")
}

# Partial-dependence curve for one driver from a fitted GAM: vary that
# driver over its observed range while holding every other driver at its
# median, predict plume_area from the model, and return the curve with a +/-2SE band. 
# Works with fit_gam() as-is (pairwise te() tensor smooths only, no univariate s() terms)
# since partial dependence is a post-hoc prediction technique.
# It doesn't need a particular smooth-term structure, just a model to predict from. 
# Used by func/figure.R::plot_gam_partial_effects().
# gam_partial_effect(fit_gam(driver_matrices[["GULF_OF_LION"]]), "wind_spd", driver_matrices[["GULF_OF_LION"]])
gam_partial_effect <- function(gam_model, driver_name, df, n_points = 50){
  drivers <- available_drivers(df)
  df_valid <- tidyr::drop_na(df, dplyr::all_of(c("plume_area", drivers)))

  # Hold every driver but driver_name at a "typical" value: median for
  # numeric drivers, most-frequent level (as a plain string) for categorical
  # ones (median() errors on a factor).
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


# Step 3: per-metric models --------------------------------------------------

fit_metric_models <- function(df, response){
  list(
    glm_additive    = fit_baseline_glm(df, response),
    glm_interaction = fit_interaction_glm(df, response),
    gam             = fit_gam(df, response)
  )
}

metric_responses <- c("plume_area", "mean_SPM_in_the_plume_area", "mass_SPM_in_the_plume_area_in_t",
                      "lon_weighted_centroid_of_the_plume_area", "lat_weighted_centroid_of_the_plume_area")


# Step 4: exploratory random forest importance --------------------------------
# The H-statistic interaction diagnostic that used to sit here (one extra
# ranger fit + iml::Interaction$new() per zone) was removed 2026-09-28: its
# output (driver_rf_interaction_hstat.csv) had no manuscript consumer and was
# the single slowest part of this step.

fit_rf_diagnostic <- function(zone_name, driver_matrices, response = "plume_area", n_repeats = 10,
                              num_threads = NULL){
  df <- driver_matrices[[zone_name]]
  drivers <- available_drivers(df)
  df_complete <- tidyr::drop_na(df, dplyr::all_of(c(response, drivers)))
  if(nrow(df_complete) < 30 || stats::var(df_complete[[response]]) < 1e-6){
    message("fit_rf_diagnostic: skipping ", zone_name, " (insufficient variation in ", response, ")")
    return(NULL)
  }

  rf_formula <- stats::as.formula(paste(response, "~", paste(drivers, collapse = " + ")))

  message("  [", Sys.time(), "] ", zone_name, ": fitting ", n_repeats, " ranger repeats on ",
          nrow(df_complete), " rows...")
  repeats_start <- Sys.time()

  # Permutation importance is known to redistribute unpredictably among
  # correlated predictors and vary between repeated fits of the same forest
  # (Strobl et al. 2007; Nicodemus et al. 2010; Wang et al. 2016), worse the
  # smaller the sample -- refit n_repeats times with different seeds and
  # report the mean and SD across repeats, rather than trusting a single fit,
  # so a driver whose ranking is unstable shows a large importance_sd instead
  # of silently looking as solid as a stable one.
  importance_repeats <- purrr::map(seq_len(n_repeats), function(i){
    ranger::ranger(rf_formula, data = df_complete[, c(response, drivers)],
                   importance = "permutation", num.trees = 500, seed = i,
                   num.threads = num_threads)$variable.importance
  })
  importance_mat <- do.call(rbind, importance_repeats)
  importance_mean <- colMeans(importance_mat)
  importance_sd <- apply(importance_mat, 2, stats::sd)

  message("  [", Sys.time(), "] ", zone_name, ": ranger repeats done (",
          round(difftime(Sys.time(), repeats_start, units = "secs"), 1), "s)")

  list(importance = importance_mean, importance_sd = importance_sd)
}


# Runner: execute all steps for one set of plume results ---------------------

run_full_analysis <- function(plume_dir, stats_dir, rf_num_threads = NULL, overwrite = TRUE){

  if(!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

  message("[", Sys.time(), "] Step 0/4: building daily driver matrices for ", length(zones), " zones (", plume_dir, ")...")

  # Step 0: build driver matrices
  driver_matrices <- purrr::map(zones, function(zone_name){
    message("  [", Sys.time(), "] ", zone_name, ": building driver matrix...")
    build_driver_matrix(zone_name, plume_dir = plume_dir)
  }) |> purrr::set_names(zones)

  # Save daily combined tables
  purrr::iwalk(driver_matrices, function(df, zone){
    readr::write_csv(df, file.path(stats_dir, paste0("daily_driver_matrix_", zone, ".csv")))
  })

  message("[", Sys.time(), "] Step 1-2/4: fitting baseline + interaction GLMs and comparing via LRT/AIC...")

  # Step 2: GLM comparison
  glm_comparison_stats <- purrr::map(zones, compare_glms, driver_matrices = driver_matrices) |>
    purrr::compact() |> dplyr::bind_rows()
  readr::write_csv(glm_comparison_stats, file.path(stats_dir, "driver_glm_comparison.csv"))

  # (The zone-level GAM and regime-GLM steps that used to sit here were
  # removed 2026-09-28 -- see the header note. Their GAM R^2/deviance
  # numbers were redundant with the plume_area cell below, which is what
  # paragraph_source_registry.csv cites.)

  message("[", Sys.time(), "] Step 3/4: fitting per-metric models for ", length(metric_responses), " response variables...")

  # Step 3: per-metric models. Each (response, zone) cell is checkpointed to
  # its own file under metric_cache_dir as soon as it's fitted -- this is
  # what makes overwrite = FALSE useful: a rerun after an interrupted run
  # reads any cell that's already on disk instead of refitting it, while
  # overwrite = TRUE (the default) always refits and just refreshes the
  # cache file for next time.
  metric_cache_dir <- file.path(stats_dir, ".checkpoints", "metric_models")
  dir.create(metric_cache_dir, recursive = TRUE, showWarnings = FALSE)

  metric_model_stats <- purrr::map(metric_responses, function(resp){
    purrr::imap_dfr(driver_matrices, function(df, zone_name){
      cache_path <- file.path(metric_cache_dir, paste0(resp, "__", zone_name, ".csv"))
      if(!overwrite && file.exists(cache_path)){
        message("  [", Sys.time(), "] ", resp, " / ", zone_name, ": using cached result (overwrite = FALSE)")
        return(readr::read_csv(cache_path, show_col_types = FALSE))
      }
      if(!(resp %in% names(df))){
        message("  [", Sys.time(), "] ", resp, " / ", zone_name, ": skipped (response not present)")
        return(NULL)
      }
      df_resp <- tidyr::drop_na(df, dplyr::all_of(c(resp, available_drivers(df))))
      if(nrow(df_resp) < 30){
        message("  [", Sys.time(), "] ", resp, " / ", zone_name, ": skipped (< 30 complete rows)")
        return(NULL)
      }
      message("  [", Sys.time(), "] ", resp, " / ", zone_name, ": fitting GLM + GAM on ", nrow(df_resp), " rows...")
      fit_start <- Sys.time()
      models <- fit_metric_models(df_resp, resp)
      message("  [", Sys.time(), "] ", resp, " / ", zone_name, ": done (",
              round(difftime(Sys.time(), fit_start, units = "secs"), 1), "s)")
      result <- tibble::tibble(zone = zone_name, response = resp,
                     aic_additive = stats::AIC(models$glm_additive),
                     aic_interaction = stats::AIC(models$glm_interaction),
                     gam_r_sq_adj = summary(models$gam)$r.sq,
                     gam_deviance_explained = summary(models$gam)$dev.expl)
      readr::write_csv(result, cache_path)
      result
    })
  }) |> purrr::set_names(metric_responses)

  purrr::iwalk(metric_model_stats, function(stats_df, resp){
    readr::write_csv(stats_df, file.path(stats_dir, paste0("driver_metric_models_", resp, ".csv")))
  })

  message("[", Sys.time(), "] Step 4/4: fitting random forest importance diagnostics...")

  # Random forest importance. Same per-cell checkpointing as step 3, one file
  # per zone.
  rf_cache_dir <- file.path(stats_dir, ".checkpoints", "rf")
  dir.create(rf_cache_dir, recursive = TRUE, showWarnings = FALSE)

  rf_results <- purrr::map(zones, function(zone_name){
    cache_path <- file.path(rf_cache_dir, paste0(zone_name, ".rds"))
    if(!overwrite && file.exists(cache_path)){
      message("  [", Sys.time(), "] ", zone_name, ": using cached RF result (overwrite = FALSE)")
      return(readRDS(cache_path))
    }
    result <- fit_rf_diagnostic(zone_name, driver_matrices, num_threads = rf_num_threads)
    if(!is.null(result)) saveRDS(result, cache_path)
    result
  }) |> purrr::set_names(zones) |> purrr::compact()

  rf_importance <- purrr::imap_dfr(rf_results, function(res, zone_name){
    tibble::tibble(zone = zone_name, driver = names(res$importance), importance = res$importance,
                   importance_sd = res$importance_sd[names(res$importance)])
  })
  readr::write_csv(rf_importance, file.path(stats_dir, "driver_rf_importance.csv"))

  message("run_full_analysis() complete. Outputs written to ", stats_dir, " at ", Sys.time())
}


# Entry point: called from code/4_time_series.py via Rscript subprocess -----
# Runs the dynamic-threshold (main results) driver analysis. The
# static-threshold (supplementary) pass was dropped (2026-09-24): its GAM
# fitting (step 5/6) ran an order of magnitude slower than the dynamic
# pass for reasons not worth chasing down, and the static results aren't
# individually manuscript-referenced.

run_driver_interactions_analysis <- function(overwrite = TRUE){

  message("== Driver interactions: dynamic threshold (main results) ==")
  run_full_analysis(
    plume_dir = "output/panache/dynamic",
    stats_dir = "output/STATS",
    overwrite = overwrite
  )

  message("func/driver_interactions.R::run_driver_interactions_analysis() complete.")
  invisible(TRUE)
}
