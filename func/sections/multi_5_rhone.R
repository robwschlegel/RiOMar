# Rhone only analyses ------------------------------------------------------

# Response to Claude's e-mail on the Grand Rhone plume-area
# trend. The e-mail raised four ideas; disposition of each below:
#   1. Concentration-detrending sensitivity test -- IMPLEMENTED,
#      rhone_detrend_test().
#   2. Has the balance of northern-tributary vs Cevennes ("cevenol")
#      floods shifted, explaining the concentration rise? -- PARTIALLY
#      implemented, rhone_flood_timing_shift(). RiOMar only has
#      combined Grand+Petit Rhone discharge at the mouth, not upstream
#      sub-basin gauges, so this cannot attribute the shift to a specific
#      tributary -- that needs the Observatoire des Sediments du Rhone
#      data / Olivier Radakovitch the e-mail points to. What's implemented
#      is a first-pass check of whether the *calendar timing* of high-flow
#      days has moved, using data already in the repo.
#   3. Does wind/wave forcing move the plume beyond what flow explains
#      (reframing "does the plume have its own seasonality" as "does it
#      respond to wind/wave directly"), including the specific Mistral vs
#      onshore/calm hypothesis -- IMPLEMENTED, rhone_wind_wave_effect().
#   4. Would a more stratified surface layer let the plume slide further
#      offshore? -- NOT implemented. The e-mail itself frames this as a
#      question for idealised numerical modelling of academic cases, not
#      something to test against the RiOMar time series, and RiOMar has no
#      surface density/stratification product to test it against even if
#      it were. Left as a discussion point for the reply, not a script.
#
# Her reply to that first round raised two more, smaller follow-ups, handled
# below as their own functions rather than folded into the numbered list
# above: does the wind effect (item 3) show a winter-vs-summer nuance
# (rhone_wind_wave_seasonal_effect()), and has the flow trend itself moved
# differently across the flow distribution rather than just at the mean
# (rhone_flow_quantile_trend()).
#
# All functions below are GULF_OF_LION / Grand Rhone only, they are not
# generalised to the other three zones.


## 1. Concentration-detrending sensitivity test -------------------------------
# The plume-area trend could be, at least partly, an artefact of the SPM
# concentration near the Rhone mouth having risen over the record (see
# output/STATS/mass_SPM_trend_summary.csv and the earlier e-mail's claim)
# rather than the plume physically growing. The e-mail's proposed check:
# imagine the concentration had NO trend (a "trend-free" counterfactual,
# built "at least as a first approximation" by flattening the trend line),
# recompute the plume area under that counterfactual, and see whether the
# area trend survives.
#
# Re-running the full pixel-level panache detection under a hypothetical,
# de-trended SPM field is out of scope here -- it needs the raw daily SPM
# maps (not exposed to R) and a re-run of the Python flood-fill algorithm
# (func/plume.py::find_SPM_threshold()) for every day. An earlier version of
# this function used panache's own daily SPM_threshold (Results.csv) as a
# stand-in for concentration, but that is panache's classification cutoff
# for that day (func/plume.py:1254, `SPM_criterion = ds_reduced >
# SPM_threshold`) -- a stricter cutoff mechanically shrinks the classified
# area regardless of any real dilution physics, which confounds exactly the
# effect this test is trying to isolate. Instead, this uses real, independent
# in-situ SPM concentration measured at Arles on the Rhone itself
# (data/INSITU_data/OSR/ARLES_CMES-2.txt, hourly, mg/L, 2005-2023), from the
# Observatoire des Sediments du Rhone (OSR) -- the exact data source Claude's
# e-mail pointed to. Only quality flag "v" (validated/measured) rows are
# kept; flag "e" rows are the OSR's own gap-fill, estimated from a fixed
# discharge-to-concentration rating curve (see
# data/INSITU_data/OSR/Report_ARLES_CMES-2.txt), which would reintroduce a
# collinearity with flow that using real in-situ data is meant to avoid.
# Because the OSR record only covers 2005-2023, the joined record below is
# shorter than the full 1998-2025 satellite plume record.
# rhone_detrend_test()
rhone_detrend_test <- function(){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")

  # Daily-mean in-situ SPM concentration at Arles, validated readings only.
  df_conc <- read_delim("data/INSITU_data/OSR/ARLES_CMES-2.txt", delim = ";", skip = 3,
                        col_names = c("datetime", "value", "quality", "min", "max"),
                        col_types = "cdccc", show_col_types = FALSE) |>
    dplyr::filter(quality == "v") |>
    dplyr::mutate(date = as.Date(datetime, format = "%d/%m/%Y %H:%M:%S")) |>
    dplyr::summarise(conc = mean(value, na.rm = TRUE), .by = "date")

  # Standard plume-area + flow object used by every other driver in this
  # file. inner_join (not left_join) restricts the record to where the
  # satellite plume record and the OSR in-situ record overlap.
  df_flow <- combine_plume_driver("flow", meta)
  df <- df_flow |>
    dplyr::inner_join(df_conc, by = "date") |>
    zoo::na.trim()  # drops leading/trailing NA across plume_area, flow, and conc together

  # Interior gaps (e.g. days with no validated hourly readings) are linearly
  # interpolated so every remaining day has a value -- the same treatment
  # trend fitting already gets elsewhere in this file (ar_weights_func()).
  df$conc <- as.numeric(zoo::na.approx(df$conc, x = df$date, na.rm = FALSE))

  # a) Linear trend of the concentration proxy, and its "trend-free"
  #    counterfactual: the same series with the fitted slope subtracted
  #    back out (mean unchanged, so it stays on the original scale).
  conc_lm <- lm(conc ~ date, data = df)
  df$conc_fit <- predict(conc_lm)
  df$conc_detrended <- df$conc - (df$conc_fit - mean(df$conc_fit))

  # b) Empirical area ~ concentration sensitivity, controlling for flow so
  #    the concentration coefficient isn't just re-capturing flow's own,
  #    already-established effect on area (see driver_plume_trend("flow")).
  #    Caveat worth keeping in mind when reading beta_conc: flow and
  #    concentration are themselves correlated (bigger floods carry more
  #    sediment), so this is a partial, not fully independent, effect.
  area_conc_lm <- lm(plume_area ~ value + conc, data = df)
  beta_conc <- coef(area_conc_lm)["conc"]

  # c) Counterfactual area: subtract, day by day, the modelled contribution
  #    of that day's concentration sitting above (or below) its trend-free
  #    counterfactual.
  df$plume_area_conc_adj <- df$plume_area - beta_conc * (df$conc - df$conc_detrended)

  # d) Refit the trend on the raw area, the concentration-adjusted area,
  #    and the concentration proxy itself, all with the same AR-weighted /
  #    HAC estimator used throughout this file (fit_wls_hac_trend()), so
  #    the comparison is apples-to-apples with every other trend reported
  #    for this zone.
  trend_conc <- fit_wls_hac_trend("ar", df$conc, df$date) |> dplyr::mutate(series = "SPM concentration (Arles, OSR in-situ)")
  trend_raw  <- fit_wls_hac_trend("ar", df$plume_area, df$date) |> dplyr::mutate(series = "plume area (raw)")
  trend_adj  <- fit_wls_hac_trend("ar", df$plume_area_conc_adj, df$date) |> dplyr::mutate(series = "plume area (concentration-adjusted)")

  stats <- dplyr::bind_rows(trend_conc, trend_raw, trend_adj) |>
    dplyr::mutate(slope_yearly = slope * 365.25, beta_conc = beta_conc) |>
    dplyr::select(series, n, slope, slope_yearly, slope_se, slope_p, beta_conc)
  write_csv(stats, "output/STATS/rhone_conc_detrend.csv")

  # Plot: a) in-situ concentration with its fitted trend and trend-free
  # counterfactual; b) raw vs. concentration-adjusted plume area with their
  # respective AR-weighted trend lines, so the effect of the adjustment
  # (if any) is visible directly.
  # NB: conc is a high-variance daily series (floods vs. base flow), so
  # connecting every day with a line (geom_line) paints a near-solid band
  # that hides everything underneath it. Points + a clean geom_abline for
  # the fitted trend stays legible; the trend-free counterfactual is, at
  # this scale, visually indistinguishable from the raw series (the whole
  # record's drift is small next to the daily noise), so it is reported in
  # `stats` rather than over-plotted here.
  # SPM concentration is strongly right-skewed (flood spikes into the
  # thousands of mg/L vs. a base-flow median around 14 mg/L), so a linear
  # y-axis crushes almost every point flat against zero -- log10 is needed
  # to actually see the day-to-day/seasonal structure. The fitted trend
  # line itself is still the linear-scale fit used for the detrending math
  # above (conc_lm); log10() only changes how it is displayed here.
  # NB: geom_abline() does NOT work here -- it applies intercept/slope in
  # the already-log10-transformed coordinate space, not the original mg/L
  # space conc_lm was fit in, so the line would be drawn many orders of
  # magnitude off the visible panel. Instead, get the two endpoint values
  # by predicting conc_lm (linear space) at the start/end dates, then plot
  # those as ordinary data with geom_line() so they go through the same
  # log10 transform as the points.
  conc_trend_line <- tibble::tibble(date = range(df$date)) |>
    dplyr::mutate(conc_fit_line = predict(conc_lm, newdata = tibble::tibble(date = date)))

  pl_conc <- ggplot(df, aes(x = date)) +
    geom_point(aes(y = conc), colour = "grey50", alpha = 0.15, size = 0.5) +
    geom_line(data = conc_trend_line, aes(y = conc_fit_line), colour = "firebrick", linewidth = 1.2) +
    scale_y_log10() +
    labs(x = NULL, y = "SPM concentration, Arles (OSR in-situ, mg/L, log scale)",
         title = "Grand Rhone: in-situ SPM concentration at Arles",
         subtitle = "Red = fitted linear trend (see `stats` for the trend-free counterfactual comparison)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  fit_area_raw <- coef(lm(plume_area ~ date, data = df))
  fit_area_adj <- coef(lm(plume_area_conc_adj ~ date, data = df))

  pl_area <- ggplot(df, aes(x = date)) +
    geom_point(aes(y = plume_area), colour = "sienna", alpha = 0.15) +
    geom_point(aes(y = plume_area_conc_adj), colour = "darkblue", alpha = 0.15) +
    geom_abline(intercept = fit_area_raw[1], slope = fit_area_raw[2], colour = "sienna", linewidth = 1.2) +
    geom_abline(intercept = fit_area_adj[1], slope = fit_area_adj[2], colour = "darkblue", linewidth = 1.2) +
    labs(x = NULL, y = "Plume area (km²)",
         title = "Grand Rhone: plume area, raw vs. concentration-adjusted",
         subtitle = "Brown = raw area (+ OLS trend); blue = area with the concentration-trend contribution removed (+ OLS trend)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  pl_combi <- ggpubr::ggarrange(pl_conc, pl_area, ncol = 1, nrow = 2)
  ggsave(filename = "figures/rhone_side_analyses/rhone_detrend_test.png", plot = pl_combi, width = 12, height = 10)

  return(list(stats = stats, data = df, plot = pl_combi))
}


## 2. Has the seasonal timing of Rhone floods shifted? ------------------------
# The e-mail's second idea asks *why* Rhone SPM concentration would be
# rising: has the balance between northern-tributary floods (winter/spring,
# snowmelt-and-rain driven) and Cevennes floods ("episodes cevenols",
# roughly Sep-Nov) shifted? A proper attribution needs the sub-basin gauge
# data the e-mail points to (Observatoire des Sediments du Rhone / Olivier
# Radakovitch) -- RiOMar only has the combined Grand+Petit Rhone discharge
# at the mouth (data/RIVER_FLOW/GULF_OF_LION), so this cannot say which
# tributary is responsible. What follows is a first-pass check of whether
# the *calendar timing* of high-flow days at the combined gauge has moved
# over the record -- worth reporting even though it can't yet attribute a
# cause.
# rhone_flood_timing_shift()
rhone_flood_timing_shift <- function(high_flow_quantile = 0.90){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")
  df_flow <- load_driver("flow", meta) |>
    dplyr::mutate(year = year(date), month = month(date), doy = yday(date))

  # Meteorological seasons (DJF/MAM/JJA/SON), extended from the Cevenol
  # (SON) focus to all four so winter/spring/summer can be compared on the
  # same footing. December is relabelled into the *following* year
  # (season_year = year + 1) so winter isn't split across two different
  # true winters -- e.g. Dec 2010 + Jan/Feb 2011 are one winter, both
  # tagged season_year 2011. For the other three seasons season_year is
  # just the calendar year. This has to happen *before* the high-flow
  # threshold below, because that threshold is also grouped by season_year
  # (see next comment) rather than raw calendar year.
  df_flow <- df_flow |>
    dplyr::mutate(season = dplyr::case_when(
                    month %in% 3:5  ~ "spring (MAM)",
                    month %in% 6:8  ~ "summer (JJA)",
                    month %in% 9:11 ~ "autumn (SON, Cevenol)",
                    TRUE            ~ "winter (DJF)"),
                  season_year = ifelse(month == 12, year + 1, year))
  season_levels <- c("winter (DJF)", "spring (MAM)", "summer (JJA)", "autumn (SON, Cevenol)")
  df_flow$season <- factor(df_flow$season, levels = season_levels)

  # "High flow" is defined relative to each season_year's own range (that
  # season_year's 90th percentile), not a single record-wide threshold -- a
  # fixed absolute threshold would mechanically flag more/fewer days in
  # years where mean flow itself has trended, which is exactly the
  # ambiguity this check is trying to avoid. Grouped by season_year, not
  # raw calendar year, so a single true winter always gets judged against
  # one consistent threshold -- grouping by calendar year instead would
  # let December and the following January/February of the same winter be
  # compared against two different years' flow distributions.
  df_flow <- df_flow |>
    dplyr::mutate(year_thresh = quantile(value, high_flow_quantile, na.rm = TRUE), .by = "season_year") |>
    dplyr::mutate(is_high_flow = value >= year_thresh)

  # Per (calendar) year: how many high-flow days, and their mean day-of-year
  # (NB: a circular statistic tied to the calendar year specifically,
  # because day-of-year itself resets at Jan 1 -- averaging doy within a
  # Dec-Feb season_year group would run straight across that reset and
  # produce a meaningless mid-year value, so this stays on raw "year").
  df_year <- df_flow |>
    dplyr::filter(is_high_flow) |>
    dplyr::summarise(n_high_flow = dplyr::n(), mean_doy_high_flow = mean(doy), .by = "year")

  # Per season-year: each season's share of that season-year's high-flow
  # days -- the four seasons' shares sum to ~1 within a season_year, so
  # this is the direct four-season generalisation of the single prop_autumn
  # check from the first pass. tidyr::complete() adds an explicit 0 (not a
  # missing row) for any season_year x season combination that had no
  # high-flow days at all that year -- "this season had none of the year's
  # extreme days" is real information (and changes the trend below
  # noticeably), not something to silently drop from the regression.
  df_season_year <- df_flow |>
    dplyr::filter(is_high_flow) |>
    dplyr::summarise(n_season = dplyr::n(), .by = c("season_year", "season")) |>
    tidyr::complete(season_year, season, fill = list(n_season = 0)) |>
    dplyr::mutate(n_total = sum(n_season), .by = "season_year") |>
    dplyr::mutate(prop_season = n_season / n_total)

  # Trend of each season's share over the years, and of the mean
  # day-of-year. One value per year/season-year here (not the
  # daily/monthly autocorrelated case fit_wls_hac_trend() is built for), so
  # a plain OLS trend is enough.
  season_trend <- function(season_name){
    d <- dplyr::filter(df_season_year, season == season_name)
    fit <- summary(lm(prop_season ~ season_year, data = d))$coefficients
    tibble::tibble(metric = paste0("proportion_of_high_flow_days_in_", season_name),
                   slope_per_year = fit["season_year", "Estimate"], p_value = fit["season_year", "Pr(>|t|)"])
  }
  stats_season <- purrr::map_dfr(season_levels, season_trend)

  fit_doy <- summary(lm(mean_doy_high_flow ~ year, data = df_year))$coefficients
  stats_doy <- tibble::tibble(metric = "mean_doy_of_high_flow_days",
                              slope_per_year = fit_doy["year", "Estimate"], p_value = fit_doy["year", "Pr(>|t|)"])

  stats <- dplyr::bind_rows(stats_doy, stats_season)
  write_csv(stats, "output/STATS/rhone_flood_seasonality.csv")

  pl <- ggplot(df_season_year, aes(x = season_year, y = prop_season)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
    facet_wrap(~season, ncol = 2) +
    labs(x = NULL, y = paste0("Proportion of top ", round((1 - high_flow_quantile) * 100), "% flow days falling in that season"),
         title = "Grand Rhone: has the seasonal timing of high-flow days shifted?",
         subtitle = "A rising trend means that season is becoming relatively more prominent for floods (Cevenol = autumn/SON panel)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  ggsave(filename = "figures/rhone_side_analyses/rhone_flood_timing_shift.png", plot = pl, width = 10, height = 8)

  return(list(stats = stats, data = df_season_year, plot = pl))
}


## 3. Does the plume respond to wind/wave forcing beyond flow? ---------------
# The e-mail's third idea reframes "does plume size have a seasonal cycle
# beyond river flow" as "does plume size respond directly to wind/wave
# forcing" -- the hypothesis being that winter's stronger wind/waves push
# the plume further offshore (larger detected area) for the same
# discharge, and that THIS, not a genuine seasonal cycle in the plume's own
# dynamics, is what a naive flow-vs-season comparison would pick up. It
# also raises a specific, testable version of that: Mistral (a strong,
# roughly N/NNW wind for the Gulf of Lion) should push the plume clear of
# the coast (larger detected area); weak or onshore/easterly wind should
# leave it hugging the coast, where it may be under-detected because of
# confounding coastal resuspension.
# rhone_wind_wave_effect()
rhone_wind_wave_effect <- function(){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")
  df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, wind_spd = value, wind_dir)
  df_wave <- load_driver("wave", meta) |> dplyr::select(date, wave_height = value)

  df <- df_flow |>
    dplyr::left_join(df_wind, by = "date") |>
    dplyr::left_join(df_wave, by = "date") |>
    tidyr::drop_na(plume_area, flow, wind_spd, wave_height)  # complete cases only, so the nested models below are fit to identical rows

  # a) Does wind/wave explain plume-area variance beyond flow alone? Nested
  #    model comparison (F-test via anova()) -- the same "does adding this
  #    term help" logic as code/6_driver_interactions.R::fit_interaction_glm(),
  #    scoped here to a plain additive comparison since the question is
  #    "does it matter at all", not "how do drivers interact".
  m_flow      <- lm(plume_area ~ flow, data = df)
  m_flow_wind <- lm(plume_area ~ flow + wind_spd, data = df)
  m_flow_wave <- lm(plume_area ~ flow + wave_height, data = df)

  anova_wind <- anova(m_flow, m_flow_wind)
  anova_wave <- anova(m_flow, m_flow_wave)

  model_comparison <- tibble::tibble(
    added_term = c("wind_spd", "wave_height"),
    f_statistic = c(anova_wind$F[2], anova_wave$F[2]),
    p_value = c(anova_wind$`Pr(>F)`[2], anova_wave$`Pr(>F)`[2]),
    r2_flow_only = summary(m_flow)$r.squared,
    r2_with_term = c(summary(m_flow_wind)$r.squared, summary(m_flow_wave)$r.squared)
  )
  write_csv(model_comparison, "output/STATS/rhone_wind_wave_models.csv")

  # b) Mistral vs onshore/easterly vs calm, holding flow's effect fixed by
  #    working with the *residual* plume area after removing the flow
  #    relationship (m_flow above) -- so what gets compared across wind
  #    categories is "area given today's flow", not raw area (which would
  #    just reflect whichever category happened to have wetter days).
  #    Sector definitions are an approximation, not a precise meteorological
  #    classification: Mistral ~ NNW-N (e.g. Guenard et al. 2006 use
  #    ~300-340 deg); the onshore/"Levant" sector ~E-SE is the rough
  #    opposite. wind_dir follows the meteorological "from" convention (see
  #    util.R::.speed_direction()), so these are compass bearings the wind
  #    is blowing FROM.
  df$area_resid <- residuals(m_flow)
  df$wind_category <- dplyr::case_when(
    df$wind_spd < 3                          ~ "calm (<3 m s⁻¹)",
    df$wind_dir >= 300 & df$wind_dir <= 350  ~ "Mistral (NNW-N)",
    df$wind_dir >= 90  & df$wind_dir <= 150  ~ "onshore/easterly",
    TRUE                                     ~ "other"
  )
  # Ordered calm -> onshore/easterly -> Mistral -> other (weakest/most
  # coast-hugging to strongest/most offshore-pushing, per the e-mail's
  # hypothesis), rather than the alphabetical default, so this order is
  # what every plot legend, boxplot axis, and summarise() output below uses.
  df$wind_category <- factor(df$wind_category,
                              levels = c("calm (<3 m s⁻¹)", "onshore/easterly", "Mistral (NNW-N)", "other"))
  wind_category_colours <- c("calm (<3 m s⁻¹)" = "grey50", "onshore/easterly" = "steelblue",
                              "Mistral (NNW-N)" = "firebrick", "other" = "goldenrod")

  category_summary <- df |>
    dplyr::summarise(n_days = dplyr::n(),
                      mean_area_resid = mean(area_resid, na.rm = TRUE),
                      sd_area_resid = sd(area_resid, na.rm = TRUE), .by = "wind_category") |>
    dplyr::arrange(wind_category)  # .by = doesn't preserve factor level order, arrange() does
  write_csv(category_summary, "output/STATS/rhone_wind_wave_categories.csv")
  category_anova <- summary(aov(area_resid ~ wind_category, data = df))

  pl_terms <- ggplot(df, aes(x = wind_spd, y = plume_area)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm", colour = "purple") +
    labs(x = "Wind speed (m s⁻¹)", y = "Plume area (km²)",
         title = "Grand Rhone: plume area vs. wind speed",
         subtitle = "Raw relationship, not controlled for flow (see model_comparison)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  # c) Same flow-controlled residual as the category summary above, but as a
  #    scatter against wind speed (not a categorical boxplot) with a
  #    separate per-category trend line -- shows whether wind speed's
  #    relationship with (flow-controlled) plume area actually differs by
  #    wind category, rather than just comparing category means.
  pl_category <- ggplot(df, aes(x = wind_spd, y = area_resid, colour = wind_category)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    scale_colour_manual(values = wind_category_colours) +
    labs(x = "Wind speed (m s⁻¹)", y = "Plume area residual after removing the flow effect (km²)",
         colour = NULL, title = "Grand Rhone: flow-controlled area vs. wind, by category") +
    theme(panel.border = element_rect(fill = NA, colour = "black"), legend.position = "none")

  # d) The e-mail's "or even wave, which may be simpler?" alternative to
  #    wind -- this plot did not exist before. Plume area vs. wave height,
  #    with wind category as colour so the wind effect is visible alongside
  #    the wave effect on the same panel, per the e-mail's suggestion.
  pl_wave <- ggplot(df, aes(x = wave_height, y = plume_area, colour = wind_category)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    scale_colour_manual(values = wind_category_colours) +
    labs(x = "Wave height (m)", y = "Plume area (km²)", colour = NULL,
         title = "Grand Rhone: plume area vs. wave height, by wind category",
         subtitle = "Raw relationship, not controlled for flow") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  # NB: ggpubr's common.legend = TRUE (sharing one legend for pl_category
  # and pl_wave) renders a malformed black panel with this ggplot2 version,
  # avoided here by just dropping pl_category's legend (identical
  # categories/colours to pl_wave's) and keeping pl_wave's own.
  pl_combi <- ggpubr::ggarrange(pl_terms, pl_category, pl_wave, ncol = 3, nrow = 1)
  ggsave(filename = "figures/rhone_side_analyses/rhone_wind_wave_effect.png", plot = pl_combi, width = 18, height = 6)

  return(list(model_comparison = model_comparison, wind_category_summary = category_summary,
              wind_category_anova = category_anova, data = df, plot = pl_combi))
}


## Follow-up: does the wind effect on plume area differ between winter and the stratified season? ----
# Follow-up to rhone_wind_wave_effect() above: her reply asked whether the
# wind effect (Mistral vs onshore/calm, item 3 above) shows "nuances" between
# winter and "la periode stratifiee". RiOMar has no near-mouth
# density/stratification product to test this against directly (GLORYS MLD
# exists but is too coarse/open-ocean for a coastal plume, per her own
# reservation about item 4 in the original e-mail) -- calendar season is
# used here as a climatological proxy for stratification state, not a
# per-day measurement of it. Winter (DJF) stands in for a well-mixed water
# column, summer (JJA) for a stratified one, using the same meteorological
# season definitions already used in rhone_flood_timing_shift() above.
# Spring/autumn (MAM/SON) are transitional -- Gulf of Lion stratification
# builds and breaks down gradually through those months -- so they are
# dropped from this comparison rather than folded into either state, which
# would blur exactly the contrast being tested.
# Caveat worth keeping alongside the result: Mistral itself is a mechanism
# that can erode near-surface stratification on the day it blows, so a
# Mistral day tagged JJA is not guaranteed to be under stratified conditions
# that day -- season is a proxy for the *typical* state, not a direct
# per-day measurement.
# rhone_wind_wave_seasonal_effect()
rhone_wind_wave_seasonal_effect <- function(){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")
  df_flow <- combine_plume_driver("flow", meta) |> dplyr::select(date, plume_area, flow = value)
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, wind_spd = value, wind_dir)

  df <- df_flow |>
    dplyr::left_join(df_wind, by = "date") |>
    tidyr::drop_na(plume_area, flow, wind_spd)

  # Flow-only baseline model fit on the full record (not per season), so the
  # residual below isolates wind's effect alone -- fitting a separate flow
  # model per season would let any season-specific area-flow relationship
  # contaminate what is meant to be a wind-only comparison.
  m_flow <- lm(plume_area ~ flow, data = df)
  df$area_resid <- residuals(m_flow)

  # Same wind categories, ordering, and colours as rhone_wind_wave_effect()
  # above, so the two analyses stay directly comparable.
  df$wind_category <- dplyr::case_when(
    df$wind_spd < 3                          ~ "calm (<3 m s⁻¹)",
    df$wind_dir >= 300 & df$wind_dir <= 350  ~ "Mistral (NNW-N)",
    df$wind_dir >= 90  & df$wind_dir <= 150  ~ "onshore/easterly",
    TRUE                                     ~ "other"
  )
  df$wind_category <- factor(df$wind_category,
                              levels = c("calm (<3 m s⁻¹)", "onshore/easterly", "Mistral (NNW-N)", "other"))
  wind_category_colours <- c("calm (<3 m s⁻¹)" = "grey50", "onshore/easterly" = "steelblue",
                              "Mistral (NNW-N)" = "firebrick", "other" = "goldenrod")

  # Meteorological DJF/JJA, same definition as rhone_flood_timing_shift()
  # above. MAM/SON are dropped from this comparison (see header comment).
  df$month <- lubridate::month(df$date)
  df$season <- dplyr::case_when(
    df$month %in% c(12, 1, 2) ~ "winter (DJF)",
    df$month %in% 6:8         ~ "summer (JJA)",
    TRUE                      ~ NA_character_
  )
  df <- dplyr::filter(df, !is.na(season))
  df$season <- factor(df$season, levels = c("winter (DJF)", "summer (JJA)"))

  # a) The actual test of "does the wind effect differ by season": the
  #    wind_category:season interaction term in a two-way ANOVA on the
  #    flow-controlled residual. Two separate within-season ANOVAs (one for
  #    DJF, one for JJA) would only support an informal eyeball comparison;
  #    the interaction term is the formal test of whether the categories'
  #    effect actually differs between the two seasons.
  fit_interaction <- aov(area_resid ~ wind_category * season, data = df)
  interaction_table <- broom::tidy(fit_interaction)
  write_csv(interaction_table, "output/STATS/rhone_wind_wave_seasonal_models.csv")

  # b) By-season, by-category summary -- the interpretable companion to (a).
  category_summary <- df |>
    dplyr::summarise(n_days = dplyr::n(),
                      mean_area_resid = mean(area_resid, na.rm = TRUE),
                      sd_area_resid = sd(area_resid, na.rm = TRUE), .by = c("season", "wind_category")) |>
    dplyr::arrange(season, wind_category)  # .by = doesn't preserve factor level order, arrange() does
  write_csv(category_summary, "output/STATS/rhone_wind_wave_seasonal_categories.csv")

  # c) Flow-controlled area by wind category, faceted by season -- the
  #    visual companion to the interaction test above.
  interaction_row <- dplyr::filter(interaction_table, term == "wind_category:season")
  pl <- ggplot(df, aes(x = wind_category, y = area_resid, fill = wind_category)) +
    geom_boxplot(outlier.alpha = 0.2) +
    scale_fill_manual(values = wind_category_colours) +
    facet_wrap(~season, ncol = 2) +
    labs(x = NULL, y = "Plume area residual after removing the flow effect (km²)",
         title = "Grand Rhone: does the wind effect on plume area differ between winter and summer?",
         subtitle = paste0("wind_category x season interaction: F = ", round(interaction_row$statistic, 2),
                            ", p = ", signif(interaction_row$p.value, 3)),
         fill = NULL) +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
  ggsave(filename = "figures/rhone_side_analyses/rhone_wind_wave_seasonal.png", plot = pl, width = 10, height = 6)

  return(list(interaction_table = interaction_table, category_summary = category_summary, data = df, plot = pl))
}


## Follow-up: has the Rhone flow trend moved differently across the flow distribution? ----
# Raised as pure curiosity in her reply, not a work request: particle
# concentration seems to rise non-linearly with flow, so does the flow
# trend itself look different depending on where in the flow distribution
# you look, rather than only at the mean? The single trend reported
# elsewhere in this file and in the e-mail (-3.3 to -3.4 m3/s/yr, p~0.6-0.65)
# is a mean/OLS-style trend and would miss a shape change, e.g. a declining
# flood tail masked by a flat median. Quantile regression (quantreg::rq())
# fits a separate linear trend at each of several percentiles of the flow
# distribution so that shape change, if present, is visible directly --
# this only characterises the flow distribution's own trend, it does not
# model the concentration-vs-flow relationship itself.
# rhone_flow_quantile_trend()
rhone_flow_quantile_trend <- function(taus = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95)){

  meta <- get_zone_meta(mouth_name = "Grand Rhone")
  df <- load_driver("flow", meta) |> tidyr::drop_na(value)
  # rq() needs a numeric regressor; Date's own numeric encoding is already
  # days-since-epoch, so the fitted slope is directly per-day (matching the
  # slope/slope_annualised convention used by fit_wls_hac_trend() elsewhere
  # in this file).
  df$date_num <- as.numeric(df$date)

  # se = "boot": quantile regression has no closed-form standard error, so a
  # bootstrap is used for the CI/p-value at every quantile.
  fit <- quantreg::rq(value ~ date_num, data = df, tau = taus)
  fit_summary <- summary(fit, se = "boot")

  stats <- purrr::map2_dfr(fit_summary, taus, function(s, tau){
    coefs <- s$coefficients
    tibble::tibble(tau = tau,
                   slope = coefs["date_num", "Value"],
                   slope_annualised = coefs["date_num", "Value"] * 365.25,
                   slope_se = coefs["date_num", "Std. Error"],
                   slope_p = coefs["date_num", "Pr(>|t|)"])
  })
  write_csv(stats, "output/STATS/rhone_flow_quantile_trend.csv")

  # a) Trend point estimate +/- 95% CI at each quantile -- shows directly
  #    whether the trend's sign/magnitude changes across the distribution.
  pl_trend <- ggplot(stats, aes(x = tau, y = slope_annualised)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_pointrange(aes(ymin = slope_annualised - 1.96 * slope_se * 365.25,
                         ymax = slope_annualised + 1.96 * slope_se * 365.25)) +
    scale_x_continuous(breaks = taus, labels = scales::percent) +
    labs(x = "Flow quantile", y = "Flow trend (m³ s⁻¹ yr⁻¹)",
         title = "Grand Rhone: does the flow trend differ across the flow distribution?",
         subtitle = "Point = quantile regression slope, error bars = 95% CI (bootstrap SE)") +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  # b) Same fits shown against the raw daily flow series, for context.
  pl_lines <- ggplot(df, aes(x = date, y = value)) +
    geom_point(alpha = 0.08, size = 0.5, colour = "grey40") +
    geom_quantile(quantiles = taus, formula = y ~ x, colour = "firebrick", linewidth = 0.7) +
    labs(x = NULL, y = "Grand Rhone flow (m³ s⁻¹)",
         title = "Grand Rhone: flow with per-quantile trend lines",
         subtitle = paste0("Quantiles shown: ", paste0(taus * 100, "%", collapse = ", "))) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))

  pl_combi <- ggpubr::ggarrange(pl_lines, pl_trend, ncol = 2, nrow = 1)
  ggsave(filename = "figures/rhone_side_analyses/rhone_flow_quantile_trend.png", plot = pl_combi, width = 14, height = 6)

  return(list(stats = stats, data = df, plot = pl_combi))
}


# NB: not run automatically on source() -- call explicitly
# rhone_detrend_test()
# rhone_flood_timing_shift()
# rhone_wind_wave_effect()
# rhone_wind_wave_seasonal_effect()
# rhone_flow_quantile_trend()

