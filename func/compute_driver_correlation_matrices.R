# func/compute_driver_correlation_matrices.R
#
# Exploratory Pearson correlation matrices, per zone (and separately, per
# zone x season), across the physical drivers (river discharge, wind speed,
# tidal range, coastal current, wave height) and three plume-response
# metrics (plume area, mean SPM, SPM mass) -- requested by Robert alongside
# the driver-interactions rerun, to sanity-check the GAM/RF results (see the
# RiOMar memory project_driver_rerun_key_findings.md: RF permutation
# importance ranks wave height above river discharge at 3/4 zones, and one
# open question is whether that's distorted by correlated predictors -- this
# gives a first, simple look at exactly that).
#
# Wind direction and wave direction are deliberately excluded from the
# matrix: both are circular (compass-angle) variables, and a Pearson
# correlation computed on raw degrees is not a valid measure of association
# for a circular quantity (it ignores the 0/360 wrap-around, so e.g. 359 deg
# and 1 deg -- nearly identical directions -- look maximally different to a
# linear correlation). Direction's relationship with plume area is instead
# assessed elsewhere via a categorical (8-octant) GAM term (manuscript
# Table 6), not a linear correlation. If a circular measure is wanted here
# too, that would need a different statistic (e.g. circular-linear
# correlation, or entering direction as sin/cos components) -- flagging
# rather than silently approximating it with a meaningless Pearson r.
#
# Seasons use the same DJF/MAM/JJA/SON meteorological convention, and the
# same December-belongs-to-the-following-winter attribution, already
# established in this file's rhone_flood_timing_shift() (Dec 2010 +
# Jan/Feb 2011 are one winter, not split across two calendar years) -- see
# that function's comment for the reasoning. season_year is computed here
# for consistency with that precedent even though the pooled (all-years)
# per-season matrices below only group by season, not season_year: any
# December is unconditionally winter regardless of which season_year label
# it gets, so season_year doesn't change which pool a day falls into here,
# it just keeps the convention identical for anyone extending this to a
# per-winter (not per-season-type) breakdown later.
#
# Run from repo root: Rscript func/compute_driver_correlation_matrices.R

source("func/driver_interactions.R")  # brings in multi.R (zones, zone_title, build_driver_matrix, get_zone_meta, etc.)

# Response metrics + drivers included in the matrix, in display order -------
VARS <- c(
  plume_area                          = "Plume area",
  mean_SPM_in_the_plume_area          = "Mean SPM",
  mass_SPM_in_the_plume_area_in_g_m   = "SPM mass",
  flow                                = "River discharge",
  wind_spd                            = "Wind speed",
  tide_range                          = "Tidal range",
  current                             = "Current speed",
  wave_height                         = "Wave height"
)
VAR_LEVELS <- unname(VARS)

# Reads the daily_driver_matrix_<zone>.csv already saved by
# run_full_analysis()'s Step 0 (the same build_driver_matrix() output, from
# the completed 2026-07-31 rerun) rather than rebuilding it from the raw
# wind/wave/current NetCDF archives -- seconds instead of ~20 minutes per
# zone, and identical data since nothing about driver loading has changed
# since that rerun.
load_saved_driver_matrix <- function(zone_name, stats_dir = "output/STATS"){
  readr::read_csv(file.path(stats_dir, paste0("daily_driver_matrix_", zone_name, ".csv")),
                  show_col_types = FALSE)
}

compute_cor_matrix <- function(df){
  df_sub <- df |> dplyr::select(dplyr::all_of(names(VARS))) |> as.data.frame()
  colnames(df_sub) <- unname(VARS)
  list(cor = cor(df_sub, use = "pairwise.complete.obs"),
       n = sum(stats::complete.cases(df_sub)))
}

plot_correlation_matrix <- function(cor_result, title){
  cor_df <- as.data.frame(as.table(cor_result$cor))
  names(cor_df) <- c("Var1", "Var2", "r")
  cor_df$Var1 <- factor(cor_df$Var1, levels = VAR_LEVELS)
  cor_df$Var2 <- factor(cor_df$Var2, levels = rev(VAR_LEVELS))

  ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2, fill = r)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 3.1) +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                  midpoint = 0, limits = c(-1, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, subtitle = paste0("n = ", cor_result$n, " complete days"),
                 x = NULL, y = NULL, fill = "r") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, colour = "grey40"),
                   panel.grid = ggplot2::element_blank())
}

# Adds meteorological season + season_year to a build_driver_matrix() output,
# same convention as rhone_flood_timing_shift() below in this file.
assign_season <- function(df){
  df |>
    dplyr::mutate(month = lubridate::month(date), year = lubridate::year(date),
                  season = dplyr::case_when(
                    month %in% 3:5  ~ "Spring (MAM)",
                    month %in% 6:8  ~ "Summer (JJA)",
                    month %in% 9:11 ~ "Autumn (SON)",
                    TRUE            ~ "Winter (DJF)"),
                  season_year = ifelse(month == 12, year + 1, year))
}
SEASON_LEVELS <- c("Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Autumn (SON)")


# 1. Overall (all-season) correlation matrix, one per zone ------------------

message("Building overall (all-season) correlation matrices...")
overall_plots <- purrr::map(zones, function(zone_name){
  df <- load_saved_driver_matrix(zone_name)
  plot_correlation_matrix(compute_cor_matrix(df), zone_title(zone_name))
})

overall_grid <- patchwork::wrap_plots(overall_plots, ncol = 2) +
  patchwork::plot_annotation(title = "Driver / plume-metric correlation matrix, per zone (all seasons pooled)",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15)))
ggplot2::ggsave("output/STATS/driver_correlation_matrix_overall.png", overall_grid, width = 13, height = 11, dpi = 200)
message("Wrote output/STATS/driver_correlation_matrix_overall.png")


# 2. Per-season correlation matrices, one per zone x season ------------------

message("Building per-season correlation matrices...")
seasonal_plots <- purrr::map(zones, function(zone_name){
  df <- load_saved_driver_matrix(zone_name) |> assign_season()
  purrr::map(SEASON_LEVELS, function(s){
    df_season <- dplyr::filter(df, season == s)
    plot_correlation_matrix(compute_cor_matrix(df_season), paste0(zone_title(zone_name), " -- ", s))
  })
}) |> purrr::flatten()

# Row = zone (north to south), column = season (Winter/Spring/Summer/Autumn)
seasonal_grid <- patchwork::wrap_plots(seasonal_plots, ncol = length(SEASON_LEVELS)) +
  patchwork::plot_annotation(title = "Driver / plume-metric correlation matrix, per zone x season",
                             theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15)))
ggplot2::ggsave("output/STATS/driver_correlation_matrix_seasonal.png", seasonal_grid, width = 24, height = 20, dpi = 200)
message("Wrote output/STATS/driver_correlation_matrix_seasonal.png")
