# func/compute_wave_category_lagged_correlation.R
#
# Supplementary Fig. daily_flow (func/figure.R::Figure_S_daily_flow(), formerly
# manuscript Figure 5 before the seasonal-analysis section repurposed that
# slot) shows, per zone: (left) a scatter of daily plume area vs. river flow with an OLS
# fit; (right) the lagged correlation between the two daily series (plume
# lagged 0-14 days behind flow), with the lag of maximum correlation marked.
# Robert asked for the same analysis for wave height, restricted separately
# to onshore-wind and offshore-wind days (the calm/onshore/offshore
# categories established for manuscript Figure 8 and reused throughout
# func/compute_driver_x11_figures.R).
#
# Filtering to one category's days *before* computing the lag would corrupt
# the lag itself: func/util.R::lagged_correlation() shifts by row position
# (dplyr::lag()), so removing the other categories' days would make a
# "lag of 3" mean "3 onshore-days ago" (spanning a variable number of actual
# calendar days) rather than 3 real days. Instead, category_lag_correlation()
# below computes the lagged driver value from the *full, continuous* daily
# series (so a lag really is a lag in days) and only restricts which day's
# *plume area* the correlation is computed against to the category's days --
# i.e. "given that today is an onshore day, how does today's plume area
# correlate with wave height N days ago".
#
# Run from repo root: Rscript func/compute_wave_category_lagged_correlation.R

source("func/multi.R")

OUT_DIR <- "figures/driver_x11_comparison/lagged_by_category"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CATEGORIES_TO_PLOT <- c("onshore", "offshore")
MAX_LAG_DAILY <- 14

# df must be date-ordered, complete (load_plume_ts()/combine_plume_driver()
# already gap-fill to a continuous daily sequence), with plume_area, value,
# and wind_category columns.
category_lag_correlation <- function(df, category, max_lag = MAX_LAG_DAILY){
  in_category <- df$wind_category == category
  purrr::map_dfr(0:max_lag, function(l){
    lagged_value <- dplyr::lag(df$value, l)
    tibble::tibble(lag = l, cor = cor(df$plume_area[in_category], lagged_value[in_category], use = "pairwise.complete.obs"))
  })
}

# Same wind-category derivation used throughout compute_driver_x11_figures.R.
add_wind_category <- function(df, meta){
  df_wind <- load_driver("wind", meta) |> dplyr::select(date, cat_wind_spd = value, cat_wind_direction = direction)
  df <- df |> dplyr::left_join(df_wind, by = "date") |> tidyr::drop_na(cat_wind_spd, cat_wind_direction)
  df$wind_category <- dplyr::case_when(
    df$cat_wind_spd < 3            ~ "calm (<3 m/s)",
    df$cat_wind_direction == "off" ~ "offshore",
    TRUE                           ~ "onshore"
  )
  df
}

panel_theme <- theme(plot.title = element_text(hjust = 0.5, size = 20),
                     axis.title = element_text(size = 15, colour = "black"),
                     axis.text = element_text(size = 12, colour = "black"),
                     plot.margin = margin(t = 8, r = 12, b = 5, l = 5),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(linetype = "solid", fill = NA))

message("Building Figure 5-style wave-height lagged correlation, per wind category...")
for(category in CATEGORIES_TO_PLOT){
  panels <- purrr::map(zones, function(zone_name){
    meta <- get_zone_meta(zone_name = zone_name)
    df <- combine_plume_driver("wave", meta) |> add_wind_category(meta)
    df_cat <- dplyr::filter(df, wind_category == category)

    cor_df <- category_lag_correlation(df, category)
    peak <- cor_df |> dplyr::slice_max(cor, n = 1)

    zone_label <- zone_title(zone_name)

    scatter_plot <- ggplot(df_cat, aes(x = value, y = plume_area)) +
      geom_point(alpha = 0.3, colour = "grey30", size = 0.8) +
      geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1.2) +
      labs(x = "Wave height (m)", y = "Plume area (km²)", title = zone_label,
          subtitle = paste0(category, " days only, n = ", nrow(df_cat))) +
      panel_theme + theme(plot.subtitle = element_text(hjust = 0.5, size = 11, colour = "grey40"))

    lag_plot <- ggplot(cor_df, aes(x = lag, y = cor)) +
      geom_line(colour = "grey30") +
      geom_point(colour = "grey30") +
      geom_point(data = peak, colour = "firebrick", size = 3) +
      geom_text(data = peak, aes(label = paste0(lag, "d")), vjust = -1.2, colour = "firebrick", size = 4) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(x = "Lag, plume after wave height (days)", y = "Correlation (r)", title = zone_label) +
      panel_theme

    list(scatter = scatter_plot, lag = lag_plot)
  })

  plotlist <- panels |> purrr::map(function(p) list(p$scatter, p$lag)) |> purrr::flatten()
  panel_labels <- paste0(letters[seq_along(plotlist)], ")")
  full_plot <- ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(panels), align = "v",
                                 labels = panel_labels, font.label = list(size = 18, face = "bold"),
                                 hjust = -0.3, vjust = 1.3)

  out_file <- file.path(OUT_DIR, paste0("wave_lagged_", category, ".png"))
  ggplot2::ggsave(out_file, full_plot, width = 14, height = 20, dpi = 200)
  message("Wrote ", out_file)
}
