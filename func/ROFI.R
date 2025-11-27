# func/ROFI.R

# The analyses used to compare ROFI to plume size


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(patchwork)
library(doParallel); doParallel::registerDoParallel(cores = 14)

# Zones
zones <- c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION")


# Load data ---------------------------------------------------------------

## Plumes ------------------------------------------------------------------

# Load plume time series
df_plume_ts <- map_dfr(zones, load_plume_ts)

# Load all daily maps into one data.frame
## NB: There are a lot of files to load, need some heavy lifting to get it done
# df_plume_surface <- map_dfr(zones, load_plume_surface)


## ROFI --------------------------------------------------------------------

# Load ROFI
df_ROFI <- map_dfr(dir("data/ROFI", full.names = TRUE), load_ROFI)


# Analyses ----------------------------------------------------------------

# Merge plume and ROFI data
df_ROFI_plume <- df_plume_ts |> 
  left_join(df_ROFI, by = c("zone", "date")) |> 
  filter(!is.na(ROFI_surface)) |> 
  dplyr::select(zone, date, plume_area, ROFI_surface) |> 
  make_pretty_title()

# Monthly means
# df_ROFI_plume_monthly <- df_ROFI_plume |> 
#   mutate(year = lubridate::year(date),
#          month = lubridate::month(date)) |> 
#   group_by(zone, year, month) |> 
#   summarise(plume_area = mean(plume_area, na.rm = TRUE),
#             ROFI_surface = mean(ROFI_surface, na.rm = TRUE),
#             .groups = "drop") |> 
#   mutate(date = as.Date(paste(year, month, "15", sep = "-"))) |> 
#   dplyr::select(zone, date, plume_area, ROFI_surface)

# Lagged correlations
ROFI_plume_cor <- df_ROFI_plume |> 
  reframe(
    lag = 0:365,
    cor = map_dbl(0:365, ~ cor(ROFI_surface, lag(plume_area, .), use = "pairwise.complete.obs")),
    .by = c("zone", "plot_title")
  )

# Plot plume area and ROFI surface as line plots with different colour for each variable, faceted by zone
p_ROFI_plume_ts <- df_ROFI_plume |> 
  mutate(plume_area_lag = lag(plume_area, 30)) |>
  pivot_longer(cols = c(plume_area, ROFI_surface), names_to = "variable", values_to = "value") |> 
  mutate(variable = recode(variable,
                           plume_area = "Plume Area (km²)",
                           ROFI_surface = "ROFI Surface (km²)")) |>
  ggplot(aes(x = date, y = value, colour = variable)) +
  geom_line() +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(title = "Daily time series of plume area and ROFI surface", y = "Value", x = "Date", colour = "Variable") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "bottom")
p_ROFI_plume_ts
# ggsave("output/ROFI/ROFI_plume_time_series.png", p_ROFI_plume_ts, width = 12, height = 8)

# Plot wind speed and panache size correlation
p_ROFI_plume_scatter <- ggplot(df_ROFI_plume, aes(x = ROFI_surface, y = plume_area)) + 
  # geom_hex(bins = 30, show.legend = FALSE) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linewidth = 1.5, linetype = "dashed", color = "orchid") +
  geom_smooth(method = "lm", linewidth = 1.5) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(y = "Plume Surface (km²)", x = "ROFI Surface (km²))") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "bottom")
p_ROFI_plume_scatter
# ggsave("output/ROFI/ROFI_plume_scatter.png", p_ROFI_plume_scatter, width = 8, height = 10)

# Plot lag results
p_ROFI_plume_cor_lag <- ggplot(ROFI_plume_cor, aes(x = lag, y = cor)) +
  geom_point() +
  labs(x = "lag plume surface area after ROFI surface area (days)", y = "correlation (r)") +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
p_ROFI_plume_cor_lag

# Combine plots and save
plot_title <- grid::textGrob("ROFI surface area vs plume surface area (km^2)", gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
ts_plot <- ggpubr::ggarrange(p_ROFI_plume_ts, labels = c("a)"))
cor_plot <- ggpubr::ggarrange(p_ROFI_plume_scatter, p_ROFI_plume_cor_lag, ncol = 1, nrow = 2, labels = c("b)", "c)"), heights = c(1, 0.3))
full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
full_plot_title <- ggpubr::ggarrange(plot_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
ggsave(filename = paste0("figures/cor_plot_ROFI_plume.png"), plot = full_plot_title, width = 36, height = 18, dpi = 600)

