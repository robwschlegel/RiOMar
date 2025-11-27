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
df_plume_surface <- map_dfr(zones, load_plume_surface)


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
df_ROFI_plume_monthly <- df_ROFI_plume |> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) |> 
  group_by(zone, year, month) |> 
  summarise(plume_area = mean(plume_area, na.rm = TRUE),
            ROFI_surface = mean(ROFI_surface, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(date = as.Date(paste(year, month, "15", sep = "-"))) |> 
  dplyr::select(zone, date, plume_area, ROFI_surface)

# Plot plume area and ROFI surface as line plots with different colour for each variable, faceted by zone
p_ROFI_plume_ts <- df_ROFI_plume |> 
  pivot_longer(cols = c(plume_area, ROFI_surface), names_to = "variable", values_to = "value") |> 
  mutate(variable = recode(variable,
                           plume_area = "Plume Area (km²)",
                           ROFI_surface = "ROFI Surface (km²)")) |>
  ggplot(aes(x = date, y = value, colour = variable)) +
  geom_line() +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(title = "Daily time series of plume area and ROFI surface", y = "Value", x = "Date", colour = "Variable") +
  theme_minimal()
p_ROFI_plume_ts
# ggsave("output/ROFI/ROFI_plume_time_series.png", p_ROFI_plume_ts, width = 12, height = 8)

# Scatter plot of plume area vs ROFI surface, with a linear regression line, faceted by zone
p_ROFI_plume_scatter <- df_ROFI_plume |> 
  ggplot(aes(x = ROFI_surface, y = plume_area)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", colour = "red") +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(title = "Scatter plot of plume area vs ROFI surface", y = "Plume Surface (km²)", x = "ROFI Surface (km²)") +
  theme_minimal()
p_ROFI_plume_scatter
# ggsave("output/ROFI/ROFI_plume_scatter.png", p_ROFI_plume_scatter, width = 8, height = 10)

# Compare panache size against wind speed
# Two separate comparisons based on upwelling or downwelling times
ROFI_plume_cor <- df_ROFI_plume |> 
  reframe(
    lag = 0:365,
    cor = map_dbl(0:365, ~ cor(ROFI_surface, lag(plume_area, .), method = "kendall")),
    .by = c("zone", "plot_title")
  )

# Plot wind speed and panache size correlation
p_plume_ROFI_cor <- ggplot(df_ROFI_plume, aes(x = ROFI_surface, y = plume_area)) + 
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm") +
  labs(y = "Plume Surface (km²)", x = "ROFI Surface (km²))") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "bottom")
# p_plume_ROFI_cor

# Plot lag results
p_plume_ROFI_cor_lag <- ggplot(ROFI_plume_cor, aes(x = lag, y = cor)) +
  geom_point() +
  labs(x = "lag plume after tidal range (days)", y = "correlation (r)") +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
# flow_plume_cor_lag_plot

# Combine plots and save
flow_plume_title <- grid::textGrob(paste0(mouth_info$mouth_name," : river flow vs plume size"), gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
ts_plot <- ggpubr::ggarrange(flow_plot, panache_plot, ncol = 1, nrow = 2, labels = c("a)", "b)"), align = "v")
cor_plot <- ggpubr::ggarrange(flow_plume_cor_plot, flow_plume_cor_lag_plot, ncol = 1, nrow = 2, labels = c("c)", "d)"), heights = c(1, 0.3))
full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 2, nrow = 1)
full_plot_title <- ggpubr::ggarrange(flow_plume_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
ggsave(filename = paste0("figures/cor_plot_flow_plume_",mouth_info$mouth_name,".png"), width = 12, height = 6, dpi = 600)

