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

# Plot the analysis one zone at a time for more clarity
# zone_name <- zones[1]
comp_ROFI_plume <- function(zone_name, df = df_ROFI_plume){
  
  # Subset the zone
  df_zone <- df |> filter(zone == zone_name) |> 
    mutate(timestep = "daily")
  
  # Create monthly means
  df_zone_monthly <- df_zone |> 
    mutate(year = lubridate::year(date),
           month = lubridate::month(date)) |> 
    summarise(plume_area = mean(plume_area, na.rm = TRUE),
              ROFI_surface = mean(ROFI_surface, na.rm = TRUE),
              .by = c("zone", "plot_title", "year", "month")) |> 
    mutate(date = as.Date(paste(year, month, "15", sep = "-")),
           timestep = "monthly") |> 
    dplyr::select(date, zone, plot_title, plume_area, ROFI_surface, timestep)
  
  # Create annual means
  df_zone_annual <- df_zone |> 
    mutate(year = lubridate::year(date)) |> 
    summarise(plume_area = mean(plume_area, na.rm = TRUE),
              ROFI_surface = mean(ROFI_surface, na.rm = TRUE),
              .by = c("zone", "plot_title", "year")) |> 
    mutate(date = as.Date(paste(year, "07", "01", sep = "-")),
           timestep = "annual") |> 
    dplyr::select(date, zone, plot_title, plume_area, ROFI_surface, timestep)
  
  # df with all timesteps
  df_timestep <- bind_rows(df_zone, df_zone_monthly, df_zone_annual) |> 
    mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual"))) |> 
    pivot_longer(cols = c(plume_area, ROFI_surface), names_to = "variable", values_to = "value") |> 
    mutate(variable = recode(variable,
                             plume_area = "Plume Surface (km²)",
                             ROFI_surface = "ROFI Surface (km²)")) 
  
  # Lagged correlations
  ROFI_plume_cor_daily <- df_zone |> 
    reframe(lag = 0:365,
            cor = map_dbl(0:365, ~ cor(ROFI_surface, lag(plume_area, .), use = "pairwise.complete.obs"))) |> 
    mutate(timestep = "daily")
  ROFI_plume_cor_monthly <- df_zone_monthly |> 
    reframe(lag = 0:12,
            cor = map_dbl(0:12, ~ cor(ROFI_surface, lag(plume_area, .), use = "pairwise.complete.obs"))) |>
    mutate(timestep = "monthly")
  ROFI_plume_cor_annual <- df_zone_annual |> 
    reframe(lag = 0:10,
            cor = map_dbl(0:10, ~ cor(ROFI_surface, lag(plume_area, .), use = "pairwise.complete.obs"))) |> 
    mutate(timestep = "annual")
  
  # Create label df
  ROFI_plume_cor_all <- bind_rows(ROFI_plume_cor_daily, ROFI_plume_cor_monthly, ROFI_plume_cor_annual)
  
  # Determine daily y label position based on zone
  if(zone_name == "BAY_OF_SEINE"){
    y_pos_daily <- c(8500, 8000)
  } else if(zone_name == "BAY_OF_BISCAY"){
    y_pos_daily <- c(15000, 14200)
  } else if(zone_name == "SOUTHERN_BRITTANY"){
    y_pos_daily <- c(10200, 9600)
  } else if(zone_name == "GULF_OF_LION"){
    y_pos_daily <- c(900, 800)
  }
  
  # Calculate trend and create labels
  # trend_labels_all <- df_timestep |> 
  #   summarise(slope = round(coef(lm(value ~ date))["date"] * 365.25), 
  #             .by = c("zone", "plot_title", "timestep", "variable")) |>
  #   mutate(variable = recode(variable,
  #                            plume_area = "Plume Surface (km²)",
  #                            ROFI_surface = "ROFI Surface (km²)")) |> 
  #   mutate(x = as.Date("2002-01-01"),
  #          y = c(9000, 8500, 
  #                7500, 7000, 
  #                5500, 5000))
  
  # Or just the daily values
  trend_labels_daily <- df_timestep |> 
    filter(timestep == "daily") |>
    summarise(slope = round(coef(lm(value ~ date))["date"] * 365.25), 
              .by = c("zone", "plot_title", "timestep", "variable")) |>
    mutate(variable = recode(variable,
                             plume_area = "Plume Surface (km²)",
                             ROFI_surface = "ROFI Surface (km²)")) |> 
    mutate(x = as.Date("2002-01-01"),
           y = y_pos_daily)
  
  # Plot plume area and ROFI surface as line plots with different colour for each variable
  p_ROFI_plume_ts <- ggplot(df_timestep, aes(x = date, y = value, colour = variable, linetype = timestep)) +
    geom_line(data = filter(df_timestep, timestep == "daily")) +
    geom_line(data = filter(df_timestep, timestep == "monthly"), linewidth = 1.5) +
    geom_line(data = filter(df_timestep, timestep == "annual"), linewidth = 1.5) +
    geom_label(data = trend_labels_daily, show.legend = FALSE,
               aes(x = x, y = y, colour = variable, size = 6, hjust = 0,
                   label = paste0(timestep," ",variable," slope = ", round(slope, 2), " km^2 yr-1", sep = ""))) +
    scale_linetype_manual(values = c("daily" = "dotted", "monthly" = "dashed", "annual" = "solid")) +
    scale_x_date(expand = c(0, 0)) +
    labs(title = paste0("Time series of ROFI and plume surface area"), 
         y = "Surface Area (km^2)", x = NULL, colour = "Variable") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  # p_ROFI_plume_ts
  
  # Plot wind speed and panache size correlation
  # TODO: Add a second scatterplot panel based on the lag day of the best correlation
  # TODO: Add points for the monthly and annual values by colour and size
  p_ROFI_plume_scatter <- ggplot(df_zone, aes(x = ROFI_surface, y = plume_area)) + 
    # geom_hex(bins = 30, show.legend = FALSE) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1.5, linetype = "dashed", color = "orchid") +
    geom_smooth(method = "lm", linewidth = 1.5) +
    labs(title = "Scatterplot of ROFI surface vs plume surface",
         subtitle = "Purple line shows 1:1 ratio, blue line shows linear regression fit",
         y = "Plume Surface (km²)", x = "ROFI Surface (km²))") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df_zone$ROFI_surface, na.RM = TRUE))) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "bottom")
  # p_ROFI_plume_scatter
  
  # Plot lag results as individual plots for easier combinations
  p_ROFI_plume_cor_lag_daily <- ggplot(ROFI_plume_cor_daily, aes(x = lag, y = cor)) +
    geom_line() + geom_point() +
    labs(title = "Daily lag of plume after ROFI",
         x = "lagged days", y = "correlation (r)") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # p_ROFI_plume_cor_lag_daily
  p_ROFI_plume_cor_lag_monthly <- ggplot(ROFI_plume_cor_monthly, aes(x = lag, y = cor)) +
    geom_line() + geom_point() +
    labs(title = "Monthly lag of plume after ROFI",,
         x = "lagged months", y = "correlation (r)") +
    scale_x_continuous(breaks = seq(0, 12, 1), expand = c(0, 0)) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # p_ROFI_plume_cor_lag_monthly
  p_ROFI_plume_cor_lag_annual <- ggplot(ROFI_plume_cor_annual, aes(x = lag, y = cor)) +
    geom_line() + geom_point() +
    labs(title = "Annual lag of plume after ROFI",
         x = "lagged years", y = "correlation (r)") +
    scale_x_continuous(breaks = seq(0, 10, 1), expand = c(0, 0)) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  # p_ROFI_plume_cor_lag_annual
  
  # Combine plots and save
  plot_title <- grid::textGrob(paste0(df_zone$plot_title[1], " - ROFI vs river plume surface area (km^2)"), gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"))
  ts_plot <- ggpubr::ggarrange(p_ROFI_plume_ts, labels = c("a)"))
  # cor_plot_top <- ggpubr::ggarrange(p_ROFI_plume_scatter, p_ROFI_plume_cor_lag_daily, ncol = 1, nrow = 2, labels = c("b)", "c)"), heights = c(1, 0.3))
  # cor_plot_bot <- ggpubr::ggarrange(p_ROFI_plume_cor_lag_monthly, p_ROFI_plume_cor_lag_annual, ncol = 2, nrow = 1, labels = c("d)", "e)"))
  cor_plot_left <- ggpubr::ggarrange(p_ROFI_plume_scatter, labels = c("b)"))
  cor_plot_right_top <- ggpubr::ggarrange(p_ROFI_plume_cor_lag_daily, labels = c("c)"))
  cor_plot_right_bottom <- ggpubr::ggarrange(p_ROFI_plume_cor_lag_monthly, p_ROFI_plume_cor_lag_annual, ncol = 2, nrow = 1, labels = c("d)", "e)"))
  cor_plot_right <- ggpubr::ggarrange(cor_plot_right_top, cor_plot_right_bottom, ncol = 1, nrow = 2)
  cor_plot <- ggpubr::ggarrange(cor_plot_left, cor_plot_right, ncol = 2, nrow = 1, widths = c(1, 0.8))
  full_plot <- ggpubr::ggarrange(ts_plot, cor_plot, ncol = 1, nrow = 2, heights = c(0.8, 1))
  full_plot_title <- ggpubr::ggarrange(plot_title, full_plot, ncol = 1, nrow = 2, heights = c(0.05, 1)) + ggpubr::bgcolor("white")
  ggsave(filename = paste0("figures/cor_plot_ROFI_plume_",zone_name,".png"), plot = full_plot_title, width = 18, height = 18, dpi = 300)
}

# Run for all zones
# NB: No ROFI data for Gulf of Lion
walk(zones[1:3], comp_ROFI_plume)

