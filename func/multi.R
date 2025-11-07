# func/multi.R
# Loads all drivers of plume size, performs stats, plots results

# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(doParallel); doParallel::registerDoParallel(cores = 4)

# Zones
zones <- c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION")


# Functions ---------------------------------------------------------------

# Load all plume and driver data and perform stl
# zone <- zones[1]
multi_stl <- function(zone){
  
  # Determine meta-data based on zone
  if(zone == "BAY_OF_SEINE"){
    gauge = "LE_HAVRE"; mouth_name = "Seine"
    mouth_lon = 0.145; mouth_lat = 49.43
  } else if(zone == "BAY_OF_BISCAY"){
    gauge = "PORT-BLOC"; mouth_name = "Gironde"
    mouth_lon = -1.05; mouth_lat = 45.59
  } else if(zone == "SOUTHERN_BRITTANY"){
    gauge = "SAINT-NAZAIRE"; mouth_name = "Loire"
    mouth_lon = -2.10; mouth_lat = 47.29
  } else if(zone == "GULF_OF_LION"){
    gauge = "MARSEILLE"; mouth_name = "Grand Rhone"
    mouth_lon = 4.83; mouth_lat = 43.41
  } else {
    stop("Zone not recognised.")
  }
  
  # Load panache time series based on river mouth name
  df_plume <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,"/PLUME_DETECTION/Time_series_of_DAILY_plume_area_and_SPM_threshold.csv")) |> 
    dplyr::select(date:confidence_index_in_perc) |>
    complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
    dplyr::rename(plume_area = area_of_the_plume_mask_in_km2) |> 
    mutate(plume_area = ifelse(plume_area > 20000, NA, plume_area)) |> 
    zoo::na.trim()
  
  # Standardise last column name so it can be combined across sites
  colnames(df_plume)[ncol(df_plume)] <- "SPM_threshold"
  
  # Load river flow data
  df_river_flow <- load_river_flow(paste0("data/RIVER_FLOW/",zone))

  # Load tide data
  df_tide <- load_tide_gauge(paste0("data/TIDES/",gauge))
  
  # Load and subset wind data to + 0.5 N and E from mouth
  lon_round <- plyr::round_any(mouth_lon, 0.5); lat_round <- plyr::round_any(mouth_lat, 0.5)
  lon_range <- c(lon_round-0.5, lon_round+0.5); lat_range <- c(lat_round-0.5, lat_round+0.5); 
  wind_files <- dir(paste0("~/pCloudDrive/data/WIND/",zone), pattern = "_daily_", full.names = TRUE)
  df_wind <- purrr::map_dfr(wind_files, load_wind_sub, lon_range, lat_range)
  
  # Determine simple upwelling/downwelling index based on coastal direction (based on river mouth name)
  # TODO: Think of a more sophisticated way to do this
  if(zone %in% c("BAY_OF_BISCAY", "SOUTHERN_BRITTANY")){
    df_wind_welling <- df_wind |> 
      mutate(welling = ifelse(u < 0, "up", "down"))
  } else if (zone == "BAY_OF_SEINE"){
    df_wind_welling <- df_wind |> 
      mutate(welling = ifelse(v > 0, "up", "down"))
  } else if (zone == "GULF_OF_LION"){
    df_wind_welling <- df_wind |> 
      mutate(welling = ifelse(v < 0, "up", "down"))
  }
  
  # Calculate wind speed and direction
  # NB: wind_dir is where the wind is coming from, not going to
  df_wind_full <- df_wind_welling |> 
    mutate(wind_spd = round(sqrt(u^2 + v^2), 2),
           wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))
  
  # Combine all dataframes for further stats
  # NB: The trailing NAs are problematic...
  df_all <- left_join(df_plume, df_river_flow, by = "date") |> 
    left_join(df_tide, by = "date") |> 
    left_join(df_wind_full, by = "date") |> 
    mutate(plume_seas = stl_single(plume_area, out_col = "seas", start_date = min(df_plume$date)),
           plume_inter = stl_single(plume_area, out_col = "inter", start_date = min(df_plume$date)),
           flow_seas = stl_single(flow, out_col = "seas", start_date = min(df_plume$date)),
           flow_inter = stl_single(flow, out_col = "inter", start_date = min(df_plume$date)),
           tide_seas = stl_single(tide_range, out_col = "seas", start_date = min(df_plume$date)),
           tide_inter = stl_single(tide_range, out_col = "inter", start_date = min(df_plume$date)),
           wind_seas = stl_single(wind_spd, out_col = "seas", start_date = min(df_plume$date)),
           wind_inter = stl_single(wind_spd, out_col = "inter", start_date = min(df_plume$date))) |> 
    mutate(zone = zone, .before = "date")
  # print(ncol(df_all))

  # Exit
  return(df_all)
}

# Plot the results 
# df_stl <- stl_all
multi_plot <- function(df_stl){
  
  # Make pretty plot titles
  df_pretty <- df_stl |> 
    mutate(plot_title = case_when(zone == "BAY_OF_SEINE" ~ "Seine estuary",
                                  zone == "BAY_OF_BISCAY" ~ "Gironde estuary",
                                  zone == "SOUTHERN_BRITTANY" ~ "Loire estuary",
                                  zone == "GULF_OF_LION" ~ "Rhône estuary"), .after = "zone") |> 
    mutate(plot_title = factor(plot_title, levels = c("Seine estuary", "Loire estuary", "Gironde estuary", "Rhône estuary")))
  unique_years <- df_pretty$date |> year() |> unique()
  
  # One year of data for seasonal plots
  # df_mean <- df_pretty
  # TODO: Reduce to monthly points and take the min max as bands for a ribbon plot
  ## Calculate the mean values per site to bump up the seasonal clims correctly
  df_seas <- df_pretty |> 
    filter(year(date) == 1999) |> 
    mutate(month = month(date, label = TRUE, abbr = TRUE),
           doy = yday(date)) |> 
    dplyr::select(zone, plot_title, month, doy, plume_seas, flow_seas, tide_seas, wind_seas) |> 
    distinct()
  
  # Daily ts of river plumes
  plume_daily <- ggplot(data = df_pretty) + 
    # geom_point(aes(x = date, y = plume_area), color = "brown") + 
    geom_path(aes(x = date, y = plume_area), color = "brown") +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    # X-axis labels
    scale_x_date(name = "", 
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "Plume area (km²)") +
    ggplot_theme()
  ggsave(filename = "figures/plume_daily.png", plot = plume_daily, width = 24, height = 20, dpi = 300)
  
  # Seasonal ts of river plumes
  plume_seas <- ggplot(data = df_seas) + 
    # geom_point(aes(x = date, y = plume_area), color = "brown") + 
    geom_path(aes(x = doy, y = plume_seas), color = "brown", linewidth = 2) +
    # Facet
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    # X-axis labels
    # scale_x_date(name = "", 
    #              breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
    #              labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "Plume area (km²)") +
    ggplot_theme()
  ggsave(filename = "figures/plume_seas.png", plot = plume_seas, width = 24, height = 20, dpi = 300)
  
  # Interannual ts of river plumes
  plume_inter <- ggplot(data = df_pretty) + 
    # geom_point(aes(x = date, y = plume_area), color = "brown") + 
    geom_path(aes(x = date, y = plume_inter), color = "brown", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    # X-axis labels
    scale_x_date(name = "", 
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "Plume area (km²)") +
    ggplot_theme()
  ggsave(filename = "figures/plume_inter.png", plot = plume_inter, width = 24, height = 20, dpi = 300)
  

  # Get scaling factors for plotting
  # TODO: Multiple scalng of second y-axes does not appear to be possible
  ## Will need to create each facet individualy
  # scaling_factor_flow <- sec_axis_adjustement_factors(var_to_scale = df_pretty$flow_inter, 
  #                                                     var_ref = df_pretty$plume_inter)
  # scaling_factor_tide <- sec_axis_adjustement_factors(var_to_scale = df_pretty$tide_inter, 
  #                                                     var_ref = df_pretty$plume_inter)
  # scaling_factor_wind <- sec_axis_adjustement_factors(var_to_scale = df_pretty$wind_inter, 
  #                                                     var_ref = df_pretty$plume_inter)
  # 
  # df_scaling_flow <- summarise(df_pretty, sec_axis_adjustement_factors(flow_inter, plume_inter), .by = plot_title)
  # df_scale_flow <- left_join(df_pretty, df_scaling_flow, by = "plot_title") |> 
  #   mutate(flow_scaled = flow_inter * diff + adjust)
  #   
  # 
  # 
  #   mutate(flow_scaled = flow_inter * scaling_factor_flow$diff + scaling_factor_flow$adjust,
  #          tide_scaled = tide_inter * scaling_factor_tide$diff + scaling_factor_tide$adjust,
  #          wind_scaled = wind_inter * scaling_factor_wind$diff + scaling_factor_wind$adjust)
  # 
  # # Interannual plume vs river flow
  # ggplot(data = df_scale_flow) + 
  #   # Plume data
  #   geom_point(aes(x = date, y = plume_inter), color = "brown") + 
  #   geom_path(aes(x = date, y = plume_inter), color = "brown") + 
  #   # Wind data
  #   # geom_point(aes(x = date, y = flow_scaled), color = "blue") + 
  #   # geom_path(aes(x = date, y = flow_scaled), color = "blue") + 
  #   # Facet
  #   facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  #   # X-axis labels
  #   scale_x_date(name = "", 
  #                breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
  #                labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
  #   # Y-axis labels
  #   scale_y_continuous(name = "Plume area (km²)",
  #                      sec.axis = sec_axis(transform = ~ {. - df_scaling_flow$adjust} / df_scaling_flow$diff, 
  #                                          name = "River flow (m³/s)")) +
  #   # Extra bits
  #   ggplot_theme() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  #         plot.subtitle = element_text(hjust = 0.5),
  #         axis.text.y.left = element_text(color = "brown"), 
  #         axis.ticks.y.left = element_line(color = "brown"),
  #         axis.line.y.left = element_line(color = "brown"),
  #         axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
  #         axis.text.y.right = element_text(color = "blue"), 
  #         axis.ticks.y.right = element_line(color = "blue"),
  #         axis.line.y.right = element_line(color = "blue"),
  #         axis.title.y.right = element_text(color = "blue", margin = unit(c(0, 0, 0, 7.5), "mm")),
  #         panel.border = element_rect(linetype = "solid", fill = NA))
  
  # Seasonal ts of river flow
  flow_seas <- ggplot(data = df_seas) + 
    geom_path(aes(x = doy, y = flow_seas), color = "blue", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    scale_y_continuous(name = "River flow (m^3/s)") +
    ggplot_theme()
  ggsave(filename = "figures/flow_seas.png", plot = flow_seas, width = 24, height = 20, dpi = 300)
  
  # Interannual ts of river plumes
  flow_inter <- ggplot(data = df_pretty) + 
    # geom_point(aes(x = date, y = plume_area), color = "brown") + 
    geom_path(aes(x = date, y = flow_inter), color = "blue", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    # X-axis labels
    scale_x_date(name = "", 
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "River flow (m^3/s)") +
    ggplot_theme()
  ggsave(filename = "figures/flow_inter.png", plot = flow_inter, width = 24, height = 20, dpi = 300)
  
  # Seasonal ts of tide
  tide_seas <- ggplot(data = df_seas) + 
    geom_path(aes(x = doy, y = tide_seas), color = "green", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    scale_y_continuous(name = "River tide (m^3/s)") +
    ggplot_theme()
  ggsave(filename = "figures/tide_seas.png", plot = tide_seas, width = 24, height = 20, dpi = 300)
  
  # Interannual ts of tide
  tide_inter <- ggplot(data = df_pretty) + 
    # geom_point(aes(x = date, y = plume_area), color = "brown") + 
    geom_path(aes(x = date, y = tide_inter), color = "green", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    # X-axis labels
    scale_x_date(name = "", 
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "River tide (m^3/s)") +
    ggplot_theme()
  ggsave(filename = "figures/tide_inter.png", plot = tide_inter, width = 24, height = 20, dpi = 300)
  
  # Seasonal ts of wind
  wind_seas <- ggplot(data = df_seas) + 
    geom_path(aes(x = doy, y = wind_seas), color = "purple", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    scale_y_continuous(name = "River wind (m^3/s)") +
    ggplot_theme()
  ggsave(filename = "figures/wind_seas.png", plot = wind_seas, width = 24, height = 20, dpi = 300)
  
  # Interannual ts of wind
  wind_inter <- ggplot(data = df_pretty) + 
    # geom_point(aes(x = date, y = plume_area), color = "brown") + 
    geom_path(aes(x = date, y = wind_inter), color = "purple", linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
    # X-axis labels
    scale_x_date(name = "", 
                 breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                 labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
    # Y-axis labels
    scale_y_continuous(name = "River wind (m^3/s)") +
    ggplot_theme()
  ggsave(filename = "figures/wind_inter.png", plot = wind_inter, width = 24, height = 20, dpi = 300)
  
  # Everything on one plot
  df_all_scaled <- df_pretty |> 
    group_by(plot_title) |> 
    mutate(plum_scaled = plume_inter/max(plume_inter, na.rm = TRUE),
           flow_scaled = flow_inter/max(flow_inter, na.rm = TRUE),
           tide_scaled = tide_inter/max(tide_inter, na.rm = TRUE),
           wind_scaled = wind_inter/max(wind_inter, na.rm = TRUE)) |> 
    dplyr::select(plot_title, date, plum_scaled:wind_scaled) |> 
    pivot_longer(plum_scaled:wind_scaled)
  
  all_plot <- ggplot(df_all_scaled, aes(x = date, y = value)) +
    geom_path(aes(colour = name), linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1) +
    ggplot_theme()
  ggsave(filename = "figures/all_plot.png", plot = all_plot, width = 20, height = 20, dpi = 300)
}


# Run ---------------------------------------------------------------------

# Compute all STL stats and save
stl_all <- plyr::ldply(zones, multi_stl, .parallel = TRUE)
save(stl_all, file = "output/STATS/stl_all.RData")

# Create plots



# Missing data ------------------------------------------------------------

# Get missing dates of
SPM_files_NA <- data.frame(file_name = dir("data/SEXTANT/SPM/", pattern = ".nc", recursive = TRUE)) |> 
  mutate(base_name = basename(file_name)) |> 
  separate(base_name, "-", extra = "drop") |> 
  dplyr::rename(date = `-`) |> 
  mutate(date = as.Date(date, format = "%Y%m%d")) |> 
  complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
  filter(is.na(file_name))
write_csv(SPM_files_NA, "output/STATS/missing_SPM.csv")
chla_files_NA <- data.frame(file_name = dir("data/SEXTANT/CHLA/", pattern = ".nc", recursive = TRUE)) |> 
  mutate(base_name = basename(file_name)) |> 
  separate(base_name, "-", extra = "drop") |> 
  dplyr::rename(date = `-`) |> 
  mutate(date = as.Date(date, format = "%Y%m%d")) |> 
  complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
  filter(is.na(file_name))
write_csv(chla_files_NA, "output/STATS/missing_chla.csv")

# Filter down to missing days
SPM_files_NA_count <- SPM_files |>  
  mutate(year = year(date),
         month = month(date, label = TRUE, abbr = TRUE)) |> 
  summarise(miss_count_month_year = n(), .by = c("year", "month"))
chla_files_NA_count <- chla_files |> 
  mutate(year = year(date),
         month = month(date, label = TRUE, abbr = TRUE)) |> 
  summarise(miss_count_month_year = n(), .by = c("year", "month"))

# Plot
ggplot(SPM_files_NA_count, aes(x = month, y = miss_count_month_year)) +
  geom_col() +
  facet_wrap(~year) +
  labs(x = NULL, y = "count", title = "Monthly count of missing SPM SEXTANT files") +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
ggsave("figures/missng_SPM.png", width = 9, height = 9, dpi = 600)
ggplot(chla_files_NA_count, aes(x = month, y = miss_count_month_year)) +
  geom_col() +
  facet_wrap(~year) +
  labs(x = NULL, y = "count", title = "Monthly count of missing chl a SEXTANT files") +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
ggsave("figures/missng_chla.png", width = 9, height = 9, dpi = 600)

