# func/multi.R
# Loads all drivers of plume size, performs stats, plots results


# Analysis ideas ----------------------------------------------------------

# Fix y-axis scaling

# Create GIFs

# Basic time series comparisons to get seasonal and interannual comparisons
## Perform seasonal smoothing with heatwaveR
## also look into fixing the tidal range time series
## ultimately the point is to reduce the data in such a way that the time series can be related to one another in some way
## what does an extreme event analysis reveal?
### how do X11 and STL differ?
## Look at seasonal Trends per month, not long term
## Get correlations of interannual and seasonal time series

# Treat each pixel like its own time series and see what is happening with the forces when the pixel is triggered as a panache
## also how high SPM is while all this is happening
## show primary wind direction when pixel is triggered
## also relationship with SPM and tide range or category
## number of times pixel is flagged related to the size of the total panache when it is flagged
## could also tally the shape of the panache whenever pixel is flagged


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(seasonal)
library(patchwork)
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
    mutate(plot_title = case_when(zone == "BAY_OF_SEINE" ~ "Bay of Seine",
                                  zone == "SOUTHERN_BRITTANY" ~ "Southern Brittany",
                                  zone == "BAY_OF_BISCAY" ~ "Bay of Biscay",
                                  zone == "GULF_OF_LION" ~ "Gulf of Lion"), .after = "zone") |> 
    mutate(plot_title = factor(plot_title, levels = c("Bay of Seine", "Southern Brittany", "Bay of Biscay", "Gulf of Lion")))
  
  # One year of data for seasonal plots
  df_mean <- df_pretty |> 
    summarise(plume_mean = mean(plume_inter, na.rm = TRUE), 
              flow_mean = mean(flow_inter, na.rm = TRUE), 
              wind_mean = mean(wind_inter, na.rm = TRUE), 
              tide_mean = mean(tide_inter, na.rm = TRUE), .by = c(zone, plot_title))
  df_seas <- df_pretty |> 
    filter(year(date) == 1999) |> 
    mutate(month = month(date, label = TRUE, abbr = TRUE),
           doy = yday(date)) |> 
    dplyr::select(zone, plot_title, month, doy, plume_seas, flow_seas, tide_seas, wind_seas) |> 
    distinct() |> 
    left_join(df_mean, by = c("zone", "plot_title")) |>
    mutate(plume_seas = plume_seas + plume_mean,
           flow_seas = flow_seas + flow_mean,
           tide_seas = tide_seas + tide_mean,
           wind_seas = wind_seas + wind_mean)
  
  # Convenience wrappers for daily, seasonal, and interannual plot
  plot_daily <- function(df, y_col, line_colour, y_label, file_stub){
    unique_years <- df$date |> year() |> unique()
    pl_daily <- ggplot(data = df) + 
      geom_path(aes_string(x = "date", y = y_col), color = line_colour) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_date(name = "", expand = c(0,0),
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
      scale_y_continuous(name = y_label) +
      labs( x = NULL) +
      ggplot_theme()
    ggsave(filename = paste0("figures/",file_stub,"_daily.png"), plot = pl_daily, width = 24, height = 24, dpi = 300)
    # pl_daily
  }
  plot_seas <- function(df, y_col, line_colour, y_label, file_stub){
    df_sub <- df[,c("plot_title", "month", y_col)]
    colnames(df_sub)[3] = "val"
    df_sub <- df_sub |> 
      summarise(val_min = min(val, na.rm = TRUE),
                val_mean = mean(val, na.rm = TRUE),
                val_max = max(val, na.rm = TRUE), .by = c("plot_title", "month")) |> 
      mutate(month_int = as.integer(month))
    pl_seas <- ggplot(data = df_sub, aes(x = month_int)) + 
      geom_ribbon(aes(ymin = val_min, ymax = val_max), fill = line_colour, alpha = 0.3) +
      geom_path(aes(y = val_mean), color = line_colour, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_y_continuous(name = y_label) +
      scale_x_continuous(expand = c(0, 0), breaks = 1:12, labels = month.abb) +
      labs(x = NULL) +
      ggplot_theme()
    ggsave(filename = paste0("figures/",file_stub,"_seas.png"), plot = pl_seas, width = 24, height = 24, dpi = 300)
    # pl_seas
  }
  plot_inter <- function(df, y_col, line_colour, y_label, file_stub){
    unique_years <- df$date |> year() |> unique()
    pl_inter <- ggplot(data = df) + 
      geom_path(aes_string(x = "date", y = y_col), color = line_colour, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_date(name = "", 
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(), 
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) +
      scale_y_continuous(name = y_label) +
      ggplot_theme()
    ggsave(filename = paste0("figures/",file_stub,"_inter.png"), plot = pl_inter, width = 24, height = 24, dpi = 300)
    # pl_inter
  }
  
  # Daily time series
  plot_daily(df_pretty, "plume_area", "brown", "Plume area (km^2)", "plume")
  
  # Seasonal time series
  plot_seas(df_seas, "plume_seas", "brown", "Plume area (km^2)", "plume")
  plot_seas(df_seas, "flow_seas", "blue", "River flow (m^3 s-1)", "flow")
  plot_seas(df_seas, "tide_seas", "darkgreen", "Tidal range (m)", "tide")
  plot_seas(df_seas, "wind_seas", "purple", "Wind speed (m s-1)", "wind")
  
  # Interannual time series
  plot_inter(df_pretty, "plume_inter", "brown", "Plume area (km^2)", "plume")
  plot_inter(df_pretty, "flow_inter", "blue", "River flow (m^3 s-1)", "flow")
  plot_inter(df_pretty, "tide_inter", "darkgreen", "Tidal range (m)", "tide")
  plot_inter(df_pretty, "wind_inter", "purple", "Wind speed (m s-1)", "wind")

  # Comparison plots
  # df <- df_pretty; var_1 <- "plume_inter"; var_2 <- "flow_inter"
  # df <- df_seas; var_1 <- "plume_seas"; var_2 <- "flow_seas"
  # colour_1 <- "brown"; colour_2 <- "blue"; label_1 <- "Plume area (km^2)"; label_2 <- "River flow (m^3 s-1)"; file_stub <- "comparison_plume_flow_inter"
  comparison_plot <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2){
    
    if(grepl("seas", var_1)){
      df_sub <- df[,c("plot_title", "month", var_1, var_2)]
      colnames(df_sub) <- c("plot_title", "month", "var_1", "var_2")
    } else {
      df_sub <- df[,c("plot_title", "date", var_1, var_2)]
      colnames(df_sub) <- c("plot_title", "date", "var_1", "var_2")
    }
    
    # Scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_sub$var_2, df_sub$var_1)
    df_scaling <- summarise(df_sub, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_sub, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    
    # Plot base
    if(grepl("seas", var_1)){
      
      # Get range for ribbon plot
      df_scale_sub <- df_scale |> 
        summarise(var_1_min = min(var_1, na.rm = TRUE),
                  var_1_mean = mean(var_1, na.rm = TRUE),
                  var_1_max = max(var_1, na.rm = TRUE),
                  var_2_min = min(var_2_scaled, na.rm = TRUE),
                  var_2_mean = mean(var_2_scaled, na.rm = TRUE),
                  var_2_max = max(var_2_scaled, na.rm = TRUE), .by = c("plot_title", "month")) |> 
        mutate(month_int = as.integer(month))
      
      # Plot them
      pl_base <- ggplot(data = df_scale_sub, aes(x = month_int)) + 
        # Var 1
        geom_ribbon(aes(ymin = var_1_min, ymax = var_1_max), fill = colour_1, alpha = 0.3) +
        geom_path(aes(y = var_1_mean), color = colour_1, linewidth = 2) +
        # Var 2
        geom_ribbon(aes(ymin = var_2_min, ymax = var_2_max), fill = colour_2, alpha = 0.3) +
        geom_path(aes(y = var_2_mean), color = colour_2, linewidth = 2) +
        facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
        scale_x_continuous(expand = c(0, 0), breaks = 1:12, labels = month.abb)
    } else {
      unique_years <- df_scale$date |> year() |> unique()
      pl_base <- ggplot(data = df_scale) +
        # Var 1 data
        geom_point(aes(x = date, y = var_1), color = colour_1) +
        geom_path(aes(x = date, y = var_1), color = colour_1) +
        # Var 2 data
        geom_point(aes(x = date, y = var_2_scaled), color = colour_2) +
        geom_path(aes(x = date, y = var_2_scaled), color = colour_2) +
        # Facet
        facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
        # X-axis labels
        scale_x_date(name = "", expand = c(0, 0),
                     breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                     labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) 
    }
    
    # Finish up the comparison plot
    pl_comp <- pl_base +
      # Y-axis labels
      scale_y_continuous(name = label_1,
                         sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
                                             name = label_2)) +
      labs( x = NULL) +
      # Extra bits
      ggplot_theme() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.y.left = element_text(color = colour_1),
            axis.ticks.y.left = element_line(color = colour_1),
            axis.line.y.left = element_line(color = colour_1),
            axis.title.y.left = element_text(color = colour_1, margin = unit(c(0, 7.5, 0, 0), "mm")),
            axis.text.y.right = element_text(color = colour_2),
            axis.ticks.y.right = element_line(color = colour_2),
            axis.line.y.right = element_line(color = colour_2),
            axis.title.y.right = element_text(color = colour_2, margin = unit(c(0, 0, 0, 7.5), "mm")),
            panel.border = element_rect(linetype = "solid", fill = NA))
    return(pl_comp)
  }
  
  # Convenience wrapper to run and save comparison plots
  # NB: this is hqrd coded to work with four plots
  comparison_plot_save <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2, file_stub){
    comp_list <- plyr::dlply(df, c("zone"), comparison_plot, var_1 = var_1, var_2 = var_2, 
                             colour_1 = colour_1, colour_2 = colour_2, label_1 = label_1, label_2 = label_2)
    comp_fig <- comp_list[[2]] + comp_list[[4]] + comp_list[[1]] + comp_list[[3]] + plot_layout(ncol = 1, axes = "collect")
    ggsave(filename = paste0("figures/",file_stub,".png"), plot = comp_fig, width = 24, height = 24, dpi = 300)
  }
  
  # Seasonal comparison plots
  comparison_plot_save(df_seas, "plume_seas", "flow_seas", "brown", "blue", "Plume area (km^2)", "River flow (m^3 s-1)", "comparison_plume_flow_seas")
  comparison_plot_save(df_seas, "plume_seas", "wind_seas", "brown", "purple", "Plume area (km^2)", "Wind speed (m s-1)", "comparison_plume_wind_seas")
  comparison_plot_save(df_seas, "plume_seas", "tide_seas", "brown", "darkgreen", "Plume area (km^2)", "Tidal range (m)", "comparison_plume_tide_seas")
  
  # Interannual comparison plots
  comparison_plot_save(df_pretty, "plume_inter", "flow_inter", "brown", "blue", "Plume area (km^2)", "River flow (m^3 s-1)", "comparison_plume_flow_inter")
  comparison_plot_save(df_pretty, "plume_inter", "tide_inter", "brown", "darkgreen", "Plume area (km^2)", "Tidal range (m)", "comparison_plume_tide_inter")
  comparison_plot_save(df_pretty, "plume_inter", "wind_inter", "brown", "purple", "Plume area (km^2)", "Wind speed (m s-1)", "comparison_plume_wind_inter")
  
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
# stl_all <- plyr::ldply(zones, multi_stl, .parallel = TRUE)
# save(stl_all, file = "output/STATS/stl_all.RData")
load("output/STATS/stl_all.RData")

# Create plots
multi_plot(stl_all)


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


# Decomposition comparison ------------------------------------------------

# X11

# STL

# heatwaveR

