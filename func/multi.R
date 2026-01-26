# func/multi.R
# Loads all drivers of plume size, performs stats, plots results


# Analysis ideas ----------------------------------------------------------

# Create GIFs

# Basic time series comparisons to get seasonal and interannual comparisons
## Perform seasonal smoothing with heatwaveR
## also look into fixing the tidal range time series
## ultimately the point is to reduce the data in such a way that the time series can be related to one another in some way
## what does an extreme event analysis reveal?
### how do X11 and STL differ?
## Look at seasonal Trends per month, not long term
## Get correlations of interannual and seasonal time series
## Dynamic Linear Models (DLM) is another option: library(dlm)

# Treat each pixel like its own time series and see what is happening with the forces when the pixel is triggered as a panache
## also how high SPM is while all this is happening
## show primary wind direction when pixel is triggered
## also relationship with SPM and tide range or category
## number of times pixel is flagged related to the size of the total panache when it is flagged
### Would need to relate wind with time lag to this as well
## could also tally the shape of the panache whenever pixel is flagged

# Other analyses
## Get mean offshore distance of centroid
## when creating GIFs, would be cool to have the centroid visualised as a 21 dot
### it could leave a trail of 1 dots behind it per day
### or fill colour could be left to show tide, wind, etc. on that day

# More ideas
## nmds of mean characteristics during plume events
## extreme event analysis
### need to establish a reasonable baseline, and go from there
### could be interesting to use the time varying X11 seasonality
### but then is it relevant to calculate event stats from the median signal?
  ### one could modify the analysis to always take from the base, being zero
### how would one establish the 90th percentile? Or rather, just always take everything over the seasonal signal
## EMD - see Vincent email
### percent contribution of each component to time series
### changes over time as well
## like a CTD cast, figure out a way to measure from when a plume peaks and then goes down below a certain threshold as a way of determining if it is an individual event
### and then get statistics from that
### also add up number of days with onshore wind, neap tide, etc.
### the spatial threshold could be based on the percentile of the total plume area over the full time series
### start searching by creating a contour plot of every 10th percentile
## also account for the lon/lat of the centroid and how that relates to the other drivers
# be able to say if the plume is increasing or decreasing in size so that the drivers on that day can be categorised under what the plume is doing

# 


# Libraries ---------------------------------------------------------------

source("func/util.R")
library(tidyverse)
library(tidync)
library(heatwaveR) # For seasonal smoothing analysis
library(seasonal) # For X11 analysis (currently not used)
library(RcppRoll) # For running means to get STL interannual signals closer to X11
library(patchwork)
library(doParallel); doParallel::registerDoParallel(cores = 14)

# Zones
zones <- c("BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION")


# STL ---------------------------------------------------------------------

# Load all plume and driver data and perform stl
# zone <- zones[4]
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
    df_wind_updown <- df_wind |> 
      mutate(direction = ifelse(u < 0, "off", "on"))
  } else if (zone == "BAY_OF_SEINE"){
    df_wind_updown <- df_wind |> 
      mutate(direction = ifelse(v > 0, "off", "on"))
  } else if (zone == "GULF_OF_LION"){
    df_wind_updown <- df_wind |> 
      mutate(direction = ifelse(v < 0, "off", "on"))
  }
  
  # Calculate wind speed and direction
  # NB: wind_dir is where the wind is coming from, not going to
  df_wind_full <- df_wind_updown |> 
    mutate(wind_spd = round(sqrt(u^2 + v^2), 2),
           wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))
  
  # Combine all dataframes for further stats
  # NB: The trailing NAs are problematic...
  df_all <- left_join(df_plume, df_river_flow, by = "date") |> 
    left_join(df_tide, by = "date") |> 
    left_join(df_wind_full, by = "date") |> 
    mutate(plume_seas = stl_single(plume_area, out_col = "seas", start_date = min(df_plume$date)),
           plume_inter = stl_single(plume_area, out_col = "inter", start_date = min(df_plume$date)),
           plume_resid = stl_single(plume_area, out_col = "remain", start_date = min(df_plume$date)),
           flow_seas = stl_single(flow, out_col = "seas", start_date = min(df_plume$date)),
           flow_inter = stl_single(flow, out_col = "inter", start_date = min(df_plume$date)),
           flow_resid = stl_single(flow, out_col = "remain", start_date = min(df_plume$date)),
           tide_seas = stl_single(tide_range, out_col = "seas", start_date = min(df_plume$date)),
           tide_inter = stl_single(tide_range, out_col = "inter", start_date = min(df_plume$date)),
           tide_resid = stl_single(tide_range, out_col = "remain", start_date = min(df_plume$date)),
           wind_seas = stl_single(wind_spd, out_col = "seas", start_date = min(df_plume$date)),
           wind_inter = stl_single(wind_spd, out_col = "inter", start_date = min(df_plume$date)),
           wind_resid = stl_single(wind_spd, out_col = "remain", start_date = min(df_plume$date))) |> 
    mutate(zone = zone, .before = "date")
  # print(ncol(df_all))

  # Exit
  return(df_all)
}

# Compute all STL stats and save
stl_all <- plyr::ldply(zones, multi_stl, .parallel = TRUE)
save(stl_all, file = "output/STATS/stl_all.RData")


# heatwaveR ---------------------------------------------------------------

# Load data
plume_clim <- map_dfr(zones, plume_clim_calc)

  
# Multi-driver comparison -------------------------------------------------

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
    colnames(df)[which(colnames(df) == y_col)] <- "value"
    df_sub <- df |> 
      dplyr::select(plot_title, date, value) |>
      mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
      # mutate(date = round_date(date, unit = "months")) |>
      filter(date >= min(df$date)) |> 
      group_by(plot_title, date) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "keep") |> 
      group_by(plot_title) |>
      mutate(running_mean = roll_mean(value, n = 48, fill = NA, align = "center")) |> 
      ungroup()
    pl_inter <- ggplot(data = df_sub) + 
      geom_path(aes(x = date, y = running_mean), color = line_colour, linewidth = 2) +
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
    
    # Plot base
    if(grepl("seas", var_1)){
      
      # Scaling factor
      scaling_factor <- sec_axis_adjustement_factors(df_sub$var_2, df_sub$var_1)
      df_scaling <- summarise(df_sub, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
      df_scale <- left_join(df_sub, df_scaling, by = "plot_title") |>
        mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
      
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
      
      # Perform rolling mean
      df_roll_mean <- df_sub |>
        mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
        # mutate(date = round_date(date, unit = "months")) |>
        filter(date >= min(df$date)) |>
        group_by(plot_title, date) |>
        summarise(var_1 = mean(var_1, na.rm = TRUE),
                  var_2 = mean(var_2, na.rm = TRUE), .groups = "keep") |>
        group_by(plot_title) |>
        mutate(var_1 = roll_mean(var_1, n = 48, fill = NA, align = "center"),
               var_2 = roll_mean(var_2, n = 48, fill = NA, align = "center")) |>
        ungroup()
      
      # Then get the scaling factor
      scaling_factor <- sec_axis_adjustement_factors(df_roll_mean$var_2, df_roll_mean$var_1)
      df_scaling <- summarise(df_roll_mean, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
      df_scale <- left_join(df_roll_mean, df_scaling, by = "plot_title") |>
        mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
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
    mutate(plum_scaled = plum_scaled/mean(plum_scaled, na.rm = TRUE),
           flow_scaled = flow_scaled/mean(flow_scaled, na.rm = TRUE),
           tide_scaled = tide_scaled/mean(tide_scaled, na.rm = TRUE),
           wind_scaled = wind_scaled/mean(wind_scaled, na.rm = TRUE)) |>
    dplyr::select(plot_title, date, plum_scaled:wind_scaled) |> 
    pivot_longer(plum_scaled:wind_scaled) |> 
    mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
    # mutate(date = round_date(date, unit = "months")) |>
    filter(date >= min(df_pretty$date)) |> 
    group_by(plot_title, name, date) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "keep") |> 
    group_by(plot_title, name) |>
    mutate(running_mean = roll_mean(value, n = 48, fill = NA, align = "center")) |> 
    ungroup()
  
  all_plot <- ggplot(df_all_scaled, aes(x = date, y = running_mean)) +
    geom_path(aes(colour = name), linewidth = 2) +
    facet_wrap(~plot_title, ncol = 1) +
    ggplot_theme()
  ggsave(filename = "figures/all_plot.png", plot = all_plot, width = 20, height = 20, dpi = 300)
}

# Load STL calculated above
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
load_X11 <- function(zone, type = "plume"){
  if(type == "plume"){
    file_stub = "/X11_ANALYSIS/area_of_the_plume_mask_in_km2/SEXTANT_merged_Standard_WEEKLY.csv"
  } else {
    file_stub = "/X11_ANALYSIS/river_flow/River_flow___WEEKLY.csv"
  }
  suppressMessages(
  df <- read_csv(paste0("output/FIXED_THRESHOLD/",zone,file_stub)) |> 
    mutate(zone = zone, .before = "dates") |> 
    dplyr::rename(date = dates) |> 
    dplyr::select(zone:Residual_signal) |> 
    rename_with(~ paste0(.x, "_X11_",type), everything()) 
  )
  colnames(df)[1:2] <- c("zone", "date")
  return(df)
}
X11_plume <- map_dfr(zones, load_X11)
X11_flow <- map_dfr(zones, load_X11, type = "flow")

# STL
load("output/STATS/stl_all.RData")
stl_sub <- dplyr::select(stl_all, zone, date, flow, tide_range, direction:wind_resid) |> 
  rename_with(~ paste0(.x, "_STL"), everything())
colnames(stl_sub)[1:7] <- c("zone", "date", "flow", "tide_range", "direction", "wind_spd", "wind_dir")

# heatwaveR
plume_clim <- map_dfr(zones, plume_clim_calc) |> 
  # mutate(plume_inter = plume_area - plume_seas) |> 
  group_by(zone) |> 
  mutate(plume_seas = plume_seas - mean(plume_seas)) |> 
  ungroup() |> 
  mutate(plume_inter = plume_area - plume_seas)

# Combine into one big df
decomp_df <- left_join(stl_sub, X11_plume, by = c("zone", "date")) |> 
  left_join(X11_flow, by = c("zone", "date")) |> 
  mutate(plot_title = case_when(zone == "BAY_OF_SEINE" ~ "Bay of Seine",
                                zone == "SOUTHERN_BRITTANY" ~ "Southern Brittany",
                                zone == "BAY_OF_BISCAY" ~ "Bay of Biscay",
                                zone == "GULF_OF_LION" ~ "Gulf of Lion"), .after = "zone") |> 
  mutate(plot_title = factor(plot_title, 
                             levels = c("Bay of Seine", "Southern Brittany", "Bay of Biscay", "Gulf of Lion"))) |> 
  mutate(year = year(date),
         doy = yday(date))
decomp_df$doy <- mapply(adjust_doy, decomp_df$year, decomp_df$doy)
decomp_df <- left_join(decomp_df, plume_clim, by = c("zone", "date", "doy"))


## Plume comparison --------------------------------------------------------

### Interannual -------------------------------------------------------------

# Extract and prep the interannual values
plume_inter <- decomp_df |> 
  dplyr::select(zone, plot_title, date, plume_inter_STL, Interannual_signal_X11_plume, plume_inter) |> 
  pivot_longer(cols = c(plume_inter_STL, Interannual_signal_X11_plume, plume_inter)) |> 
  mutate(name = case_when(name == "plume_inter_STL" ~ "plume STL",
                          name == "Interannual_signal_X11_plume" ~ "plume X11",
                          name == "plume_inter" ~ "plume smooth")) |> 
  filter(!is.na(value))

# Line plot of plume size - interannual
line_plume_inter <- ggplot(plume_inter, aes(x = date, y = value)) +
  geom_path(aes(colour = name), linewidth = 2) +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = NULL, y = "Plume area (km^2)") +
  # scale_colour_manual(key_glyph = "point") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of plume size - interannual
scatter_plume_inter <- plume_inter |> 
  pivot_wider(names_from = name, values_from = value) |> 
  mutate(month = month(date, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `plume X11`, y = `plume STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = month)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = NULL) +
  # coord_cartesian(ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_plume_inter <- line_plume_inter + scatter_plume_inter + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_plume_inter_comp.png", plot = multi_plume_inter, height = 20, width = 30)


### Seasonal ----------------------------------------------------------------

# Extract and prep the monthly values
plume_seas <- decomp_df |> 
  dplyr::select(zone, plot_title, date, doy, plume_seas_STL, Seasonal_signal_X11_plume, plume_seas) |> 
  pivot_longer(cols = c(plume_seas_STL, Seasonal_signal_X11_plume, plume_seas)) |> 
  mutate(name = case_when(name == "plume_seas_STL" ~ "plume STL",
                          name == "Seasonal_signal_X11_plume" ~ "plume X11",
                          name == "plume_seas" ~ "plume smooth")) |> 
  filter(!is.na(value))
  
# Line plot of plume size - seas
line_plume_seas <- ggplot(plume_seas, aes(x = date, y = value)) +
  geom_path(aes(colour = name), alpha = 0.8, linewidth = 2) +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = NULL, y = "Plume area (km^2)") +
  # scale_colour_manual(key_glyph = "point") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of plume size - seas
scatter_plume_seas <- plume_seas |> 
  pivot_wider(names_from = name, values_from = value) |> 
  mutate(month = month(date, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `plume X11`, y = `plume STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = month)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = NULL) +
  # coord_cartesian(ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_plume_seas <- line_plume_seas + scatter_plume_seas + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_plume_seas_comp.png", plot = multi_plume_seas, height = 20, width = 30)


### DOY ---------------------------------------------------------------------

# get the average doy values
plume_doy <- plume_seas |> 
  summarise(val_min = min(value, na.rm = TRUE),
            val_mean = mean(value, na.rm = TRUE),
            val_max = max(value, na.rm = TRUE), 
            .by = c("zone", "plot_title", "name", "doy")) |> 
  arrange(doy)

# Ribbon plot of plume - doy
line_plume_doy <- ggplot(plume_doy, aes(x = doy, y = val_mean)) +
  geom_ribbon(aes(fill = name, ymin = val_min, ymax = val_max), alpha = 0.2, show.legend = FALSE) +
  geom_path(aes(colour = name), linewidth = 2)  +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = "day-of-year", y = "Plume area (km^2)") +
  # scale_colour_manual(key_glyph = "point") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of plume size - doy
scatter_plume_doy <- plume_doy |> 
  dplyr::select(-val_min, -val_max) |> 
  pivot_wider(names_from = name, values_from = val_mean) |> 
  # mutate(month = month(doy_adj, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `plume X11`, y = `plume STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = doy), size = 7) +
  scale_colour_viridis_c() +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = "day-of-year") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_plume_doy <- line_plume_doy + scatter_plume_doy + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_plume_doy_comp.png", plot = multi_plume_doy, height = 20, width = 30)


### Residual ----------------------------------------------------------------

# Extract and prep the interannual values
plume_resid <- decomp_df |> 
  dplyr::select(zone, plot_title, date, plume_resid_STL, Residual_signal_X11_plume) |> 
  pivot_longer(cols = c(plume_resid_STL, Residual_signal_X11_plume)) |> 
  mutate(name = case_when(name == "plume_resid_STL" ~ "plume STL",
                          name == "Residual_signal_X11_plume" ~ "plume X11")) |> 
  filter(!is.na(value))

# Line plot of plume size - residual
line_plume_resid <- ggplot(plume_resid, aes(x = date, y = value)) +
  geom_path(aes(colour = name), linewidth = 2) +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = NULL, y = "Plume area (km^2)") +
  # scale_colour_manual(key_glyph = "point") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of plume size - residual
scatter_plume_resid <- plume_resid |> 
  pivot_wider(names_from = name, values_from = value) |> 
  mutate(month = month(date, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `plume X11`, y = `plume STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = month)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = NULL) +
  # coord_cartesian(ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_plume_resid <- line_plume_resid + scatter_plume_resid + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_plume_resid_comp.png", plot = multi_plume_resid, height = 20, width = 30)


### Trends ------------------------------------------------------------------

# Quickly show the trend line and stats for all sites
plume_trend_daily_all <- decomp_df |> 
  dplyr::select(zone, plot_title, date, plume_area)
## Monthly
plume_trend_monthly_all <- plume_trend_daily_all |> 
  mutate(date = round_date(date, "month") + days(14)) |> 
  summarise(plume_area = mean(plume_area, na.rm = TRUE), .by = c("zone", "plot_title", "date"))
## Annual
plume_trend_annual_all <- plume_trend_daily_all |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
  summarise(plume_area = mean(plume_area, na.rm = TRUE), .by = c("zone", "plot_title", "date"))

# Calculate trend 
trend_daily_all <- plume_trend_daily_all |> 
  summarise(slope = coef(lm(plume_area ~ date))["date"] * 365.25, .by = c("zone", "plot_title"))
trend_monthly_all <- plume_trend_monthly_all |> 
  summarise(slope = coef(lm(plume_area ~ date))["date"] * 365.25, .by = c("zone", "plot_title"))
trend_annual_all <- plume_trend_annual_all |> 
  summarise(slope = coef(lm(plume_area ~ date))["date"] * 365.25, .by = c("zone", "plot_title"))
trend_labels_all <- rbind(trend_daily_all, trend_monthly_all, trend_annual_all) |> 
  mutate(time_step = rep(c("daily", "monthly", "annual"), each = nrow(trend_daily_all)), .before = "slope") |> 
  arrange(plot_title, time_step) |>
  mutate(x = as.Date("1998-01-01"),
         y = c(2300, 2100, 1900, 
               11000, 10000, 9000, 
               15000, 13800, 12600, 
               5500, 5000, 4500))

# Plot all sites as facets
line_trend_base_all <- ggplot(plume_trend_daily_all, aes(x = date, y = plume_area)) +
  geom_point(alpha = 0.5, aes(colour = "daily")) +
  geom_smooth(method = "lm", colour = "black", linewidth = 2) +
  geom_point(data = plume_trend_monthly_all, size = 4, alpha = 0.5, aes(colour = "monthly")) +
  geom_smooth(data = plume_trend_monthly_all, method = "lm", colour = "darkblue", linewidth = 2) +
  geom_point(data = plume_trend_annual_all, size = 7, alpha = 0.5, aes(colour = "annual")) +
  geom_smooth(data = plume_trend_annual_all, method = "lm", colour = "blue", linewidth = 2) +
  geom_label(data = trend_labels_all, show.legend = FALSE,
             aes(x = x, y = y, colour = time_step, size = 6, hjust = 0,
                 label = paste0(time_step," data slope = ", round(slope, 2), " km^2 yr-1", sep = ""))) +
  labs(x = NULL, y = "Plume area (km^2)", title = "Trend for daily, monthly, and annual mean plume data") +
  scale_color_manual(name = "Time step",
                     values = c("daily" = "black", "monthly" = "darkblue", "annual" = "blue"),
                     breaks = c("daily", "monthly", "annual")) +
  facet_wrap(~plot_title, ncol = 2, scales = "free_y") +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_trend_comparison_base_all.png", plot = line_trend_base_all, height = 9, width = 12)


# Here we look at only Gulf of Lion to keep the output simpler

# First we start with a linear trend of the raw time series, and another plot with the linear analysis per month
## Daily data
plume_trend_daily <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |>
  dplyr::select(zone, plot_title, date, plume_area)
## Monthly
plume_trend_monthly <- plume_trend_daily |> 
  mutate(date = round_date(date, "month") + days(14)) |> 
  summarise(plume_area = mean(plume_area, na.rm = TRUE), .by = c("zone", "plot_title", "date"))
## Annual
plume_trend_annual <- plume_trend_daily |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
  summarise(plume_area = mean(plume_area, na.rm = TRUE), .by = c("zone", "plot_title", "date"))
# Fit the linear models and extract slopes
trend_daily <- coef(lm(plume_area ~ date, data = plume_trend_daily))["date"] * 365.25
trend_monthly <- coef(lm(plume_area ~ date, data = plume_trend_monthly))["date"] * 365.25
trend_annual <- coef(lm(plume_area ~ date, data = plume_trend_annual))["date"] * 365.25

# Create a little dataframe for plotting the labels
trend_labels <- data.frame(
  time_step = c("daily", "monthly", "annual"),
  slope = c(trend_daily, trend_monthly, trend_annual),
  x = as.Date(c("1998-01-01", "1998-01-01", "1998-01-01")),
  y = c(5500, 5000, 4500)
)

# Plot the three dataframes as time series with linear trends on the same panel
line_trend_base <- ggplot(plume_trend_daily, aes(x = date, y = plume_area)) +
  geom_point(alpha = 0.5, aes(colour = "daily")) +
  geom_smooth(method = "lm", colour = "black", linewidth = 2) +
  geom_point(data = plume_trend_monthly, size = 4, alpha = 0.5, aes(colour = "monthly")) +
  geom_smooth(data = plume_trend_monthly, method = "lm", colour = "darkblue", linewidth = 2) +
  geom_point(data = plume_trend_annual, size = 7, alpha = 0.5, aes(colour = "annual")) +
  geom_smooth(data = plume_trend_annual, method = "lm", colour = "blue", linewidth = 2) +
  geom_label(data = trend_labels, aes(x = x, y = y, label = paste0(time_step," data slope = ", round(slope, 2), " km^2 yr-1", sep = "")),
             colour = c("black", "darkblue", "blue"), size = 6, hjust = 0) +
  labs(x = NULL, y = "Plume area (km^2)", title = "Trend for daily, monthly, and annual mean plume data") +
  scale_color_manual(name = "Time step",
                     values = c("daily" = "black", "monthly" = "darkblue", "annual" = "blue"),
                     breaks = c("daily", "monthly", "annual")) +
  # ggplot_theme() +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_trend_comparison_base.png", plot = line_trend_base, height = 6, width = 12)

# The same plot comparing the seasonal components
plume_seas <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  dplyr::select(zone, plot_title, date, plume_seas_STL, Seasonal_signal_X11_plume, plume_seas) |> 
  pivot_longer(cols = c(plume_seas_STL, Seasonal_signal_X11_plume, plume_seas)) |> 
  mutate(name = case_when(name == "plume_seas_STL" ~ "STL seasonal",
                          name == "Seasonal_signal_X11_plume" ~ "X11 seasonal",
                          name == "plume_seas" ~ "smoothed seasonal")) |> 
  filter(!is.na(value),
         date <= as.Date("2025-02-14"),
         date >= as.Date("1999-02-15"))

# Create small dataframe of linear model slopes by name
trend_seas_STL <- coef(lm(value ~ date, data = filter(plume_seas, name == "STL seasonal")))["date"] * 365.25
trend_seas_X11 <- coef(lm(value ~ date, data = filter(plume_seas, name == "X11 seasonal")))["date"] * 365.25
trend_seas_smooth <- coef(lm(value ~ date, data = filter(plume_seas, name == "smoothed seasonal")))["date"] * 365.25
trend_seas_labels <- data.frame(
  name = c("STL seasonal", "X11 seasonal", "smoothed seasonal"),
  slope = c(trend_seas_STL, trend_seas_X11, trend_seas_smooth),
  x = as.Date(c("1998-01-01", "1998-01-01", "1998-01-01")),
  y = c(1500, 1300, 1100)
)

line_trend_seas <- ggplot(plume_seas, aes(x = date, y = value)) +
  geom_path(alpha = 0.5, aes(colour = name)) +
  geom_smooth(method = "lm", linewidth = 2, aes(colour = name)) +
  labs(x = NULL, y = "Plume area (km^2)", title = "Trend for plume seasonal components") +
  geom_label(data = trend_seas_labels, aes(x = x, y = y, label = paste0(name," data slope = ", round(slope, 2), " km^2 yr-1", sep = "")),
             colour = c("turquoise4", "chartreuse4", "indianred4"), size = 6, hjust = 0) +
  scale_color_manual(name = "Decomposition method",
                     values = c("STL seasonal" = "turquoise4", "X11 seasonal" = "chartreuse4", "smoothed seasonal" = "indianred4"),
                     breaks = c("STL seasonal", "X11 seasonal", "smoothed seasonal")) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_trend_comparison_seas.png", plot = line_trend_seas, height = 6, width = 12)

# The same plot comparing the residual components
plume_resid <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  dplyr::select(zone, plot_title, date, plume_resid_STL, Residual_signal_X11_plume) |> 
  pivot_longer(cols = c(plume_resid_STL, Residual_signal_X11_plume)) |> 
  mutate(name = case_when(name == "plume_resid_STL" ~ "STL residual",
                          name == "Residual_signal_X11_plume" ~ "X11 residual")) |> 
  filter(!is.na(value)) |> 
  filter(date <= as.Date("2024-12-31"),
         date >= as.Date("1999-01-01"))

# Create small dataframe of linear model slopes by name
trend_resid_STL <- coef(lm(value ~ date, data = filter(plume_resid, name == "STL residual")))["date"] * 365.25
trend_resid_X11 <- coef(lm(value ~ date, data = filter(plume_resid, name == "X11 residual")))["date"] * 365.25
trend_resid_labels <- data.frame(
  name = c("STL residual", "X11 residual"),
  slope = c(trend_resid_STL, trend_resid_X11),
  x = as.Date(c("1998-01-01", "1998-01-01")),
  y = c(4500, 4000)
)

line_trend_resid <- ggplot(plume_resid, aes(x = date, y = value)) +
  geom_path(alpha = 0.5, aes(colour = name)) +
  geom_smooth(method = "lm", linewidth = 2, aes(colour = name)) +
  labs(x = NULL, y = "Plume area (km^2)", title = "Trend for plume residual components") +
  geom_label(data = trend_resid_labels, aes(x = x, y = y, label = paste0(name," data slope = ", round(slope, 2), " km^2 yr-1", sep = "")),
             colour = c("turquoise4", "chartreuse4"), size = 6, hjust = 0) +
  scale_color_manual(name = "Decomposition method",
                     values = c("STL residual" = "turquoise4", "X11 residual" = "chartreuse4"),
                     breaks = c("STL residual", "X11 residual")) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_trend_comparison_resid.png", plot = line_trend_resid, height = 6, width = 12)

# The same plot comparing the interannual components
plume_inter <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  dplyr::select(zone, plot_title, date, plume_inter_STL, Interannual_signal_X11_plume) |> 
  pivot_longer(cols = c(plume_inter_STL, Interannual_signal_X11_plume)) |> 
  mutate(name = case_when(name == "plume_inter_STL" ~ "STL interannual",
                          name == "Interannual_signal_X11_plume" ~ "X11 interannual")) |> 
  filter(!is.na(value)) |> 
  filter(date <= as.Date("2024-12-31"),
         date >= as.Date("1999-01-01"))

# Create small dataframe of linear model slopes by name
trend_inter_STL <- coef(lm(value ~ date, data = filter(plume_inter, name == "STL interannual")))["date"] * 365.25
trend_inter_X11 <- coef(lm(value ~ date, data = filter(plume_inter, name == "X11 interannual")))["date"] * 365.25
trend_inter_labels <- data.frame(
  name = c("STL interannual", "X11 interannual"),
  slope = c(trend_inter_STL, trend_inter_X11),
  x = as.Date(c("1998-01-01", "1998-01-01")),
  y = c(2800, 2500)
)

# Create monthly means for more approximate comparison
plume_inter_monthly <- plume_inter |> 
  mutate(date = round_date(date, "month") + days(14)) |> 
  summarise(value = mean(value, na.rm = TRUE), .by = c("zone", "plot_title", "date", "name")) |> 
  mutate(timestep = "monthly")

# Create annual means
plume_inter_annual <- plume_inter |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
  summarise(value = mean(value, na.rm = TRUE), .by = c("zone", "plot_title", "date", "name")) |> 
  mutate(timestep = "annual")

# Combine to plot all at once
plume_inter_timesteps <- mutate(plume_inter, timestep = "daily") |> 
  rbind(plume_inter_monthly) |> 
  rbind(plume_inter_annual) |> 
  mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))

# Plot interannual trends
line_trend_inter <- ggplot(plume_inter_timesteps, aes(x = date, y = value)) +
  geom_line(data = filter(plume_inter_timesteps, timestep == "daily"), 
            aes(colour = name, linetype = timestep), alpha = 0.7) +
  geom_line(data = filter(plume_inter_timesteps, timestep == "monthly"), 
            aes(colour = name, linetype = timestep), alpha = 0.8, linewidth = 1.5) +
  geom_line(data = filter(plume_inter_timesteps, timestep == "annual"), 
            aes(colour = name, linetype = timestep), alpha = 0.9, linewidth = 2.0) +
  geom_smooth(method = "lm", linewidth = 2, aes(colour = name)) +
  labs(x = NULL, y = "Plume area (km^2)", title = "Trend for plume interannual components") +
  geom_label(data = trend_inter_labels, aes(x = x, y = y, label = paste0(name," data slope = ", round(slope, 2), " km^2 yr-1", sep = "")),
             colour = c("turquoise4", "chartreuse4"), size = 6, hjust = 0) +
  scale_color_manual(name = "Decomposition",
                     values = c("STL interannual" = "turquoise4", "X11 interannual" = "chartreuse4"),
                     breaks = c("STL interannual", "X11 interannual")) +
  scale_linetype_manual(name = "Time step",
                     values = c("daily" = "dotted", "monthly" = "dashed", "annual" = "solid"),
                     breaks = c("daily", "monthly", "annual"),
                     labels = c("daily/weekly", "monthly", "annual")) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_trend_comparison_inter.png", plot = line_trend_inter, height = 6, width = 12)


### Proportion --------------------------------------------------------------

# Plot the interannual time series against the raw data
plume_ts <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  filter(date >= as.Date("1999-01-01"), date <= as.Date("2024-12-31")) |> 
  dplyr::select(zone, plot_title, date, plume_area,
                plume_seas_STL, plume_inter_STL, plume_resid_STL,
                Interannual_signal_X11_plume, Seasonal_signal_X11_plume, Residual_signal_X11_plume) |> 
  mutate(total_STL = plume_seas_STL + plume_inter_STL + plume_resid_STL,
         total_X11 = Seasonal_signal_X11_plume + Interannual_signal_X11_plume + Residual_signal_X11_plume) |>
  pivot_longer(cols = plume_area:total_X11) |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
  summarise(value = mean(value, na.rm = TRUE), .by = c("zone", "plot_title", "date", "name")) |> 
  mutate(decomp_group = case_when(name == "plume_area" ~ "base",
                                  name %in% c("total_STL", "plume_seas_STL", "plume_inter_STL", "plume_resid_STL") ~ "STL",
                                  name %in% c("total_X11", "Seasonal_signal_X11_plume", "Interannual_signal_X11_plume", "Residual_signal_X11_plume") ~ "X11"),
         component_group = case_when(name %in% c("plume_area", "total_STL", "total_X11") ~ "total",
                                     name %in% c("plume_seas_STL", "Seasonal_signal_X11_plume") ~ "seasonal",
                                     name %in% c("plume_inter_STL", "Interannual_signal_X11_plume") ~ "interannual",
                                     name %in% c("plume_resid_STL", "Residual_signal_X11_plume") ~ "residual")) |> 
  mutate(decomp_group = factor(decomp_group, levels = c("base", "X11", "STL")),
         component_group = factor(component_group, levels = c("total", "interannual", "seasonal", "residual")))

# Plot the time series' of values
line_plume_comp_vs_raw <-ggplot(plume_ts, aes(x = date, y = value, colour = component_group, linetype = decomp_group)) +
  geom_path(linewidth = 1.0) +
  labs(colour = NULL, x = NULL, y = "Plume area (km^2)",
       title = "Decomposition components compared to base plume area (average mean values)") +
  scale_colour_brewer(palette = "Dark2") +
  scale_linetype_manual(name = "Decomposition",
                 values = c("solid", "dashed", "dotted")) +
                 # values = c("total" = "solid", "seasonal" = "dotdash", "interannual" = "dashed", "residual" = "dotted")) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", expand = c(0, 0)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_components_vs_raw.png", plot = line_plume_comp_vs_raw, height = 6, width = 12)

# Calculate proportion values
plume_prop <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  filter(date >= as.Date("1999-01-01"), date <= as.Date("2024-12-31")) |> 
  dplyr::select(zone, plot_title, date, plume_area,
                plume_seas_STL, plume_inter_STL, plume_resid_STL,
                Interannual_signal_X11_plume, Seasonal_signal_X11_plume, Residual_signal_X11_plume) |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |>
  group_by(zone, plot_title, date) |> 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  mutate(total_STL = plume_inter_STL + plume_seas_STL + plume_resid_STL,
         prop_inter_STL = plume_inter_STL / plume_area,
         prop_seas_STL = plume_seas_STL / plume_area,
         prop_resid_STL = plume_resid_STL / plume_area,
         total_X11 = Interannual_signal_X11_plume + Seasonal_signal_X11_plume + Residual_signal_X11_plume,
         prop_inter_X11 = Interannual_signal_X11_plume / plume_area,
         prop_seas_X11 = Seasonal_signal_X11_plume / plume_area,
         prop_resid_X11 = Residual_signal_X11_plume / plume_area) |> 
  dplyr::select(zone, plot_title, date, plume_area, plume_inter_STL, Interannual_signal_X11_plume,
                total_STL, prop_inter_STL, prop_seas_STL, prop_resid_STL,
                total_X11, prop_inter_X11, prop_seas_X11, prop_resid_X11) |>
  pivot_longer(cols = plume_area:prop_resid_X11) |> 
  mutate(decomp_group = case_when(name == "plume_area" ~ "base",
                                  name %in% c("total_STL", "plume_inter_STL",
                                              "prop_seas_STL", "prop_inter_STL", "prop_resid_STL") ~ "STL",
                                  name %in% c("total_X11", "Interannual_signal_X11_plume",
                                              "prop_seas_X11", "prop_inter_X11", "prop_resid_X11") ~ "X11"),
         component_group = case_when(name %in% c("plume_area", "total_STL", "total_X11") ~ "total",
                                     name %in% c("prop_seas_STL", "prop_seas_X11") ~ "seasonal",
                                     name %in% c("prop_inter_STL", "prop_inter_X11") ~ "interannual",
                                     name %in% c("prop_resid_STL", "prop_resid_X11") ~ "residual"),
         linear_plot = case_when(name %in% c("plume_area", "plume_inter_STL", "Interannual_signal_X11_plume") ~ "yes")) |> 
  mutate(decomp_group = factor(decomp_group, levels = c("base", "X11", "STL")),
         component_group = factor(component_group, levels = c("total", "interannual", "seasonal", "residual")))

# Get scaling factors for dual axis plot
scaling_factor_plume_prop <- sec_axis_adjustement_factors(plume_prop$value[plume_prop$component_group == "total"], 
                                                          plume_prop$value[plume_prop$component_group != "total"])

# Create scaled columns
plume_prop <- plume_prop |> 
  mutate(scaled_value = case_when(component_group != "total" ~ value,
                                  linear_plot == "yes" ~ value * scaling_factor_plume_prop$diff + scaling_factor_plume_prop$adjust))
                                  # component_group == "total" ~ value * scaling_factor_prop$diff + scaling_factor_prop$adjust))

# Double y-axis plot that shows the total values on the second Y axis, and the proportion values on the first Y axis
line_plume_prop <- ggplot(filter(plume_prop, component_group != "total"), aes(x = date, y = scaled_value)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "black", linewidth = 1.5) +
  geom_path(data = filter(plume_prop, linear_plot == "yes"),#component_group == "total"),
            aes(y = scaled_value, colour = decomp_group), linewidth = 1.5) +
  geom_col(aes(fill = component_group, colour = decomp_group), linewidth = 1.0, position = "dodge", alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_manual(values = c("brown", "yellow", "purple")) +
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    name = "Proportion of total plume area",
    sec.axis = sec_axis(transform = ~ {. - scaling_factor_plume_prop$adjust} / scaling_factor_plume_prop$diff,
                        name = "Average annual plume area (km^2)", breaks = c(600, 800, 1000), labels = c("600", "800", "1000"))
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", expand = c(0, 0)) +
  labs(x = NULL, colour = "Decomposition", fill = "Component",
       title = "Proportion of annual mean base plume area by decomposition method",
       subtitle = "Second axis shows annual mean base plume area and the STL and X11 interannual values") +
  # scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_proportion_comparison.png", plot = line_plume_prop, height = 6, width = 12)


## River flow comparison ---------------------------------------------------

#### Interannual -------------------------------------------------------------

# Extract and prep the interannual values
flow_inter <- decomp_df |> 
  dplyr::select(zone, plot_title, date, flow_inter_STL, Interannual_signal_X11_flow) |> 
  pivot_longer(cols = c(flow_inter_STL, Interannual_signal_X11_flow)) |> 
  mutate(name = case_when(name == "flow_inter_STL" ~ "flow STL",
                          name == "Interannual_signal_X11_flow" ~ "flow X11")) |> 
  filter(!is.na(value))

# Line plot of flow size - interannual
line_flow_inter <- ggplot(flow_inter, aes(x = date, y = value)) +
  geom_path(aes(colour = name), linewidth = 2) +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = NULL, y = "River flow (m^3 s-1)") + 
  scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of flow size - interannual
scatter_flow_inter <- flow_inter |> 
  pivot_wider(names_from = name, values_from = value) |> 
  mutate(month = month(date, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `flow X11`, y = `flow STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = month)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_flow_inter <- line_flow_inter + scatter_flow_inter + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_flow_inter_comp.png", plot = multi_flow_inter, height = 20, width = 30)


### Seasonal ----------------------------------------------------------------

# Extract and prep the monthly values
flow_seas <- decomp_df |> 
  dplyr::select(zone, plot_title, date, flow_seas_STL, Seasonal_signal_X11_flow) |> 
  pivot_longer(cols = c(flow_seas_STL, Seasonal_signal_X11_flow)) |> 
  mutate(name = case_when(name == "flow_seas_STL" ~ "flow STL",
                          name == "Seasonal_signal_X11_flow" ~ "flow X11")) |> 
  filter(!is.na(value))

# Line plot of flow size - seas
line_flow_seas <- ggplot(flow_seas, aes(x = date, y = value)) +
  geom_path(aes(colour = name)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = NULL, y = "River flow (m^3 s-1)") +
  scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of flow size - seas
scatter_flow_seas <- flow_seas |> 
  pivot_wider(names_from = name, values_from = value) |> 
  mutate(month = month(date, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `flow X11`, y = `flow STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = month)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_flow_seas <- line_flow_seas + scatter_flow_seas + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_flow_seas_comp.png", plot = multi_flow_seas, height = 20, width = 30)


### DOY ---------------------------------------------------------------------

# get the average doy values
# TODO: Improve the doy workflow. Get the source code from heatwaveR
flow_doy <- flow_seas |> 
  mutate(year = year(date),
         doy = yday(date)) #|> 
# mutate(doy_adj = adjust_doy(year, doy))
flow_doy$doy_adj <- mapply(adjust_doy, flow_doy$year, flow_doy$doy)
flow_doy <- flow_doy |> 
  summarise(val_min = min(value, na.rm = TRUE),
            val_mean = mean(value, na.rm = TRUE),
            val_max = max(value, na.rm = TRUE), 
            .by = c("zone", "plot_title", "name", "doy_adj")) |> 
  arrange(doy_adj)

# Ribbon plot of flow - doy
line_flow_doy <- ggplot(flow_doy, aes(x = doy_adj, y = val_mean)) +
  geom_ribbon(aes(fill = name, ymin = val_min, ymax = val_max), alpha = 0.2, show.legend = FALSE) +
  geom_path(aes(colour = name), linewidth = 2)  +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = "day-of-year", y = "River flow (m^3 s-1)") +
  scale_colour_brewer(palette = "Dark2", aesthetics = c("colour", "fill")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of flow size - doy
scatter_flow_doy <- flow_doy |> 
  dplyr::select(-val_min, -val_max) |> 
  pivot_wider(names_from = name, values_from = val_mean) |> 
  # mutate(month = month(doy_adj, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `flow X11`, y = `flow STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = doy_adj), size = 7) +
  scale_colour_viridis_c() +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = "day-of-year") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_flow_doy <- line_flow_doy + scatter_flow_doy + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_flow_doy_comp.png", plot = multi_flow_doy, height = 20, width = 30)


### Residual ----------------------------------------------------------------

# Extract and prep the interannual values
flow_resid <- decomp_df |> 
  dplyr::select(zone, plot_title, date, flow_resid_STL, Residual_signal_X11_flow) |> 
  pivot_longer(cols = c(flow_resid_STL, Residual_signal_X11_flow)) |> 
  mutate(name = case_when(name == "flow_resid_STL" ~ "flow STL",
                          name == "Residual_signal_X11_flow" ~ "flow X11")) |> 
  filter(!is.na(value))

# Line plot of flow size - residual
line_flow_resid <- ggplot(flow_resid, aes(x = date, y = value)) +
  geom_path(aes(colour = name), linewidth = 2) +
  facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
  labs(colour = NULL, x = NULL, y = "flow area (km^2)") +
  scale_colour_brewer(palette = "Dark2", aesthetics = c("colour", "fill")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(colour = "black", linewidth = 1),
        panel.grid.major.x = element_line(colour = "black", linewidth = 2))

# Scatterplot of flow size - seas
scatter_flow_resid <- flow_resid |> 
  pivot_wider(names_from = name, values_from = value) |> 
  mutate(month = month(date, label = TRUE, abbr = FALSE)) |> 
  ggplot(aes(x = `flow X11`, y = `flow STL`)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 3, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", colour = "black", linewidth = 3) +
  geom_point(aes(colour = month)) +
  facet_wrap(~plot_title, ncol = 1, scales = "free") +
  labs(colour = NULL) +
  # coord_cartesian(ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggplot_theme() +
  theme(legend.position = "bottom")

# Combine and save
multi_flow_resid <- line_flow_resid + scatter_flow_resid + patchwork::plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/STL_X11_flow_resid_comp.png", plot = multi_flow_resid, height = 20, width = 30)


### Trends ------------------------------------------------------------------

# First we start with a linear trend of the raw time series, and another plot with the linear analysis per month
## Daily data
flow_trend_daily <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |>
  dplyr::select(zone, plot_title, date, flow_seas_STL, flow_inter_STL, flow_resid_STL) |> 
  mutate(flow = flow_seas_STL + flow_inter_STL + flow_resid_STL) |> 
  filter(date <= as.Date("2023-12-31"),
         date >= as.Date("1999-01-01"))
## Monthly
flow_trend_monthly <- flow_trend_daily |> 
  mutate(date = round_date(date, "month") + days(14)) |> 
  summarise(flow = mean(flow, na.rm = TRUE), .by = c("zone", "plot_title", "date"))
## Annual
flow_trend_annual <- flow_trend_daily |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
  summarise(flow = mean(flow, na.rm = TRUE), .by = c("zone", "plot_title", "date"))
# Fit the linear models and extract slopes
trend_daily <- coef(lm(flow ~ date, data = flow_trend_daily))["date"] * 365.25
trend_monthly <- coef(lm(flow ~ date, data = flow_trend_monthly))["date"] * 365.25
trend_annual <- coef(lm(flow ~ date, data = flow_trend_annual))["date"] * 365.25

# Create a little dataframe for plotting the labels
trend_labels <- data.frame(
  time_step = c("daily", "monthly", "annual"),
  slope = c(trend_daily, trend_monthly, trend_annual),
  x = as.Date(c("1998-01-01", "1998-01-01", "1998-01-01")),
  y = c(9000, 8000, 7000)
)

# Plot the three dataframes as time series with linear trends on the same panel
line_trend_base <- ggplot(flow_trend_daily, aes(x = date, y = flow)) +
  geom_point(alpha = 0.5, aes(colour = "daily")) +
  geom_smooth(method = "lm", colour = "black", linewidth = 2) +
  geom_point(data = flow_trend_monthly, size = 4, alpha = 0.5, aes(colour = "monthly")) +
  geom_smooth(data = flow_trend_monthly, method = "lm", colour = "darkblue", linewidth = 2) +
  geom_point(data = flow_trend_annual, size = 7, alpha = 0.5, aes(colour = "annual")) +
  geom_smooth(data = flow_trend_annual, method = "lm", colour = "blue", linewidth = 2) +
  geom_label(data = trend_labels, aes(x = x, y = y, label = paste0(time_step," data slope = ", round(slope, 2), " m^3 s-1 y-1", sep = "")),
             colour = c("black", "darkblue", "blue"), size = 6, hjust = 0) +
  labs(x = NULL, y = "River flow (m^3 s-1)", title = "Trend for daily, monthly, and annual mean river flow data") +
  scale_color_manual(name = "Time step",
                     values = c("daily" = "black", "monthly" = "darkblue", "annual" = "blue"),
                     breaks = c("daily", "monthly", "annual")) +
  # ggplot_theme() +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/flow_trend_comparison_base.png", plot = line_trend_base, height = 6, width = 12)

# The same plot comparing the interannual components
flow_inter <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  dplyr::select(zone, plot_title, date, flow_inter_STL, Interannual_signal_X11_flow) |> 
  pivot_longer(cols = c(flow_inter_STL, Interannual_signal_X11_flow)) |> 
  mutate(name = case_when(name == "flow_inter_STL" ~ "STL interannual",
                          name == "Interannual_signal_X11_flow" ~ "X11 interannual")) |> 
  filter(!is.na(value)) |> 
  filter(date <= as.Date("2023-12-31"),
         date >= as.Date("1999-01-01"))

# Create small dataframe of linear model slopes by name
trend_inter_STL <- coef(lm(value ~ date, data = filter(flow_inter, name == "STL interannual")))["date"] * 365.25
trend_inter_X11 <- coef(lm(value ~ date, data = filter(flow_inter, name == "X11 interannual")))["date"] * 365.25
trend_inter_labels <- data.frame(
  name = c("STL interannual", "X11 interannual"),
  slope = c(trend_inter_STL, trend_inter_X11),
  x = as.Date(c("1998-01-01", "1998-01-01")),
  y = c(6000, 5500)
)

# Create monthly means for more approximate comparison
flow_inter_monthly <- flow_inter |> 
  mutate(date = round_date(date, "month") + days(14)) |> 
  summarise(value = mean(value, na.rm = TRUE), .by = c("zone", "plot_title", "date", "name")) |> 
  mutate(timestep = "monthly")

# Create annual means
flow_inter_annual <- flow_inter |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |> 
  summarise(value = mean(value, na.rm = TRUE), .by = c("zone", "plot_title", "date", "name")) |> 
  mutate(timestep = "annual")

# Combine to plot all at once
flow_inter_timesteps <- mutate(flow_inter, timestep = "daily") |> 
  rbind(flow_inter_monthly) |> 
  rbind(flow_inter_annual) |> 
  mutate(timestep = factor(timestep, levels = c("daily", "monthly", "annual")))

line_trend_inter <- ggplot(flow_inter_timesteps, aes(x = date, y = value)) +
  geom_line(data = filter(flow_inter_timesteps, timestep == "daily"), 
            aes(colour = name, linetype = timestep), alpha = 0.7) +
  geom_line(data = filter(flow_inter_timesteps, timestep == "monthly"), 
            aes(colour = name, linetype = timestep), alpha = 0.8, linewidth = 1.5) +
  geom_line(data = filter(flow_inter_timesteps, timestep == "annual"), 
            aes(colour = name, linetype = timestep), alpha = 0.9, linewidth = 2.0) +
  geom_smooth(method = "lm", linewidth = 2, aes(colour = name)) +
  labs(x = NULL, y = "River flow (m^3 s-1)", title = "Trend for river flow interannual components") +
  geom_label(data = trend_inter_labels, aes(x = x, y = y, label = paste0(name," data slope = ", round(slope, 2), " m^3 s-1 y-1", sep = "")),
             colour = c("turquoise4", "chartreuse4"), size = 6, hjust = 0) +
  scale_color_manual(name = "Decomposition method",
                     values = c("STL interannual" = "turquoise4", "X11 interannual" = "chartreuse4"),
                     breaks = c("STL interannual", "X11 interannual")) +
  scale_linetype_manual(name = "Time step",
                        values = c("daily" = "dotted", "monthly" = "dashed", "annual" = "solid"),
                        breaks = c("daily", "monthly", "annual"),
                        labels = c("daily/weekly", "monthly", "annual")) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/flow_trend_comparison_inter.png", plot = line_trend_inter, height = 6, width = 12)


### Proportion --------------------------------------------------------------

# Calculate proportion values
flow_prop <- decomp_df |> 
  filter(zone == "GULF_OF_LION") |> 
  filter(date >= as.Date("1999-01-01"), date <= as.Date("2024-12-31")) |> 
  dplyr::select(zone, plot_title, date, flow,
                flow_seas_STL, flow_inter_STL, flow_resid_STL,
                Interannual_signal_X11_flow, Seasonal_signal_X11_flow, Residual_signal_X11_flow) |> 
  mutate(date = as.Date(paste0(year(date),"-07-01"))) |>
  group_by(zone, plot_title, date) |> 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  mutate(total_STL = flow_seas_STL + flow_inter_STL + flow_resid_STL,
         prop_seas_STL = flow_seas_STL / flow,
         prop_inter_STL = flow_inter_STL / flow,
         prop_resid_STL = flow_resid_STL / flow,
         total_X11 = Seasonal_signal_X11_flow + Interannual_signal_X11_flow + Residual_signal_X11_flow,
         prop_seas_X11 = Seasonal_signal_X11_flow / flow,
         prop_inter_X11 = Interannual_signal_X11_flow / flow,
         prop_resid_X11 = Residual_signal_X11_flow / flow) |> 
  dplyr::select(zone, plot_title, date, flow, flow_inter_STL, Interannual_signal_X11_flow,
                total_STL, prop_seas_STL, prop_inter_STL, prop_resid_STL,
                total_X11, prop_seas_X11, prop_inter_X11, prop_resid_X11) |>
  pivot_longer(cols = flow:prop_resid_X11) |> 
  mutate(decomp_group = case_when(name == "flow" ~ "base",
                                  name %in% c("total_STL", "flow_inter_STL",
                                              "prop_seas_STL", "prop_inter_STL", "prop_resid_STL") ~ "STL",
                                  name %in% c("total_X11", "Interannual_signal_X11_flow",
                                              "prop_seas_X11", "prop_inter_X11", "prop_resid_X11") ~ "X11"),
         component_group = case_when(name %in% c("flow", "total_STL", "total_X11") ~ "total",
                                     name %in% c("prop_seas_STL", "prop_seas_X11") ~ "seasonal",
                                     name %in% c("prop_inter_STL", "prop_inter_X11") ~ "interannual",
                                     name %in% c("prop_resid_STL", "prop_resid_X11") ~ "residual"),
         linear_plot = case_when(name %in% c("flow", "flow_inter_STL", "Interannual_signal_X11_flow") ~ "yes")) |> 
  mutate(decomp_group = factor(decomp_group, levels = c("base", "STL", "X11")),
         component_group = factor(component_group, levels = c("total", "interannual", "seasonal", "residual")))

# Get scaling factors for dual axis plot
scaling_factor_flow_prop <- sec_axis_adjustement_factors(flow_prop$value[flow_prop$component_group == "total"], 
                                                          flow_prop$value[flow_prop$component_group != "total"])

# Create scaled columns
flow_prop <- flow_prop |> 
  mutate(scaled_value = case_when(component_group != "total" ~ value,
                                  linear_plot == "yes" ~ value * scaling_factor_flow_prop$diff + scaling_factor_flow_prop$adjust)) |> 
  filter(date >= as.Date("1999-01-01"), date <= as.Date("2023-12-31"))

# Double y-axis plot that shows the total values on the second Y axis, and the proportion values on the first Y axis
line_flow_prop <- ggplot(filter(flow_prop, component_group != "total"), aes(x = date, y = scaled_value)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "black", linewidth = 1.5) +
  geom_path(data = filter(flow_prop, linear_plot == "yes"),#component_group == "total"),
            aes(y = scaled_value, colour = decomp_group), linewidth = 1.5) +
  geom_col(aes(fill = component_group, colour = decomp_group), linewidth = 1.0, position = "dodge", alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_manual(values = c("brown", "yellow", "purple")) +
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
                     name = "Proportion of total river flow",
                     sec.axis = sec_axis(transform = ~ {. - scaling_factor_flow_prop$adjust} / scaling_factor_flow_prop$diff,
                                         name = "Average annual river flow (m^3 s-1)", 
                                         breaks = c(1000, 1500, 2000), labels = c("1000", "1500", "2000"))
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", expand = c(0, 0)) +
  labs(x = NULL, colour = "Decomposition", fill = "Component",
       title = "Proportion of annual mean base river flow by decomposition method",
       subtitle = "Second axis shows annual mean base river flow and the STL and X11 interannual values") +
  # scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/flow_proportion_comparison.png", plot = line_flow_prop, height = 6, width = 12)

# Combine flow Prop and plume_prop dataframes
plume_flow_prop <- rbind(mutate(plume_prop, var_name = "plume"),
                         mutate(flow_prop, var_name = "flow")) |> 
  filter(date >= as.Date("1999-01-01"), date <= as.Date("2023-12-31"))

# Plot the proportion of each component of X11 for flow next to plume area
bar_plume_flow_prop <-ggplot(filter(plume_flow_prop, component_group != "total", decomp_group == "X11"), aes(x = date, y = value)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "black", linewidth = 1.5) +
  geom_col(aes(fill = component_group, colour = var_name), linewidth = 1.0, position = "dodge", alpha = 0.5) +
  scale_colour_brewer("Variable", palette = "Set2") +
  scale_fill_manual("Component", values = c("brown", "yellow", "purple")) +
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of mean area/flow",
       title = "Annual proportion of mean plume area and river flow by X11 decomposition components") +
  # scale_colour_brewer(palette = "Dark2") +
  guides(colour = guide_legend(override.aes = list(linewidth = 5))) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/plume_flow_proportion_comparison.png", plot = bar_plume_flow_prop, height = 6, width = 12)



# EMD example -------------------------------------------------------------

# The needed package
library(EMD)
library(ggplot2)

# Generate some dummy data
set.seed(13)
t <- seq(0, 10, by = 0.1)
signal <- sin(2 * pi * t) + 0.5 * sin(2 * pi * 5 * t) + 0.2 * rnorm(length(t))

# Perform EMD
emd_result <- emd(signal, t)

# Extract IMFs
imfs <- emd_result$imf

# Create a data frame for the original signal
df_signal <- data.frame(Time = t, Signal = signal, Type = "Original")

# Create a data frame for the IMFs
df_imfs <- data.frame(Time = rep(t, times = ncol(imfs)),
                      Signal = as.vector(t(imfs)),
                      Type = rep(paste0("IMF ", 1:ncol(imfs)), each = length(t)))

# Combine the data frames
df_combined <- rbind(df_signal, df_imfs)

# Plot the original signal and IMFs
ggplot(df_combined, aes(x = Time, y = Signal, color = Type)) +
  geom_line() +
  facet_wrap(~ Type, ncol = 1, scales = "free_y") +
  labs(title = "Empirical Mode Decomposition (EMD)",
       x = "Time",
       y = "Amplitude") +
  theme_minimal() +
  theme(legend.position = "none")


# BEAST example -----------------------------------------------------------

library(Rbeast)

# TODO: Create example


# POC/DOC -----------------------------------------------------------------

# Load MOOSE data
rhone_moose <- read_csv("~/Downloads/Water_sample_analyses_-_MOOSE_-_Rhone_river/SEDOO-MOOSE-Rhone  Biogenic data-2005-2022.csv")

# Melt and columns of interest
rhone_moose_long <- rhone_moose |> 
  dplyr::rename(debit_m3s = `Débit  moyen – m3/s`,
                SPM_mgL = `Matière en suspension – mg/litre`,
                DOC_µmCL = `carbone organique dissous – µmoles(C)/litre`,
                POC_µmCL = `Carbone organique particulaire – µmoles(C)/litre`) |>
  mutate(debit_m3s = as.numeric(debit_m3s),
         SPM_mgL = as.numeric(SPM_mgL),
         DOC_µmCL = as.numeric(DOC_µmCL),
         POC_µmCL = as.numeric(POC_µmCL)) |>
  dplyr::select(date, debit_m3s, SPM_mgL, DOC_µmCL, POC_µmCL) |>
  pivot_longer(cols = c(debit_m3s, SPM_mgL, DOC_µmCL, POC_µmCL), names_to = "variable", values_to = "value") |> 
  mutate(date = as.Date(date, format = "%d/%m/%Y")) |> 
  mutate(variable = factor(variable, levels = c("debit_m3s", "SPM_mgL", "DOC_µmCL", "POC_µmCL"),
                           labels = c("River discharge (m3 s-1)", "Suspended particulate matter (mg L-1)",
                                      "Dissolved organic carbon (µmol C L-1)", "Particulate organic carbon (µmol C L-1)"))) |> 
  mutate(date = case_when(year(date) < 2005 ~ date + years(2000),
                          TRUE ~ date))

# Get linear model stats
rhone_moose_lm_stats <- rhone_moose_long |> 
  summarise(var_slope = coef(lm(value ~ date))["date"] * 365.25,
            # var_perc = round((value/mean(rhone_moose_long$value, na.rm = TRUE))*100, 2),
            var_p = round(summary(lm(value ~ date))[["coefficients"]][2,4], 4), .by = "variable")

# Quick line plot faceted by variable
line_rhone_moose <- ggplot(rhone_moose_long, aes(x = date, y = value)) +
  geom_line() +
  geom_smooth(method = "lm") +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = NULL, title = "Rhône river (MOOSE)") +
  ggplot_theme() +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
ggsave(filename = "figures/rhone_moose_biogenic_timeseries.png",
       plot = line_rhone_moose, height = 12, width = 12)
