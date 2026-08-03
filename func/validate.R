# func/validate.R

# This script contains the code necessary to run validation on a number of
# satellite products and variables against multiple sources of in situ data

# Each section can run independently of the others


# Setup -------------------------------------------------------------------

# The packages needed for this script
library(tidyverse)
library(furrr)
library(raster) # Used here to select specific pixels in a .nc file
library(ncdf4)
library(sf) # Used for complex shape files
library(scales) # For better plot labels
library(ggExtra) # For histogram border plots
library(gt) # For fancy tables
library(geosphere) # For pixel distances
# library(doParallel); registerDoParallel(cores = detectCores()-6) # NB: Trying to move away from this

# The shared functions
source("func/util.R")


# Satellite times ---------------------------------------------------------

# Get all satellite daily start and end times when overhead France
# This could be used to correctly filter the in situ data by time when possible
## SEXTANT
### NB: Handled differently because the values are always the same and downloading all files takes an hour...
times_SEXTANT <- data.frame(date_start = ymd(unlist(lapply(str_split(basename(dir("~/pCloudDrive/data/SEXTANT/SPM",
                                                                                  recursive = TRUE, pattern = "SPIM")), "-"), "[[", 1)))) |>
  mutate(date_end = date_start + 1,
         start_time = ymd_hms(paste(date_start, "12:00:00"), tz = "GMT"),
         end_time = ymd_hms(paste(date_end, "12:00:00"), tz = "GMT")) |>
  dplyr::select(start_time, end_time)
save(times_SEXTANT, file = "metadata/times_SEXTANT.RData")
load("metadata/times_SEXTANT.RData")

# Set cores for following calculations
plan(multicore, workers = parallel::detectCores() - 2)

## MODIS
times_MODIS <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily",
                               pattern = "SPM", full.names = TRUE), get_start_end_time)
save(times_MODIS, file = "metadata/times_MODIS.RData")
load("metadata/times_MODIS.RData")

### Double check that times are the same for any region
# Yes, start and end times are exactly the same
times_MODIS_check <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/GULF_OF_LION/daily", 
                                     pattern = "SPM", full.names = TRUE), get_start_end_time)
times_MODIS_test <- bind_cols(times_MODIS, times_MODIS_check) |>
  mutate(start_diff = start_time...1 - start_time...3, end_diff = end_time...2 - end_time...4)

### SST-NIGHT are meant to have a different time
times_MODIS_SST_NIGHT <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily",
                                         pattern = "SST-NIGHT", full.names = TRUE), get_start_end_time)
save(times_MODIS_SST_NIGHT, file = "metadata/times_MODIS_SST_NIGHT.RData")
load("metadata/times_MODIS_SST_NIGHT.RData")
### SST day time
## These are calculated on rs-pro2 where the expert data are located
load("metadata/times_MODIS_SST.RData")

## MERIS
times_MERIS <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/BAY_OF_SEINE/daily",
                               pattern = "SPM", full.names = TRUE), get_start_end_time)
save(times_MERIS, file = "metadata/times_MERIS.RData")
load("metadata/times_MERIS.RData")

## OLCI-A
times_OLCI_A <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-A/BAY_OF_SEINE/daily",
                                pattern = "SPM", full.names = TRUE), get_start_end_time)
save(times_OLCI_A, file = "metadata/times_OLCI_A.RData")
load("metadata/times_OLCI_A.RData")

## OLCI-B
times_OLCI_B <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-B/BAY_OF_SEINE/daily",
                                pattern = "SPM", full.names = TRUE), get_start_end_time)
save(times_OLCI_B, file = "metadata/times_OLCI_B.RData")
load("metadata/times_OLCI_B.RData")

## ALL
times_ALL <- bind_rows(mutate(times_SEXTANT, sat_name = "SEXTANT"),
                       mutate(times_MODIS, sat_name = "MODIS"),
                       mutate(times_MODIS_SST, sat_name = "MODIS (SST)"),
                       mutate(times_MODIS_SST_NIGHT, sat_name = "MODIS (SST-NIGHT)"),
                       mutate(times_MERIS, sat_name = "MERIS"),
                       mutate(times_OLCI_A, sat_name = "OLCI-A"),
                       mutate(times_OLCI_B, sat_name = "OLCI-B")) |>
  mutate(date = as.Date(start_time),
         start_hms = hms::as_hms(start_time),
         end_hms = hms::as_hms(end_time)) |>
  complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA), nesting(sat_name))
save(times_ALL, file = "metadata/times_ALL.RData") # Saving in R format to maintain time stamps
load("metadata/times_ALL.RData")

# Plot the times of day during which the satellites are available
time_plot <- ggplot(data = times_ALL, aes(x = date, y = start_hms)) +
  geom_point(data = filter(times_ALL, sat_name == "SEXTANT"), aes(colour = sat_name), size = 3) +
  geom_ribbon(aes(ymin = start_hms, ymax = end_hms, fill = sat_name)) +
  scale_fill_brewer("Sensor", palette = "Dark2", aesthetics = c("colour", "fill")) +
  labs(y = "Overhead time (UTC)", x = NULL,
       title = "SEXTANT and ODATIS-MR sensor availability and overhead times") +
  theme_bw()
# time_plot
ggsave("figures/validation/all_sensors_time.png", time_plot, width = 8, height = 4)


# Load in situ ------------------------------------------------------------

# REPHY
## NB: The script used to process the data, and necessary files, were provided by Victor Pochic
# "data/INSITU_data/REPHY/REPHY_Dataset_20250408.R"
REPHY <- read.csv2("data/INSITU_data/REPHY/Table1_REPHY_hydro_20250408.csv", fileEncoding = "ISO-8859-1") |> 
  mutate(lon = as.numeric(lon), lat = as.numeric(lat), source = "REPHY") |> 
  dplyr::rename(site = Code_point_Libelle, date = Date, time = Heure, variable = Code.parametre, value = Valeur_mesure)

# SOMLIT
# "data/INSITU_data/SOMLIT/SOMLIT_prep.R"
SOMLIT <- read_csv("data/INSITU_data/SOMLIT/Somlit_clean.csv") |> 
  mutate(source = "SOMLIT")

# Filter out sites that are within the RiOMar regions
in_situ_site_list <- bind_rows(dplyr::select(REPHY, source, lon, lat, site),
                               dplyr::select(SOMLIT, source, lon, lat, site)) |>
  distinct() |>
  summarise(lon = mean(lon), lat = mean(lat), .by = c("source", "site")) |>
  # Lanveoc (Rade de Brest, ~48.29N) sits just outside the Southern Brittany
  # bbox (lat_max = 48.00), so the bbox check below would otherwise drop it
  # entirely -- mirrors func/validate.py's region_mapping fix for the same
  # station.
  mutate(zone = case_when(site == "Lanvéoc" ~ "SOUTHERN_BRITTANY",
                          lon >= zones_bbox$lon_min[1] & lon <= zones_bbox$lon_max[1] &
                            lat >= zones_bbox$lat_min[1] & lat <= zones_bbox$lat_max[1] ~ "GULF_OF_LION",
                          lon >= zones_bbox$lon_min[2] & lon <= zones_bbox$lon_max[2] &
                            lat >= zones_bbox$lat_min[2] & lat <= zones_bbox$lat_max[2] ~ "BAY_OF_SEINE",
                          lon >= zones_bbox$lon_min[3] & lon <= zones_bbox$lon_max[3] &
                            lat >= zones_bbox$lat_min[3] & lat <= zones_bbox$lat_max[3] ~ "BAY_OF_BISCAY",
                          lon >= zones_bbox$lon_min[4] & lon <= zones_bbox$lon_max[4] &
                            lat >= zones_bbox$lat_min[4] & lat <= zones_bbox$lat_max[4] ~ "SOUTHERN_BRITTANY")) |>
  mutate(zone_pretty = factor(zone,
                              levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                              labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")), .after = "zone") |>
  mutate(source = factor(source, levels = c("SOMLIT", "REPHY")))
# Written for func/figure.py::save_insitu_stations_for_plot() (Figure 1's
# station overlay) and re-read further below in this script's own
# Validation tables section -- both need this on disk.
write_csv(in_situ_site_list, "metadata/in_situ_site_list.csv")

# Filter in situ stations to just those within a zone
zone_sites <- in_situ_site_list |> filter(!is.na(zone))

# Clean up and combine all in situ data into one dataframe
## REPHY
REPHY_clean <- right_join(REPHY, zone_sites, by = c("source", "site", "lon", "lat")) |> 
  filter(Qualite.resultat == "Bon") |> # Only keep 'good' measurements
  filter(as.numeric(Profondeur.metre) <= 10) |>  # Only keep values within 10 meters of surface
  dplyr::select(source, site, lon, lat, date, time, variable, value) |> 
  mutate(variable = case_when(variable == "SALI" ~ "SAL",
                              variable == "TURB" ~ "TUR",
                              variable == "CHLOROA" ~ "CHLA", TRUE ~ variable),
         date = as.Date(date),
         time = case_when(time == "" ~ NA, TRUE ~ time)) |> 
  mutate(time = hms::as_hms(time)) |> # Doesn't seem to like case_when() ...
  filter(value >= 0)

## SOMLIT
SOMLIT_clean <- right_join(SOMLIT, zone_sites, by = c("source", "site", "lon", "lat")) |> 
  filter(prof_num <= 10) |> # Filter data deeper than 10 meters
  # Filter bad flags
  mutate(TEMP = case_when(temp_QC %in% c(2, 6, 7) ~ temp),
         SAL = case_when(sal_QC %in% c(2, 6, 7) ~ sal),
         POC = case_when(POC_QC %in% c(2, 6, 7) ~ POC),
         SPM = case_when(SPM_QC %in% c(2, 6, 7) ~ SPM),
         CHLA = case_when(CHLA_QC %in% c(2, 6, 7) ~ CHLA)) |> 
  dplyr::select(source, site, lon, lat, date, time, TEMP, SAL, POC, SPM, CHLA) |> 
  pivot_longer(TEMP:CHLA, values_to = "value", names_to = "variable") |>
  mutate(time = hms::as_hms(time)) |> 
  filter(value >= 0)

# Combine with the satellite times to filter variables appropriately
# Ultimately decided to just filter from 10:00 to 15:00 to match Devreker et al. 2026
# "~/pCloudDrive/Documents/OMTAB/RiOMar/articles/Devreker_et_al_Archimer_2026.pdf"

## Combine
zone_data_in_situ <- bind_rows(REPHY_clean, SOMLIT_clean) |> 
  # Get only variables of interest
  filter(variable %in% c('TEMP', 'SAL', 'POC', 'SPM', 'CHLA', 'TUR')) |> 
  # Create additional HMS filter for night-time TEMP
  # Filter from times 10:00 to 15:00; 00:00 to 04:00 for TEMP
  mutate(time_class = case_when(time >= hms("10:00:00") & time <= hms("14:00:00") ~ "day",
                                variable == "TEMP" & time >= hms("00:00:00") | time <= hms("09:00:00") ~ "night")) |> 
  filter(time_class %in% c("day", "night")) |>
  filter(!(variable != "TEMP" & time_class == "night")) |>
  # Create daily means
  summarise(value = mean(value, na.rm = TRUE), .by = c("source", "site", "lon", "lat", "date", "time_class", "variable")) |> 
  left_join(zone_sites, by = join_by(source, site, lon, lat)) |> 
  dplyr::select(zone, zone_pretty, source, everything())
write_csv(zone_data_in_situ, "data/INSITU_data/zone_data_in_situ.csv")


# Map In situ -------------------------------------------------------------

# Load high-res shape files
# Tuto here : https://www.etiennebacher.com/posts/2021-12-27-mapping-french-rivers-network/
if(!exists("borders_FR")) borders_FR <- read_sf("data/FRANCE_shapefile/gadm41_FRA_0.shp")
if(!exists("rivers_FR")) rivers_FR <- st_intersection(read_sf("data/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"), borders_FR)

# Labels for number of sites per zone
zone_site_counts <- zone_sites |>
  count(zone, zone_pretty, source, name = "site_count") |>
  complete(nesting(zone, zone_pretty), source = levels(zone_sites[["source"]]), fill = list(site_count = 0)) |>
  arrange(zone_pretty, source) |>
  summarise(zone_label = paste0(source, " = ", site_count, collapse = "\n"),
            .by = c("zone", "zone_pretty"))

# Map all in situ stations, highlighting the zones and the stations used
zone_site_labels <- zone_site_counts |>
  left_join(dplyr::select(zones_bbox, zone, lon_min, lon_max, lat_min, lat_max), by = "zone") |>
  mutate(label_lon = case_when(zone == "GULF_OF_LION" ~ lon_min + 0.3,
                               zone == "BAY_OF_SEINE" ~ lon_max + 0.1,
                               zone == "BAY_OF_BISCAY" ~ lon_max + 0.5,
                               zone == "SOUTHERN_BRITTANY" ~ lon_max + 0.5),
         label_lat = case_when(zone == "GULF_OF_LION" ~ lat_max + 0.8,
                               zone == "BAY_OF_SEINE" ~ lat_min,
                               zone == "BAY_OF_BISCAY" ~ lat_min + 1.3,
                               zone == "SOUTHERN_BRITTANY" ~ lat_min + 1))

in_situ_station_map <- ggplot() +
  geom_sf(data = borders_FR, color = "black", fill = "sienna4", inherit.aes = FALSE) +
  geom_sf(data = rivers_FR, color = "lightblue", inherit.aes = FALSE, linewidth = 0.2) +
  geom_rect(data = zones_bbox, fill = NA, linewidth = 2, show.legend = FALSE,
            aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max, colour = zone_pretty)) +
  geom_label(data = zone_site_labels, show.legend = FALSE,
             aes(x = label_lon, y = label_lat, label = zone_label, colour = zone_pretty),
             hjust = 0, vjust = 1, size = 7, linewidth = 3.0,
             fill = "white", alpha = 1, inherit.aes = FALSE) +
  geom_label(data = zone_site_labels,
             aes(x = label_lon, y = label_lat, label = zone_label),
             hjust = 0, vjust = 1, size = 7, linewidth = 0.0,
             fill = "white", alpha = 1, inherit.aes = FALSE) +
  coord_sf(xlim = c(-5, 10), ylim = c(41.5, 51)) +
  geom_point(data = in_situ_site_list,
             aes(x = lon, y = lat, shape = source), color = "black", size = 4) +
  geom_point(data = filter(in_situ_site_list, is.na(zone)),
             aes(x = lon, y = lat, shape = source), color = "red", size = 3) +
  geom_point(data = filter(in_situ_site_list, !is.na(zone)),
             aes(x = lon, y = lat, shape = source, color = zone_pretty), size = 3) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(labels = scales::unit_format(unit = "°E")) +
  scale_y_continuous(labels = scales::unit_format(unit = "°N")) +
  scale_color_manual(name = "zone", values = colours_of_stations(), drop = FALSE) +
  cowplot::theme_map() +
  ggplot_theme() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.73, 0.8),
        legend.box.background = element_rect(colour = "black", fill = "white"),
        legend.box.margin = margin(5, 5, 5, 5),
        axis.text = element_text(size = 20))
# in_situ_station_map
ggsave("figures/validation/map_in_situ_stations.png", height = 14, width = 15.5, bg = "white")


# Prep satellite pixels ---------------------------------------------------

# Create the data.frames of pixel matchups
## SEXTANT
write_pixel_coords("SEXTANT", "SPM")
write_pixel_coords("SEXTANT", "CHLA")

# NB: The ODATIS-MR extractions are now made on the expert data, which are housed on Rspro2
# Therefore the following code is not run on the local machine, but on Rspro2, 
# and the resulting files are then transferred to the local machine for validation analyses

# MODIS
# MERIS
# OLCI-A
# OLCI-B


# Extract satellite data --------------------------------------------------

# SEXTANT
## Roughly 1 hour per variable
extract_pixels_all("SEXTANT")

# NB: The ODATIS-MR extractions are now made on the expert data, which are housed on Rspro2
# Therefore the following code is not run on the local machine, but on Rspro2, 
# and the resulting files are then transferred to the local machine for validation analyses

# MODIS
# MERIS
# OLCI-A
# OLCI-B


# Satellite grid-window match-up (Table 4) --------------------------------

# Concentric N x N pixel-window statistics (mean/median/std/n) around one
# grid index, for grid sizes 1, 3, 5, 7. This is a literal square lat/lon-
# index window (as opposed to get_pixels()'s nearest-by-distance selection
# above), matching the retired func/validate.py::compute_grid_stats()
# convention so Table 4 stays numerically consistent with its values before
# this port.
sextant_window_stats <- function(window, centre_lon_idx, centre_lat_idx, grid_sizes = c(1, 3, 5, 7)){

  stats <- list()
  for(size in grid_sizes){
    h <- size %/% 2
    lon_range <- (centre_lon_idx - h):(centre_lon_idx + h)
    lat_range <- (centre_lat_idx - h):(centre_lat_idx + h)
    lon_range <- lon_range[lon_range >= 1 & lon_range <= dim(window)[1]]
    lat_range <- lat_range[lat_range >= 1 & lat_range <= dim(window)[2]]
    vals <- as.vector(window[lon_range, lat_range])
    stats[[paste0(size, "x", size, "_mean")]] <- mean(vals, na.rm = TRUE)
    stats[[paste0(size, "x", size, "_median")]] <- median(vals, na.rm = TRUE)
    stats[[paste0(size, "x", size, "_std")]] <- sd(vals, na.rm = TRUE)
    stats[[paste0(size, "x", size, "_n")]] <- sum(!is.na(vals))
  }
  as.data.frame(stats, check.names = FALSE)
}

# Build the satellite validation summary (Table 4) by matching each SOMLIT
# SPM / REPHY turbidity in-situ record to the SEXTANT SPM pixel window
# around its station, for every date where both exist. Writes
# output/MATCH_UP_DATA/FRANCE/summary.csv, replacing the retired
# func/validate.py::Match_up_with_insitu_measurements() path.
#
# NB: zone_data_in_situ (built above) is already a QC'd, depth-filtered,
# day/night-filtered daily mean per site/date/variable -- reused directly
# here rather than re-deriving in-situ filtering logic. This means Table 4
# now uses one row per site/day (a daily mean) rather than one row per
# individual in-situ sample as the retired Python path did; the actual
# comparison (grid-window satellite vs in-situ value) is unaffected, and
# per-record satellite/in-situ overhead times are no longer tracked (left
# NA), since zone_data_in_situ collapses same-day records before this point.
compile_sextant_validation_summary <- function(){

  sextant_dir <- path.expand("~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY")

  # SEXTANT's lat/lon grid is fixed across every daily file -- build the
  # index lookup once from a reference file.
  ref_nc <- nc_open(dir(sextant_dir, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)[1])
  grid_lon <- ncvar_get(ref_nc, "lon")
  grid_lat <- ncvar_get(ref_nc, "lat")
  nc_close(ref_nc)

  in_situ_matches <- zone_data_in_situ |>
    filter(variable %in% c("SPM", "TUR")) |>
    mutate(Insitu_variable = if_else(variable == "TUR", "TURB", variable)) |>
    filter(date >= as.Date("1998-01-01"))

  site_grid_index <- in_situ_matches |>
    distinct(source, site, zone, lon, lat) |>
    rowwise() |>
    mutate(lon_idx = which.min(abs(grid_lon - lon)),
           lat_idx = which.min(abs(grid_lat - lat))) |>
    ungroup()

  half <- 3 # max(grid_sizes) %/% 2 for grid_sizes = c(1, 3, 5, 7)
  summary_rows <- list()

  for(the_date in sort(unique(in_situ_matches$date))){

    the_date <- as.Date(the_date, origin = "1970-01-01")
    date_str <- format(the_date, "%Y%m%d")
    nc_path <- file.path(sextant_dir, format(the_date, "%Y"), format(the_date, "%m"), format(the_date, "%d"),
                         paste0(date_str, "-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"))

    if(!file.exists(nc_path)) next

    records_of_the_day <- in_situ_matches |> filter(date == the_date)
    sites_of_the_day <- site_grid_index |> semi_join(records_of_the_day, by = c("source", "site"))

    nc <- nc_open(nc_path)

    for(i in seq_len(nrow(sites_of_the_day))){

      site_row <- sites_of_the_day[i, ]

      lon_lo <- max(1, site_row$lon_idx - half); lon_hi <- min(length(grid_lon), site_row$lon_idx + half)
      lat_lo <- max(1, site_row$lat_idx - half); lat_hi <- min(length(grid_lat), site_row$lat_idx + half)

      window <- ncvar_get(nc, "analysed_spim",
                          start = c(lon_lo, lat_lo, 1),
                          count = c(lon_hi - lon_lo + 1, lat_hi - lat_lo + 1, 1))

      window_stats <- sextant_window_stats(window,
                                           centre_lon_idx = site_row$lon_idx - lon_lo + 1,
                                           centre_lat_idx = site_row$lat_idx - lat_lo + 1)

      site_records <- records_of_the_day |> filter(source == site_row$source, site == site_row$site)

      for(j in seq_len(nrow(site_records))){

        rec <- site_records[j, ]

        summary_rows[[length(summary_rows) + 1]] <- cbind(
          data.frame(SOURCE = rec$source, SITE = rec$site, REGION = gsub("_", " ", rec$zone),
                    LATITUDE = site_row$lat, LONGITUDE = site_row$lon,
                    Data_source = "SEXTANT", sensor_name = "merged", atmospheric_correction = "Standard",
                    Satellite_variable = "SPM", Temporal_resolution = "DAILY", Satellite_algorithm = "SPIM",
                    DATE = format(the_date, "%Y-%m-%d"), Satellite_time = NA_character_, Insitu_time = NA_character_,
                    Insitu_variable = rec$Insitu_variable, Insitu_value = rec$value, check.names = FALSE),
          window_stats)
      }
    }

    nc_close(nc)
  }

  summary_df <- bind_rows(summary_rows)

  dir.create("output/MATCH_UP_DATA/FRANCE", recursive = TRUE, showWarnings = FALSE)
  write_csv(summary_df, "output/MATCH_UP_DATA/FRANCE/summary.csv")

  summary_df
}

compile_sextant_validation_summary()


# Validation stats + figures  --------------------------------------------

# SEXTANT
validate_sensor("SEXTANT", "all")
validate_sensor("SEXTANT", "small")

# MODIS
validate_sensor("MODIS", "all")
validate_sensor("MODIS", "small")

# MERIS
validate_sensor("MERIS", "all")
validate_sensor("MERIS", "small")

# OLCI-A
validate_sensor("OLCI-A", "all")
validate_sensor("OLCI-A", "small")

# OLCI-B
validate_sensor("OLCI-B", "all")
validate_sensor("OLCI-B", "small")


# Validation tables -------------------------------------------------------

# Load in situ site list
in_situ_site_list <- read_csv("metadata/in_situ_site_list.csv")

# Filter in situ stations to just those within a zone
zone_sites <- in_situ_site_list |> filter(!is.na(zone))

# Get stats files
files_stats <- dir("output/MATCH_UP_DATA/FRANCE/STATISTICS", pattern = "_stats_", full.names = TRUE)
files_stats_lm <- files_stats[grepl("_lm_", files_stats)]
files_stats_area <- files_stats[!files_stats %in% files_stats_lm]
files_stats_all <- files_stats_area[grepl("_all", files_stats_area)]
files_stats_small <-  files_stats_area[grepl("_small", files_stats_area)]

# Basic tables
zone_all_stats <- map_dfr(files_stats_all, read_csv, show_col_types = FALSE)
zone_small_stats <- map_dfr(files_stats_small, read_csv, show_col_types = FALSE)
write_csv(zone_all_stats, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_all.csv")

# Get just Southern Brittany stats
zone_SB_stats <- zone_all_stats |> 
  filter(zone == "SOUTHERN_BRITTANY")
write_csv(zone_SB_stats, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_all_SB.csv")
zone_SB_stats_global <- zone_SB_stats |> 
  filter(source == "ALL", site == "ALL", season == "ALL") |> 
  arrange(Error) |> 
  dplyr::select(zone, sensor, correction, processing, grid_size, source, site, season, variable, variable_sat, n, Slope, Slope_log, Bias, Error, MAPE, MSA, RMSE)
write_csv(zone_SB_stats_global, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_global_SB.csv")
zone_SB_sites <- zone_sites |> filter(zone == "SOUTHERN_BRITTANY")
write_csv(zone_SB_sites, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_sites_SB.csv")

# Get global stats
zone_all_stats_global <- zone_all_stats |> 
  filter(zone == "GLOBAL", source == "ALL", site == "ALL", season == "ALL") |> 
  # filter(variable != "TEMP") |> # Remove this as it's not ocean colour related
  # dplyr::select(sensor, variable, variable_sat, n, Slope_log, Bias, Error) |> 
  arrange(Error) |> 
  mutate(Error = round(Error),
         Bias = round(Bias))
write_csv(zone_all_stats_global, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_global_all.csv")
zone_small_stats_global <- zone_small_stats |> 
  filter(zone == "GLOBAL", source == "ALL", site == "ALL", season == "ALL") |> 
  # filter(variable != "TEMP") |> # Remove this as it's not ocean colour related
  # dplyr::select(sensor, variable, variable_sat, n, Slope_log, Bias, Error) |> 
  arrange(Error) |> 
  mutate(Error = round(Error),
         Bias = round(Bias))
write_csv(zone_small_stats_global, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_global_small.csv")

# Get high level stats
zone_all_stats_top <- zone_all_stats |> 
  filter(site == "ALL", season == "ALL")

# Fancy tables
validation_tables("output/MATCH_UP_DATA/FRANCE/STATISTICS/", "SEXTANT")


# Time series analyses ----------------------------------------------------

# Load MODIS SST and in situ temperatures
zone_median_SST <- map_dfr(dir("output/MATCH_UP_DATA/FRANCE", pattern = "zone_median_MODIS_", 
                               full.names = TRUE), data.table::fread) |> 
  mutate(date = as.Date(date)) |> 
  filter(variable == "SST")
# write_csv(zone_median_SST, "output/MATCH_UP_DATA/FRANCE/zone_MODIS_SST.csv")
# zone_median_SST <- read_csv("output/MATCH_UP_DATA/FRANCE/zone_MODIS_SST.csv")

# Load in situ data
zone_data_in_situ <- read_csv("data/INSITU_data/zone_data_in_situ.csv")

# List number of unique sites per source
zone_data_in_situ |> dplyr::select(source, site) |> distinct() |> 
  summarise(site_count = n(), .by = "source")

# List date range for each source of in situ data
zone_data_in_situ |> dplyr::select(source, date) |> 
  summarise(min_date = min(date), max_date = max(date), .by = "source")

# List the variables used in the in situ data
zone_data_in_situ |> dplyr::select(source, variable) |> distinct()

# Join
zone_all_TEMP <- zone_data_in_situ |> 
  filter(variable == "TEMP") |> 
  group_by(site) |> 
  mutate(count = n()) |> 
  ungroup() |> 
  filter(count > 90) |> 
  mutate(variable_sat = "SST") |> 
  left_join(zone_median_SST, by = c("zone", "source", "site", "date", "variable_sat" = "variable")) |> 
  filter(value > 0, median > 0) |> 
  dplyr::rename(value_in_situ = value, value_satellite = median) |> 
  mutate(season = case_when(
    month(date) %in% c(12, 1, 2) ~ "Winter", month(date) %in% 3:5  ~ "Spring",
    month(date) %in% 6:8  ~ "Summer", month(date) %in% 9:11 ~ "Autumn"), .after = "date") |> 
  dplyr::select(zone, dplyr::everything())

# Basic plot
plot_temp_compare <- zone_all_TEMP |> 
  mutate(zone_pretty = factor(zone_pretty,
                              levels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion"))) |> 
  filter(date >= "2004-07-04") |> 
  ggplot(aes(x = date, y = value_in_situ)) +
  geom_point(aes(group = site, shape = source), show.legend = FALSE) +
  geom_point(aes(y = value_satellite, group = site, shape = source), colour = "red", show.legend = FALSE) +
  geom_smooth(method = "lm", colour = "black", se = FALSE) +
  geom_smooth(aes(y = value_satellite), method = "lm", colour = "red", se = FALSE) +
  facet_wrap(~zone_pretty) +
  labs(title = "Comparison of in situ TEMP (black) and MODIS SST (red)",
       y = "Temperature (°C)", x = NULL) +
  theme(plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 23),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 23),
        legend.title = element_text(size = 23),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA))
plot_temp_compare
ggsave("figures/validation/comparison_TEMP_SST.png", height = 10, width = 16)

