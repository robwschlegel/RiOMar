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

# The shared functions
source("func/util.R")


# Satellite times ---------------------------------------------------------

# Get all satellite daily start and end times when overhead France
# This could be used to correctly filter the in situ data by time when possible
## SEXTANT
### NB: Handled differently because the values are always the same and downloading all files takes an hour...
# times_SEXTANT <- data.frame(date_start = ymd(unlist(lapply(str_split(basename(dir("~/pCloudDrive/data/SEXTANT/SPM",
#                                                                                   recursive = TRUE, pattern = "SPIM")), "-"), "[[", 1)))) |>
#   mutate(date_end = date_start + 1,
#          start_time = ymd_hms(paste(date_start, "12:00:00"), tz = "GMT"),
#          end_time = ymd_hms(paste(date_end, "12:00:00"), tz = "GMT")) |>
#   dplyr::select(start_time, end_time)
# save(times_SEXTANT, file = "metadata/times_SEXTANT.RData")
load("metadata/times_SEXTANT.RData")

## MODIS
# times_MODIS <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily",
#                                pattern = "SPM", full.names = TRUE), get_start_end_time)
# save(times_MODIS, file = "metadata/times_MODIS.RData")
load("metadata/times_MODIS.RData")

### Double check that times are the same for any region
# Yes, start and end times are exactly the same
# times_MODIS_check <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/GULF_OF_LION/daily", 
#                                      pattern = "SPM", full.names = TRUE), get_start_end_time)
# times_MODIS_test <- bind_cols(times_MODIS, times_MODIS_check) |>
#   mutate(start_diff = start_time...1 - start_time...3, end_diff = end_time...2 - end_time...4)

### SST-NIGHT are meant to have a different time
# times_MODIS_SST_NIGHT <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily",
#                                          pattern = "SST-NIGHT", full.names = TRUE), get_start_end_time)
# save(times_MODIS_SST_NIGHT, file = "metadata/times_MODIS_SST_NIGHT.RData")
load("metadata/times_MODIS_SST_NIGHT.RData")
### SST day time
## These are calculated on rs-pro2 where the expert data are located
load("metadata/times_MODIS_SST.RData")

## MERIS
# times_MERIS <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/BAY_OF_SEINE/daily",
#                                pattern = "SPM", full.names = TRUE), get_start_end_time)
# save(times_MERIS, file = "metadata/times_MERIS.RData")
load("metadata/times_MERIS.RData")

## OLCI-A
# times_OLCI_A <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-A/BAY_OF_SEINE/daily",
#                                 pattern = "SPM", full.names = TRUE), get_start_end_time)
# save(times_OLCI_A, file = "metadata/times_OLCI_A.RData")
load("metadata/times_OLCI_A.RData")

## OLCI-B
# times_OLCI_B <- future_map_dfr(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-B/BAY_OF_SEINE/daily",
#                                 pattern = "SPM", full.names = TRUE), get_start_end_time)
# save(times_OLCI_B, file = "metadata/times_OLCI_B.RData")
load("metadata/times_OLCI_B.RData")

## ALL
# times_ALL <- bind_rows(mutate(times_SEXTANT, sat_name = "SEXTANT"),
#                        mutate(times_MODIS, sat_name = "MODIS"),
#                        mutate(times_MODIS_SST, sat_name = "MODIS (SST)"),
#                        mutate(times_MODIS_SST_NIGHT, sat_name = "MODIS (SST-NIGHT)"),
#                        mutate(times_MERIS, sat_name = "MERIS"),
#                        mutate(times_OLCI_A, sat_name = "OLCI-A"),
#                        mutate(times_OLCI_B, sat_name = "OLCI-B")) |>
#   mutate(date = as.Date(start_time),
#          start_hms = hms::as_hms(start_time),
#          end_hms = hms::as_hms(end_time)) |>
#   complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA), nesting(sat_name))
# save(times_ALL, file = "metadata/times_ALL.RData") # Saving in R format to maintain time stamps
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
  mutate(zone = case_when(lon >= zones_bbox$lon_min[zones_bbox$zone == "GULF_OF_LION"] & lon <= zones_bbox$lon_max[zones_bbox$zone == "GULF_OF_LION"] &
                            lat >= zones_bbox$lat_min[zones_bbox$zone == "GULF_OF_LION"] & lat <= zones_bbox$lat_max[zones_bbox$zone == "GULF_OF_LION"] ~ "GULF_OF_LION",
                          lon >= zones_bbox$lon_min[zones_bbox$zone == "BAY_OF_SEINE"] & lon <= zones_bbox$lon_max[zones_bbox$zone == "BAY_OF_SEINE"] &
                            lat >= zones_bbox$lat_min[zones_bbox$zone == "BAY_OF_SEINE"] & lat <= zones_bbox$lat_max[zones_bbox$zone == "BAY_OF_SEINE"] ~ "BAY_OF_SEINE",
                          lon >= zones_bbox$lon_min[zones_bbox$zone == "BAY_OF_BISCAY"] & lon <= zones_bbox$lon_max[zones_bbox$zone == "BAY_OF_BISCAY"] &
                            lat >= zones_bbox$lat_min[zones_bbox$zone == "BAY_OF_BISCAY"] & lat <= zones_bbox$lat_max[zones_bbox$zone == "BAY_OF_BISCAY"] ~ "BAY_OF_BISCAY",
                          lon >= zones_bbox$lon_min[zones_bbox$zone == "SOUTHERN_BRITTANY"] & lon <= zones_bbox$lon_max[zones_bbox$zone == "SOUTHERN_BRITTANY"] &
                            lat >= zones_bbox$lat_min[zones_bbox$zone == "SOUTHERN_BRITTANY"] & lat <= zones_bbox$lat_max[zones_bbox$zone == "SOUTHERN_BRITTANY"] ~ "SOUTHERN_BRITTANY")) |>
  mutate(zone_pretty = factor(zone,
                              levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                              labels = c("Bay of Seine", "S. Brittany", "Gironde shelf", "Rhône shelf")), .after = "zone") |>
  mutate(source = factor(source, levels = c("SOMLIT", "REPHY")))
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
# Tutorial here : https://www.etiennebacher.com/posts/2021-12-27-mapping-french-rivers-network/
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
  # Non-zone stations
  geom_point(data = filter(in_situ_site_list, is.na(zone)), show.legend = FALSE,
             aes(x = lon, y = lat, shape = source), color = "red", size = 3) +
  # REPHY points
  geom_point(data = filter(in_situ_site_list, !is.na(zone), source == "REPHY"), 
             aes(x = lon, y = lat, shape = source, color = zone_pretty), size = 3) +
  # SOMLIT points
  geom_point(data = filter(in_situ_site_list, !is.na(zone), source == "SOMLIT"),
             aes(x = lon, y = lat, shape = source), color = "black", size = 6) +
  geom_point(data = filter(in_situ_site_list, !is.na(zone), source == "SOMLIT"), 
             aes(x = lon, y = lat, shape = source, color = zone_pretty), size = 5) +
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

