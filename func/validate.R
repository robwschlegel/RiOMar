# func/validate.R

# This script contains the code necessary to run validation on a number of
# satellite products and variables against multiple sources of in situ data

# Each section can run independently of the others

# Setup -------------------------------------------------------------------

# The packages needed for this script
library(tidyverse)
library(raster) # Used here to select specific pixels in a .nc file
library(ncdf4)
library(sf) # Used for complex shape files
library(scales) # For better plot labels
library(ggExtra) # For histogram border plots
library(gt) # For fancy tables
library(geosphere) # For pixel distances
library(doParallel); registerDoParallel(cores = detectCores()-2)

# The shared functions
source("func/util.R")


# Satellite times ---------------------------------------------------------

# Get all L3 satellite daily start and end times when overhead France
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

## MODIS
times_MODIS <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily",
                               pattern = "SPM", full.names = TRUE), get_start_end_time, .parallel = TRUE)
save(times_MODIS, file = "metadata/times_MODIS.RData")
# Double check that times are the same for any region
# Yes, start and end times are exactly the same
times_MODIS_check <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/GULF_OF_LION/daily", pattern = "SPM", full.names = TRUE),
                                 get_start_end_time, .parallel = TRUE)
times_MODIS_test <- bind_cols(times_MODIS, times_MODIS_check) |>
  mutate(start_diff = start_time...1 - start_time...3, end_diff = end_time...2 - end_time...4)

# ## MERIS
times_MERIS <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/BAY_OF_SEINE/daily",
                               pattern = "SPM", full.names = TRUE), get_start_end_time, .parallel = TRUE)
save(times_MERIS, file = "metadata/times_MERIS.RData")

## OLCI-A
times_OLCI_A <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-A/BAY_OF_SEINE/daily",
                                pattern = "SPM", full.names = TRUE), get_start_end_time, .parallel = TRUE)
save(times_OLCI_A, file = "metadata/times_OLCI_A.RData")

## OLCI-B
times_OLCI_B <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-B/BAY_OF_SEINE/daily",
                                pattern = "SPM", full.names = TRUE), get_start_end_time, .parallel = TRUE)
save(times_OLCI_B, file = "metadata/times_OLCI_B.RData")

## ALL
times_ALL <- bind_rows(mutate(times_SEXTANT, sat_name = "SEXTANT"),
                       mutate(times_MODIS, sat_name = "MODIS"),
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
  labs(y = "Passover time (UTC)", x = NULL,
       title = "SEXTANT and ODATIS-MR sensor availability and overhead times") +
  theme_bw()
# time_plot
ggsave("figures/all_sensors_time.png", width = 8, height = 4)


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
  mutate(zone = case_when(lon >= zones_bbox$lon_min[1] & lon <= zones_bbox$lon_max[1] &
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
write_csv(in_situ_site_list, "metadata/in_situ_site_list.csv")
# in_situ_site_list <- read_csv("metadata/in_situ_site_list.csv")

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
  # Filter from times 10:00 to 15:00
  filter(time >= hms("10:00:00"), time <= hms("15:00:00")) |> 
  # Create daily means
  summarise(value = mean(value, na.rm = TRUE), .by = c("source", "site", "lon", "lat", "date", "variable"))
write_csv(zone_data_in_situ, "data/INSITU_data/zone_data_in_situ.csv")


# Map In situ -------------------------------------------------------------

# Load high-res shape files
# Tuto here : https://www.etiennebacher.com/posts/2021-12-27-mapping-french-rivers-network/
if(!exists("borders_FR")) borders_FR <- read_sf("data/FRANCE_shapefile/gadm41_FRA_0.shp")
if(!exists("rivers_FR")) rivers_FR <- st_intersection(read_sf("data/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"), borders_FR)

# Map all in situ stations, highlighting the zones and the stations used
# TODO: Add text labels showing the number of REPHY and SOMLIT stations per zone
in_situ_station_map <- ggplot() +
  geom_sf(data = borders_FR, color = "black", fill = "sienna4", inherit.aes = FALSE) +
  geom_sf(data = rivers_FR, color = "lightblue", inherit.aes = FALSE, linewidth = 0.2) +
  geom_rect(data = zones_bbox, fill = NA, linewidth = 2, show.legend = FALSE,
            aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max, colour = zone_pretty)) +
  # coord_map("moll") +
  coord_sf(xlim = c(-5, 10), ylim = c(41.5, 51)) +
  geom_point(data = in_situ_site_list,
             aes(x = lon, y = lat, shape = source), color = "black", size = 4) +
  geom_point(data = filter(in_situ_site_list, is.na(zone)),
             aes(x = lon, y = lat, shape = source), color = "red", size = 3) +
  geom_point(data = filter(in_situ_site_list, !is.na(zone)),
             aes(x = lon, y = lat, shape = source, color = zone_pretty), size = 3) +
  # ggrepel::geom_text_repel(data =  filter(in_situ_site_list, source == "SOMLIT"), aes(x = lon, y = lat, label = site),
  #                 color = "red", size = 5, max.overlaps = 20) +
  # labs(title = paste(nrow(SOMLIT_stations), "SOMLIT stations"), x = NULL, y = NULL) +
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
ggsave("figures/map_in_situ_stations.png", height = 14, width = 15.5, bg = "white")


# Prep satellite pixels ---------------------------------------------------

# Create the data.frames of pixel matchups
## SEXTANT
write_pixels(zone_sites, "SEXTANT")

## MODIS
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "MODIS")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "MODIS")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "MODIS")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "MODIS")

## MERIS
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "MERIS")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "MERIS")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "MERIS")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "MERIS")

## OLCI-A
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "OLCI-A")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "OLCI-A")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "OLCI-A")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "OLCI-A")

## OLCI-B
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "OLCI-B")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "OLCI-B")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "OLCI-B")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "OLCI-B")


# Extract satellite data --------------------------------------------------

# SEXTANT
extract_pixels_all("SEXTANT")

# MODIS
# TODO: This can be optimized with a map or walk function
# But it isn't bad to have more control over the process...
# Though, one could walk through this entire section, all sensors and zones, in one function call...
# Could also use m_ply()
extract_pixels_all("MODIS", "BAY_OF_SEINE"); gc()
extract_pixels_all("MODIS", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("MODIS", "BAY_OF_BISCAY"); gc()
extract_pixels_all("MODIS", "GULF_OF_LION"); gc()

# MERIS
extract_pixels_all("MERIS", "BAY_OF_SEINE"); gc()
extract_pixels_all("MERIS", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("MERIS", "BAY_OF_BISCAY"); gc()
extract_pixels_all("MERIS", "GULF_OF_LION"); gc()

# OLCI-A
extract_pixels_all("OLCI-A", "BAY_OF_SEINE"); gc()
extract_pixels_all("OLCI-A", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("OLCI-A", "BAY_OF_BISCAY"); gc()
extract_pixels_all("OLCI-A", "GULF_OF_LION"); gc()

# OLCI-B
extract_pixels_all("OLCI-B", "BAY_OF_SEINE"); gc()
extract_pixels_all("OLCI-B", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("OLCI-B", "BAY_OF_BISCAY"); gc()
extract_pixels_all("OLCI-B", "GULF_OF_LION"); gc()


# Validation stats + figures  --------------------------------------------

validate_sensor("SEXTANT")
validate_sensor("MODIS")
validate_sensor("MERIS")
validate_sensor("OLCI-A")
validate_sensor("OLCI-B")


# Validation tables -------------------------------------------------------

validation_tables("output/MATCH_UP_DATA/FRANCE/STATISTICS/", "SEXTANT")

