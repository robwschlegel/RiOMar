# func/run_validate_stats_rerun.R
#
# Standalone rerun of func/validate.R's "Validation stats + figures" and
# "Validation tables" sections only -- NOT "Extract satellite data" (~1hr/
# variable, assumed already run and unchanged on disk -- only the zone
# display-name rename below needs picking up, not new match-up data; verify
# output/MATCH_UP_DATA/FRANCE/zone_median_SEXTANT* exists before running
# this, and if it doesn't, stop and rerun validate.R's full extraction step
# instead) and NOT "Map In situ" (needs library(sf), which conflicts with
# panache's already-loaded conda GDAL/GEOS/PROJ stack when run in-process
# via rpy2 -- see func/run_figure_1.R's identical comment). This script
# deliberately never loads sf, so it's safe to source via rpy2 too, but is
# still meant to be run as its own process for consistency with that
# established workaround pattern.
#
# Picks up the func/util.R zone-label rename (Rhone shelf/Gironde shelf,
# 2026-09-04) in the on-disk validation STATISTICS CSVs and figures/
# validation/*.png outputs, which were last regenerated before that rename.
#
# Regenerates data/INSITU_data/zone_data_in_situ.csv first (validate.R's
# "Load in situ" section's logic, inlined below -- cheap, no sf, and
# validate_sensor() reads this file directly from disk with no fallback) so
# this script has no in-memory dependency on validate.R's earlier
# Setup/Satellite times/Load in situ sections ever having run in this R
# session.
#
# Usage: Rscript func/run_validate_stats_rerun.R  (no arguments)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0) {
  stop("Usage: Rscript func/run_validate_stats_rerun.R  (no arguments)")
}

library(tidyverse)
library(furrr)
library(scales)
library(ggExtra)
library(gt)
source("func/util.R")

# --- Load in situ (validate.R's "Load in situ" section, unchanged logic, no sf) ---

REPHY <- read.csv2("data/INSITU_data/REPHY/Table1_REPHY_hydro_20250408.csv", fileEncoding = "ISO-8859-1") |>
  mutate(lon = as.numeric(lon), lat = as.numeric(lat), source = "REPHY") |>
  dplyr::rename(site = Code_point_Libelle, date = Date, time = Heure, variable = Code.parametre, value = Valeur_mesure)
SOMLIT <- read_csv("data/INSITU_data/SOMLIT/Somlit_clean.csv", show_col_types = FALSE) |> mutate(source = "SOMLIT")

in_situ_site_list <- bind_rows(dplyr::select(REPHY, source, lon, lat, site),
                               dplyr::select(SOMLIT, source, lon, lat, site)) |>
  distinct() |>
  summarise(lon = mean(lon), lat = mean(lat), .by = c("source", "site")) |>
  mutate(zone = case_when(
    lon >= zones_bbox$lon_min[zones_bbox$zone == "GULF_OF_LION"] & lon <= zones_bbox$lon_max[zones_bbox$zone == "GULF_OF_LION"] &
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

zone_sites <- in_situ_site_list |> filter(!is.na(zone))

REPHY_clean <- right_join(REPHY, zone_sites, by = c("source", "site", "lon", "lat")) |>
  filter(Qualite.resultat == "Bon") |> filter(as.numeric(Profondeur.metre) <= 10) |>
  dplyr::select(source, site, lon, lat, date, time, variable, value) |>
  mutate(variable = case_when(variable == "SALI" ~ "SAL", variable == "TURB" ~ "TUR",
                              variable == "CHLOROA" ~ "CHLA", TRUE ~ variable),
         date = as.Date(date), time = case_when(time == "" ~ NA, TRUE ~ time)) |>
  mutate(time = hms::as_hms(time)) |> filter(value >= 0)

SOMLIT_clean <- right_join(SOMLIT, zone_sites, by = c("source", "site", "lon", "lat")) |>
  filter(prof_num <= 10) |>
  mutate(TEMP = case_when(temp_QC %in% c(2, 6, 7) ~ temp), SAL = case_when(sal_QC %in% c(2, 6, 7) ~ sal),
         POC = case_when(POC_QC %in% c(2, 6, 7) ~ POC), SPM = case_when(SPM_QC %in% c(2, 6, 7) ~ SPM),
         CHLA = case_when(CHLA_QC %in% c(2, 6, 7) ~ CHLA)) |>
  dplyr::select(source, site, lon, lat, date, time, TEMP, SAL, POC, SPM, CHLA) |>
  pivot_longer(TEMP:CHLA, values_to = "value", names_to = "variable") |>
  mutate(time = hms::as_hms(time)) |> filter(value >= 0)

zone_data_in_situ <- bind_rows(REPHY_clean, SOMLIT_clean) |>
  filter(variable %in% c('TEMP', 'SAL', 'POC', 'SPM', 'CHLA', 'TUR')) |>
  mutate(time_class = case_when(time >= hms("10:00:00") & time <= hms("14:00:00") ~ "day", # Doesn't seem to like case_when() ...
                                variable == "TEMP" & time >= hms("00:00:00") | time <= hms("09:00:00") ~ "night")) |>
  filter(time_class %in% c("day", "night")) |> filter(!(variable != "TEMP" & time_class == "night")) |>
  summarise(value = mean(value, na.rm = TRUE), .by = c("source", "site", "lon", "lat", "date", "time_class", "variable")) |>
  left_join(zone_sites, by = join_by(source, site, lon, lat)) |>
  dplyr::select(zone, zone_pretty, source, everything())
write_csv(zone_data_in_situ, "data/INSITU_data/zone_data_in_situ.csv")

# --- Validation stats + figures (validate.R's own section, unchanged) -----

for (sat in c("SEXTANT", "MODIS", "MERIS", "OLCI-A", "OLCI-B")) {
  validate_sensor(sat, "all")
  validate_sensor(sat, "small")
}

# --- Validation tables (validate.R's own section, unchanged) --------------

in_situ_site_list <- read_csv("metadata/in_situ_site_list.csv", show_col_types = FALSE)
zone_sites <- in_situ_site_list |> filter(!is.na(zone))

files_stats <- dir("output/MATCH_UP_DATA/FRANCE/STATISTICS", pattern = "_stats_", full.names = TRUE)
files_stats_lm <- files_stats[grepl("_lm_", files_stats)]
files_stats_area <- files_stats[!files_stats %in% files_stats_lm]
files_stats_all <- files_stats_area[grepl("_all", files_stats_area)]
files_stats_small <- files_stats_area[grepl("_small", files_stats_area)]

zone_all_stats <- map_dfr(files_stats_all, read_csv, show_col_types = FALSE)
zone_small_stats <- map_dfr(files_stats_small, read_csv, show_col_types = FALSE)
write_csv(zone_all_stats, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_all.csv")

zone_SB_stats <- zone_all_stats |> filter(zone == "SOUTHERN_BRITTANY")
write_csv(zone_SB_stats, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_all_SB.csv")
zone_SB_stats_global <- zone_SB_stats |> filter(source == "ALL", site == "ALL", season == "ALL") |>
  arrange(Error) |> dplyr::select(zone, sensor, correction, processing, grid_size, source, site, season, variable, variable_sat, n, Slope, Slope_log, Bias, Error, MAPE, MSA, RMSE)
write_csv(zone_SB_stats_global, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_global_SB.csv")
zone_SB_sites <- zone_sites |> filter(zone == "SOUTHERN_BRITTANY")
write_csv(zone_SB_sites, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_sites_SB.csv")

zone_all_stats_global <- zone_all_stats |> filter(zone == "GLOBAL", source == "ALL", site == "ALL", season == "ALL") |>
  arrange(Error) |> mutate(Error = round(Error), Bias = round(Bias))
write_csv(zone_all_stats_global, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_global_all.csv")
zone_small_stats_global <- zone_small_stats |> filter(zone == "GLOBAL", source == "ALL", site == "ALL", season == "ALL") |>
  arrange(Error) |> mutate(Error = round(Error), Bias = round(Bias))
write_csv(zone_small_stats_global, "output/MATCH_UP_DATA/FRANCE/STATISTICS/table_global_small.csv")

validation_tables("output/MATCH_UP_DATA/FRANCE/STATISTICS/", "SEXTANT")

message("Done: validation stats + tables refreshed with current zone labels.")
