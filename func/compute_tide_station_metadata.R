# One-off: derive the tide-gauge station list (station name, coordinates,
# date range) from the SHOM SensorML metadata files under
# data/TIDES/<station>/*.sml (name, coordinates) and each station's own
# *.txt archive filenames (start/end year of the record actually downloaded
# and usable here). Originally fed a dedicated tide-gauge table in
# manuscript.tex; that table was removed in favour of one sentence of prose
# (Section~\ref{sec:tide_gauges}), so this script is kept only for the
# station-metadata derivation logic, not as an active table generator.
#
# NB: the .sml file's own date_prem_obs field is NOT used for the start
# year -- for Saint-Nazaire it claims 1821, a much earlier historical
# origin than what SHOM has actually made available, so the
# archived .txt file years are the reliable source for what's usable.
#
# Run from repo root: Rscript func/compute_tide_station_metadata.R
library(dplyr)
library(xml2)

sml_files <- dir("data/TIDES", pattern = "\\.sml$", recursive = TRUE, full.names = TRUE)

station_metadata <- purrr::map_dfr(sml_files, function(sml_path) {

  doc <- read_xml(sml_path)
  ns <- xml_ns(doc)

  field_value <- function(field_name) {
    xml_text(xml_find_first(doc, paste0(".//swe:field[@name='", field_name, "']/*/swe:value"), ns))
  }

  station_dir <- dirname(sml_path)
  years <- as.integer(gsub(".*_(\\d{4})\\.txt$", "\\1", dir(station_dir, pattern = "_\\d{4}\\.txt$")))

  data.frame(
    station = field_value("spm"),
    longitude = as.numeric(field_value("longitude")),
    latitude = as.numeric(field_value("latitude")),
    archive_start_year = min(years, na.rm = TRUE),
    archive_end_year = max(years, na.rm = TRUE)
  )
}) |>
  arrange(station)

print(station_metadata, row.names = FALSE)
