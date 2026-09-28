# func/util.R
# The storage point for many functions re-used by other scripts


# Meta-data ---------------------------------------------------------------

# Project-wide settings (data root, zones...) from metadata/riomar_config.yml
source("func/config.R")

# Canonical zone order -- north to south
ZONE_ORDER <- c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION")

order_zones <- function(zone_vector){
  zone_vector[order(match(zone_vector, ZONE_ORDER))]
}

# Generated from panache.utils.define_parameters()
# func/util.py::export_panache_zone_metadata()
# Rerun that function whenever panache's zone parameters change
# do not hand-edit: metadata/panache_zone_metadata.csv.
panache_zone_metadata <- read_csv("metadata/panache_zone_metadata.csv", show_col_types = FALSE)

river_mouths <- panache_zone_metadata |>
  dplyr::transmute(row_name = dplyr::row_number(), mouth_name, mouth_lon, mouth_lat)

# Zone bounding boxes -- matched to panache library
zones_bbox <- panache_zone_metadata |>
  dplyr::select(zone, lon_min, lon_max, lat_min, lat_max) |>
  mutate(zone_pretty = factor(zone,
                              levels = ZONE_ORDER,
                              labels = c("Bay of Seine", "S. Brittany",
                                         "Gironde shelf", "Rhône shelf")), .after = "zone") |>
  dplyr::arrange(match(zone, ZONE_ORDER))

# Create FRANCE bounding box with same structure as zones_bbox
france_bbox <- data.frame(zone = "FRANCE",
                          lon_min = c(-7.8),
                          lon_max = c(10.3),
                          lat_min  = c(41.2),
                          lat_max = c(51.5))


# Manuscript figure/table registry ------------------------------------------
# Single source of truth for "what manuscript slot is this, what number does
# it currently have, which R function renders it" -- see
# metadata/figure_table_registry.csv. Renumbering a figure/table is a
# one-row edit to that CSV; figure.R functions look up their output folder
# via get_registry_row() instead of hardcoding a "FIGURE_N"/"TABLE_N" string.
figure_table_registry <- read_csv("metadata/figure_table_registry.csv", show_col_types = FALSE)

# Single source of truth for "which script/data file produced the numbers in
# this manuscript paragraph" -- see metadata/paragraph_source_registry.csv.
# Checked by metadata/make_figures_tables.R's make_all_figures_tables(); read here
# rather than there so the checklist script can source() this file the same
# way it already does for figure_table_registry above.
paragraph_source_registry <- read_csv("metadata/paragraph_source_registry.csv", show_col_types = FALSE)

get_registry_row <- function(slot_key){
  row <- dplyr::filter(figure_table_registry, .data$slot_key == !!slot_key)
  if(nrow(row) != 1) stop("get_registry_row(): expected exactly 1 row for slot_key '",
                          slot_key, "', found ", nrow(row))
  row
}

# "FIGURE_S1" -> "Figure_S1", matching every figure's filename convention
# (Figure_1.png, Figure_S1.png, ...) from its registry output_subdir. Use
# this for save helpers (e.g. save_plot_as_png()) that append their own
# extension; use registry_filename() below when the full filename is needed.
registry_basename <- function(output_subdir){
  sub("^FIGURE_", "Figure_", output_subdir)
}

# "FIGURE_S1" -> "Figure_S1.png", matching every figure's filename convention
# (Figure_1.png, Figure_S1.png, ...) from its registry output_subdir.
registry_filename <- function(output_subdir, ext = "png"){
  paste0(registry_basename(output_subdir), ".", ext)
}


