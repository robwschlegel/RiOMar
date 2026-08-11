# func/export_driver_x11_inputs.R
#
# Writes the daily plume-area + driver series func/compute_driver_x11_correlation_table.py
# and func/compute_driver_x11_figures.py need, so the Python-side weekly
# Census-X11 decomposition (func/X11.py) can run without a Python port of
# func/multi.R's driver-loading/joining pipeline (combine_plume_driver(),
# load_driver(), get_zone_meta(), etc.) -- those stay R-only and unmodified,
# called here exactly as they already are everywhere else in the pipeline.
#
# One CSV per zone x driver: output/STATS/driver_x11_inputs/<zone>_<driver>.csv
# (columns: date, plume_area, value).
#
# Run from repo root: Rscript func/export_driver_x11_inputs.R

source("func/multi.R")

OUT_DIR <- "output/STATS/driver_x11_inputs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

other_drivers <- c("wind", "tide", "wave", "current")

for(zone_name in zones){
  meta <- get_zone_meta(zone_name = zone_name)
  for(driver_name in other_drivers){
    df <- combine_plume_driver(driver_name, meta) |>
      dplyr::select(date, plume_area, value)
    out_file <- file.path(OUT_DIR, paste0(zone_name, "_", driver_name, ".csv"))
    readr::write_csv(df, out_file)
    message("Wrote ", out_file, " (", nrow(df), " rows)")
  }
}
