# func/surface.R
# Analyses specifically of surface data
#
# DEPRECATED (2026-07) -- merged into func/multi.R as part of consolidating
# all the individual driver-vs-plume analysis scripts (flow.R, wind.R,
# tide.R, ROFI.R, surface.R) into one deduplicated module. See the
# "REFACTOR (2026-07)" comment block at the top of func/multi.R for what
# changed and why.
#
# Equivalents in func/multi.R (kept in their own "Surface / pixel-level
# multi-driver maps" section since these operate per-pixel, not per-day, so
# they are conceptually distinct from the time-series functions above them):
#   surface_plot(zone)             -> surface_plot(zone) (same name, same behaviour; now uses zone_meta/load_driver() instead of an inline zone lookup)
#   surface_plot_daily_maps(zone)  -> surface_plot_daily_maps(zone) (unchanged)
#
# The pre-refactor version of this file is preserved in git history
# (`git log -- func/surface.R`) if anything here needs to be recovered verbatim.

source("func/multi.R")
