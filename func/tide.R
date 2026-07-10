# func/tide.R
# Comparisons of tides against plume size
#
# DEPRECATED (2026-07) -- merged into func/multi.R as part of consolidating
# all the individual driver-vs-plume analysis scripts (flow.R, wind.R,
# tide.R, ROFI.R, surface.R) into one deduplicated module. See the
# "REFACTOR (2026-07)" comment block at the top of func/multi.R for what
# changed and why.
#
# Equivalents in func/multi.R:
#   tide_calc(mouth_row) -> combine_plume_driver("tide", meta) + plot_driver_plume_comparison(df, "tide", meta$mouth_name) + plot_driver_plume_dual_axis(df, "tide", meta$zone)
#   plyr::d_ply(river_mouths, "row_name", tide_calc) -> run_driver_suite("tide")
#
# The pre-refactor version of this file is preserved in git history
# (`git log -- func/tide.R`) if anything here needs to be recovered verbatim.

source("func/multi.R")
