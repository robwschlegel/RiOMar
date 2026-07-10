# func/flow.R
# Comparisons of river flow against plume size
#
# DEPRECATED (2026-07) -- merged into func/multi.R as part of consolidating
# all the individual driver-vs-plume analysis scripts (flow.R, wind.R,
# tide.R, ROFI.R, surface.R) into one deduplicated module. See the
# "REFACTOR (2026-07)" comment block at the top of func/multi.R for what
# changed and why.
#
# Equivalents in func/multi.R:
#   flow_comp(mouth_row)             -> combine_plume_driver("flow", meta) + plot_driver_plume_comparison(df, "flow", meta$mouth_name)
#   flow_trend(mouth_row)            -> superseded by flow_plume_trend_plus() below; not carried forward separately
#   flow_plume_trend_plus(mouth_row) -> driver_plume_trend(df, "flow", meta$mouth_name)
#   plyr::d_ply(river_mouths, "row_name", flow_comp) etc. -> run_driver_suite("flow")
#
# The pre-refactor version of this file is preserved in git history
# (`git log -- func/flow.R`) if anything here needs to be recovered verbatim.

source("func/multi.R")
