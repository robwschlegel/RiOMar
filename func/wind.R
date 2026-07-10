# func/wind.R
# Comparisons of wind forcing against plume size
#
# DEPRECATED (2026-07) -- merged into func/multi.R as part of consolidating
# all the individual driver-vs-plume analysis scripts (flow.R, wind.R,
# tide.R, ROFI.R, surface.R) into one deduplicated module. See the
# "REFACTOR (2026-07)" comment block at the top of func/multi.R for what
# changed and why.
#
# Equivalents in func/multi.R:
#   spatial_wind_calc(mouth_row) -> combine_plume_driver("wind", meta) + plot_driver_plume_comparison(df, "wind", meta$mouth_name) + plot_driver_plume_dual_axis(df, "wind", meta$zone)
#   on-/off-shore classification -> wind_add_direction()
#   plyr::d_ply(river_mouths, "row_name", spatial_wind_calc) -> run_driver_suite("wind")
#
# NB: the unrelated "Etang de Berre" ERA5 exploratory section that was at the
# bottom of this file (a different embayment, not one of the 4 RiOMar zones)
# was NOT part of the driver-vs-plume comparison pipeline and was not carried
# into func/multi.R. Recover it from git history if still needed.
#
# The pre-refactor version of this file is preserved in git history
# (`git log -- func/wind.R`) if anything here needs to be recovered verbatim.

source("func/multi.R")
