# func/ROFI.R
# The analyses used to compare ROFI to plume size
#
# DEPRECATED (2026-07) -- merged into func/multi.R as part of consolidating
# all the individual driver-vs-plume analysis scripts (flow.R, wind.R,
# tide.R, ROFI.R, surface.R) into one deduplicated module. See the
# "REFACTOR (2026-07)" comment block at the top of func/multi.R for what
# changed and why.
#
# Equivalents in func/multi.R:
#   comp_ROFI_plume(mouth_row)                 -> combine_plume_driver("current", meta) + plot_driver_plume_comparison(df, "current", meta$mouth_name)
#   daily/monthly/annual (lag) correlation      -> driver_plume_timesteps() / driver_plume_correlation() (now used for every driver, not just currents)
#   walk(zones[1:3], comp_ROFI_plume) (no Gulf of Lion ROFI data) -> run_driver_suite("current") (same zones[1:3] restriction, applied automatically)
#
# The pre-refactor version of this file is preserved in git history
# (`git log -- func/ROFI.R`) if anything here needs to be recovered verbatim.

source("func/multi.R")
