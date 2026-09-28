# func/multi.R
# Loader for the multi-driver analysis functions. Split 2026-09-28 into topic files under
# func/sections/ -- each part is an unchanged, contiguous slice of the former
# single multi.R, and they are sourced here in their original order, so
# source("func/multi.R") behaves exactly as before. Edit the parts;
# `python tools/check_split.py` verifies the parts still reassemble.

source("func/sections/multi_1_setup.R", local = TRUE)  # Analysis notes, libraries, zone/gauge metadata (was lines 1-126)
source("func/sections/multi_2_drivers.R", local = TRUE)  # Driver loading, plume+driver joins, multi-timestep correlation (was lines 127-258)
source("func/sections/multi_3_plotting.R", local = TRUE)  # Driver comparison plotting (was lines 259-831)
source("func/sections/multi_4_surface_missing.R", local = TRUE)  # Surface maps and missing-data summaries (runs on source) (was lines 832-907)
source("func/sections/multi_5_rhone.R", local = TRUE)  # Rhone-only analyses (was lines 908-1499)
