# func/util.R
# Loader for the shared helper functions. Split 2026-09-28 into topic files under
# func/sections/ -- each part is an unchanged, contiguous slice of the former
# single util.R, and they are sourced here in their original order, so
# source("func/util.R") behaves exactly as before. Edit the parts;
# `python tools/check_split.py` verifies the parts still reassemble.

source("func/sections/util_1_metadata.R", local = TRUE)  # Zone metadata and manuscript figure/table/paragraph registries (was lines 1-79)
source("func/sections/util_2_tide_qc.R", local = TRUE)  # Tide gauge sub-daily QC (was lines 80-272)
source("func/sections/util_3_pixels.R", local = TRUE)  # Satellite pixel extraction (was lines 273-672)
source("func/sections/util_4_loading.R", local = TRUE)  # Data loading (plume, river flow, drivers) (was lines 673-934)
source("func/sections/util_5_statistics.R", local = TRUE)  # Statistics helpers (was lines 935-1047)
source("func/sections/util_6_plotting.R", local = TRUE)  # Plotting helpers (was lines 1048-1807)
source("func/sections/util_7_tables.R", local = TRUE)  # Table helpers (was lines 1808-1873)
