# func/figure.R
# Loader for the publication figure functions. Split 2026-09-28 into topic files under
# func/sections/ -- each part is an unchanged, contiguous slice of the former
# single figure.R, and they are sourced here in their original order, so
# source("func/figure.R") behaves exactly as before. Edit the parts;
# `python tools/check_split.py` verifies the parts still reassemble.

source("func/sections/figure_1_setup_utils.R", local = TRUE)  # Libraries and plotting utilities (was lines 1-261)
source("func/sections/figure_2_maps_methodology.R", local = TRUE)  # Maps (Figure 1), validation and methodology panels (was lines 262-702)
source("func/sections/figure_3_timeseries.R", local = TRUE)  # Plume time series, seasonal heatmap, X11 interannual vs flow (was lines 703-1019)
source("func/sections/figure_4_drivers.R", local = TRUE)  # Driver rose and GAM partial effects (was lines 1020-1165)
source("func/sections/figure_5_x11_thresholds.R", local = TRUE)  # X11 components and dynamic-vs-static threshold comparisons (was lines 1166-1471)
