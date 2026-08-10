#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to create the figures used in the publication of this workflow.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl
import rpy2.robjects as robjects

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, figure
from figure import (Figure_1, Figure_2, Figure_3, Figure_3_panels, Figure_3_zone_maps, Figure_5_seasonal_analysis,
                    Figure_7_driver_rose, Figure_8_gam_partial, Figure_S_daily_flow,
                    Figure_4_S1_timeseries, Figure_X11_weekly_results, Figure_S3_seasonal_boxplots)

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')


# =============================================================================
# ### Main figures (dynamic threshold)
# =============================================================================

Figure_1(where_are_saved_satellite_data = "../pCloudDrive/data",
         where_to_save_the_figure = "figures")

Figure_2(where_to_save_the_figure = "figures")

Figure_3_panels(where_are_saved_regional_maps = "output",
                where_to_save_the_figure = "figures")

# Writes its regional-zone-maps panel into FIGURE_3
Figure_3_zone_maps(where_are_saved_regional_maps = "output",
                   where_to_save_the_figure = "figures")

# Figure 3: plume-detection methodology composite (Figure_3_panels()' A-D
# panels + Figure_3_zone_maps()' zone-maps panel); must run after both
Figure_3(where_to_save_the_figure = "figures")

# Figure 4 (plume area + mean SPM time series) + Figure S1 (dynamic vs.
# static threshold comparison)
Figure_4_S1_timeseries(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                       where_are_saved_plume_results_with_fixed_threshold = "output/panache/static",
                       where_to_save_the_figure = "figures")

# Figure 5 (sec:results_seasonal): monthly boxplots of all plume properties
# and drivers, per zone, dynamic threshold; also writes the shared
# dynamic+static data that Figure_S3_seasonal_boxplots() below re-reads, so
# it must run before that call. Per-month trends themselves are computed
# separately by func/compute_seasonal_trend.R (feeds manuscript Table 7),
# not by this figure.
Figure_5_seasonal_analysis(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                           where_are_saved_plume_results_with_static_threshold = "output/panache/static",
                           where_to_save_the_figure = "figures")

# Figure 6 (X11 interannual signal, dynamic threshold, vs. river flow) +
# Figure S2 (X11 seasonal + residual, dynamic threshold, vs. river flow) +
# Figure S5/S6/S6_residual (X11 interannual/seasonal/residual, dynamic vs.
# static threshold comparison of plume area itself)
Figure_X11_weekly_results(where_are_saved_X11_results_dynamic = "output/panache/dynamic",
                          where_are_saved_X11_results_static = "output/panache/static",
                          where_to_save_the_figure = "figures")

# manuscript Figures 7 and 9: wind/wave roses and GAM partial-dependence
# curves. Figure 8 (flow-controlled wave/wind-category scatter) removed
# 2026-08-05 (per Robert): its information is now covered by Figure 7.
Figure_7_driver_rose(where_to_save_the_figure = "figures")
Figure_8_gam_partial(where_to_save_the_figure = "figures")


# =============================================================================
# ### Supplementary figures
# =============================================================================

# Figure S3: monthly boxplots of plume properties and drivers, dynamic vs.
# static threshold
Figure_S3_seasonal_boxplots(where_to_save_the_figure = "figures")

# Fig. daily_flow ("Sx. Lagged daily correlations"): daily plume area vs.
# river flow scatter + lagged correlation, per zone. Restored 2026-08-07
# under a new name/folder (previously Figure_5_driver_comparison(), writing
# to FIGURE_5/ -- now repurposed for the new main-text Figure 5 above).
Figure_S_daily_flow(where_to_save_the_figure = "figures")

# Figure S4a/b: flow/plume/ROFI succession + lagged correlation. Previously
# run by hand -- see func/ROFI.R::run_rofi_plume_succession().
rofi_R_path = os.path.join(func_dir, 'ROFI.R')
robjects.r['source'](rofi_R_path)
robjects.r['run_rofi_plume_succession']()

