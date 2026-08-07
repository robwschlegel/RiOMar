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
from figure import (Figure_1, Figure_2, Figure_3, Figure_3_panels, Figure_3_zone_maps, Figure_5_driver_comparison,
                    Figure_7_driver_rose, Figure_9_gam_partial,
                    Figure_4_S1_timeseries, Figure_X11_weekly_results, Figure_S3_seasonal_boxplots)

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')


# =============================================================================
# ### Main figures (dynamic threshold)
# =============================================================================

Figure_1(where_are_saved_satellite_data = "../pCloudDrive/data",
         where_are_saved_regional_maps = "output",
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

# Figure 5: Seasonal analysis comprison of plume properties and drivers
# My current thinking is to provide proportional boxplots of all plume properties and drivers per zone per month
# Also somehow want to add in the trends per plume propoerty and driver for each monthly group

# Figure_5_driver_comparison(where_to_save_the_figure = "figures") # Deprecated code

# Figure 6 (X11 interannual signal, dynamic threshold) + Figure S3 (X11
# seasonal + residual, dynamic threshold) + Figure S6/S7 (static-threshold
# equivalents, supplementary)
Figure_X11_weekly_results(where_are_saved_X11_results_dynamic = "output/panache/dynamic",
                          where_are_saved_X11_results_static = "output/panache/static",
                          where_to_save_the_figure = "figures")

# manuscript Figures 7 and 9: wind/wave roses and GAM partial-dependence
# curves. Figure 8 (flow-controlled wave/wind-category scatter) removed
# 2026-08-05 (per Robert): its information is now covered by Figure 7.
Figure_7_driver_rose(where_to_save_the_figure = "figures")
Figure_9_gam_partial(where_to_save_the_figure = "figures")


# =============================================================================
# ### Supplementary figures
# =============================================================================

# Figure S3: seasonal (JJA vs. NDJ) boxplots of plume metrics, dynamic vs.
# static threshold
Figure_S3_seasonal_boxplots(where_to_save_the_figure = "figures")

# Figure S4a/b: flow/plume/ROFI succession + lagged correlation. Previously
# run by hand -- see func/ROFI.R::run_rofi_plume_succession().
rofi_R_path = os.path.join(func_dir, 'ROFI.R')
robjects.r['source'](rofi_R_path)
robjects.r['run_rofi_plume_succession']()

