#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to create the figures used in the publication of this workflow.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, figure
from figure import (Figure_1, Figure_2, Figure_3, Figure_4, Figure_5, Figure_5_driver_comparison,
                    Figure_7_driver_rose, Figure_8_driver_category, Figure_9_gam_partial,
                    Figure_6_7, Figure_X11_weekly_results, Figure_S4_seasonal_boxplots)

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')


# =============================================================================
# ### Main figures (dynamic threshold)
# =============================================================================

Figure_1(where_are_saved_satellite_data = "../pCloudDrive/data",
         where_to_save_the_figure = "figures")

Figure_2(where_to_save_the_figure = "figures")

Figure_4(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

# Figure_5() writes its regional-zone-maps panel into FIGURE_3/ (folded into
# the Figure 3 composite, not a standalone figure -- see figure.py's docstring)
Figure_5(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

# Figure 3: plume-detection methodology composite (FIGURE_4 panels A-D +
# Figure_5()'s zone-maps panel); must run after both Figure_4() and Figure_5()
Figure_3(where_to_save_the_figure = "figures")

# manuscript Figure 5: plume-area-vs-flow scatter + lagged correlation (distinct
# from the Figure_5() regional-zone-maps call above -- see figure.py's docstring)
Figure_5_driver_comparison(where_to_save_the_figure = "figures")

# Figure 6-7: compare dynamic vs static threshold plume detection
Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
           where_are_saved_plume_results_with_fixed_threshold = "output/panache/static",
           where_to_save_the_figure = "figures")

# manuscript Figures 7-9: wind/wave roses, flow-controlled wave/wind-category
# scatter, and GAM partial-dependence curves -- replace the original X11-
# decomposed wind/wave/tide time series concept (direction matters more than
# magnitude for wind/wave; tide has no meaningful long-term trend to show)
Figure_7_driver_rose(where_to_save_the_figure = "figures")
Figure_8_driver_category(where_to_save_the_figure = "figures")
Figure_9_gam_partial(where_to_save_the_figure = "figures")

# Figure 6 (X11 interannual signal, dynamic threshold) + Figure S3 (X11
# seasonal + residual, dynamic threshold) + Figure S6/S7 (static-threshold
# equivalents, supplementary) -- replaces the deprecated Figure_8_9_10(),
# whose static-threshold output was computed but never actually cited
# anywhere in the manuscript (2026-08-01 audit).
Figure_X11_weekly_results(where_are_saved_X11_results_dynamic = "output/panache/dynamic",
                          where_are_saved_X11_results_static = "output/panache/static",
                          where_to_save_the_figure = "figures")


# =============================================================================
# ### Supplementary figures
# =============================================================================

# Figure S4: seasonal (JJA vs. NDJ) boxplots of plume metrics, dynamic vs.
# static threshold -- migrated from manuscript/make_figures_tables.R into
# the real pipeline (2026-08-01).
Figure_S4_seasonal_boxplots(where_to_save_the_figure = "figures")
