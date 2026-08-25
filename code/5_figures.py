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

# Writes its regional-zone-maps panel into the plume_methodology_panel
# slot's output folder (see manuscript/figure_table_registry.csv)
Figure_3_zone_maps(where_are_saved_regional_maps = "output",
                   where_to_save_the_figure = "figures")

# plume_methodology_panel: plume-detection methodology composite
# (Figure_3_panels()' A-D panels + Figure_3_zone_maps()' zone-maps panel);
# must run after both
Figure_3(where_to_save_the_figure = "figures")

# plume_area_timeseries (plume area + mean SPM time series) +
# thresholds_comparison (dynamic vs. static threshold comparison)
Figure_4_S1_timeseries(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                       where_are_saved_plume_results_with_fixed_threshold = "output/panache/static",
                       where_to_save_the_figure = "figures")

# seasonal_boxplot_heatmap (sec:results_seasonal): monthly boxplots of all
# plume properties and drivers, per zone, dynamic threshold; also writes the
# shared dynamic+static data that Figure_S3_seasonal_boxplots() below
# re-reads, so it must run before that call. Per-month trends themselves are
# computed separately by func/compute_seasonal_trend.R (feeds the
# monthly_trends_table_main slot), not by this figure.
Figure_5_seasonal_analysis(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                           where_are_saved_plume_results_with_static_threshold = "output/panache/static",
                           where_to_save_the_figure = "figures")

# x11_interannual_river_flow (X11 interannual signal, dynamic threshold, vs.
# river flow) + x11_components_dynamic (X11 seasonal + residual, dynamic
# threshold, vs. river flow) + x11_interannual_dynamic_vs_static /
# x11_seasonal_dynamic_vs_static / x11_residual_dynamic_vs_static (dynamic
# vs. static threshold comparison of plume area itself)
Figure_X11_weekly_results(where_are_saved_X11_results_dynamic = "output/panache/dynamic",
                          where_are_saved_X11_results_static = "output/panache/static",
                          where_to_save_the_figure = "figures")

# driver_rose_diagram: wind/wave roses
Figure_7_driver_rose(where_to_save_the_figure = "figures")

# gam_partial_effects: GAM partial-dependence curves
Figure_8_gam_partial(where_to_save_the_figure = "figures")


# =============================================================================
# ### Supplementary figures
# =============================================================================
# See manuscript/figure_table_registry.csv for every slot's current number.

# seasonal_boxplots_dynamic_vs_static: monthly boxplots of plume properties
# and drivers, dynamic vs. static threshold
Figure_S3_seasonal_boxplots(where_to_save_the_figure = "figures")

# daily_flow_lagged_correlation ("Sx. Lagged daily correlations"): daily
# plume area vs. river flow scatter + lagged correlation, per zone.
Figure_S_daily_flow(where_to_save_the_figure = "figures")

# rofi_succession + rofi_lagged_correlation: flow/plume/ROFI succession +
# lagged correlation. Previously run by hand -- see
# func/ROFI.R::run_rofi_plume_succession().
rofi_R_path = os.path.join(func_dir, 'ROFI.R')
robjects.r['source'](rofi_R_path)
robjects.r['run_rofi_plume_succession']()

