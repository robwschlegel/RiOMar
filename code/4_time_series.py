#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to run the time series analyses.


# =============================================================================
#### Modules
# =============================================================================

import os
import sys
import matplotlib as mpl
import rpy2.robjects as robjects

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from X11 import Apply_X11_method_on_time_series, Apply_X11_method_on_time_series_per_river

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')

# The zones for mapping
zones_list = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY']

# Basic arguments to be used throughout the script
sextant_spm_all = {'Data_sources':['SEXTANT'],
                   'Sensor_names':["merged"],
                   'Satellite_variables':['SPM'],
                   'Atmospheric_corrections':['Standard'],
                   'Temporal_resolution':['DAILY'],
                   'start_day':'1998/01/01',
                   'end_day':'2025/12/31'}


# =============================================================================
# ### X11 analyses — dynamic threshold (main results)
# =============================================================================

# NB: X11 can only be used on weekly or monthly data, not daily
Apply_X11_method_on_time_series(sextant_spm_all,
                                Zones = zones_list,
                                plume_time_step = "WEEKLY",
                                plume_dir_in = "output/panache/dynamic",
                                X11_dir_out = "output/panache/dynamic",
                                include_river_flow = True)

# Per-river variant (metadata/river_discharge_mapping.csv): now that panache
# v5.0.0+ tracks individual river-mouth plume time series, also run the same
# X11 plume-vs-flow comparison once per river instead of only per zone.
Apply_X11_method_on_time_series_per_river(sextant_spm_all,
                                          plume_time_step = "WEEKLY",
                                          plume_dir_in = "output/panache/dynamic",
                                          X11_dir_out = "output/panache/dynamic")


# =============================================================================
# ### X11 analyses — static threshold (supplementary)
# =============================================================================

Apply_X11_method_on_time_series(sextant_spm_all,
                                Zones = zones_list,
                                plume_time_step = "WEEKLY",
                                plume_dir_in = "output/panache/static",
                                X11_dir_out = "output/panache/static",
                                include_river_flow = True)

Apply_X11_method_on_time_series_per_river(sextant_spm_all,
                                          plume_time_step = "WEEKLY",
                                          plume_dir_in = "output/panache/static",
                                          X11_dir_out = "output/panache/static")


# =============================================================================
# ### Multi-driver interaction analysis (GLM / GAM / RF, both thresholds)
# =============================================================================

# Source the R script
driver_interactions_R_path = os.path.join(func_dir, 'driver_interactions.R')
robjects.r['source'](driver_interactions_R_path)

r_function = robjects.r['run_driver_interactions_analysis']

# Call the R function (runs both the dynamic-threshold main analysis and the
# static-threshold supplementary analysis; see func/driver_interactions.R)
r_function()


# =============================================================================
# ### Monthly multi-driver interaction analysis (sec:seasonal_methods)
# =============================================================================

# Re-runs the same six-step GLM/GAM/RF sequence above independently within
# each calendar month's data subset, dynamic threshold only. Feeds the
# Supplementary monthly driver-dominance table (manuscript.tex); see
# func/driver_interactions.R::run_monthly_driver_interactions_analysis().
r_function_monthly = robjects.r['run_monthly_driver_interactions_analysis']
r_function_monthly()


# =============================================================================
# ### Plume shape (compactness), both thresholds
# =============================================================================

# Derives PlumeShape.csv per zone/threshold from panache's PlumeMasks.nc (see
# func/compute_plume_shape.py); read by func/figure.R's compactness panels
# and func/compute_shape_alongcoast_trend.R.
import compute_plume_shape  # noqa: F401


# =============================================================================
# ### Manuscript stats scripts
# =============================================================================

# These are one-off scripts that generate stats used in the manuscript; 
# wiring them in here keeps their outputs (output/STATS/*.csv) in sync 
# with the panache/X11/ driver-interactions data above. 
# Order matters: generate_monthly_trend_pct_heatmap.R and 
# generate_table_s_monthly_trends.R read compute_seasonal_trend.R's output; 
# generate_table_s_octant_trends.R reads compute_direction_octant_trend.R's output.
stats_scripts = [
    'compute_area_trend.R',
    'compute_mass_spm_trend.R',
    'compute_shape_alongcoast_trend.R',
    'compute_driver_correlation_trend.R',
    'compute_seasonal_trend.R',
    'generate_monthly_trend_pct_heatmap.R',
    'generate_table_s_monthly_trends.R',
    'compute_direction_octant_trend.R',
    'generate_table_s_octant_trends.R',
    'compute_driver_correlation_matrices.R',
]
for script in stats_scripts:
    robjects.r['source'](os.path.join(func_dir, script))

