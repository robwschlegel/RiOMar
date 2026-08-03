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

from X11 import Apply_X11_method_on_time_series

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


# =============================================================================
# ### X11 analyses — static threshold (supplementary)
# =============================================================================

Apply_X11_method_on_time_series(sextant_spm_all,
                                Zones = zones_list,
                                plume_time_step = "WEEKLY",
                                plume_dir_in = "output/panache/static",
                                X11_dir_out = "output/panache/static",
                                include_river_flow = True)


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

