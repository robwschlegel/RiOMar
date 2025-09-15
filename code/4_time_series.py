#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to run the time series analyses used in the RiOMar project. 
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, X11
from X11 import Apply_X11_method_on_time_series

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')

# The zones for mapping
zones_list = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY']

# Basic arguments to be used throughout the script
core_arguments = {'Data_sources':['SEXTANT'],
                  'Sensor_names':["merged"],
                  'Satellite_variables':['SPM'],
                  'Atmospheric_corrections':['polymer'],
                  'Temporal_resolution':['DAILY'],
                  'start_day':'1998/01/01',
                  'end_day':'2025/12/31'}


# =============================================================================
# ### X11 analyses
# =============================================================================

# Full analysis in one go
Apply_X11_method_on_time_series(core_arguments,
                                Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                # Zones = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                nb_of_cores_to_use = 4,
                                on_which_temporal_resolution_the_plumes_have_been_detected = 'DAILY',
                                where_are_saved_plume_time_series = "output/FIXED_THRESHOLD",
                                where_to_save_X11_results = "output/FIXED_THRESHOLD",
                                include_river_flow = True)

