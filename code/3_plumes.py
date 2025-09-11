#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to analyse the river plumes used in the RiOMar project. 
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, plume
from plume import apply_plume_mask, make_and_plot_time_series_of_plume_areas

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
# ### Detect plumes
# =============================================================================

# Basic plume detection
for zone in zones_list:
    apply_plume_mask(core_arguments,
                     Zones = [zone],
                     detect_plumes_on_which_temporal_resolution_data = 'DAILY',
                     nb_of_cores_to_use = 14,
                     use_dynamic_threshold = False,
                     where_are_saved_regional_maps = "output/REGIONAL_MAPS",
                     where_to_save_plume_results = "output/FIXED_THRESHOLD")
    
# TODO: Dynamic thresholding currently fails for some zones (e.g. BAY_OF_SEINE)


# =============================================================================
# ### Create time series of plume surface
# =============================================================================

# All in one go
for zone in zones_list:
    make_and_plot_time_series_of_plume_areas(core_arguments,
                                             Zones = [zone],
                                        #  Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                         # Zones = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                         nb_of_cores_to_use = 14,
                                         on_which_temporal_resolution_the_plumes_have_been_detected = 'WEEKLY',
                                         where_are_saved_plume_results = "output/FIXED_THRESHOLD",
                                         where_to_save_plume_time_series = "output/FIXED_THRESHOLD")
    
# TODO: Add dynamic thresholding once it works properly

