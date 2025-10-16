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

# The time steps to investigate
time_steps = ['DAILY', 'WEEKLY', 'MONTHLY', 'ANNUAL']

# Basic arguments to be used throughout the script
sextant_spm_all = {'Data_sources':['SEXTANT'],
                   'Sensor_names':["merged"],
                   'Satellite_variables':['SPM'],
                   'Atmospheric_corrections':['Standard'],
                   'Temporal_resolution':['DAILY'],
                   'start_day':'1998/01/01',
                   'end_day':'2025/12/31'}


# =============================================================================
# ### Detect plumes
# =============================================================================

# Basic plume detection for all zones and daily results
for zone in zones_list:
    apply_plume_mask(sextant_spm_all,
                     Zones = [zone],
                     time_step = 'DAILY',
                     nb_cores = 14,
                     dynamic_thresh = False,
                     regional_map_dir = "output/REGIONAL_MAPS",
                     plume_dir = "output/FIXED_THRESHOLD")

# For weekly results
for zone in zones_list:
    apply_plume_mask(sextant_spm_all,
                     Zones = [zone],
                     time_step = 'WEEKLY',
                     nb_cores = 14,
                     dynamic_thresh = False,
                     regional_map_dir = "output/REGIONAL_MAPS",
                     plume_dir = "output/FIXED_THRESHOLD")


# For weekly results with dynamic threshold
for zone in zones_list:
    apply_plume_mask(sextant_spm_all,
                     Zones = [zone],
                     time_step = 'WEEKLY',
                     nb_cores = 14,
                     dynamic_thresh = True,
                     regional_map_dir = "output/REGIONAL_MAPS",
                     plume_dir = "output/DYNAMIC_THRESHOLD")
    

# =============================================================================
# ### Create time series of plume surface
# =============================================================================

# All in one go for daily results
# NB: X11 will not run on the daily results - need to use STL decomposition instead
for zone in zones_list:
    make_and_plot_time_series_of_plume_areas(sextant_spm_all,
                                             Zones = [zone],
                                             nb_cores = 14,
                                             time_step = 'DAILY',
                                             plume_dir_in = "output/FIXED_THRESHOLD",
                                             plume_dir_out = "output/FIXED_THRESHOLD")

# Rather use weekly results for the plots
for zone in zones_list:
    make_and_plot_time_series_of_plume_areas(sextant_spm_all,
                                             Zones = [zone],
                                             nb_cores = 14,
                                             time_step = 'WEEKLY',
                                             plume_dir_in = "output/FIXED_THRESHOLD",
                                             plume_dir_out = "output/FIXED_THRESHOLD")

# Weekly results with dynamic threshold
for zone in zones_list:
    make_and_plot_time_series_of_plume_areas(sextant_spm_all,
                                             Zones = [zone],
                                             nb_cores = 14,
                                             time_step = 'WEEKLY',
                                             plume_dir_in = "output/DYNAMIC_THRESHOLD",
                                             plume_dir_out = "output/DYNAMIC_THRESHOLD")
    