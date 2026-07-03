#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to analyse the river plumes used in the RiOMar project. 
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os
import sys
import subprocess
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

# Run panache module directly via terminal
# NB: Takes roughly 20 minutes per zone
# TODO: At the moment this does not run via Python. Needs to be called directly within a terminal.
subprocess.run("panache metadata/zone_config_GULF_OF_LION.json", shell=True)
subprocess.run("panache metadata/zone_config_BAY_OF_BISCAY.json", shell=True)
subprocess.run("panache metadata/zone_config_SOUTHERN_BRITTANY.json", shell=True)
subprocess.run("panache metadata/zone_config_BAY_OF_SEINE.json", shell=True)

# Plume detection for all zones and daily results
# for zone in zones_list:
#     apply_plume_mask(sextant_spm_all,
#                      Zones = [zone],
#                      time_step = 'DAILY',
#                      nb_cores = 14,
#                      dynamic_thresh = False,
#                      regional_map_dir = "output/REGIONAL_MAPS",
#                      plume_dir = "output/FIXED_THRESHOLD")
    
# Daily plumes with a dynamic threshold
# for zone in zones_list:
#     apply_plume_mask(sextant_spm_all,
#                      Zones = [zone],
#                      time_step = 'DAILY',
#                      nb_cores = 14,
#                      dynamic_thresh = True,
#                      regional_map_dir = "output/REGIONAL_MAPS",
#                      plume_dir = "output/DYNAMIC_THRESHOLD")

# For weekly results
# for zone in zones_list:
#     apply_plume_mask(sextant_spm_all,
#                      Zones = [zone],
#                      time_step = 'WEEKLY',
#                      nb_cores = 14,
#                      dynamic_thresh = False,
#                      regional_map_dir = "output/REGIONAL_MAPS",
#                      plume_dir = "output/FIXED_THRESHOLD")


# For weekly results with dynamic threshold
# for zone in zones_list:
#     apply_plume_mask(sextant_spm_all,
#                      Zones = [zone],
#                      time_step = 'WEEKLY',
#                      nb_cores = 14,
#                      dynamic_thresh = True,
#                      regional_map_dir = "output/REGIONAL_MAPS",
#                      plume_dir = "output/DYNAMIC_THRESHOLD")
    

# =============================================================================
# ### Create time series of plume surface
# =============================================================================

# All in one go for daily results
# NB: X11 will not run on the daily results - need to use STL decomposition instead
# for zone in zones_list:
#     make_and_plot_time_series_of_plume_areas(sextant_spm_all,
#                                              Zones = [zone],
#                                              nb_cores = 14,
#                                              time_step = 'DAILY',
#                                              plume_dir_in = "output/FIXED_THRESHOLD",
#                                              plume_dir_out = "output/FIXED_THRESHOLD")
    
# Daily dynamic threshold time series
# for zone in zones_list:
#     make_and_plot_time_series_of_plume_areas(sextant_spm_all,
#                                              Zones = [zone],
#                                              nb_cores = 14,
#                                              time_step = 'DAILY',
#                                              plume_dir_in = "output/DYNAMIC_THRESHOLD",
#                                              plume_dir_out = "output/DYNAMIC_THRESHOLD")


# Rather use weekly results for the plots
# for zone in zones_list:
#     make_and_plot_time_series_of_plume_areas(sextant_spm_all,
#                                              Zones = [zone],
#                                              nb_cores = 14,
#                                              time_step = 'WEEKLY',
#                                              plume_dir_in = "output/FIXED_THRESHOLD",
#                                              plume_dir_out = "output/FIXED_THRESHOLD")

# Weekly results with dynamic threshold
# for zone in zones_list:
#     make_and_plot_time_series_of_plume_areas(sextant_spm_all,
#                                              Zones = [zone],
#                                              nb_cores = 14,
#                                              time_step = 'WEEKLY',
#                                              plume_dir_in = "output/DYNAMIC_THRESHOLD",
#                                              plume_dir_out = "output/DYNAMIC_THRESHOLD")
    