#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to create the regional maps used in the RiOMar project. 
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, regmap
from regmap import create_regional_maps, QC_of_regional_maps

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')

# The zones for mapping
zones_list = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY']


# =============================================================================
# ### Make regional maps
# =============================================================================

# All years of data
sextant_chla_all = {'Data_sources':['SEXTANT'],
                    'Sensor_names':["merged"],
                    'Satellite_variables':['CHLA'],
                    'Atmospheric_corrections':['polymer'],
                    'Temporal_resolution':['DAILY'],
                    'start_day':'1998/01/01',
                    'end_day':'2025/12/31'}
for zone in zones_list:
    create_regional_maps(sextant_chla_all,
                         Zones = [zone],
                         overwrite_existing_regional_maps = False,
                         save_map_plots_of_which_time_frequency = {'DAILY' : False, 'WEEKLY' : False, 'MONTHLY' : True, 'ANNUAL' : True},
                         nb_of_cores_to_use = 14,
                         where_are_saved_satellite_data = "data",
                         where_to_save_regional_maps = "output/REGIONAL_MAPS")


# =============================================================================
# ### QC of created maps
# =============================================================================

# All years and zones at once
QC_of_regional_maps(sextant_chla_all,
                    # Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY'],
                    Zones = zones_list,
                    nb_of_cores_to_use = 14,
                    where_are_saved_regional_maps = "output/REGIONAL_MAPS")

