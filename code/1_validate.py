#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script is meant to house all of the code needed to download the full dataset
# used in the RiOMar project. It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, validate
from validate import Match_up_with_insitu_measurements

# Set matplotlib backend to prevent plots from displaying
# mpl.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM!)
mpl.use('agg') # Prevent showing plot in the Plot panel (this saves RAM)


# =============================================================================
# ### Data validation
# =============================================================================

# Match up in situ and satellite data
sextant_chla_all = {'Data_sources':['SEXTANT'],
                    'Sensor_names':["merged"],
                    'Satellite_variables':['CHLA'],
                    'Atmospheric_corrections':['polymer'],
                    'Temporal_resolution':['DAILY'],
                    'start_day':'1998/01/01',
                    'end_day':'2025/12/31'}
Match_up_with_insitu_measurements(sextant_chla_all,
                                  zones = ['FRANCE'],
                                  # zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                  redo_the_MU_database = True, 
                                  nb_of_cores_to_use = 14,
                                  where_are_saved_satellite_data = 'data',
                                  where_to_save_Match_Up_data = 'output')

