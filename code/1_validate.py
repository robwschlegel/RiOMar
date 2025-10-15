#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to validate the full dataset used in the RiOMar project.
# It can be designed to be called by the Makefile.


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
mpl.use('agg')


# =============================================================================
# ### Satellite - in situ Match up for Chl a and SPM data
# =============================================================================

# Match up in situ and satellite Chl a data
sextant_chla_all = {'Data_sources':['SEXTANT'],
                    'Sensor_names':["merged"],
                    'Satellite_variables':['CHLA'],
                    'Atmospheric_corrections':['Standard'],
                    'Temporal_resolution':['DAILY'],
                    'start_day':'1998/01/01',
                    'end_day':'2025/12/31'}
Match_up_with_insitu_measurements(sextant_chla_all,
                                  zones = ['FRANCE'],
                                  # zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                  redo_the_MU_database = False, # Change to True is running for the first time 
                                  nb_of_cores_to_use = 14,
                                  where_are_saved_satellite_data = 'data',
                                  where_to_save_Match_Up_data = 'output')

# Match up in situ and satellite SPM data
sextant_spm_all = {'Data_sources':['SEXTANT'],
                    'Sensor_names':["merged"],
                    'Satellite_variables':['SPM'],
                    'Atmospheric_corrections':['Standard'],
                    'Temporal_resolution':['DAILY'],
                    'start_day':'1998/01/01',
                    'end_day':'2025/12/31'}
Match_up_with_insitu_measurements(sextant_spm_all,
                                  zones = ['FRANCE'],
                                  # zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                  redo_the_MU_database = True, # Change to True is running for the first time 
                                  nb_of_cores_to_use = 14,
                                  where_are_saved_satellite_data = 'data',
                                  where_to_save_Match_Up_data = 'output')


# =============================================================================
# ### Satellite - in situ Match up for drivers
# =============================================================================

# I'm imagining to put here some sort of initial comparison of the driver data against some
# measurement of the rivers and their plumes. But this may be better to leave until 
# step 4) time series analysis.

