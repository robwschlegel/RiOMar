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

import util, dl
from dl import Download_satellite_data, Plot_and_Save_the_map

# Set matplotlib backend to prevent plots from displaying
# mpl.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM!)
mpl.use('agg') # Prevent showing plot in the Plot panel (this saves RAM)


# =============================================================================
#### Download Chl a data
# =============================================================================

# Download all CHl A data from 1998 to 2025
# NB: This takes a few hours and usesd ~280 GB of disk spaces
sextant_chla_all = {'Data_sources':['SEXTANT'],
                    'Sensor_names':["merged"],
                    'Satellite_variables':['CHLA'],
                    'Atmospheric_corrections':['polymer'],
                    'Temporal_resolution':['DAILY'],
                    'start_day':'1998/01/01',
                    'end_day':'2025/12/31'}
Download_satellite_data(sextant_chla_all,
                        nb_of_cores_to_use = 14,
                        overwrite_existing_satellite_files = False,
                        where_to_save_satellite_data = 'data')


# =============================================================================
#### Download wind data
# =============================================================================

# First test at downloading one year in a single go
sextant_wind_1999 = {'Data_sources':['SEXTANT'],
                     'Sensor_names':["merged"],
                     'Satellite_variables':['WIND'],
                     'Atmospheric_corrections':['polymer'],
                     'Temporal_resolution':['DAILY'],
                     'start_day':'1999/01/01',
                     'end_day':'1999/01/31'} 
Download_satellite_data(sextant_wind_1999,
                        nb_of_cores_to_use = 14,
                        overwrite_existing_satellite_files = False,
                        where_to_save_satellite_data = 'data')


# =============================================================================
#### Create and save maps of the downloaded data
# =============================================================================

# Plot everything in one go (this takes a while)
for year in range(1998, 2026):
    Plot_and_Save_the_map(
        sextant_chla_all,
        nb_of_cores_to_use = 14,
        where_are_saved_satellite_data = 'data',
        start_day_of_maps_to_plot = f'{year}/01/01',
        end_day_of_maps_to_plot = f'{year}/12/31'
    )
    print(f'Year {year} done!')

