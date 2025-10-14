#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to download the full dataset used in the RiOMar project
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, dl
from util import daily_integral
from dl import Download_satellite_data, Plot_and_Save_the_map, download_cmems_subset

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg') # Prevent showing plot in the Plot panel (this saves RAM)
# Or, to show plots on the Plot panel (be careful as it consumes RAM!)
# mpl.use('module://matplotlib_inline.backend_inline')

# The zones for downloading
zones_list = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY']


# =============================================================================
#### Download Chl a data
# =============================================================================

# Download all Chl a data from 1998 to 2025
# NB: This takes a few hours and usesd ~280 GB of disk spaces
# NB: Several days of data are missing
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
#### Download SPM data
# =============================================================================

# Download all SPM data from 1998 to 2025
# NB: This takes a few hours and usesd ~280 GB of disk spaces
# NB: Several days of data are missing
sextant_spm_all = {'Data_sources':['SEXTANT'],
                    'Sensor_names':["merged"],
                    'Satellite_variables':['SPM'],
                    'Atmospheric_corrections':['polymer'],
                    'Temporal_resolution':['DAILY'],
                    'start_day':'1998/01/01',
                    'end_day':'2025/12/31'}
Download_satellite_data(sextant_spm_all,
                        nb_of_cores_to_use = 14,
                        overwrite_existing_satellite_files = False,
                        where_to_save_satellite_data = 'data')


# =============================================================================
#### Download wind data
# =============================================================================

# The historic WIND data (1998-2007) at 0.25° resolution
for zone in zones_list:
    download_cmems_subset(
        zone,
        'cmems_obs-wind_glo_phy_my_l4_0.25deg_PT1H',
        ['eastward_wind', 'northward_wind'],
        '1998-01-01T00:00:00', '2008-01-01T00:00:00', # Remember to get one extra hour to cover the integral up to 2007-12-31 23:00
        f'~/pCloudDrive/data/WIND/{zone}'
    )
# The recent WIND data (2008-2024) at 0.125° resolution
for zone in zones_list:
    download_cmems_subset(
        zone,
        'cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H',
        ['eastward_wind', 'northward_wind'],
        '2008-01-01T00:00:00', '2025-01-01T00:00:00',
        f'~/pCloudDrive/data/WIND/{zone}'
    )
# The near-real-time data (2024-2025 so far)
for zone in zones_list:
    download_cmems_subset(
        zone,
        'cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H',
        ['eastward_wind', 'northward_wind'],
        '2025-01-01T00:00:00', '2025-09-01T00:00:00',
        f'~/pCloudDrive/data/WIND/{zone}'
    )


# =============================================================================
#### Download other surface variables
# =============================================================================

# The historic GLORYS data (1998-2021)
for zone in zones_list:
    download_cmems_subset(
        zone,
        'cmems_mod_glo_phy_my_0.083deg_P1D-m',
        ['uo', 'vo', 'zos', 'thetao', 'bottomT', 'mlotst', 'so'],
        '1998-01-01T00:00:00', '2021-06-30T00:00:00',
        f'~/pCloudDrive/data/GLORYS/{zone}'
    )
# The recent GLORYS data (2021-2025 so far)
for zone in zones_list:
    download_cmems_subset(
        zone,
        'cmems_mod_glo_phy_myint_0.083deg_P1D-m',
        ['uo', 'vo', 'zos', 'thetao', 'bottomT', 'mlotst', 'so'],
        '2021-07-01T00:00:00', '2025-07-31T00:00:00',
        f'~/pCloudDrive/data/GLORYS/{zone}'
    )

# =============================================================================
#### Create daily integrals of hourly wind data
# =============================================================================

for zone in zones_list:
    daily_integral(f'/home/calanus/pCloudDrive/data/WIND/{zone}',
                   overwrite = False) # Change to True if running for the first time
    

# =============================================================================
#### Create maps of the Chla and SPM data
# =============================================================================

# Plot everything in one go (this takes a while)
for year in range(1998, 2025):
    Plot_and_Save_the_map(
        sextant_chla_all,
        nb_of_cores_to_use = 14,
        where_are_saved_satellite_data = 'data',
        start_day_of_maps_to_plot = f'{year}/01/01',
        end_day_of_maps_to_plot = f'{year}/12/31'
    )
    print(f'Year {year} done!')

# Plot everything in one go (this takes a while)
for year in range(1998, 2025):
    Plot_and_Save_the_map(
        sextant_spm_all,
        nb_of_cores_to_use = 14,
        where_are_saved_satellite_data = 'data',
        start_day_of_maps_to_plot = f'{year}/01/01',
        end_day_of_maps_to_plot = f'{year}/12/31'
    )
    print(f'Year {year} done!')

