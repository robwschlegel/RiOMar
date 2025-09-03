#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
# ### Import libraries and functions
# =============================================================================

# NB: This method does not appear to add the functions etc. to the environment
# import runpy
# runpy.run_path(path_name='../code/0_data_management.py')

exec(open("code/0_data_management.py").read())
exec(open("code/1_data_validation.py").read())
exec(open("code/2_regional_maps.py").read())


# =============================================================================
# ### Prep backend and core arguments
# =============================================================================

# Set matplotlib backend to prevent plots from displaying
# mpl.use('module://matplotlib_inline.backend_inline') # To show plots on the Plot panel (be careful as it consumes RAM!)
mpl.use('agg') # Prevent showing plot in the Plot panel (this saves RAM)

core_arguments = {'Data_sources':['SEXTANT'],
                  'Sensor_names':["merged"],
                  'Satellite_variables':['CHLA'],
                  'Atmospheric_corrections':['polymer'],
                  'Temporal_resolution':['DAILY'],
                  'start_day':'2018/08/01',
                  'end_day':'2018/08/31'}


# =============================================================================
# ### Import libraries and functions
# =============================================================================

## 0 Data management

# Test 1
## Download some target data
Download_satellite_data(core_arguments,
                        nb_of_cores_to_use = 1,
                        overwrite_existing_satellite_files = False,
                        where_to_save_satellite_data = 'data')

# Test 2
## Create and save maps of the downloaded data
Plot_and_Save_the_map(core_arguments,
                      nb_of_cores_to_use = 6,
                      where_are_saved_satellite_data = 'data',
                      start_day_of_maps_to_plot = '2018/08/01',
                      end_day_of_maps_to_plot = '2018/08/06')

## 1 Data validation

# Test 3
## Match up in situ and satellite data
Match_up_with_insitu_measurements(core_arguments, 
                                  # zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                  zones = ['FRANCE'],
                                  redo_the_MU_database = False, 
                                  nb_of_cores_to_use = 6,
                                  where_are_saved_satellite_data = "data",
                                  where_to_save_Match_Up_data = "output/MATCH_UP_DATA")

## 2 Regional Maps

# Test 4
## Make regional map
create_regional_maps(core_arguments,
                     Zones = ['SOUTHERN_BRITTANY'],
                     overwrite_existing_regional_maps = True,
                     save_map_plots_of_which_time_frequency = {'DAILY' : False, 'WEEKLY' : False, 'MONTHLY' : True, 'ANNUAL' : True},
                     nb_of_cores_to_use = 4,
                     where_are_saved_satellite_data = "data",
                     where_to_save_regional_maps = "figures")
