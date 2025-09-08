#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
# ### Import libraries and functions
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

# Check directories
# for p in sys.path :
#     print( p )

import util, dl, validate, regmap, plume, X11, figure
from dl import Download_satellite_data, Plot_and_Save_the_map
from validate import Match_up_with_insitu_measurements
from regmap import create_regional_maps, QC_of_regional_maps
from plume import apply_plume_mask, make_and_plot_time_series_of_plume_areas
from X11 import Apply_X11_method_on_time_series
from figure import Figure_1, Figure_2, Figure_4, Figure_5, Figure_6_7, Figure_8_9_10


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
# ### Test 1: data downloading
# =============================================================================

# Download some target data
Download_satellite_data(core_arguments,
                        nb_of_cores_to_use = 1,
                        overwrite_existing_satellite_files = False,
                        where_to_save_satellite_data = 'data')

# Create and save maps of the downloaded data
Plot_and_Save_the_map(core_arguments,
                      nb_of_cores_to_use = 6,
                      where_are_saved_satellite_data = 'data',
                      start_day_of_maps_to_plot = '2018/08/01',
                      end_day_of_maps_to_plot = '2018/08/06')


# =============================================================================
# ### Test 2: data validation
# =============================================================================

# Match up in situ and satellite data
Match_up_with_insitu_measurements(core_arguments,
                                  # zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                  zones = ['FRANCE'],
                                  redo_the_MU_database = False, 
                                  nb_of_cores_to_use = 6,
                                  where_are_saved_satellite_data = "data",
                                  where_to_save_Match_Up_data = "output")


# =============================================================================
# ### Test 3: regional Maps
# =============================================================================

# Make regional map
create_regional_maps(core_arguments,
                     Zones = ['SOUTHERN_BRITTANY'],
                     overwrite_existing_regional_maps = True,
                     save_map_plots_of_which_time_frequency = {'DAILY' : False, 'WEEKLY' : False, 'MONTHLY' : True, 'ANNUAL' : True},
                     nb_of_cores_to_use = 4,
                     where_are_saved_satellite_data = "data",
                     where_to_save_regional_maps = "figures")

# QC of these regional maps
QC_of_regional_maps(core_arguments,
                    Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY'],
                    nb_of_cores_to_use = 4,
                    where_are_saved_regional_maps = "output")


# =============================================================================
# ### Test 4: detection of plumes
# =============================================================================

# Basic plume detection
apply_plume_mask(core_arguments,
                 # Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                 Zones = ['SOUTHERN_BRITTANY'],
                 detect_plumes_on_which_temporal_resolution_data = 'WEEKLY',
                 nb_of_cores_to_use = 4,
                 use_dynamic_threshold = True,
                 where_are_saved_regional_maps = "figures",
                 where_to_save_plume_results = "output/DYNAMIC_THRESHOLD")


# Plot time-series of plume surface
make_and_plot_time_series_of_plume_areas(core_arguments,
                                         Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                         # Zones = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                         nb_of_cores_to_use = 4,
                                         on_which_temporal_resolution_the_plumes_have_been_detected = 'WEEKLY',
                                         where_are_saved_plume_results = "output/DYNAMIC_THRESHOLD",
                                         where_to_save_plume_time_series = "output/DYNAMIC_THRESHOLD")


# =============================================================================
# ### Test 5: X11 analyses
# =============================================================================

# Do X11 analysis on plume time-series
Apply_X11_method_on_time_series(core_arguments,
                                Zones = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                # Zones = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY'],
                                nb_of_cores_to_use = 4,
                                on_which_temporal_resolution_the_plumes_have_been_detected = 'WEEKLY',
                                where_are_saved_plume_time_series = "output/DYNAMIC_THRESHOLD",
                                where_to_save_X11_results = "output/DYNAMIC_THRESHOLD",
                                include_river_flow = True)


# =============================================================================
# ### Test 6: figures for the article
# =============================================================================


Figure_1(where_are_saved_satellite_data = "data/OCEAN_COLOR",
         where_to_save_the_figure = "figures")

Figure_2(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_4(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_5(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold = "output/DYNAMIC_THRESHOLD",
           where_are_saved_plume_results_with_fixed_threshold = "output/FIXED_THRESHOLD",
           where_to_save_the_figure = "figures")

Figure_8_9_10(where_are_saved_X11_results = "output/DYNAMIC_THRESHOLD",
              where_to_save_the_figure = "figures")
