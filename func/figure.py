#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import os, sys, pickle, re, glob
import numpy as np
import xarray as xr
import pandas as pd
import rpy2.robjects as robjects
from functools import reduce

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (load_csv_files_in_the_package_folder, path_to_fill_to_where_to_save_satellite_files,
                  align_bathymetry_to_resolution, define_parameters)
from validate import get_insitu_measurements
from plume import (reduce_resolution,
                   preprocess_annual_dataset_and_compute_land_mask, create_polygon_mask,
                   Create_the_plume_mask, Pipeline_to_delineate_the_plume)


# =============================================================================
#### Utility functions
# =============================================================================


def save_files_for_Figure_1(where_are_saved_satellite_data, where_to_save_the_figure, 
                            date_of_the_map, coordinates_of_the_map) :

    folder_where_to_save_Figure_1_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURES', 'FIGURE_1', 'DATA')
    os.makedirs(folder_where_to_save_Figure_1_data, exist_ok = True)
    
    path_to_nc_file = (path_to_fill_to_where_to_save_satellite_files(where_are_saved_satellite_data)
                       .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                                'SEXTANT/SPM/merged/Standard/DAILY')
                       .replace('[YEAR]/[MONTH]/[DAY]',
                                date_of_the_map))
        
    with xr.open_dataset( glob.glob(path_to_nc_file + "/*.nc")[0] ) as ds :
        SPM_map = (ds['analysed_spim']
                   .sel(lat=slice(coordinates_of_the_map['lat_min'], coordinates_of_the_map['lat_max']), 
                        lon=slice(coordinates_of_the_map['lon_min'], coordinates_of_the_map['lon_max'])) 
                   .to_dataframe()
                   .reset_index()
                   .drop(columns=["time"]))
        
        SPM_map.to_csv(folder_where_to_save_Figure_1_data + "/SPM_map.csv")
          
    extract_insitu_stations_and_save_the_file_for_plot(folder_where_to_save_Figure_1_data)
    
    
def load_the_regional_maps_and_save_them_for_plotting(where_are_saved_regional_maps, where_to_save_the_figure, dates_for_each_zone) :
        
    folder_where_to_save_Figure_2_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURES', 'FIGURE_2', 'DATA')
    os.makedirs(folder_where_to_save_Figure_2_data, exist_ok = True)
    
    path_to_regional_maps = {key : (path_to_fill_to_where_to_save_satellite_files( os.path.join(where_are_saved_regional_maps, 'REGIONAL_MAPS', key) )
                                       .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                                                'SEXTANT/SPM/merged/Standard/MAPS/DAILY')
                                       .replace('[YEAR]/[MONTH]/[DAY]', f'{date[:4]}/{date}.pkl')) 
                               for key, date in dates_for_each_zone.items()}
    
    for key, path_to_map in path_to_regional_maps.items() : 
    
        coordinates_of_the_map = define_parameters(key)    
    
        with open(path_to_map, 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data']  
            
            SPM_map = (ds
                       .sel(lat=slice(coordinates_of_the_map['lat_range_of_the_map_to_plot'][0], 
                                      coordinates_of_the_map['lat_range_of_the_map_to_plot'][1]), 
                            lon=slice(coordinates_of_the_map['lon_range_of_the_map_to_plot'][0],
                                      coordinates_of_the_map['lon_range_of_the_map_to_plot'][1])) 
                       .to_dataframe()
                       .reset_index())
            
            SPM_map.to_csv(folder_where_to_save_Figure_2_data + f"/{key}.csv")
            
    extract_insitu_stations_and_save_the_file_for_plot(folder_where_to_save_Figure_2_data)


def extract_insitu_stations_and_save_the_file_for_plot(folder_where_to_save_Figure_data) :     

    coordinates_of_the_RIOMARS = { zone_name : {'lat_min' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][0],
                                                'lat_max' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][1],
                                                'lon_min' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][0],
                                                'lon_max' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][1]}
                                  for zone_name in ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY', 'BAY_OF_SEINE'] }
    
    pd.DataFrame.from_dict(coordinates_of_the_RIOMARS).to_csv(folder_where_to_save_Figure_data + "/RIOMAR_limits.csv")
    
    _, insitu_stations = get_insitu_measurements()
    station_LATITUDES = insitu_stations.LATITUDE.to_numpy(dtype=float)
    station_LONGITUDES = insitu_stations.LONGITUDE.to_numpy(dtype=float)
    
    index_of_stations_in_the_RIOMARS = [ np.where((station_LATITUDES >= coords['lat_min']) & 
                                                  (station_LATITUDES <= coords['lat_max']) & 
                                                  (station_LONGITUDES >= coords['lon_min']) & 
                                                  (station_LONGITUDES <= coords['lon_max']) )[0] 
                                        for coords in coordinates_of_the_RIOMARS.values() ]
    
    index_of_stations_in_the_RIOMARS = np.unique(np.concatenate(index_of_stations_in_the_RIOMARS))
    
    stations_in_the_RIOMARS = insitu_stations.iloc[index_of_stations_in_the_RIOMARS]
    
    stations_in_the_RIOMARS.to_csv(folder_where_to_save_Figure_data + "/Stations_position.csv")
    
    
def dates_for_each_zone() : 
    return {'GULF_OF_LION' : '2014-01-04',
            'BAY_OF_BISCAY' : '2009-04-22',
            'SOUTHERN_BRITTANY' : '2016-05-23',# '2022-01-21',
            'BAY_OF_SEINE' : '2018-02-25'}


# =============================================================================
#### Main functions
# =============================================================================


def Figure_1(where_are_saved_satellite_data, where_to_save_the_figure):
    save_files_for_Figure_1(where_are_saved_satellite_data,
                            where_to_save_the_figure,
                            date_of_the_map="2011/02/02",
                            coordinates_of_the_map={"lat_min": 42, "lat_max": 51.5, "lon_min": -6, "lon_max": 8})

    # Source the R scrip
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_1']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def Figure_2(where_are_saved_regional_maps, where_to_save_the_figure, include_station_points=True):
    # dates_for_each_zone = {'GULF_OF_LION' : '2014-01-04',
    #                        'BAY_OF_BISCAY' : '2009-04-22',
    #                        'SOUTHERN_BRITTANY' : '2022-01-21',
    #                        'BAY_OF_SEINE' : '2018-02-25'}

    the_dates_for_each_zone = dates_for_each_zone()

    load_the_regional_maps_and_save_them_for_plotting(where_are_saved_regional_maps,
                                                      where_to_save_the_figure,
                                                      the_dates_for_each_zone)

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_2']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               include_station_points=robjects.BoolVector([include_station_points]))


def Figure_4(where_are_saved_regional_maps, where_to_save_the_figure):
    
    # Static date for each zone to illustrate the plume detection steps
    Zone = 'BAY_OF_SEINE'
    plume_name = 'Seine'
    Date = '2018-02-25'

    parameters = define_parameters(Zone)

    path_to_the_satellite_file_to_use = os.path.join(where_are_saved_regional_maps, 'REGIONAL_MAPS', Zone, 'SEXTANT', 'SPM',
                                                     'merged',
                                                     'Standard', 'MAPS', 'DAILY', Date[:4], f'{Date}.pkl')

    # Open and load the file (binary file assumed to contain data)
    with open(path_to_the_satellite_file_to_use, 'rb') as f:
        ds = pickle.load(f)['Basin_map']['map_data']

    # Reduce the resolution of the dataset to the specified latitude and longitude resolutions
    ds_reduced = (reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])
                  if parameters['lat_new_resolution'] is not None
                  else ds)

    bathymetry_data_aligned_to_reduced_map = align_bathymetry_to_resolution(ds_reduced,
                                                                            f'{where_are_saved_regional_maps}/REGIONAL_MAPS/{Zone}/Bathy_data.pkl')

    (_, land_mask) = preprocess_annual_dataset_and_compute_land_mask(
        (path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + "/REGIONAL_MAPS/" + Zone)
         .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                  'SEXTANT/SPM/merged/Standard/MAPS/MULTIYEAR')
         .replace('[YEAR]/[MONTH]/[DAY]', 'Averaged_over_multi-years.pkl')
         ), parameters)

    inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)

    where_to_save_the_figure_4 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_4")

    the_plume = Create_the_plume_mask(ds_reduced,
                                      bathymetry_data_aligned_to_reduced_map,
                                      land_mask,
                                      parameters,
                                      plume_name)
    the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
                        name_of_the_plot='A')

    # TODO: Change this so it reacts to fixed or dynamic. Or just alway use fixed.
    the_plume.determine_SPM_threshold(dynamic_determination_of_SPM_threshold=True)
    the_plume.SPM_threshold = 5.5
    the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
                        name_of_the_plot='B')

    the_plume.do_a_raw_plume_detection()
    the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
                        name_of_the_plot='C')

    the_plume.include_cloudy_regions_to_plume_area()

    the_plume.remove_the_areas_with_sediment_resuspension(
        maximal_bathymetry=parameters['maximal_bathymetric_for_zone_with_resuspension'][plume_name],
        minimal_distance_from_estuary=parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][
            plume_name])
    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
    #                    name_of_the_plot='before_shallow_water_removal')
    ##

    the_plume.remove_shallow_waters()
    the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
                        name_of_the_plot='D')

    the_plume.remove_close_river_mouth(the_plume.parameters['pixel_starting_points_close_river_mouth'])

    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
    #                    name_of_the_plot='before_')
    ##

    the_plume.dilate_the_main_plume_area_to_merge_close_plume_areas()

    the_plume.remove_small_shapes_that_do_not_meet_a_minimum_size_criterion()

    the_plume.set_pixels_to_False_if_outside_of_the_searching_area(inside_polygon_mask)

    the_plume.identify_the_main_plume_shape_based_on_the_plume_core_location()

    the_plume.remove_shallow_waters()

    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
    #                    name_of_the_plot='before_shrink_widen')
    ##

    if not np.isin(plume_name, ['Seine']):
        the_plume.remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase()

    the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
                        name_of_the_plot='E')


def Figure_5(where_are_saved_regional_maps, where_to_save_the_figure):
    the_dates_for_each_zone = dates_for_each_zone()

    plume_fixed_thresholds = {
        'Seine': 5.5,
        # 'Gironde' : 5,
        'Grand Rhone': 3,
        'Petit Rhone': 3,
        'Loire': 5,
        'Sevre': 11}
    # 'SOUTHERN_BRITTANY' : '2008-01-26',# '2022-01-21',
    # 'BAY_OF_SEINE' : '2018-02-25'}

    where_to_save_the_figure_5 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_5")
    os.makedirs(where_to_save_the_figure_5, exist_ok=True)

    for Zone, Date in the_dates_for_each_zone.items():

        parameters = define_parameters(Zone)

        path_to_the_satellite_file_to_use = os.path.join(where_are_saved_regional_maps, 'REGIONAL_MAPS', Zone, 'SEXTANT',
                                                         'SPM', 'merged',
                                                         'Standard', 'MAPS', 'DAILY', Date[:4], f'{Date}.pkl')

        # Open and load the file (binary file assumed to contain data)
        with open(path_to_the_satellite_file_to_use, 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data']

            # Reduce the resolution of the dataset to the specified latitude and longitude resolutions
        ds_reduced = (reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])
                      if parameters['lat_new_resolution'] is not None
                      else ds)

        bathymetry_data_aligned_to_reduced_map = align_bathymetry_to_resolution(ds_reduced,
                                                                                f'{where_are_saved_regional_maps}/REGIONAL_MAPS/{Zone}/Bathy_data.pkl')

        (_, land_mask) = preprocess_annual_dataset_and_compute_land_mask(
            (path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + "/REGIONAL_MAPS/" + Zone)
             .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                      'SEXTANT/SPM/merged/Standard/MAPS/MULTIYEAR')
             .replace('[YEAR]/[MONTH]/[DAY]', 'Averaged_over_multi-years.pkl')
             ), parameters)

        inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)

        all_mask_area = []
        all_river_mouth_to_remove = []
        thresholds = {key: None for key in parameters['starting_points']}
        # Loop through each plume starting point to process plume detection
        for plume_name, starting_point in parameters['starting_points'].items():

            if plume_name in plume_fixed_thresholds:
                parameters['fixed_threshold'][plume_name] = plume_fixed_thresholds[plume_name]
                use_dynamic_threshold = False
            else:
                use_dynamic_threshold = True

            the_plume = Pipeline_to_delineate_the_plume(ds_reduced,
                                                        bathymetry_data_aligned_to_reduced_map,
                                                        land_mask,
                                                        parameters,
                                                        plume_name,
                                                        inside_polygon_mask,
                                                        use_dynamic_threshold=use_dynamic_threshold)

            thresholds[plume_name] = the_plume.SPM_threshold
            all_mask_area.append(the_plume.plume_mask)
            if "close_river_mouth_mask" in vars(the_plume):
                all_river_mouth_to_remove.append(the_plume.close_river_mouth_mask)

        # Combine all detected plume areas using logical OR
        final_mask_area = reduce(np.logical_or, all_mask_area)
        final_close_river_mouth_area = reduce(np.logical_or, all_river_mouth_to_remove)

        coordinates_of_the_map = define_parameters(Zone)

        final_mask_area = (final_mask_area
                           .sel(lat=slice(coordinates_of_the_map['lat_range_of_the_map_to_plot'][0],
                                          coordinates_of_the_map['lat_range_of_the_map_to_plot'][1]),
                                lon=slice(coordinates_of_the_map['lon_range_of_the_map_to_plot'][0],
                                          coordinates_of_the_map['lon_range_of_the_map_to_plot'][1])))

        ds_reduced = (ds_reduced
                      .sel(lat=slice(coordinates_of_the_map['lat_range_of_the_map_to_plot'][0],
                                     coordinates_of_the_map['lat_range_of_the_map_to_plot'][1]),
                           lon=slice(coordinates_of_the_map['lon_range_of_the_map_to_plot'][0],
                                     coordinates_of_the_map['lon_range_of_the_map_to_plot'][1])))

        SPM_map = ds_reduced.to_dataframe().reset_index()
        SPM_map['plume'] = final_mask_area.values.flatten()

        if Zone == 'GULF_OF_LION':
            SPM_map.plume[np.where((SPM_map.plume == False) &
                                   (SPM_map.analysed_spim > 10) &
                                   (SPM_map.lat > 43.2) &
                                   (SPM_map.lat < 43.375) &
                                   (SPM_map.lon < 5) &
                                   (SPM_map.lon > 4.7))[0]] = True
            SPM_map.plume[np.where((SPM_map.plume == True) &
                                   (SPM_map.analysed_spim < 10) &
                                   (SPM_map.lon < 4.7))[0]] = False
            # SPM_map.plume[np.where((SPM_map.plume == False) &
            #                        (SPM_map.analysed_spim > thresholds['Grand Rhone']) &
            #                        (SPM_map.lat > 43.2) &
            #                        (SPM_map.lat < 43.35) &
            #                        (SPM_map.lon < 5) &
            #                        (SPM_map.lon > 4.6))[0]] = True

        SPM_map.to_csv(where_to_save_the_figure_5 + f"/DATA/{Zone}.csv")

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_5']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure_5]))

    Figure_2(where_are_saved_regional_maps="output",
             where_to_save_the_figure="figures",
             include_station_points=False)


def Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold,
               where_are_saved_plume_results_with_fixed_threshold,
               where_to_save_the_figure):
    where_to_save_the_figures_6_7 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURES_6_7")
    os.makedirs(where_to_save_the_figures_6_7 + '/DATA', exist_ok=True)

    ts_files_with_dynamic_threshold = glob.glob(os.path.join(where_are_saved_plume_results_with_dynamic_threshold,
                                                             '*', 'SEXTANT', 'SPM', '*', 'Standard', 'PLUME_DETECTION',
                                                             'WEEKLY',
                                                             'Time_series_of_plume_area_and_SPM_threshold.csv'))

    ts_data_with_dynamic_threshold = []
    for ts_file in ts_files_with_dynamic_threshold:
        ts_data_with_dynamic_threshold.append(pd.read_csv(ts_file))

    ts_data_with_dynamic_threshold = pd.concat(ts_data_with_dynamic_threshold)
    ts_data_with_dynamic_threshold['Dynamic_threshold'] = True

    ts_files_with_fixed_threshold = glob.glob(os.path.join(where_are_saved_plume_results_with_fixed_threshold,
                                                           '*', 'SEXTANT', 'SPM', 'merged', 'Standard',
                                                           'PLUME_DETECTION',
                                                           'WEEKLY', 'Time_series_of_plume_area_and_SPM_threshold.csv'))

    ts_data_with_fixed_threshold = []
    for ts_file in ts_files_with_fixed_threshold:
        ts_data_with_fixed_threshold.append(pd.read_csv(ts_file))

    ts_data_with_fixed_threshold = pd.concat(ts_data_with_fixed_threshold)
    ts_data_with_fixed_threshold['Dynamic_threshold'] = False

    pd.concat([ts_data_with_dynamic_threshold, ts_data_with_fixed_threshold]).to_csv(
        os.path.join(where_to_save_the_figures_6_7, 'DATA', 'ts_data.csv'))

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figures_6_7']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figures_6_7]))


def Figure_8_9_10(where_are_saved_X11_results, where_to_save_the_figure):
    where_to_save_the_figures_8_9_10 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURES_8_9_10")
    os.makedirs(where_to_save_the_figures_8_9_10 + '/DATA', exist_ok=True)

    ts_plume_files = glob.glob(
        os.path.join(where_are_saved_X11_results, '*', 'X11_ANALYSIS', 'area_of_the_plume_mask_in_km2',
                     'SEXTANT_merged_Standard_WEEKLY.csv'))

    ts_river_files = glob.glob(
        os.path.join(where_are_saved_X11_results, '*', 'X11_ANALYSIS', 'river_flow', 'River_flow___WEEKLY.csv'))

    regions = ["BAY_OF_BISCAY", "GULF_OF_LION", "BAY_OF_SEINE", "SOUTHERN_BRITTANY"]

    ts_plume_data = []
    for ts_file in ts_plume_files:
        region_found = next((region for region in regions if region in ts_file), None)
        ts_data = pd.read_csv(ts_file)
        ts_data['Zone'] = region_found
        ts_plume_data.append(ts_data)

    ts_river_data = []
    for ts_file in ts_river_files:
        region_found = next((region for region in regions if region in ts_file), None)
        ts_data = pd.read_csv(ts_file)
        ts_data['Zone'] = region_found
        ts_river_data.append(ts_data)

    pd.concat(ts_plume_data).to_csv(os.path.join(where_to_save_the_figures_8_9_10, 'DATA', 'ts_plume_data.csv'))
    pd.concat(ts_river_data).to_csv(os.path.join(where_to_save_the_figures_8_9_10, 'DATA', 'ts_river_data.csv'))

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figures_8_9_10']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figures_8_9_10]))

