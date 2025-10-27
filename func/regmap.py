#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import os, re, gc, pickle, datetime, glob, calendar, math, sys, multiprocess
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from multiprocess import Pool
from scipy.ndimage import binary_dilation
from scipy.stats import lognorm, kstest

multiprocess.set_start_method('spawn', force = True) # MacOS friendly multiprocessing

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (load_file, align_bathymetry_to_resolution, degrees_to_km, find_sat_data_files, expand_grid,
                            path_to_fill_to_where_to_save_satellite_files, fill_the_sat_paths,
                            extract_and_format_date_from_path, define_parameters,
                            unique_years_between_two_dates, get_all_cases_to_process_for_regional_maps_or_plumes_or_X11)


# =============================================================================
#### Utility functions
# =============================================================================


def get_coordinates_of_the_site(Study_area) : 
        
    """
    Get geographic coordinates and spatial extents for a given study area.

    Parameters
    ----------
    Study_area : str
        Name of the study area (e.g., "GULF_OF_LION", "BAY_OF_BISCAY").

    Returns
    -------
    dict
        A dictionary containing basin limits, embouchure, and bloom coordinates and extents.
    """
    
    # Define spatial parameters based on the study area.
    if Study_area == "GULF_OF_LION" : 

        Basin_limits = [41.2, 43.6, 2.75, 9] # (Lat min , Lat max, Lon min, Lon max)
        
        Embouchure_central_point = [43.2, 4.6] # (Lat, Lon) # [43.24, 4.75] 
        Embouchure_lat_extend = 0.5 # 0.55
        Embouchure_lon_extend = 0.75 # 0.75
        
        Bloom_central_point = [42.05, 5.2] # (Lat, Lon)
        Bloom_lat_extend = 1.7 # 1.25
        Bloom_lon_extend = 3.25 # 2.25
        
    elif Study_area == "BAY_OF_BISCAY" : 

        Basin_limits = [43, 49, -7, -0.5] # (Lat min , Lat max, Lon min, Lon max)
        
        Embouchure_central_point = [45.6, -1.5] # (Lat, Lon) # [43.24, 4.75] 
        Embouchure_lat_extend = 1 # 0.55
        Embouchure_lon_extend = 1 # 0.75
        
        # Bloom_central_point = [45.5, -2.5] # (Lat, Lon)
        # Bloom_lat_extend = 2 # 1.25
        # Bloom_lon_extend = 1 # 2.25
        
        Bloom_central_point = [np.nan, np.nan] # (Lat, Lon)
        Bloom_lat_extend = 0 # 1.25
        Bloom_lon_extend = 0 # 2.25
        
    elif Study_area == "SOUTHERN_BRITTANY" : 

         Basin_limits = [46, 48.5, -5, -1.5] # (Lat min , Lat max, Lon min, Lon max)
         
         Embouchure_central_point = [47.125, -2.75] # (Lat, Lon) # [43.24, 4.75] 
         Embouchure_lat_extend = 1.25 # 0.55
         Embouchure_lon_extend = 1.5 # 0.75
         
         Bloom_central_point = [np.nan, np.nan] # (Lat, Lon)
         Bloom_lat_extend = 0 # 1.25
         Bloom_lon_extend = 0 # 2.25
        
    elif Study_area == "FRANCE" : 

        Basin_limits = [41, 52, -8, 11] # (Lat min , Lat max, Lon min, Lon max)
        
        Embouchure_central_point = [np.nan, np.nan] # (Lat, Lon) # [43.24, 4.75] 
        Embouchure_lat_extend = 0 # 0.55
        Embouchure_lon_extend = 0 # 0.75
        
        Bloom_central_point = [np.nan, np.nan] # (Lat, Lon)
        Bloom_lat_extend = 0 # 1.25
        Bloom_lon_extend = 0 # 2.25
        
    elif Study_area == "BAY_OF_SEINE" : 

        Basin_limits = [49.16, 51.38, -1.61, 2.6] # (Lat min , Lat max, Lon min, Lon max)
        
        Embouchure_central_point = [49.5, 0] # (Lat, Lon) # [43.24, 4.75] 
        Embouchure_lat_extend = 0.3 # 0.55
        Embouchure_lon_extend = 0.5 # 0.75
        
        Bloom_central_point = [49.8, 0.15] # (Lat, Lon)
        Bloom_lat_extend = 1 # 1.25
        Bloom_lon_extend = 1.5 # 2.25
        
    elif Study_area == "EASTERN_CHANNEL" : 

        Basin_limits = [49.16, 51.5, -1.61, 3] # (Lat min , Lat max, Lon min, Lon max)
        
        Embouchure_central_point = [50.45, 1.5] # (Lat, Lon) # [43.24, 4.75] 
        Embouchure_lat_extend = 0.9 # 0.55
        Embouchure_lon_extend = 0.5 # 0.75
        
        Bloom_central_point = [np.nan, np.nan] # (Lat, Lon)
        Bloom_lat_extend = 0 # 1.25
        Bloom_lon_extend = 0 # 2.25
        
    elif Study_area == "ETANG_DE_BERRE" : 

        Basin_limits = [43.37, 43.57, 4.97, 5.26] # (Lat min , Lat max, Lon min, Lon max)
        
        Embouchure_central_point = [np.nan, np.nan] # (Lat, Lon) # [43.24, 4.75] 
        Embouchure_lat_extend = 0 # 0.55
        Embouchure_lon_extend = 0 # 0.75
        
        Bloom_central_point = [np.nan, np.nan] # (Lat, Lon)
        Bloom_lat_extend = 0 # 1.25
        Bloom_lon_extend = 0 # 2.25

    else :
        print(f"The zone {Study_area} is not available. Please select one of the following zones : 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'GULF_OF_LION', 'EASTERN_CHANNEL', 'SOUTHERN_BRITTANY'.")
        return None
        
    return {'Basin_limits' : Basin_limits,
            
            'Embouchure_central_point' : Embouchure_central_point,
            'Embouchure_lat_extend' : Embouchure_lat_extend,
            'Embouchure_lon_extend' : Embouchure_lon_extend,
            
            'Bloom_central_point' : Bloom_central_point,
            'Bloom_lat_extend' : Bloom_lat_extend,
            'Bloom_lon_extend' : Bloom_lon_extend}


def extract_key_data(map_ini, info,
                     zone_limits = [], 
                     position_of_the_central_point = [], 
                     lat_extend_of_the_area = [],
                     lon_extend_of_the_area = []) : 
        
    """
    Extract key statistical data from a specific zone in the satellite map.

    Parameters
    ----------
    map_ini : xarray.Dataset
        The satellite data map.
    info : object
        Metadata containing details about the data source and variables.
    zone_limits : list, optional
        Latitude and longitude boundaries for the zone (default is []).
    position_of_the_central_point : list, optional
        Central point coordinates for the zone (default is []).
    lat_extend_of_the_area : list, optional
        Latitude range for the zone (default is []).
    lon_extend_of_the_area : list, optional
        Longitude range for the zone (default is []).

    Returns
    -------
    dict
        Dictionary containing statistical summaries and zone boundaries.

    Notes
    -----
    The function handles both fixed zone boundaries and dynamic extraction based
    on a central point with a given extent.
    """
    
    # Handle fixed zone boundaries
    if zone_limits != [] : 
        
      zone_limits_to_use = zone_limits
      
      map_of_the_zone = map_ini.sel(lat=slice(zone_limits_to_use[0], zone_limits_to_use[1]), 
                                    lon=slice(zone_limits_to_use[2], zone_limits_to_use[3]))
            
    # Handle dynamic central point extraction           
    if position_of_the_central_point != [] : 
        
        if lat_extend_of_the_area == 0 : 
            return None
        
        exact_position_of_the_central_point = map_ini.sel(lat= position_of_the_central_point[0], lon=position_of_the_central_point[1], method='nearest')
        i_lat_of_the_central_point = np.where(map_ini.lat == exact_position_of_the_central_point.lat)[0][0]
        i_lon_of_the_central_point = np.where(map_ini.lon == exact_position_of_the_central_point.lon)[0][0]
        
        lat_resolution = abs( np.mean( np.diff(map_ini.lat) ) )
        lon_resolution = abs( np.mean( np.diff(map_ini.lon) ) )
        
        half_size_of_the_area_in_pixel_LAT = int( ((lat_extend_of_the_area) / (lat_resolution)) / 2 ) 
        half_size_of_the_area_in_pixel_LON = int( ((lon_extend_of_the_area) / (lon_resolution)) / 2 ) 
        
        map_of_the_zone = map_ini.isel(lat=slice(i_lat_of_the_central_point - half_size_of_the_area_in_pixel_LAT, 
                                                 i_lat_of_the_central_point + half_size_of_the_area_in_pixel_LAT +1), 
                                       lon=slice(i_lon_of_the_central_point - half_size_of_the_area_in_pixel_LON, 
                                                 i_lon_of_the_central_point + half_size_of_the_area_in_pixel_LON +1))
        
    # Assign coordinates to the zone data
    map_of_the_zone = map_of_the_zone.assign_coords(day = info.day, week = info.week, month = info.month, year = info.year)
     
    # Determine the variable to extract based on data source and sensor
    if np.isin( info.Data_source, ["SEXTANT", 'SAMUEL'] ) : 
            
          if info.sensor_name == 'merged' :   
              
              map_of_the_zone = map_of_the_zone.drop_vars('time')
              map_of_the_zone = map_of_the_zone.squeeze('time')   
           
              variable_to_extract = [x for x in list(map_of_the_zone.keys()) 
                                     if x.startswith('analysed') and not x.endswith('anomaly')][0]
              
          else : 
              
              if info.Satellite_variable == 'CHLA' : 
                  variable_to_extract = 'chlorophyll_a'
                  
              if 'SPM' in info.Satellite_variable : 
                  variable_to_extract = 'suspended_matters'
                                
    elif info.Data_source == "ODATIS" :       
          
          variable_to_extract = [x for x in list(map_of_the_zone.keys()) 
                                 if x.endswith('mean')][0]
                
    map_values_of_the_zone = map_of_the_zone[variable_to_extract]
    
    # map_flags_of_the_zone = map_of_the_zone[variable_to_extract.replace('mean', 'flags')].astype('int')
    if any( [x in info.Satellite_variable for x in ['SPM', 'CHL']] )  : 
        map_values_of_the_zone = map_values_of_the_zone.where(map_values_of_the_zone >= 0, np.nan)
    
    # Compute statistical summaries for the zone
    map_of_the_zone_log_scale = np.log(map_values_of_the_zone + 1e-6)                                       
    mean_of_the_zone = np.exp( map_of_the_zone_log_scale.mean() ) - 1e-6
    median_of_the_zone = map_values_of_the_zone.median() - 1e-6
    sd_of_the_zone = np.exp( map_of_the_zone_log_scale.std() )
    n_of_the_zone = np.isfinite(map_of_the_zone_log_scale.values).sum()
    n_tot_of_the_zone = map_of_the_zone_log_scale.shape[0] * map_of_the_zone_log_scale.shape[1]
    
    map_values_of_the_zone = map_values_of_the_zone.assign_coords(date_for_plot = info.date_for_plot)
    
    return {'map_data' : map_values_of_the_zone,
            'mean' : float(mean_of_the_zone),
            'median' : float(median_of_the_zone),
            'sd' : float(sd_of_the_zone),
            'n' : n_of_the_zone,
            'n_tot' : n_tot_of_the_zone,
            'zone_limits' : [float(x) for x in [np.min(map_values_of_the_zone.lat), 
                                                np.max(map_values_of_the_zone.lat), 
                                                np.min(map_values_of_the_zone.lon), 
                                                np.max(map_values_of_the_zone.lon)]]}


def load_file_and_extract_key_data(nc_file, info, where_to_save_data_extended, do_the_plot = True) : 
        
    """
    Load a NetCDF file and extract key data for analysis.

    Parameters
    ----------
    nc_file : str
        Path to the NetCDF file.
    info : object
        Metadata containing details about the data source and variables.
    where_to_save_data_extended : str
        Directory to save any output plots or processed data.
    do_the_plot : bool, optional
        Whether to generate and save plots (default is True).

    Returns
    -------
    list or str
        A list containing the extracted data and key maps, or an error message
        if processing fails.

    Notes
    -----
    The function performs various checks for file validity and extracts specific
    data for predefined zones.
    """
    
    # print(nc_file)
    
    # Attempt to open the NetCDF file
    try : 
        with xr.open_dataset(nc_file) as ds :
            map_ini = ds
    except Exception as e : 
        print(f"Error on {nc_file} - Impossible to open the file - {e}")
        return "Impossible to open the file"
    
    # Validate the file structure
    if ( (len(map_ini) == 0) or (len(np.unique(map_ini.lat)) == 1) or (len(np.unique(map_ini.lon)) == 1) ) : 
        return "Format is not correct"
                
    date = extract_and_format_date_from_path(nc_file)
      
    # Update info object with date components
    info['day'] = pd.to_datetime(date).strftime("%Y-%m-%d") 
    info['month'] = pd.to_datetime(date).strftime("%m") 
    info['year'] = pd.to_datetime(date).strftime("%Y") 
    nb_of_day = pd.to_datetime(date).strftime("%d")     
    info['week'] = f'{info.month}_{determine_the_week_based_on_the_day(nb_of_day)}'
    info['date_for_plot'] = info.day
    info['period_name'] = 'daily'
    
    # Sort data by latitude and longitude
    map_ini = map_ini.sortby('lat')
    map_ini = map_ini.sortby('lon')
    
    # Get coordinates for the site
    site_coordinates = get_coordinates_of_the_site(info.Zone)

    # Extract data for the Basin zone
    Basin_data = extract_key_data(map_ini, info, zone_limits = site_coordinates['Basin_limits'])
             
    if Basin_data['n'] == 0 : 
        return "All Basin data are NAN"
    
    # Placeholder for additional site data
    # Embouchure site
    # Embouchure_data = extract_key_data(map_ini, info, 
    #                                    position_of_the_central_point = coordinates['Embouchure_central_point'], 
    #                                    lat_extend_of_the_area = coordinates['Embouchure_lat_extend'], 
    #                                    lon_extend_of_the_area = coordinates['Embouchure_lon_extend'])
    Embouchure_data = None
    
    # Bloom site
    # Bloom_data = extract_key_data(map_ini, info, 
    #                                 position_of_the_central_point = coordinates['Bloom_central_point'], 
    #                                 lat_extend_of_the_area = coordinates['Bloom_lat_extend'], 
    #                                 lon_extend_of_the_area = coordinates['Bloom_lon_extend'])
    Bloom_data = None
    
    # Generate and save plot if requested
    if do_the_plot :
    
        do_and_save_the_plot(info, 
                             Basin_data, Embouchure_data, Bloom_data, 
                             where_to_save_data_extended)
                           
    
    # Prepare the data to return                        
    key_data_to_return = ['mean', 'median', 'sd', 'n', 'n_tot']
    data_to_return = {'Basin_data' : {key: Basin_data[key] for key in key_data_to_return if Basin_data is not None},
                      'Embouchure_data' : {key: Embouchure_data[key] for key in key_data_to_return if Embouchure_data is not None},
                      'Bloom_data' : {key: Bloom_data[key] for key in key_data_to_return if Bloom_data is not None}}
    
    data_to_return = [pd.DataFrame.from_dict(data_to_return[f'{site}_data'], orient='index').T.rename(columns=lambda x: f'{site}_{x}') 
                      for site in ['Basin', 'Embouchure', 'Bloom']]
    data_to_return = pd.concat(data_to_return, axis = 1)
    data_to_return['day'] = info.day
    
    map_ini.close()
        
    gc.collect()
    
    return [data_to_return, Basin_data, Embouchure_data, Bloom_data]


def Process_each_week(Year_month_week_pattern, 
                      where_to_save_data_extended, all_days_of_the_year, 
                      suffix_ranges, info, 
                      map_files, files_have_been_processed,
                      save_map_plots_of_which_time_frequency) :

    """
    Process satellite data files for each week and generate summary outputs.

    Parameters
    ----------
    Year_month_week_pattern : str
        Pattern representing the year, month, and week to process.
    where_to_save_data_extended : str
        Directory where processed data will be saved.
    all_days_of_the_year : list of str
        List of all dates in the year.
    suffix_ranges : list
        Suffix ranges used to identify relevant files.
    info : object
        Contains metadata such as data source, sensor name, and variable details.
    map_files : list of str
        List of file paths for the current year's satellite maps.
    files_have_been_processed : DataFrame
        DataFrame tracking processing statuses for files.
    save_map_plots_of_which_time_frequency : dict
        Dictionary indicating which time frequencies to save plots for.

    Returns
    -------
    list
        A list containing:
        - weekly_results_ts : DataFrame
            Time series of the weekly results.
        - files_have_been_processed_of_the_week : DataFrame
            Updated DataFrame tracking processing statuses.

    Notes
    -----
    This function handles weekly satellite data processing, including file
    filtering, data extraction, and summary statistics computation.
    """
    
    # For debugging
    # Year_month_week_pattern = self.Year_month_week_patterns[0]
    # where_to_save_data_extended = self.where_to_save_data_extended
    # all_days_of_the_year = self.all_days_of_the_year
    # suffix_ranges = self.suffix_ranges
    # info = self.info
    # map_files = self.map_files
    # files_have_been_processed = self.files_have_been_processed

    # Get file patterns for the given week
    Year_month_week_days = get_file_patterns(Year_month_week_pattern, suffix_ranges)
    
    # # For SEXTANT merged data, include the previous day in the processing
    # if (info.Data_source == 'SEXTANT') and (info.sensor_name == 'merged') :
    #     previous_day = str( previous_day_pattern(Year_month_week_pattern, suffix_ranges) )
    #     Year_month_week_days = [previous_day] + Year_month_week_days
        
    # Filter map files matching the patterns
    map_files_of_the_week = sorted( [x for x in map_files if any([pattern in x for pattern in Year_month_week_days])] )
        
    # Load files and extract key data
    weekly_results = [load_file_and_extract_key_data(nc_file, info, where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'DAILY'), 
                                                     save_map_plots_of_which_time_frequency['DAILY']) 
                      for nc_file in map_files_of_the_week]

    # # Handle SEXTANT merged data specific computations
    # if (info.Data_source == 'SEXTANT') and (info.sensor_name == 'merged') :
        
    #     compute_diff_maps_and_save_them(weekly_results, where_to_save_data_extended, info.Satellite_variable)

    #     Year_month_week_days = Year_month_week_days[1:]
    #     if previous_day in map_files_of_the_week[0] :
    #         map_files_of_the_week = map_files_of_the_week[1:]
    #         weekly_results = weekly_results[1:]
        
    # Save daily maps to pickle files
    [pickle.dump({'Basin_map' : x[1]}, open(f"{where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'DAILY')}/{x[0].day[0]}.pkl", 'wb')) for x in weekly_results if isinstance(x, list)]
               
    # Index the files for the current week                 
    map_files_of_the_week_index = [ all_days_of_the_year.index(date) for date in [extract_and_format_date_from_path(x) for x in map_files_of_the_week ] ]

    # Update the processing status DataFrame
    files_have_been_processed_of_the_week = files_have_been_processed.iloc[map_files_of_the_week_index].copy()
    files_have_been_processed_of_the_week.Impossible_to_open_the_file = np.array([bool(x == 'Impossible to open the file') for x in weekly_results] )
    files_have_been_processed_of_the_week.Format_not_correct = [bool(x == "Format is not correct") for x in weekly_results] 
    files_have_been_processed_of_the_week.All_values_are_NAN = [bool(x == 'All Basin data are NAN') for x in weekly_results] 

    # Filter valid results
    weekly_results = [x for x in weekly_results if isinstance(x, list)]
    
    if (len(weekly_results) == 0) : 
        return [None, files_have_been_processed_of_the_week]
        # continue
       
    # Combine weekly results into a single DataFrame
    weekly_results_ts = pd.concat( [x[0] for x in weekly_results] ).reset_index(drop = True)
    
    # Compute the date of the weekly map
    date_of_the_weekly_map = np.mean([int(x.replace('/', '')) for x in Year_month_week_days]).astype(int).astype(str)
    
    # Compute and save the mean map for the week
    # NB: There is an issue in here with arguments not being ordered correctly or called explictly
    get_the_mean_map_and_save_it(where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'WEEKLY'), 
                                [{'Basin_map' : x[1], 'Embouchure_map' : x[2],'Bloom_map' : x[3]} for x in weekly_results],
                                info,
                                Year_month_week_pattern[4:], Year_month_week_pattern[4:],
                                "WEEKLY", save_map_plots_of_which_time_frequency['WEEKLY'], date_of_the_weekly_map)

    del weekly_results
    gc.collect()
    
    return [weekly_results_ts, files_have_been_processed_of_the_week]

        
def get_file_patterns(pattern, suffix_ranges):
    
    """
    Generate a list of file patterns by replacing suffixes in a given pattern with ranges of values.
    
    Parameters
    ----------
    pattern : str
        The base pattern containing a suffix to replace.
    suffix_ranges : dict
        A dictionary where keys are suffixes and values are ranges of integers.
    
    Returns
    -------
    list or str
        A list of patterns with suffixes replaced by range values, or the original pattern if no match is found.
    """
    
    for suffix, rng in suffix_ranges.items():
        if suffix in pattern:
            # Replace the suffix in the pattern with all values in the range.
            to_return = [pattern.replace(suffix, f"{x:02d}") for x in rng]
            
            to_return_formatted = [f"{x[:4]}/{x[4:6]}/{x[6:]}" if len(x) == 8 else None for x in to_return]
            
            return to_return_formatted
            
            # Return the original pattern if no suffix matches.
    return pattern    


def previous_day_pattern(current_pattern, suffix_ranges):
    
    """
    Compute the pattern for the previous day based on the current pattern and suffix ranges.

    Parameters
    ----------
    current_pattern : str
        The pattern for the current day (formatted as 'YYYYMM_WW').
    suffix_ranges : dict
        A dictionary of suffixes and their corresponding ranges.

    Returns
    -------
    int
        The pattern of the previous day based on suffix ranges.
    """
   
    # Extract year, month, and week from the current pattern.
    year = int(current_pattern[:4])
    month = int(current_pattern[4:6])
    week = int(current_pattern[7:])
    
    # Calculate the start date of the week in the current pattern.
    first_day_of_month = datetime.date(year, month, 1)
    start_of_week = first_day_of_month + datetime.timedelta(weeks=week-1, days=(7 - first_day_of_month.weekday()) % 7)
    
    # Calculate the start date of the previous week.
    start_of_previous_week = start_of_week - datetime.timedelta(weeks=1)
    
    # Extract year, month, and week number for the previous week.
    previous_year = start_of_previous_week.year
    previous_month = start_of_previous_week.month
    previous_week = ((start_of_previous_week - datetime.date(previous_year, previous_month, 1)).days // 8) + 1
    
    # Format the pattern for the previous week.
    previous_pattern = f"{previous_year:04d}{previous_month:02d}_{previous_week:02d}"
    
    # Get the maximum pattern value for the previous week from the suffix ranges.
    pattern_of_the_previous_day = np.max( [ int(x) for x in get_file_patterns(previous_pattern, suffix_ranges) ] )
    
    return pattern_of_the_previous_day


def determine_the_week_based_on_the_day(num_of_day) : 
    
    """
    Determine the week of the month based on the day number.

    Parameters
    ----------
    num_of_day : int
        Day of the month (1-31).

    Returns
    -------
    str
        The week of the month ('01', '02', '03', '04').
    """
    
    # Assign the week based on the day number.
    if int(num_of_day) in range(0,9): 
        week = '01'
    if int(num_of_day) in range(9,17): 
        week = '02'
    if int(num_of_day) in range(17,25): 
        week = '03'
    if int(num_of_day) in range(25,32): 
        week = '04'
        
    return week


def do_and_save_the_plot(info, Basin_data, Embouchure_data, Bloom_data, 
                         where_to_save_data_extended, 
                         climatological_map = False,
                         climatological_subfolder = "") : 
    
    """
    Generate and save a plot based on satellite data.

    Parameters
    ----------
    info : object
        Contains metadata about the satellite data (e.g., variable type, source, etc.).
    Basin_data : dict
        Contains map data for the basin region.
    Embouchure_data : dict or None
        Contains map data for the embouchure region or None if not applicable.
    Bloom_data : dict or None
        Contains map data for bloom regions or None if not applicable.
    where_to_save_data_extended : str
        Path where the plot will be saved.
    climatological_map : bool, optional
        Whether the plot represents climatological data. Default is False.
    climatological_subfolder : str, optional
        Subfolder name for saving climatological data. Default is an empty string.

    Notes
    -----
    This function uses matplotlib to create plots, applies custom color maps, 
    and overlays rectangles based on the provided regions.
    """
        
    # Initialize variable settings based on the satellite variable
    if info.Satellite_variable == "CHLA" : 
        unit_satellite_variable = 'mg.m-3'
        colorbar_limits = [5e-1, 5e0]
        colorbar_palette = 'viridis'
        rectangle_color = "red"
        
    elif "SPM" in info.Satellite_variable : 
        unit_satellite_variable = 'g.m-3'
        colorbar_limits = [1e0, 1e2]
        colorbar_palette = 'Spectral_r'
        rectangle_color = "black"
        
    elif "SST" in info.Satellite_variable : 
        unit_satellite_variable = '°C'
        colorbar_limits = [0.1, 30]
        colorbar_palette = 'Spectral_r'
        rectangle_color = "black"
        
    else : 
        unit_satellite_variable = 'unknown unit'
        colorbar_limits = [1e-1, 1e2]
        colorbar_palette = 'viridis'
        rectangle_color = "green"
      
    # Adjust colorbar limits for climatological maps
    if climatological_map : 
        colorbar_limits = np.array(colorbar_limits) / 3
            
    # Set figure size and adjust limits for specific regions
    if any( [x in where_to_save_data_extended for x in ["BAY_OF_SEINE", "EASTERN_CHANNEL", "GULF_OF_LION", "ETANG_DE_BERRE"]] ) : 
        fig_size = (20,10)
        colorbar_limits = colorbar_limits * 2
    
    elif any( [x in where_to_save_data_extended for x in ["BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "FRANCE"]] ) : 
        fig_size = (20,15)

    else :
        fig_size = (15,10)
                            
    size_of_axis_title = 20
    fig, ax = plt.subplots(figsize = fig_size)
    
    # Plot the basin map
    if 'map_data' in Basin_data : 
        Zone_map = Basin_data['map_data'].plot(vmin = colorbar_limits[0], vmax = colorbar_limits[1], norm=colors.LogNorm(), cmap = colorbar_palette)  
    else : 
        Zone_map = Basin_data.plot(vmin = colorbar_limits[0], vmax = colorbar_limits[1], norm=colors.LogNorm(), cmap = colorbar_palette)  
        
    # Overlay embouchure region rectangle
    if Embouchure_data is not None : 
        Embouchure_rectangle = patches.Rectangle((Embouchure_data['zone_limits'][2], Embouchure_data['zone_limits'][0]),
                                                  Embouchure_data['zone_limits'][3] - Embouchure_data['zone_limits'][2], 
                                                  Embouchure_data['zone_limits'][1] - Embouchure_data['zone_limits'][0], 
                                                  linewidth=1.5, edgecolor=rectangle_color, facecolor='none')
        ax.add_patch(Embouchure_rectangle)
    
    # Overlay bloom region rectangle
    if Bloom_data is not None : 
        Bloom_rectangle = patches.Rectangle((Bloom_data['zone_limits'][2], Bloom_data['zone_limits'][0]),
                                             Bloom_data['zone_limits'][3] - Bloom_data['zone_limits'][2], 
                                             Bloom_data['zone_limits'][1] - Bloom_data['zone_limits'][0], 
                                             linewidth=1.5, edgecolor=rectangle_color, facecolor='none')
        
        ax.add_patch(Bloom_rectangle)
        
    # Create figure title
    if 'date_for_plot' in Basin_data["map_data"].coords : 
        figure_title = f'{info.Satellite_variable} {info.day.replace("the year", str(Basin_data["map_data"].date_for_plot.values.astype("datetime64[Y]").astype(int) + 1970)).replace(f" {info.period_name[:2]}", " " + calendar.month_name[int(info.period_name[:2]) if info.period_name[:2].isdigit() else 0])}'
    else : 
        figure_title = f'{info.Satellite_variable} {info.day}'
    figure_title_suffix = f'{info.Data_source} / {info.sensor_name} / {info.atmospheric_correction}'
    figure_complete_title = f'{figure_title} ({figure_title_suffix})'
    ax.set_title(figure_complete_title, fontsize = size_of_axis_title*1.25)    
    
    # Set axis labels and ticks
    ax.set_xlabel('', fontsize=size_of_axis_title)
    ax.set_ylabel('', fontsize=size_of_axis_title)
    ax.tick_params(axis='both', which='major', labelsize=size_of_axis_title*0.75)   
            
    # Configure colorbar    
    cbar = ax.collections[0].colorbar        
    cbar.set_label(f'{info.Satellite_variable} ({unit_satellite_variable})', fontsize=size_of_axis_title)
    cbar.ax.tick_params(labelsize=size_of_axis_title*0.75)
    
    plt.tight_layout()
    
    plt.savefig(f'{where_to_save_data_extended}/{info.day.replace(" ", "_")}.png')
            
    plt.close(fig)  


def compute_diff_maps_and_save_them(weekly_results, where_to_save_data_extended, satellite_variable):
        
    """
    Computes relative difference maps from weekly results and saves them as pickle files.

    Parameters
    ----------
    weekly_results : list
        A list of dictionaries containing weekly map data.
    where_to_save_data_extended : str
        Directory path where the output difference maps will be saved.
    satellite_variable : str
        Name of the satellite variable being processed.

    Returns
    -------
    None
    """
    
    def rel_diff(da1, da2):
        
        """
        Calculates the relative difference between two data arrays.

        Parameters
        ----------
        da1, da2 : xarray.DataArray
            Data arrays to compute the relative difference.

        Returns
        -------
        xarray.DataArray
            The relative difference normalized by the time difference.
        """
        
        # Parse dates from the data arrays
        date_da1 = datetime.datetime.strptime( np.array(da1['day']).tolist() , '%Y-%m-%d')
        date_da2 = datetime.datetime.strptime( np.array(da2['day']).tolist() , '%Y-%m-%d')
        date_diff = (date_da1 - date_da2).days
        
        # Align the data arrays and compute relative difference
        aligned1, aligned2 = xr.align(da1, da2, join='outer', copy=False)
        result = (aligned1 - aligned2) / aligned2 / date_diff
        result.coords.update(aligned1.coords)
        result = result.rename(result.name + "_relative_difference")
        
        return result
    
    def compute_the_rel_diff_maps(index_name = "") :
            
        """
        Computes the relative difference maps for a specific index.

        Parameters
        ----------
        index_name : str, optional
            Key in the weekly results dictionary to process.

        Returns
        -------
        dict or None
            Dictionary with relative difference maps, or None if no maps are found.
        """
        
        # Extract maps for the given index
        maps = [ x[index_name]['map_data'] for x in maps_to_process if isinstance(x[index_name], str) == False ]        

        if len(maps) == 0 : 
            return
             
        # Compute relative difference maps                                           
        rel_diff_maps = [rel_diff(maps[i], maps[i-1]) for i in range(1, len(maps))]
               
        del maps # Free memory
        gc.collect()
        
        return {'rel_diff_maps' : rel_diff_maps}
    
    # Prepare data for processing
    maps_to_process = [{'Basin_map' : x[1], 'Embouchure_map' : x[2],'Bloom_map' : x[3]} for x in weekly_results]
    
    # Compute relative difference for basin maps
    Basin_rel_diff_map = compute_the_rel_diff_maps(index_name = "Basin_map")
    if Basin_rel_diff_map is None : 
        return
    
    # Define the save directory
    folder_where_to_save_files = where_to_save_data_extended.replace('DAILY', 'DIFF')
    os.makedirs(folder_where_to_save_files, exist_ok=True)
    
    # Save relative difference maps to pickle files
    [pickle.dump(map_obj, open(f"{folder_where_to_save_files}/Diff_{np.array(map_obj.coords['day']).tolist()}.pkl", 'wb')) 
     for map_obj in Basin_rel_diff_map['rel_diff_maps']]

    del Basin_rel_diff_map # Free memory
    
    
def get_the_mean_map_and_save_it(where_to_save_data_extended, maps_of_the_period,
                                 info, period_name, months_to_use = [str(x+1).zfill(2) for x in range(12)],
                                 climatological_subfolder = "", do_the_plot = True, date_for_plot = '') : 

    """
    Computes the mean map for a given period and saves it as a pickle file.

    Parameters
    ----------
    where_to_save_data_extended : str
        Path to save the mean maps.
    maps_of_the_period : list
        List of maps for the specified period.
    info : dict
        Information about the satellite variable and processing details.
    period_name : str
        Name of the period for averaging.
    months_to_use : list of str, optional
        List of month indices to include in the averaging (default: all months).
    climatological_subfolder : str, optional
        Subfolder for climatological data.
    do_the_plot : bool, optional
        Whether to plot the mean map (default: True).
    date_for_plot : str, optional
        Date string for plotting purposes.

    Returns
    -------
    None
    """
        
    min_number_of_maps = {'WEEKLY' : 1, # Minimal number of daily maps per week
                          'MONTHLY' : 2, # Minimal number of weekly maps per month
                          'ANNUAL' : 12} # Minimal number of monthly maps per year
    
    def compute_the_mean_map(index_name = "") :
        
        """
        Computes the mean map for a specific index.

        Parameters
        ----------
        index_name : str, optional
            Key in the maps of the period dictionary to process.

        Returns
        -------
        dict or None
            Dictionary with mean map data, or None if no maps are found.
        """
          
        # Check if the first map in the period has data for the index
        if maps_of_the_period[0][index_name] is None : 
            return
                
        zone_limits = maps_of_the_period[0][index_name]['zone_limits']
        maps = [ x[index_name]['map_data'] for x in maps_of_the_period ]
            
        the_averaged_period = f'Averaged over {period_name}'
                    
        if len(maps) == 0 : 
            return
                        
        # Log-transform maps if applicable
        if any( [x in info.Satellite_variable for x in ['SPM', 'CHL', 'SST']] ) :
            maps = [np.log(map_data + 1e-6) for map_data in maps]
            
        # Concatenate and compute the mean along the appropriate dimension
        if 'day' in list( maps[0].coords ) :
            
            combined_maps = xr.concat(maps, dim='day', coords='different', compat='equals', join='outer')
            mean_map = combined_maps.mean(dim='day', skipna=True)
            
        elif 'week' in list( maps[0].coords ) : 
            
            combined_maps = xr.concat(maps, dim='week', coords='different', compat='equals', join='outer')
            mean_map = combined_maps.mean(dim='week', skipna=True)
            
        elif 'month' in list( maps[0].coords ) : 
            
            combined_maps = xr.concat(maps, dim='month', coords='different', compat='equals', join='outer')
            mean_map = combined_maps.mean(dim='month', skipna=True)
            
        elif 'year' in list( maps[0].coords ) : 
            
            combined_maps = xr.concat(maps, dim='year', coords='different', compat='equals', join='outer')
            mean_map = combined_maps.mean(dim='year', skipna=True)

        else : 
            print("No time dimension found in the maps")
            return
        
        # Reverse log-transform if needed
        if any( [x in info.Satellite_variable for x in ['SPM', 'CHL', 'SST']] ) :
            
            mean_map = np.exp( mean_map ) - 1e-6
             
        # Assign additional metadata
        mean_map = mean_map.assign_coords(period = the_averaged_period,
                                          date_for_plot = datetime.datetime.strptime( date_for_plot , '%Y%m%d') if date_for_plot != '' else date_for_plot )
        
        del maps, combined_maps # Free memory
        gc.collect()
        
        return {'map_data' : mean_map,
                'zone_limits' : zone_limits}

    # Ensure there are maps to process
    if len(maps_of_the_period) < (min_number_of_maps[climatological_subfolder] if climatological_subfolder in min_number_of_maps.keys() else 1) : 
        return
    
    Embouchure_mean_map = None # Embouchure_mean_map = compute_the_mean_map(index_name = "Embouchure_map")

    Bloom_mean_map = None # Bloom_mean_map = compute_the_mean_map(index_name = "Bloom_map")
           
    Basin_mean_map = compute_the_mean_map(index_name = "Basin_map")
        
    # Add metadata to the info dictionary
    info['day'] = str(Basin_mean_map['map_data']['period'].values)
    info['period_name'] = period_name
        
    # Get coordinates for the site
    # site_coordinates = get_coordinates_of_the_site(info.Zone)

    # Extract data for the Basin zone
    # Basin_data = extract_key_data(map_ini, info, zone_limits = site_coordinates['Basin_limits'])
             
    # if Basin_data['n'] == 0 : 
        # return "All Basin data are NAN"

    # Save the mean maps and plot them if required
    if (Basin_mean_map is not None) or (Bloom_mean_map is not None) or (Embouchure_mean_map is not None) :

        if do_the_plot : 
            
            do_and_save_the_plot(info, 
                                 Basin_mean_map, Embouchure_mean_map, Bloom_mean_map, 
                                 where_to_save_data_extended, climatological_map=True, 
                                 climatological_subfolder = climatological_subfolder)
                    
        with open(f'{where_to_save_data_extended}/{info.day.replace(" ", "_")}.pkl', 'wb') as f:
            pickle.dump({'Embouchure_map' : Embouchure_mean_map,
                          'Bloom_map' : Bloom_mean_map,
                          'Basin_map' : Basin_mean_map}, f)
         
    del Embouchure_mean_map, Bloom_mean_map, Basin_mean_map # Free memory
    gc.collect()


def load_the_climatological_files(where_to_save_data_extended, month_nb = '', return_file_names = False):
    
    """
    Loads climatological files for the specified satellite variable and subfolder.

    Parameters
    ----------
    where_to_save_data_extended : str
        Base path where the climatological files are saved.
    satellite_variable : str
        Name of the satellite variable.
    climatological_subfolder : str
        Subfolder name for climatological data.
    month_nb : str, optional
        Month number for filtering files (default: all months).

    Returns
    -------
    list
        A list of loaded climatological data.
    """
    
    # Construct file search pattern
    if "WEEKLY" in where_to_save_data_extended : 
        pattern = os.path.join(where_to_save_data_extended, f'*{month_nb:02d}_[0-9][0-9].pkl')
    elif "MONTHLY" in where_to_save_data_extended :
        pattern = os.path.join(where_to_save_data_extended, '*[0-9][0-9].pkl')
    else :
        pattern = os.path.join(where_to_save_data_extended, '*.pkl')
        
    file_paths = glob.glob(pattern)
    
    if return_file_names : 
        
        return file_paths        
    
    else : 
        
        # Get file paths and load them using multiprocessing
        with Pool() as pool:
            data = pool.map(load_file, file_paths)
        
        return data


def Compute_the_metrics_of_one_map(file, exclude_coastal_areas, coastal_waters_mask) : 
    
    # map_data = load_file(file)['map_data']
    map_data = load_file(file)['Basin_map']['map_data']

    # 0. Mask out coastal pixels 
    if exclude_coastal_areas : 
        # map_data.values[ coastal_waters_mask.values ] = float('NaN')
        map_data_values_to_use = map_data.values[ ~ coastal_waters_mask.values ]
    else : 
        map_data_values_to_use = map_data.values
    
    # 1. Mask out NaN values
    index_finite_values = np.isfinite( map_data_values_to_use )
    # TODO: This appears to receive 0 or NaN values somehow and throughs a lot of warnings
    # I'm thinking that the masking above isn't working as intended
    # It may be that the coastal_waters_mask is not being applied correctly
    # Or that the map_data_values_to_use is not being populated correctly
    # Or that the map_data itself is not being loaded correctly
    cloud_percentage = 100 * (~index_finite_values).sum() / map_data_values_to_use.size
    filtered_values = map_data_values_to_use[ index_finite_values ]
    
    if len(filtered_values) < 10 : 
        return None
    
    # 2. Compute the 99th percentile to exclude outliers
    percentile_99 = np.nanpercentile(filtered_values, 99)
  
    # 3. Log-normal fit to the filtered values
    # Exclude zeros for log-normal fit
    positive_values = filtered_values[filtered_values > 0]
    shape, loc, scale = lognorm.fit( positive_values , floc=0)  # floc=0 for a proper log-normal fit
  
    # import seaborn as sns
    # sns.distplot(positive_values)
    # np.percentile(positive_values, 99.9)
    # np.where( positive_values > 40 )
    
    n_outliers_higher_than_50 = sum( positive_values > 50 )
    
    # 4. Kolmogorov-Smirnov test for goodness-of-fit
    ks_stat, p_value = kstest(positive_values, lognorm(s=shape, loc=loc, scale=scale).cdf)
    
    # 5. Compute skewness and kurtosis
    
    # Skewness:
    # Describes the asymmetry of the distribution.
    # A positive skew suggests a long tail on the right (a few pixels have very high values).
    # A near-zero skew indicates a symmetric distribution.

    # Kurtosis:
    # Measures the "tailedness" of the Chla distribution.
    # High kurtosis implies extreme values are more frequent (peaky distribution with heavy tails).
    # Low kurtosis indicates a flatter distribution.
    
    # data_skewness = skew(filtered_values)
    # data_kurtosis = kurtosis(filtered_values, fisher=False) 

    # 6. Compute the average
    mean_value = np.mean(filtered_values)

    # Results
    metrics = {
        "file" : file,
        "date": map_data.day.item(),
        "cloud_percentage": cloud_percentage,
        "mean_value": mean_value,
        "99th_Percentile": percentile_99,
        "Lognorm_shape": shape, 
        # "Lognorm_loc": loc, 
        # "Lognorm_scale": scale,
        # "KS_Test_p-value": p_value,
        # "Skewness": data_skewness,
        # "Kurtosis": data_kurtosis,
        "n_outliers" : n_outliers_higher_than_50
    }

    return metrics


def plot_and_save_the_QC_metrics(QC_df, metrics_to_plot, path_to_save_QC_files, info, max_cloud_cover = 100) : 

    # Plot the time series
    subplot_titles = {  "mean_value" : "Average value",
                        "99th_Percentile" : 'Are the high values within a realistic range ?',
                        "Lognorm_shape" : 'Are value distributions similar across days ?',
                        "Skewness" : 'Are there outliers ?',
                        "n_outliers" : 'How many outliers are there ? (An outlier is considered as a value > 50)'}

    # TODO: It may be that this is causing the entire dataframe to be removed,
    # so the use of QC_df_for_plot["date"] below is throwing an error because it can't subset an empty dataframe
    QC_df_for_plot = QC_df[ QC_df.cloud_percentage <= max_cloud_cover ]

    fig, axes = plt.subplots(nrows=len(metrics_to_plot), ncols=1, figsize=(14, 17), sharex=True)        
    for ax, metric in zip(axes, metrics_to_plot):
        ax.plot(QC_df_for_plot["date"], QC_df_for_plot[metric], label=metric.replace("_", " "), linewidth=2)
        ax.set_title(subplot_titles[metric], fontsize=16)
        ax.legend(loc="upper left")
    
    fig.suptitle(f"QC of {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable} in {info.Zone.replace('_', ' ')}\n",
                 fontsize=20)
    
    # Adjust layout
    plt.tight_layout()
    plt.xlabel("")
    
    plt.savefig(f'{path_to_save_QC_files}/QC_{info.Satellite_variable}_daily_maps.png')
    QC_df.to_csv(f'{path_to_save_QC_files}/QC_{info.Satellite_variable}_daily_maps.csv', index = False)
    
    plt.show()
    plt.close()
    
    return fig


# =============================================================================
#### Classes 
# =============================================================================


class Create_and_save_the_maps : 
    
    """
    A class for creating and saving satellite maps at different temporal resolutions (weekly, monthly, annual).
    """
        
    def __init__(self, where_to_save_regional_maps, where_are_saved_satellite_data, info) :
        
        """
        Initialize the Create_and_save_the_maps object.
 
        Parameters
        ----------
        working_directory : str
            Working directory where data will be saved.
        where_are_saved_satellite_data : str
            Path to the satellite data files.
        info : object
            Metadata about the satellite data (e.g., zone, data source, sensor name, year).
        """
                
        # Define where to save processed data and time series data
        where_to_save_data = fill_the_sat_paths(info, 
                                                path_to_fill_to_where_to_save_satellite_files(where_to_save_regional_maps + '/' + info.Zone).replace('[TIME_FREQUENCY]', ''),
                                                local_path = True).replace('/*/*/*', '')
        
        where_to_save_timeseries_data = f'{where_to_save_data}/TIME_SERIES/'
        time_series_file_name = f'{where_to_save_timeseries_data}/{info.Satellite_variable}_{info.Data_source}_{info.sensor_name}_{info.atmospheric_correction}.csv'

        # Find satellite data files for the current and previous years
        if isinstance(info['Year'], str) and (info['Year'] == 'MULTIYEAR') : 
            
            all_days_of_the_year = pd.date_range(start=info.date_min, end=info.date_max, freq="YE").strftime("%Y%m%d").tolist()

            map_files = find_sat_data_files(info, 
                                            fill_the_sat_paths(info.replace({info.Temporal_resolution : "DAILY"}),
                                                               (path_to_fill_to_where_to_save_satellite_files(f'{where_to_save_regional_maps}/{info.Zone}').replace('/[MONTH]/[DAY]', '').replace('[TIME_FREQUENCY]', 'MAPS/[TIME_FREQUENCY]')), 
                                                               local_path = True, 
                                                               dates = all_days_of_the_year))
            
            self.info = info
            self.where_to_save_data = where_to_save_data
            self.where_to_save_timeseries_data = where_to_save_timeseries_data
            self.time_series_file_name = time_series_file_name
            self.map_files = map_files
            self.all_days_of_the_year = all_days_of_the_year
            
            return
            
        all_days_of_the_year = pd.date_range(start = info.date_min, end = info.date_max, freq = "D").strftime("%Y%m%d").tolist()
        
        map_files = find_sat_data_files(info, 
                                        fill_the_sat_paths(info.replace({info.Temporal_resolution : "DAILY"}), 
                                                           path_to_fill_to_where_to_save_satellite_files(where_are_saved_satellite_data), 
                                                           local_path = True, 
                                                           dates = all_days_of_the_year))            
            
        # Check if satellite files are available
        if not map_files : 
            print("No satellite file found")
            return
        
        # Extract days from filenames and generate all days of the year
        map_file_days = [extract_and_format_date_from_path(x) for x in map_files]
                
        # Create a DataFrame to track file processing status
        files_have_been_processed = pd.DataFrame({'day' : all_days_of_the_year,
                                                  'Does_the_file_exist' : np.isin(all_days_of_the_year, map_file_days),
                                                  'Impossible_to_open_the_file' : '',
                                                  'Format_not_correct' : '',
                                                  'All_values_are_NAN' : ''})
        
        # Generate weekly patterns for the year
        Year_month_patterns = np.unique([x[:6] for x in map_file_days])
        Year_month_week_patterns = [f"{item}_{suffix:02d}" for item in Year_month_patterns for suffix in range(1, 5)]
        suffix_ranges = { '_01': range(1, 9), '_02': range(9, 17), '_03': range(17, 25), '_04': range(25, 32) }
        
        # TODO: This does not correctly detect if the files have been generated. Rather that the final time series file exists
        # Meaning that, if this file has been paritally generated, the code will think that all files have been processed
        are_the_maps_already_produced = os.path.isfile( time_series_file_name )

        self.info = info
        self.where_to_save_data = where_to_save_data
        self.where_to_save_timeseries_data = where_to_save_timeseries_data
        self.time_series_file_name = time_series_file_name
        self.map_files = map_files
        self.all_days_of_the_year = all_days_of_the_year
        self.files_have_been_processed = files_have_been_processed
        self.Year_month_week_patterns = Year_month_week_patterns
        self.suffix_ranges = suffix_ranges        
        self.are_the_maps_already_produced = are_the_maps_already_produced
        
    def _1_create_weekly_maps(self, nb_of_cores_to_use, save_map_plots_of_which_time_frequency) :
        
        """
        Create weekly maps by processing daily satellite files.

        Parameters
        ----------
        nb_of_cores_to_use : int
            Number of processor cores to use for parallel processing.
        """
        
        # Ensure the directory for daily maps exists
        self.where_to_save_data_extended = f'{self.where_to_save_data}/MAPS/[TIME_FREQUENCY]/{self.info.Year}'
        os.makedirs(self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'WEEKLY'), exist_ok=True)
        os.makedirs(self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'DAILY'), exist_ok=True)
                 
        # where_to_save_data_extended = self.where_to_save_data_extended
        # all_days_of_the_year = self.all_days_of_the_year 
        # suffix_ranges = self.suffix_ranges
        # info = self.info
        # map_files = self.map_files
        # files_have_been_processed = self.files_have_been_processed
        # Year_month_week_patterns = self.Year_month_week_patterns
        
        # Use multiprocess to process each week
        # pool = multiprocess.Pool(nb_of_cores_to_use)
        with multiprocess.Pool(nb_of_cores_to_use) as pool:

            results = pool.starmap(Process_each_week, [(Year_month_week_pattern, 
                                                        self.where_to_save_data_extended, self.all_days_of_the_year, 
                                                        self.suffix_ranges, self.info, 
                                                        self.map_files, 
                                                        self.files_have_been_processed,
                                                        save_map_plots_of_which_time_frequency)
                                                       for Year_month_week_pattern in self.Year_month_week_patterns ])
                
        # to_remove = []
        # for Year_month_week_pattern in self.Year_month_week_patterns:
        #     print(Year_month_week_pattern)
        #     to_remove.append(Process_each_week(   Year_month_week_pattern, 
        #                                                    self.where_to_save_data_extended, self.all_days_of_the_year, 
        #                                                    self.suffix_ranges, self.info, 
        #                                                    self.map_files, 
        #                                                    self.files_have_been_processed)) 
        
        # Update the processing status
        files_have_been_processed_results = pd.concat( [x[1] for x in results] )
        self.files_have_been_processed.iloc[files_have_been_processed_results.index] = files_have_been_processed_results
        
        # Save causes of non-processing to a CSV file
        os.makedirs(f'{self.where_to_save_timeseries_data}', exist_ok=True)
        self.files_have_been_processed.to_csv(self.time_series_file_name.replace('.csv', '_Causes_of_non_processing.csv'), index=False) 
        
        # Check if all weekly maps are None
        if all(x[0] is None for x in results) : 
            del files_have_been_processed_results, results
            gc.collect()
            print("Warning on creating weekly maps - All weekly maps are None")
            return
        
        # Concatenate and save time series results
        results_ts = pd.concat( [x[0] for x in results if x[0] is not None] ).sort_values(by=['day']).reset_index(drop = True)
        results_ts.to_csv(self.time_series_file_name, index=False) 
        
        # Clean up memory
        del results_ts, files_have_been_processed_results, results
        gc.collect()

    def _2_create_monthly_maps(self, nb_of_cores_to_use, save_map_plots_of_which_time_frequency) :

        """
        Create monthly maps by averaging weekly maps.

        Parameters
        ----------
        nb_of_cores_to_use : int
            Number of processor cores to use for parallel processing.
        """
        
        folder_where_to_save_maps = self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'MONTHLY')
        os.makedirs(folder_where_to_save_maps, exist_ok = True)

        # Use multiprocess to process each month
        # pool = multiprocess.Pool(nb_of_cores_to_use)
        with multiprocess.Pool(nb_of_cores_to_use) as pool:

            pool.starmap(get_the_mean_map_and_save_it, 
                        [(folder_where_to_save_maps, 
                          load_the_climatological_files(self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', "WEEKLY"), month_nb),
                          self.info,
                          f'{month_nb:02d}', f'{month_nb:02d}', 'MONTHLY', save_map_plots_of_which_time_frequency['MONTHLY'],
                          f'{self.info.Year}{month_nb:02d}15') for month_nb in (np.arange(12)+1 )])
        
    def _3_create_annual_maps(self, save_map_plots_of_which_time_frequency) :
        
        """
        Create an annual map by averaging monthly maps.
        """
        
        folder_where_to_save_maps = self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'ANNUAL')
        os.makedirs(folder_where_to_save_maps, exist_ok=True)
        
        get_the_mean_map_and_save_it(folder_where_to_save_maps,  
                                     load_the_climatological_files(self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', "MONTHLY")),
                                     self.info,
                                     period_name = 'the year', climatological_subfolder = 'ANNUAL', 
                                     do_the_plot = save_map_plots_of_which_time_frequency['ANNUAL'],
                                     date_for_plot = f'{self.info.Year}0701')
        
    def _4_create_the_multiyear_map(self) :
        
        """
        Create the multi-year map by averaging yearly maps.
        """
        
        self.where_to_save_data_extended = f'{self.where_to_save_data}/MAPS/[TIME_FREQUENCY]/'
        folder_where_to_save_maps = self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', 'MULTIYEAR')
        os.makedirs(folder_where_to_save_maps, exist_ok=True)
        
        get_the_mean_map_and_save_it(folder_where_to_save_maps,  
                                     maps_of_the_period = load_the_climatological_files(self.where_to_save_data_extended.replace('[TIME_FREQUENCY]', "ANNUAL")+"*/"),
                                     info = self.info,
                                     period_name = 'multi-years', 
                                     climatological_subfolder = 'MULTIYEAR', 
                                     do_the_plot = True) 


class QC_maps :
    
    def __init__(self, info, where_are_saved_regional_maps) :
        
        were_are_data_stored = fill_the_sat_paths(info, 
                                               path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + '/' + info.Zone).replace('[TIME_FREQUENCY]', ''),
                                               local_path = True).replace('/*/*/*', 'MAPS')
                
        # Find satellite data files for the current year
        map_files = load_the_climatological_files(where_to_save_data_extended = were_are_data_stored + '/DAILY/*/', 
                                                  return_file_names = True)
        
        multiannual_file = load_the_climatological_files(where_to_save_data_extended = were_are_data_stored + '/MULTIYEAR/', 
                                                         return_file_names = True)
        
        where_to_save_QC_data = fill_the_sat_paths(info, 
                                               path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + '/' + info.Zone).replace('[TIME_FREQUENCY]', 'QC'),
                                               local_path = True).replace('/*/*/*', '')
        os.makedirs(where_to_save_QC_data, exist_ok=True)
        
        self.QC_data = None
        self.QC_plot = None
        self.coastal_waters_mask = None
        self.info = info
        self.where_are_saved_regional_maps = where_are_saved_regional_maps
        self.map_files = map_files
        self.multiannual_file = multiannual_file
        self.were_are_data_stored = were_are_data_stored
        self.where_to_save_QC_data = where_to_save_QC_data
        
    def compute_mask_for_coastal_waters(self, minimal_bathymetry_in_m = 0, minimal_distance_from_land_in_km = 0) : 
        
        path_to_the_mask = f'{self.were_are_data_stored.split(self.info.Data_source)[0]}/{self.info.Data_source}/Coastal_mask_min_bathy_is_{minimal_bathymetry_in_m}m_min_dist_from_land_is_{minimal_distance_from_land_in_km}km_for_{self.info.atmospheric_correction}.nc'
        
        if os.path.isfile( path_to_the_mask ) :
            
            with xr.load_dataarray( path_to_the_mask ) as ds:
                self.coastal_waters_mask = ds
            
            return
        
        if not self.multiannual_file : 
            
            print("Need the MultiYear map to do the QC - To create the MultiYear map, regionalize at least one year of maps")
            return
        
        the_annual_map = load_file(self.multiannual_file[0])['Basin_map']['map_data']
                
        bathymetry_data_aligned_to_map_resolution = align_bathymetry_to_resolution(the_annual_map, 
                                                                                f'{self.where_are_saved_regional_maps}/{self.info.Zone}/Bathy_data.pkl')
        
        bathymetric_mask = bathymetry_data_aligned_to_map_resolution > -minimal_bathymetry_in_m
        
        land_mask = the_annual_map.isnull()
                
        lat_metric_resolution, lon_metric_resolution = degrees_to_km(   abs(the_annual_map.lat[0] - the_annual_map.lat[1]).item(), 
                                                                        abs(the_annual_map.lon[0] - the_annual_map.lon[1]).item(), 
                                                                        np.median(the_annual_map.lat) )
        
        # Dilate the main plume shape to extend its boundary and include neighboring pixels
        land_mask_diluted = binary_dilation(land_mask.values, 
                                            structure=np.full((math.ceil(2 * minimal_distance_from_land_in_km / lat_metric_resolution), 
                                                               math.ceil(2 * minimal_distance_from_land_in_km / lon_metric_resolution)), 
                                                              True))
        
        land_mask_diluted = xr.DataArray(land_mask_diluted, coords = land_mask.coords, dims = land_mask.dims)
                
        coastal_water_mask = bathymetric_mask | land_mask_diluted
        coastal_water_mask.to_netcdf( path_to_the_mask )
        
        fig, ax = plt.subplots(figsize=(20, 12))
        
        the_annual_map.plot(vmax = 5)
        
        ax.contourf(land_mask.lon, land_mask.lat, land_mask.values, levels=[0.5, 1], colors=['black']) 
        ax.contourf(land_mask.lon, land_mask.lat, ((land_mask_diluted | bathymetric_mask) & ~land_mask).values, levels=[0.5, 1], colors=['lightgray'], alpha = 0.5) 
        ax.set_title(f"Coastal water mask (bathy < {minimal_bathymetry_in_m}m or min dist from land < {minimal_distance_from_land_in_km}km)", fontsize = 20)
        ax.set_xlabel('')
        ax.set_ylabel('')
        
        plt.tight_layout()
        
        # Save the figure to the specified file
        plt.savefig(f'{path_to_the_mask.replace(".nc", ".png")}')
                
        # Close the figure to free up memory
        plt.close(fig)

        
        self.coastal_waters_mask = coastal_water_mask
        
    def compute_QC_metrics(self, nb_of_cores_to_use, exclude_coastal_areas = False) : 
        
        # Use multiprocess to process each week
        # pool = multiprocess.Pool(nb_of_cores_to_use)
        with multiprocess.Pool(nb_of_cores_to_use) as pool:

            QC_metrics = pool.starmap(Compute_the_metrics_of_one_map, 
                                   [( file_name, exclude_coastal_areas, self.coastal_waters_mask) 
                                                       for file_name in self.map_files ])
                
        # QC_metrics = []
        # for file_name in self.map_files :  
        #     print(file_name)
        #     QC_metrics.append(Compute_the_metrics_of_one_map( file_name, exclude_coastal_areas, self.coastal_waters_mask))
        
        # TODO: There is a bug here caused by 'date' not existing in some of the dictionaries in QC_metrics
        # So for now I've wrapped this last bit in a logic gate. We shall see if this causes problems later...
        QC_metrics_df = pd.DataFrame( [x for x in QC_metrics if x is not None] )
        if len(QC_metrics_df) > 0:
            QC_metrics_df['date'] = pd.to_datetime( QC_metrics_df['date'] )
            QC_metrics_df = QC_metrics_df.sort_values(by='date')
            self.QC_metrics = QC_metrics_df
            
            self.QC_plot = plot_and_save_the_QC_metrics(QC_df = QC_metrics_df, 
                                metrics_to_plot = ["mean_value", "99th_Percentile", "Lognorm_shape", "n_outliers"], 
                                path_to_save_QC_files = self.where_to_save_QC_data,
                                info = self.info,
                                max_cloud_cover = 80)
        
    def combine_QC_metrics(self) : 
                
        QC_files = glob.glob(f'{self.where_are_saved_regional_maps}/{self.info.Zone}/{self.info.Data_source}/{self.info.sensor_name}/{self.info.atmospheric_correction}/{self.info.Year}/MAPS/{self.info.Satellite_variable}/QC_daily_maps.csv')
        
        Global_QC_metrics = pd.concat( [pd.read_csv(file_name) for file_name in QC_files] )
        
        Global_QC_metrics['date'] = pd.to_datetime( Global_QC_metrics['date'] )
        Global_QC_metrics = Global_QC_metrics.sort_values(by='date')
        self.Global_QC_metrics = Global_QC_metrics
        
        self.Global_QC_plot = plot_and_save_the_QC_metrics(QC_df = Global_QC_metrics, 
                                                     metrics_to_plot = ["mean_value", "99th_Percentile", "Lognorm_shape", "n_outliers"], 
                                                     path_to_save_QC_files = self.where_to_save_QC_data,
                                                     info = self.info,
                                                     max_cloud_cover = 80)


# =============================================================================
#### Main functions
# =============================================================================


def create_regional_maps(core_arguments, Zones, overwrite_existing_regional_maps, 
                         save_map_plots_of_which_time_frequency, nb_of_cores_to_use,
                         where_are_saved_satellite_data, where_to_save_regional_maps) :
            
    core_arguments.update({'Years' : unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
                           'Zones' : Zones})
    
    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)
    
    for i, info in cases_to_process.iterrows() : 
                
        # info = cases_to_process.iloc[i]
        info = pd.concat([info, pd.Series([core_arguments['start_day'] if info.Year == core_arguments['Years'][0] else f'{info.Year}/01/01',
                                           core_arguments['end_day'] if info.Year == core_arguments['Years'][-1] else f'{info.Year}/12/31'],
                                           index = ['date_min', 'date_max'])])
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Satellite_variable})')
        
        maps_creation = Create_and_save_the_maps(where_to_save_regional_maps, where_are_saved_satellite_data, info) 
        
        if ('map_files' not in vars(maps_creation)) or len(maps_creation.map_files) == 0 : 
            print("Switch to the next iterate")
            continue
        
        if maps_creation.are_the_maps_already_produced and (overwrite_existing_regional_maps == False) : 
            print("Maps already exist - Switch to the next iterate")
            continue       
        
        maps_creation._1_create_weekly_maps(nb_of_cores_to_use, save_map_plots_of_which_time_frequency)
        
        maps_creation._2_create_monthly_maps(nb_of_cores_to_use, save_map_plots_of_which_time_frequency)
        
        maps_creation._3_create_annual_maps(save_map_plots_of_which_time_frequency)
        
    global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates()
        
    for i, info in global_cases_to_process.iterrows() :  
        
        # info = global_cases_to_process.iloc[i].copy()
        info['Year'] = 'MULTIYEAR'
        info = pd.concat([info, pd.Series([core_arguments['start_day'], core_arguments['end_day']], index = ['date_min', 'date_max'])])
        
        maps_creation = Create_and_save_the_maps(where_to_save_regional_maps, where_are_saved_satellite_data, info) 
        
        maps_creation._4_create_the_multiyear_map()

        
def QC_of_regional_maps(core_arguments, Zones, nb_of_cores_to_use, where_are_saved_regional_maps) : 
    
    core_arguments.update({'Years' : unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
                           'Zones' : Zones})
    
    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)
    
    for i, info in cases_to_process.iterrows() :  
                
        # info = cases_to_process.iloc[i].copy()
        
        print(f'{i} over {cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Satellite_variable})')
        
        QC_data_and_plots = QC_maps(info, where_are_saved_regional_maps)
        
        if len(QC_data_and_plots.map_files) == 0 : 
            print("No satellite regional maps - Switch to the next iterate")
            continue
                
        try : 
            QC_data_and_plots.compute_mask_for_coastal_waters(minimal_bathymetry_in_m = 1000,
                                                              minimal_distance_from_land_in_km = 20)
        except Exception : 
            print('End of the process - Switch to the next iterate')
            continue
        
        QC_data_and_plots.compute_QC_metrics(nb_of_cores_to_use, exclude_coastal_areas = True)
        
    global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates()
        
    for i in range(global_cases_to_process.shape[0]) : 
        
        info = global_cases_to_process.iloc[i].copy()
        info['Year'] = '*'
        
        # print(f'{i} over {global_cases_to_process.shape[0]-1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable})')
       
        QC_data_and_plots = QC_maps(info, where_are_saved_regional_maps)
               
        if len(QC_data_and_plots.map_files) == 0 : 
            continue
       
        try : 
            QC_data_and_plots.combine_QC_metrics()
        except Exception : 
            continue

