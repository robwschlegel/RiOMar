#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import os, sys, subprocess, glob, importlib.resources, tempfile, shutil, re
from pathlib import Path
import pandas as pd
import numpy as np
import geopandas as gpd
import netCDF4 as nc
from itertools import chain
from functools import reduce
from collections.abc import Mapping, Iterable
from panache.utils import expand_grid

proj_dir = os.path.dirname( os.path.abspath('__file__') )


# =============================================================================
#### Utility functions
# =============================================================================


def exit_program():
    print("Exiting the program...")
    sys.exit(0)


def get_registry_row(slot_key):
    """
    Look up one row of manuscript/figure_table_registry.csv by slot_key.

    Single source of truth for which output folder ("FIGURE_N"/"TABLE_N") and
    which R function currently render a given manuscript figure/table. Every
    figure.py entry point looks up its output_subdir/r_function here instead
    of hardcoding a "FIGURE_N" string, so renumbering a figure/table is a
    one-row edit to that CSV.
    """
    registry_path = os.path.join(proj_dir, 'manuscript', 'figure_table_registry.csv')
    registry = pd.read_csv(registry_path)
    row = registry[registry['slot_key'] == slot_key]
    if len(row) != 1:
        raise ValueError(f"get_registry_row(): expected exactly 1 row for slot_key "
                         f"'{slot_key}', found {len(row)}")
    return row.iloc[0]


def registry_basename(output_subdir):
    """"FIGURE_S1" -> "Figure_S1", matching every figure's filename
    convention (Figure_1.png, Figure_S1.png, ...) from its registry
    output_subdir. Use registry_filename() below when the extension is
    needed too."""
    return output_subdir.replace('FIGURE_', 'Figure_', 1)


def registry_filename(output_subdir, ext='png'):
    """"FIGURE_S1" -> "Figure_S1.png" -- see registry_basename() above."""
    return f"{registry_basename(output_subdir)}.{ext}"


def degrees_to_km(lat_deg, lon_deg, latitude):
    """
    Convert distances in degrees to kilometers.
    
    Parameters:
    - lat_deg: Distance in degrees of latitude
    - lon_deg: Distance in degrees of longitude
    - latitude: Latitude where the conversion is needed
    
    Returns:
    - lat_km: Distance in kilometers for latitude
    - lon_km: Distance in kilometers for longitude
    """
    # Conversion factor for latitude
    lat_km = lat_deg * 111.32
    
    # Conversion factor for longitude, adjusted by latitude
    lon_km = lon_deg * 111.32 * np.cos(np.radians(latitude))
    
    return lat_km, lon_km


def km_to_degrees(lat_km, lon_km, latitude):
    """
    Convert a distance in meters to degrees of latitude and longitude.

    Args:
        meters (float): Distance in meters.
        latitude (float): Latitude in degrees where the conversion is applied.

    Returns:
        (float, float): (Latitude degrees, Longitude degrees)
    """
    # 1 degree latitude ≈ 111.32 km (constant)
    lat_deg = lat_km / 111.320  

    # 1 degree longitude ≈ 111.32 * cos(latitude) km
    lon_deg = lon_km / (111.320 * np.cos(np.radians(latitude)))

    return lat_deg, lon_deg


def generic_variable_names() : 
    
    return ['CHLA', 'SPM', 'SST']


def find_sat_data_files(info, path_to_sat_data) : 

    if isinstance(path_to_sat_data, str) : 
        path_to_sat_data = [path_to_sat_data]

    if ('Year' in info) and isinstance(info['Year'], str) and (info['Year'] == 'MULTIYEAR') : 
        file_pattern = '/*.pkl'
    elif np.isin(info.Satellite_variable, generic_variable_names()) : 
        file_pattern = '/*.nc'
    else : 
        file_pattern = f'/*{info.Satellite_variable}*.nc'    

    map_files = []
    
    if '[YEAR]' in path_to_sat_data[0] : 
        
        for Year in info['Year'] :    
    
            # path_to_files = (path_to_sat_data 
            #                      .replace('[YEAR]', str(Year))
            #                      .replace('[MONTH]', '*')
            #                      .replace('[DAY]', '*'))
                              
            files = glob.glob(path_to_sat_data + file_pattern)
            
            map_files.extend( files )
            
    else : 
        
        map_files = list(chain.from_iterable(glob.glob(path + file_pattern) for path in path_to_sat_data))
        
    return map_files


def path_to_fill_to_where_to_save_satellite_files(where_to_save_files) :
    
    path = f'{where_to_save_files}/[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]/[YEAR]/[MONTH]/[DAY]'

    return path


def fill_the_sat_paths(info, path_to_fill, local_path = False, dates = []) : 
    
    if len(dates) > 0 and isinstance(dates[0], str) : 
        dates = pd.to_datetime(dates)
    
    path_to_fill = ( path_to_fill  
                        .replace('[DATA_SOURCE]', info.Data_source)
                        .replace('[ATMOSPHERIC_CORRECTION]', info.atmospheric_correction)
                        .replace('[SENSOR]', info.sensor_name if 'sensor_name' in info else info.Satellite_sensor) )
                      
    if 'Temporal_resolution' in info.keys() or 'Temporal_resolution' in info.keys() : 
        
        Temporal_resolution = info.Temporal_resolution if 'Temporal_resolution' in info.keys() else info.Temporal_resolution
        if local_path == False : 
            Temporal_resolution = (Temporal_resolution
                                    .replace('DAILY', 'day')
                                    .replace('MONTHLY', 'month')
                                    .replace('WEEKLY', '8-day'))        
            
    elif isinstance(info.Year, str) and (info.Year == 'MULTIYEAR') :
        
        Temporal_resolution = 'ANNUAL'
        
    else :
        
        Temporal_resolution = 'DAILY'
            
    path_to_fill = path_to_fill.replace('[TIME_FREQUENCY]', Temporal_resolution)
        
    if local_path : 
        
        Folder_name_for_the_variable = ('CHLA' if 'CHL' in info.Satellite_variable.upper()
                                        else 'SPM' if 'SPM' in info.Satellite_variable.upper()
                                        else 'SST' if 'SST' in info.Satellite_variable.upper()
                                        else info.Satellite_variable)
        
        path_to_fill = (path_to_fill.replace('[PARAMETER]', Folder_name_for_the_variable))
        
    else : 
        
        path_to_fill = (path_to_fill.replace('[PARAMETER]', info.Satellite_variable_name_on_remote_folder))

    if len(dates) > 0 : 
        
        paths_to_sat_files = [ ( path_to_fill
                                      .replace('[YEAR]', str(date.year))
                                      .replace('[MONTH]', str(date.month).zfill(2))
                                      .replace('[DAY]', str(date.day).zfill(2))
                                      .replace("[DOY]", date.strftime("%j")) ) 
                                for date in dates ]
        
    else : 
        
        paths_to_sat_files = (path_to_fill
                                .replace('[YEAR]', '*')
                                .replace('[MONTH]', '*')
                                .replace('[DAY]', '*')
                                .replace("[DOY]", '*')) 
    
    return paths_to_sat_files


def get_all_cases_to_process(core_arguments) : 
        
    cases_to_process = expand_grid(Data_source = core_arguments['Data_sources'], 
                                  sensor_name = core_arguments['Sensor_names'], 
                                  atmospheric_correction = core_arguments['Atmospheric_corrections'],
                                  Satellite_variable = core_arguments['Satellite_variables'],
                                  Temporal_resolution = core_arguments['Temporal_resolution'] if 'Temporal_resolution' in core_arguments else ['DAILY'])
    
    cases_to_process['atmospheric_correction'] = cases_to_process.apply(lambda row: 'Standard' 
                                                                        if np.isin(row['Data_source'], ['SEXTANT', 'EUMETSAT']) 
                                                                        else row['atmospheric_correction'], axis=1)
    cases_to_process['Satellite_variable'] = cases_to_process.apply(lambda row: 'SPM' 
                                                                        if row['Data_source'] == 'SEXTANT' and 'SPM' in row['Satellite_variable'] 
                                                                        else row['Satellite_variable'], axis=1)
    
    cases_to_process = cases_to_process.drop_duplicates().reset_index(drop = True)  
    
    return cases_to_process


def get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments) : 

    all_possibilities = expand_grid(Zone = core_arguments['Zones'],
                                    Data_source = core_arguments['Data_sources'], 
                                    sensor_name = core_arguments['Sensor_names'], 
                                    atmospheric_correction = core_arguments['Atmospheric_corrections'],
                                    Year = core_arguments['Years'],
                                    Satellite_variable = core_arguments['Satellite_variables'],
                                    Temporal_resolution = (core_arguments['Temporal_resolution'] if 'Temporal_resolution' in core_arguments 
                                                       else core_arguments['Temporal_resolution'] if 'Temporal_resolution' in core_arguments 
                                                       else ''))
    all_possibilities['atmospheric_correction'] = all_possibilities.apply(lambda row: 'Standard' 
                                                                        if row['Data_source'] == 'SEXTANT' 
                                                                        else row['atmospheric_correction'], axis=1)
    all_possibilities = all_possibilities.drop_duplicates()
    
    return all_possibilities


def create_arborescence(paths):
    arborescence = {}
    for path in paths:
        keys = path.split("/")
        current_level = arborescence
        for key in keys:
            if key not in current_level:
                current_level[key] = {}
            current_level = current_level[key]  # Move deeper
    return arborescence


def parameter_name_from_file_name(file_name) : 
                        
    regular_expression = r'(?:(SPM-[G|R]|SPIM|suspended_matters|TSM_NN|CHL|CHL1|CHL-OC5|CHL-GONS|chlorophyll_a|POC|NRRS[0-9]*|RRS[0-9]*|DOC|CDOM|BBP|T-FNU|SST(?:-NIGHT|)))'
    match = re.search(regular_expression, file_name)
    
    if match : 
        
        return match.group(0) 
            
    else : 
        
        print('!!! Impossible to find the parameter name from the file name (see the function parameter_name_from_file_name) !!!')                        


def add_array_to_dict(dictionary, path, array):
    """
    Adds an array to a specific position in a nested dictionary.
    
    Parameters:
        dictionary (dict): The main dictionary to update.
        path (str): The path in the form "A/B/C/D/filename.nc".
        array (numpy.ndarray): The array to be stored.
    
    Returns:
        None (modifies the dictionary in place).
    """
    keys = path.split("/")
    filename = keys[-1]  # Extract the filename
    
    # Extract parameter name (e.g., SPM-G)
    param = parameter_name_from_file_name(filename)

    # Navigate the dictionary hierarchy
    current_level = dictionary
    for key in keys[:-1]:  # Exclude filename from navigation
        if key not in current_level:
            current_level[key] = {}
        current_level = current_level[key]

    # Add the array under the parameter name
    current_level[param] = array


def get_nested_dict_value(dictionary, path):
    
    keys = path.split("/")

    item = reduce(lambda dictionary, key: dictionary[key], keys[1:], dictionary)
    
    return item
        
    
def merge_dicts(dicts):
    """
    Merges multiple nested dictionaries into a single dictionary.
    
    Parameters:
        dicts (list): List of dictionaries to merge.
    
    Returns:
        dict: A merged dictionary.
    """
    def recursive_merge(d1, d2):
        """Recursively merges two dictionaries."""
        for key, value in d2.items():
            if key in d1 and isinstance(d1[key], dict) and isinstance(value, dict):
                recursive_merge(d1[key], value)  # Merge sub-dictionaries
            else:
                d1[key] = value  # Overwrite or add new keys
        return d1

    merged_dict = {}
    for d in dicts:
        merged_dict = recursive_merge(merged_dict, d)

    return merged_dict


def get_empty_paths(dictionary, prefix=""):
    
    paths = []
    
    if isinstance(dictionary, dict) and len(dictionary) > 0:  # If it's a non-empty dictionary
        for key, val in dictionary.items():
            paths.extend(get_empty_paths(val, f"{prefix}/{key}"))  # Recursive call
    elif len(dictionary) == 0:  # Exclude NaN values
        paths.append(prefix)  # Add the path if it's a valid value
        # paths.append( [f'{prefix}/{x}' for x in list(dictionary.keys())] )  # Add the path if it's a valid value
    
    return paths


def get_non_empty_paths(dictionary, prefix=""):
    
    paths = []
    
    if isinstance(dictionary, dict) and len(dictionary) > 0 :  # If it's a non-empty dictionary
        for key, val in dictionary.items():
            paths.extend(get_non_empty_paths(val, prefix = f"{prefix}/{key}"))  # Recursive call

    elif isinstance(dictionary, dict) == False:  # Exclude NaN values
        paths.append(prefix)  # Add the path if it's a valid value
    else :  # Exclude NaN values
        for key, val in dictionary.items():
            paths.append(f"{prefix}")  # Recursive call
            break
    
    # paths = np.unique( [x.replace('Sat_values', '').replace('Time', '') for x in paths] )
    
    return paths


def extract_nested_dataframes(data):
    """Efficiently extract all DataFrames from a nested dictionary using an iterative approach."""
    stack = [data]  # Use a stack to avoid deep recursion

    while stack:
        current = stack.pop()

        if isinstance(current, pd.DataFrame):
            yield current  # Yield instead of appending to a list (memory-efficient)
        elif isinstance(current, Mapping):  # Check if it's a dictionary
            stack.extend(current.values())  # Add dictionary values to the stack
        elif isinstance(current, Iterable) and not isinstance(current, (str, bytes)):  
            stack.extend(current)  # Add list/tuple elements to the stack


def load_shapefile_data() :
    
        shp_folder = os.path.join( proj_dir, 'data', 'FRANCE_shapefile' )  # Directly get the package folder path
    
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Extract all necessary shapefile components
            for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg']:
                # shp_file = shp_folder / f'gadm41_FRA_0{ext}'
                shp_file = os.path.join(shp_folder, f'gadm41_FRA_0{ext}')
                if os.path.exists(shp_file): # Ensure the file exists before copying
                    shutil.copy(shp_file, os.path.join(tmp_dir, f'gadm41_FRA_0{ext}'))
    
            # Read the shapefile from the temporary directory
            shapefile_path = os.path.join(tmp_dir, 'gadm41_FRA_0.shp')
            return gpd.read_file(shapefile_path)
            
    # try : 
    #     france_shapefile = pygadm.Items(name="FRANCE", content_level=0)
    #     return france_shapefile
    # except Exception as e :
    #     print(f"The France shapefile can't be accessed through pygadm : {e}")
    #     print("The France shapefiles can be manually downloaded for free : e.g. https://gadm.org/download_country.html ")
 
    
def date_from_path(path):
    match = re.search(r'/(\d{4})/(\d{2})/(\d{2})/', path)
    return ''.join(match.groups()) if match else None   


def load_csv_files(SOMLIT = False, REPHY = False,
                   RIVER_FLOW = False, Zone_of_river_flow = None,
                   RIVER_FLOW_time_resolution = '', river_files = None):
    
    if SOMLIT : 
        SOMLIT_dir = os.path.join(proj_dir, 'data', 'INSITU_data', 'SOMLIT')
        SOMLIT_data = os.path.join(SOMLIT_dir, 'Somlit.csv')
        return (pd.read_csv(SOMLIT_data, sep = ";", header = 2).iloc[1:]
                                .rename(columns = {'gpsLat*':'LATITUDE', 
                                                   'gpsLong*':'LONGITUDE',
                                                   'nomSite*':"Site"}))
        
    if REPHY : 
        REPHY_dir = os.path.join(proj_dir, 'data', 'INSITU_data', 'REPHY')
        REPHY_data = os.path.join(REPHY_dir, 'Table1_REPHY_hydro_RIOMAR.csv.gz')
        return pd.read_csv(REPHY_data, sep = ";", header = 0, encoding = "ISO-8859-1", compression = {'method' : 'gzip'})
        
    if RIVER_FLOW : 
        
        where_are_river_data = os.path.join(proj_dir, 'data', 'RIVER_FLOW', Zone_of_river_flow)

        # river_files restricts the sum to specific rivers (e.g. one
        # individual river mouth, or a small basin-sharing group like
        # Gironde = Garonne + Dordogne) instead of every CSV in the zone
        # directory -- see metadata/river_discharge_mapping.csv.
        if river_files is not None:
            files_to_load = [os.path.join(where_are_river_data, f"{slug}.csv") for slug in river_files]
        else:
            files_to_load = glob.glob(os.path.join(where_are_river_data, '*.csv'))

        # Convert file paths to Path objects for consistent handling
        files_to_read = [Path(f) for f in files_to_load if Path(f).suffix in ('.csv')]

        # Load each file into a dictionary of DataFrames
        data_dict = {}
        for file in files_to_read:
            df = pd.read_csv(file)
            data_dict[file.name] = df  # Store in dictionary
                    
        for key, the_df in data_dict.items() :
            # Select by name rather than blindly overwriting all column
            # names positionally -- func/river_flow_prep.R's output now has
            # a third 'flag' column (data provenance: measured/estimated),
            # so `the_df.columns = ['date', 'debit']` raised a length-mismatch
            # ValueError the moment a 3-column file was read (fixed 2026-08-02).
            the_df['date'] = pd.to_datetime(the_df['date'], format = "%Y-%m-%d")
            data_dict[key] = the_df.loc[:,['date', 'debit']]
    
        # n_rivers must count ROWS per date ('size'), not non-NaN values
        # ('count') -- 'count' silently dropped the whole zone's flow on any
        # date where even one tributary's debit was NaN (fixed 2026-08-10,
        # found via the Bay of Seine: the Orne is 44.8% NaN, which was
        # discarding valid Seine flow on the same dates and showing up as
        # large gaps in Figure 6's river-flow line). 'sum' already skips NaN
        # by default, so this now matches func/multi.R::load_river_flow()'s
        # dplyr::n()/sum(na.rm=TRUE) semantics: a zone's flow is the sum of
        # whichever of its rivers reported a numeric value that day.
        final_df = pd.concat(data_dict.values()).groupby("date", as_index=False).agg(Values=('debit', 'sum'), n_rivers=('debit', 'size'))
        final_df = final_df[final_df.n_rivers == len(files_to_read)]
        
        bin_centers = [4, 12, 20, 28]
        def assign_bin(day):
            return min(bin_centers, key=lambda x: abs(x - day))
        
        # Apply the function to create a 'bin' column
        final_df_reduced = final_df.copy()
        final_df_reduced['bin'] = final_df_reduced['date'].dt.day.apply(assign_bin)
        final_df_reduced = final_df_reduced.loc[:,['date', 'bin', 'Values']]
        final_df_reduced['Values'] = pd.to_numeric(final_df_reduced['Values'])
        
        if RIVER_FLOW_time_resolution == 'DAILY' :
            final_df = final_df

        elif RIVER_FLOW_time_resolution == 'WEEKLY' : 
            final_df_binned = final_df_reduced.groupby([final_df_reduced['date'].dt.to_period('M'), 'bin']).agg({'Values': 'mean'}).reset_index()
            final_df_binned['date'] = pd.to_datetime( final_df_binned['date'].astype(str) + "-" + final_df_binned['bin'].astype(str), format = "%Y-%m-%d" )
            final_df = final_df_binned
            
        elif RIVER_FLOW_time_resolution == 'MONTHLY' :     
            final_df_binned = final_df_reduced.groupby([final_df_reduced['date'].dt.to_period('M')]).agg({'Values': 'mean'}).reset_index()
            final_df_binned['date'] = pd.to_datetime( final_df_binned['date'].astype(str) + "-" + "15", format = "%Y-%m-%d" )
            final_df = final_df_binned

        elif RIVER_FLOW_time_resolution == 'ANNUAL' :     
            final_df_binned = final_df_reduced.groupby([final_df_reduced['date'].dt.to_period('Y')]).agg({'Values': 'mean'}).reset_index()
            final_df_binned['date'] = pd.to_datetime( final_df_binned['date'].astype(str) + "-06-30", format = "%Y-%m-%d" )
            final_df = final_df_binned
            
        else : 
            print("The RIVER_FLOW_time_resolution must be one of the following: 'DAILY', 'WEEKLY', 'MONTHLY', 'ANNUAL'")
            exit_program()

        return final_df
         
        
def flatten_and_dedupe(lst_comprehension, return_the_unique_values) : 
    
    the_values = list(chain(*lst_comprehension))
    
    if return_the_unique_values : 
        the_values = np.unique(the_values)
    
    return the_values


ZONE_DISPLAY_NAMES = {
    'BAY_OF_SEINE': 'Bay of Seine',
    'BAY_OF_BISCAY': 'Bay of Biscay',
    'SOUTHERN_BRITTANY': 'Southern Brittany',
    'GULF_OF_LION': 'Gulf of Lion',
}


def zone_title(zone_name):
    """
    Human-facing zone display name (e.g. "BAY_OF_SEINE" -> "Bay of Seine").

    The single shared source of this mapping across the Python side of the
    codebase -- use this instead of hand-writing zone titles or leaving zone
    codes in ALL_CAPS on any human-facing plot/table/caption. See
    func/multi.R::zone_title() for the R-side equivalent.
    """
    return ZONE_DISPLAY_NAMES[zone_name]


# Canonical zone order, north to south, for consistent panel ordering across
# every multi-zone figure/table (top to bottom, or left to right): Bay of
# Seine, Southern Brittany, Bay of Biscay, Gulf of Lion. See
# func/util.R::ZONE_ORDER/order_zones() for the R-side equivalent.
ZONE_ORDER = ['BAY_OF_SEINE', 'SOUTHERN_BRITTANY', 'BAY_OF_BISCAY', 'GULF_OF_LION']


def order_zones(zone_iterable):
    """
    Sort any iterable of zone codes into the canonical north-to-south order.

    Use this instead of a plain sorted()/dict key order, which would
    otherwise come out alphabetical or reflect whatever order zones happened
    to be listed in at each call site.
    """
    return sorted(zone_iterable, key=ZONE_ORDER.index)


# The one canonical named river per zone, used only to pick which entry of
# panache's per-river starting_points dict export_panache_zone_metadata()
# reads -- func/multi.R::zone_meta hardcodes this same one-river-per-zone
# assumption (a case_when on mouth_name), so this must stay one row per zone,
# not one row per river (e.g. Gulf of Lion also has Petit Rhone in panache).
CANONICAL_RIVER_PER_ZONE = {
    'BAY_OF_SEINE': 'Seine',
    'BAY_OF_BISCAY': 'Gironde',
    'SOUTHERN_BRITTANY': 'Loire',
    'GULF_OF_LION': 'Grand Rhone',
}


def export_panache_zone_metadata(output_path=None):
    """
    Regenerate metadata/panache_zone_metadata.csv from panache.utils.define_parameters(),
    the single source read by func/util.R to build zones_bbox/river_mouths. Rerun this
    whenever panache's zone parameters change (see ~/panache/NEWS.md).
    """

    from panache.utils import define_parameters

    output_path = output_path or os.path.join(proj_dir, 'metadata', 'panache_zone_metadata.csv')

    rows = []
    for zone, mouth_name in CANONICAL_RIVER_PER_ZONE.items():
        parameters = define_parameters(zone)
        mouth_lat, mouth_lon = parameters['starting_points'][mouth_name]
        lat_min, lat_max = parameters['lat_range_of_plume_area']
        lon_min, lon_max = parameters['lon_range_of_plume_area']
        rows.append({'zone': zone, 'mouth_name': mouth_name,
                     'mouth_lon': mouth_lon, 'mouth_lat': mouth_lat,
                     'lon_min': lon_min, 'lon_max': lon_max,
                     'lat_min': lat_min, 'lat_max': lat_max})

    pd.DataFrame(rows).to_csv(output_path, index=False)


def daily_integral(file_dir, overwrite=False):
    
    """
    Creates daily integral NetCDF files from their hourly versions.
    Remember that this requires the first hour of the next day to compute the daily integral.
    file_dir : str
        Directory where the input files are located.
    overwrite : bool, optional
        Whether to overwrite existing output files. Default is False.    
    """

    # Get list of files in the directory
    dir_files = os.listdir(file_dir)

    # Split out files with 'daily' in their names
    daily_files = [f for f in dir_files if 'daily' in f]
    hourly_files = [f for f in dir_files if 'daily' not in f]

    # Remove all_files with the existing daily files if overwrite is False
    if not overwrite :
        hourly_files_check = [f for f in hourly_files if f.split('_')[0] + '_daily_' + f.split('_')[1] + '_' + f.split('_')[2] not in daily_files]
    else : 
        hourly_files_check = hourly_files
        # Remove any existing daily files if overwrite is True
        for f in daily_files:
            os.remove(os.path.join(file_dir, f))

    # Stop if no hourly files need to be processed
    if len(hourly_files_check) == 0:
        print("No new hourly files to process.")
        return

    # Process each file to create daily integral versions
    for hourly_file in hourly_files_check:

        daily_file = hourly_file.split('_')[0] + '_daily_' + hourly_file.split('_')[1] + '_' + hourly_file.split('_')[2]
        hourly_path = os.path.join(file_dir, hourly_file)
        daily_path = os.path.join(file_dir, daily_file)

        with nc.Dataset(hourly_path) as ds:
            has_wave_dir = 'VMDR' in ds.variables

        if not has_wave_dir:
            # Calculate daily means
            subprocess.run(['cdo', 'daymean', hourly_path, daily_path], check=True)
            continue

        # VMDR (mean wave direction) is a raw angle, not a vector component --
        # a plain daymean is biased on any day where the true mean direction
        # crosses the 0/360 wrap. Daymean every other variable as usual, but
        # average VMDR as a circular mean via its sin/cos components, then
        # merge the reconstructed angle back in.
        with tempfile.TemporaryDirectory() as tmp_dir:
            other_vars_file = os.path.join(tmp_dir, 'other_vars.nc')
            components_file = os.path.join(tmp_dir, 'dir_components.nc')
            vmdr_file = os.path.join(tmp_dir, 'vmdr.nc')

            subprocess.run(['cdo', 'daymean', '-delname,VMDR', hourly_path, other_vars_file], check=True)
            subprocess.run(['cdo', 'daymean', '-expr,wave_dir_sin=sin(rad(VMDR));wave_dir_cos=cos(rad(VMDR))',
                            hourly_path, components_file], check=True)
            subprocess.run(['cdo', 'setattribute,VMDR@units=degree',
                            '-expr,VMDR=mod(deg(atan2(wave_dir_sin,wave_dir_cos))+360,360)',
                            components_file, vmdr_file], check=True)
            subprocess.run(['cdo', 'merge', other_vars_file, vmdr_file, daily_path], check=True)

