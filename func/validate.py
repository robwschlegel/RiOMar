#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import re, datetime, os, sys, pickle, gc, glob, multiprocess
import pandas as pd
import numpy as np
import xarray as xr
import rpy2.robjects as robjects
from joblib import dump, load

multiprocess.set_start_method('spawn', force = True) # MacOS friendly multiprocessing

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (add_array_to_dict, path_to_fill_to_where_to_save_satellite_files, fill_the_sat_paths,  # noqa: E402
                    create_arborescence, find_sat_data_files,merge_dicts, get_empty_paths,
                    return_the_parameter_name_based_on_file_name, get_non_empty_paths,
                    extract_the_time_from_the_satellite_file, access_item_in_a_dictionnary, 
                    extract_dataframes_iterative, load_csv_files_in_the_package_folder,
                    get_the_values_from_a_list_comprehension, define_parameters, get_all_cases_to_process)


# =============================================================================
#### Utility functions
# =============================================================================


def get_insitu_measurements(zones = ['FRANCE']) : 
        
    # =============================================================================
    #     SOMLIT
    # =============================================================================
        
    QC_values_to_keep = [2,6,7]
    Depth_threshold = 10    

    SOMLIT_data = load_csv_files_in_the_package_folder(SOMLIT = True)
    SOMLIT_data['SOURCE'] = 'SOMLIT'
    SOMLIT_data["TIME"] = pd.to_datetime(SOMLIT_data["DATE"] + " " + SOMLIT_data["HEURE"], format="%Y-%m-%d %H:%M:%S")

    # Keep data with good QC
    QC_columns = SOMLIT_data.filter(regex='^q').columns.tolist()
    var_names = [x[1:] for x in QC_columns]
    
    for var_name in var_names : 
        
        mask = ~SOMLIT_data[f'q{var_name}'].isin(QC_values_to_keep)
        
        SOMLIT_data.loc[mask, var_name] = np.nan
    
    index_to_keep = np.where( SOMLIT_data['PROF_NUM'].astype(float) <= Depth_threshold )[0]
    SOMLIT_data_filtered = (SOMLIT_data
                            .iloc[index_to_keep]
                            .rename(columns = {'T':'TEMP', 
                                               'COP':'POC',
                                               'MES': 'SPM',
                                               'S': 'SAL',
                                               'Site': 'SITE'})
                            .loc[:,['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE', 'TIME', 'TEMP', 'SAL', 'CHLA', 'POC', 'SPM']])
    
    region_mapping = {
        'Point B': 'Mer ligurienne - Corse',
        'Frioul': 'Golfe du Lion', 
        'Sete': 'Golfe du Lion', 
        'Sola': 'Golfe du Lion',
        'Comprian': 'Sud Golfe de Gascogne', 
        'Eyrac': 'Sud Golfe de Gascogne', 
        'Bouee 13': 'Sud Golfe de Gascogne',
        'pk 30': 'Pays de la Loire - Pertuis', 
        'pk 52': 'Pays de la Loire - Pertuis', 
        'pk 86': 'Pays de la Loire - Pertuis', 
        'Antioche': 'Pays de la Loire - Pertuis',
        'Portzic': 'Bretagne Sud',
        'Estacade': 'Manche occidentale', 
        'Astan': 'Manche occidentale', 
        'Bizeux': 'Manche occidentale', 
        'Le Buron': 'Manche occidentale', 
        'Cézembre': 'Manche occidentale',
        'Smile': 'Baie de Seine', 
        'Luc-sur-Mer': 'Baie de Seine',
        'Point C': 'Manche orientale - Mer du Nord', 
        'Point L': 'Manche orientale - Mer du Nord'
    }
    
    zone_mapping = {
        'Mer ligurienne - Corse' : 'GULF OF LION',
        'Golfe du Lion' : 'GULF OF LION',
        'Sud Golfe de Gascogne' : 'GULF OF BISCAY',         
        'Pays de la Loire - Pertuis': 'GULF OF BISCAY',
        'Bretagne Sud' : 'SOUTHERN BRITTANY',
        'Manche occidentale' : 'SOUTHERN BRITTANY', 
        'Baie de Seine' : 'BAY OF SEINE', 
        'Manche orientale - Mer du Nord' : 'BAY OF SEINE'
    }
    
    SOMLIT_data_filtered['REGION'] = SOMLIT_data_filtered['SITE'].map(region_mapping).fillna('Region to fill').map(zone_mapping).fillna('Region to fill')
    
    # SOMLIT_stations = SOMLIT_data_filtered.loc[:,["ID_SITE", "Site", "DATE", "HEURE", "Latitude", "Longitude"]].drop_duplicates().reset_index(drop = True)

    # =============================================================================
    #     REPHY
    # =============================================================================

    var_names = ['CHLOROA', 'TURB', 'TURB-FNU', 'TEMP', 'SALI']
    QC_values_to_keep = ['Bon']
    Depth_threshold = 10

    REPHY_data = load_csv_files_in_the_package_folder(REPHY = True)         
    REPHY_data["TIME"] = pd.to_datetime(REPHY_data["Date"] + " " + REPHY_data["Heure"], format="%Y-%m-%d %H:%M:%S")
    REPHY_data["Valeur_mesure"] = REPHY_data["Valeur_mesure"].str.replace(",", ".").astype(float)
    REPHY_data['SOURCE'] = 'REPHY'
    
    index_to_keep = np.where( np.isin( REPHY_data['Code.parametre'], var_names ) &
                              np.isin( REPHY_data['Qualite.resultat'], QC_values_to_keep ) &
                              (REPHY_data['Profondeur.metre'] <= Depth_threshold) )[0]

    REPHY_data_filtered = ( REPHY_data
                               .iloc[index_to_keep]
                               .rename(columns = {'lat':'LATITUDE', 
                                                  'lon':'LONGITUDE',
                                                  'Region': 'REGION',
                                                  'Code_point_Libelle':"SITE",
                                                  'Qualite.resultat':'QC',
                                                  'Code.parametre':'VARIABLE',
                                                  'Valeur_mesure':'VALUE'})
                               .loc[:,['SOURCE', 'SITE', "REGION", 'LATITUDE', 'LONGITUDE', 'TIME', 'VARIABLE', 'VALUE', 'QC']] 
                               .pivot_table(index=['SOURCE', 'SITE', "REGION", 'LATITUDE', 'LONGITUDE', 'TIME'], 
                                            columns="VARIABLE", values="VALUE", aggfunc="mean")
                               .rename(columns = {'CHLOROA':'CHLA', 
                                                  'SALI':'SAL'})
                               .reset_index(drop = False))
    
    REPHY_data_filtered['TURB_all'] = REPHY_data_filtered[["TURB-FNU", "TURB"]].mean(axis=1)
    REPHY_data_filtered['TURB'] = REPHY_data_filtered['TURB']
    REPHY_data_filtered = REPHY_data_filtered.drop(['TURB_all', 'TURB-FNU'], axis = 1)
        
    REPHY_data_filtered['REGION'] = REPHY_data_filtered['REGION'].map(zone_mapping).fillna('Region to fill')

    # =============================================================================
    #  MERGE ALL
    # =============================================================================
    
    INSITU_measurements = pd.concat([SOMLIT_data_filtered, REPHY_data_filtered], ignore_index = True)
    INSITU_measurements['DATE'] = INSITU_measurements['TIME'].dt.date.astype(str)
    INSITU_measurements['TIME'] = INSITU_measurements['TIME'].dt.time.astype(str)
    
    BGC_columns = ['TEMP', 'SAL', 'POC', 'SPM', 'CHLA', 'TURB']
    INSITU_measurements[BGC_columns] = INSITU_measurements[BGC_columns].astype(float)
    
    # INSITU_stations = INSITU_measurements.loc[:,["SOURCE", "SITE", "LATITUDE", "LONGITUDE"]].drop_duplicates().reset_index(drop = True)

    if len(zones) > 1 or zones != ["FRANCE"] : 
        
        coordinates_of_the_zones = {zone: define_parameters(zone)['lat_range_of_the_map_to_plot'] + 
                                   define_parameters(zone)['lon_range_of_the_map_to_plot'] 
                            for zone in zones}
    
        measurement_LATITUDES = INSITU_measurements.LATITUDE.to_numpy(float)
        measurement_LONGITUDES = INSITU_measurements.LONGITUDE.to_numpy(float)
        
        index_of_measurements_in_the_zones = np.hstack([np.where((measurement_LATITUDES >= lat_min) & 
                                                                 (measurement_LATITUDES <= lat_max) & 
                                                                 (measurement_LONGITUDES >= lon_min) & 
                                                                 (measurement_LONGITUDES <= lon_max))[0] 
                                                        for lat_min, lat_max, lon_min, lon_max in coordinates_of_the_zones.values()])
                
        INSITU_measurements = INSITU_measurements.iloc[ np.unique(index_of_measurements_in_the_zones) ]
                
    INSITU_stations = (INSITU_measurements
                       .groupby(["SOURCE", "SITE", "REGION", "LATITUDE", "LONGITUDE"], as_index = False)
                       .agg({"DATE": list}))

    return INSITU_measurements, INSITU_stations


def find_indices_of_the_grid(lat_station, lon_station, lat_map, lon_map, grid_size):
    
    half_grid_size = int( (grid_size - 1) / 2 )
    
    lat_idx = np.abs(lat_map - lat_station).argmin()
    lon_idx = np.abs(lon_map - lon_station).argmin()
    
    lat_idx_all = np.arange(lat_idx-half_grid_size, lat_idx+half_grid_size+1)
    lon_idx_all = np.arange(lon_idx-half_grid_size, lon_idx+half_grid_size+1)
    
    return {'lat_index' : lat_idx_all, 'lon_index' : lon_idx_all}


def do_the_match_up_for_one_date(satellite_date, satellite_files, MU_data, MU_database_of_the_case, 
                                 MU_stations, info,
                                  grid_size, where_are_saved_satellite_data) : 

    date_pattern = re.compile(satellite_date)
    
    satellite_files_of_the_day = [file for file in satellite_files if date_pattern.search(file)]
    
    satellite_date_formatted = pd.to_datetime(satellite_date).strftime('%Y-%m-%d')
        
    for satellite_file in satellite_files_of_the_day : 
    
        try : 
            with xr.open_dataset(satellite_file, decode_times=True) as f : 
                map_ini = f
        except Exception as e:
            
            add_array_to_dict(dictionary = MU_database_of_the_case,
                              path = satellite_file.replace(f'{where_are_saved_satellite_data}/', ''),
                              array = f"❌ Impossible to load the nc file {satellite_file} - Error : {e}")
            
            continue
        
        if ( (len(map_ini) == 0) or 
             (len(np.unique(map_ini.lat)) == 1) or 
             (len(np.unique(map_ini.lon)) == 1) ) : 
            
            add_array_to_dict(dictionary = MU_database_of_the_case,
                              path = satellite_file.replace(f'{where_are_saved_satellite_data}/', ''),
                              array = f"❌ Wrong file format for {satellite_file}")
            # return 
            continue
        
        the_var_name = next(iter(map_ini.keys()))
        
        time_value = extract_the_time_from_the_satellite_file(map_ini)
        
        MU_values = [
            map_ini.isel(lat=station_info['lat_index'], lon=station_info['lon_index'])[the_var_name].values.flatten()
            if satellite_date_formatted in station_info['DATE'] else np.array(["No in-situ measurement"])
            for _, station_info in MU_data['STATIONS'].iterrows()
        ]
                
        MU_values = pd.DataFrame(MU_values)
        
        if MU_values.shape[1] > 1 : 
                   
            MU_values.columns = [f"{i}_{j}" for i in range(1, grid_size+1) for j in range(1, grid_size+1)]     
            MU_values.index = pd.MultiIndex.from_frame(
                                    MU_stations.assign(**pd.concat([info,
                                                                    pd.Series({'Satellite_algorithm' : return_the_parameter_name_based_on_file_name(satellite_file),
                                                                               'DATE' : satellite_date_formatted,
                                                                               'Satellite_time' : time_value})])))
                         
        MU_values = MU_values[ MU_values.iloc[:,0] != "No in-situ measurement" ]
        
        if MU_values.shape[0] > 0 :
            
            add_array_to_dict(dictionary = MU_database_of_the_case,
                              path = satellite_file.replace(f'{where_are_saved_satellite_data}/', ''),
                              array = MU_values)
        
        del MU_values, map_ini
        gc.collect()
            
    return MU_database_of_the_case


def find_the_path_to_satellite_MU_in_the_dict(path_to_satellite_file, where_are_saved_satellite_data) : 
    
    filename = path_to_satellite_file.split("/")[-1]
    
    # Extract parameter name (e.g., SPM-G)
    param = return_the_parameter_name_based_on_file_name(filename)
    path_in_the_dict = path_to_satellite_file.replace(filename, param).replace(where_are_saved_satellite_data, '')
    
    return path_in_the_dict


def compute_grid_stats(df, grid_size):
    """Compute mean and std of values within a square grid of given size centered in the DataFrame.
    
    Args:
        df (pd.DataFrame): The input DataFrame with grid coordinates as column names.
        grid_size (int): The size of the square grid (must be odd).
    
    Returns:
        pd.DataFrame: A DataFrame with 'mean' and 'std' columns for each row.
    """
    if grid_size % 2 == 0:
        raise ValueError("Grid size must be odd to have a centered point.")

    # Extract unique x and y coordinates
    coords = np.array([list(map(int, col.split('_'))) for col in df.columns])
    max_x, max_y = coords.max(axis=0)  # Detect max x and y values
    center_x, center_y = max_x // 2 + 1, max_y // 2 + 1  # Find center dynamically

    # Compute the range of coordinates to include
    half_size = grid_size // 2
    x_range = range(center_x - half_size, center_x + half_size + 1)
    y_range = range(center_y - half_size, center_y + half_size + 1)

    # Generate expected column names
    cols_to_include = {f"{x}_{y}" for x in x_range for y in y_range} & set(df.columns)

    if not cols_to_include:
        raise ValueError("No matching columns found for the specified grid size.")

    # Compute mean and std using vectorized operations
    selected_data = df[list(cols_to_include)]
    return pd.DataFrame({f'{grid_size}x{grid_size}_mean': selected_data.mean(axis=1, skipna = True), 
                         f'{grid_size}x{grid_size}_median': selected_data.median(axis=1, skipna = True), 
                         f'{grid_size}x{grid_size}_std': selected_data.std(axis=1, skipna = True),
                         f'{grid_size}x{grid_size}_n': selected_data.count(axis=1, numeric_only=True)})


def get_the_corresponding_insitu_parameter_name(satellite_variable_name) : 
    
    if 'SPM' in satellite_variable_name : 
        return ['SPM', 'TURB']
    
    if 'CHL' in satellite_variable_name : 
        return ['CHLA']
    
    if 'RRS' in satellite_variable_name : 
        return ['CHLA']
    
    if 'SST' in satellite_variable_name : 
        return ['TEMP']


def Summarize_the_matchup_database(path, MU_database) : 
    
    MU_values = access_item_in_a_dictionnary(MU_database, path)
    
    MU_summary_df = {size: compute_grid_stats(MU_values, size) for size in np.arange(1,9,2)}
    MU_summary_df = pd.DataFrame({key: value for res in MU_summary_df.values() for key, value in res.items()})
        
    var_name_in_INSITU_dtb = get_the_corresponding_insitu_parameter_name(path)
    
    INSITU_values = []
    for var_name in var_name_in_INSITU_dtb : 
    
        INSITU_values_of_var_name = ( MU_database['INSITU']
                         .xs( MU_summary_df.index.get_level_values('DATE')[0], level = 'DATE')
                         .rename(columns={var_name: 'Insitu_value'})
                         .rename_axis(index={'TIME': 'Insitu_time'})[['Insitu_value']] )
        INSITU_values_of_var_name.insert(0, 'Insitu_variable', var_name)
        INSITU_values.append(INSITU_values_of_var_name)
        
    INSITU_values = pd.concat(INSITU_values)
    
    merged_df = (MU_summary_df.reset_index()
                 .merge(INSITU_values.reset_index(), 
                        on=list(set(MU_summary_df.index.names) & set(INSITU_values.index.names)), how='outer')
                 .set_index(list(set(MU_summary_df.index.names) | set(INSITU_values.index.names) | set(INSITU_values.columns)))
                 .reorder_levels( list(MU_summary_df.index.names) + ['Insitu_time', 'Insitu_variable', 'Insitu_value'] ) )
    
    MU_summary = create_arborescence( [path[1:]] )

    add_array_to_dict(MU_summary, path[1:], merged_df)
    
    return MU_summary


def get_insitu_dates_for_each_parameter(INSITU_data) :

    filtered_cols = INSITU_data.columns.difference(['SOURCE', 'SITE', 'REGION', "LATITUDE", "LONGITUDE", "DATE", "TIME"])  # Select relevant columns    
    dates_with_finite_values = {col: np.unique( INSITU_data.loc[np.isfinite(INSITU_data[col]), "DATE"] ) for col in filtered_cols}
    
    return dates_with_finite_values


def get_MU_criteria_for_each_product(info) : 
        
    if info.Data_source == 'SEXTANT' : 
            
        MU_criteria = {'min_n' : 1,
                       'max_hour_diff_between_insitu_and_satellite_measurement' : 3 if info.sensor_name != "merged" else np.nan,
                       'max_CV' : np.nan,
                       'grid_size' : 1}
        
    else :
        
        MU_criteria = {'min_n' : 5,
                       'max_hour_diff_between_insitu_and_satellite_measurement' : 3,
                       'max_CV' : 30,
                       'grid_size' : 3}
    
    return MU_criteria


def scatterplot_and_save_statistics(MU_summary_df, info, where_are_saved_Match_Up_data, zones) : 

    index_to_keep = np.where((MU_summary_df.Data_source == info.Data_source) & 
                             (MU_summary_df.sensor_name == info.sensor_name) & 
                             (MU_summary_df.atmospheric_correction == info.atmospheric_correction) & 
                             (MU_summary_df.Satellite_variable == info.Satellite_variable))
    
    MU_summary_df_of_the_case = MU_summary_df.iloc[index_to_keep]
    
    MU_criteria = get_MU_criteria_for_each_product(info)
    
    # Source the R script
    validate_R_path = os.path.join(func_dir, 'validate.R')
    robjects.r['source'](validate_R_path)
    r_function = robjects.r['Save_validation_scatterplots_and_stats']
    # print("R script sourced successfully")
    
    for satellite_algorithm in np.unique(MU_summary_df_of_the_case.Satellite_algorithm) : 
        
        MU_summary_df_of_the_sat_algo = MU_summary_df_of_the_case[ MU_summary_df_of_the_case.Satellite_algorithm == satellite_algorithm ]
                
        # MU_summary_df_of_the_sat_algo = MU_summary_df_of_the_sat_algo.iloc[:100]
        
        # Call the R function
        r_function(
            satellite_median = robjects.FloatVector(MU_summary_df_of_the_sat_algo[f'{MU_criteria["grid_size"]}x{MU_criteria["grid_size"]}_mean'].to_list()),
            satellite_n = robjects.IntVector(MU_summary_df_of_the_sat_algo[f'{MU_criteria["grid_size"]}x{MU_criteria["grid_size"]}_n'].to_list()),
            satellite_sd = robjects.FloatVector(MU_summary_df_of_the_sat_algo[f'{MU_criteria["grid_size"]}x{MU_criteria["grid_size"]}_std'].to_list()),
            satellite_times = robjects.StrVector(MU_summary_df_of_the_sat_algo.Satellite_time.astype(str).to_list()),
            insitu_variable = robjects.StrVector(MU_summary_df_of_the_sat_algo.Insitu_variable.to_list()),
            insitu_value = robjects.FloatVector(MU_summary_df_of_the_sat_algo.Insitu_value.to_list()),
            insitu_time = robjects.StrVector(MU_summary_df_of_the_sat_algo.Insitu_time.to_list()),
            insitu_Data_source = robjects.StrVector(MU_summary_df_of_the_sat_algo.SOURCE.to_list()),
            site_name = robjects.StrVector(MU_summary_df_of_the_sat_algo.SITE.to_list()),
            region_name = robjects.StrVector(MU_summary_df_of_the_sat_algo.REGION.to_list()),
            zones = robjects.StrVector(zones),
            min_n = robjects.IntVector([MU_criteria["min_n"]]),
            max_CV = robjects.IntVector([MU_criteria["max_CV"]]),
            max_hour_diff = robjects.IntVector([MU_criteria["max_hour_diff_between_insitu_and_satellite_measurement"]]),
            grid_size = robjects.IntVector([MU_criteria["grid_size"]]),
            date = robjects.StrVector(MU_summary_df_of_the_sat_algo.DATE.to_list()),
            satellite_source = robjects.StrVector([info.Data_source]),
            satellite_sensor = robjects.StrVector([info.sensor_name]),
            satellite_atm_corr = robjects.StrVector([info.atmospheric_correction]),
            satellite_algorithm = robjects.StrVector([satellite_algorithm]),
            where_to_save_MU_results = robjects.StrVector([where_are_saved_Match_Up_data])
        )
        
    r_function = robjects.r['Save_validation_scatterplots_and_stats']



# =============================================================================
#### Classes 
# =============================================================================


class MU_database_processing : 
    
    def __init__(self, where_to_save_Match_Up_data, cases_to_process, zones = "FRANCE", nb_of_cores_to_use = 1, redo_the_MU_database = False) :
        
        grid_size = 9 # 9x9
        
        path_to_MU_of_the_zone = os.path.join(where_to_save_Match_Up_data, "MATCH_UP_DATA", "_&_".join(zones))
        os.makedirs( path_to_MU_of_the_zone , exist_ok=True)
                     
        path_to_the_MU_database = os.path.join(path_to_MU_of_the_zone, "MU_database.joblib")
        
        if os.path.isfile(path_to_the_MU_database) and redo_the_MU_database == False : 
            MU_database = load( path_to_the_MU_database )        
        else :
            MU_database = []           
            
        self.MU_database = MU_database
        self.where_to_save_Match_Up_data = path_to_MU_of_the_zone
        self.path_to_the_MU_database = path_to_the_MU_database
        self.grid_size = grid_size
        self.cases_to_process = cases_to_process
        self.zones = zones
        self.nb_of_cores_to_use = nb_of_cores_to_use
            
    def Create_the_MU_database(self, where_are_saved_satellite_data) : 
        
        INSITU_data, INSITU_stations = get_insitu_measurements(self.zones)

        MU_data = {'STATIONS' : INSITU_stations,
                   'DATES' : get_insitu_dates_for_each_parameter(INSITU_data),
                   'INSITU' : INSITU_data.set_index(list(INSITU_stations.columns) + ['TIME'])}
                   
        MU_stations = MU_data['STATIONS'].loc[:,list(INSITU_stations.columns)[:-1]]
        
        MU_databases = [MU_data]
        filled_base_path = path_to_fill_to_where_to_save_satellite_files(where_are_saved_satellite_data)
        
        # Parallel processing with context manager
        # pool = multiprocess.Pool(self.nb_of_cores_to_use)
        with multiprocess.Pool(self.nb_of_cores_to_use) as pool:
            
            for i, info in self.cases_to_process.iterrows() : 
                        
                # info = cases_to_process.iloc[i].copy()
                                          
                print(f'{i} over {self.cases_to_process.shape[0]-1} ({info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable})')
                        
                filled_destination_path = fill_the_sat_paths(info = info, 
                                                             path_to_fill = filled_base_path, 
                                                             local_path = True,
                                                             dates = get_the_values_from_a_list_comprehension( 
                                                                         [MU_data['DATES'][the_param] 
                                                                          for the_param in get_the_corresponding_insitu_parameter_name(info.Satellite_variable)],
                                                                         True))
                
                MU_database_of_the_case = create_arborescence( [x.replace(f'{where_are_saved_satellite_data}/', '') for x in filled_destination_path] )
                
                satellite_files = find_sat_data_files(info, path_to_sat_data = filled_destination_path)
                if not satellite_files : 
                    continue            
                
                satellite_dates = np.unique( [match_obj.group(1) for file in satellite_files if (match_obj := re.search(r"(\d{4}/\d{2}/\d{2})", file))] )
                
                # Extract satellite coordinates only once
                with xr.open_dataset(satellite_files[0]) as ds:
                    satellite_coordinates = {dim: arr.values for dim, arr in ds.coords.items()}
                        
                MU_data['STATIONS'][['lat_index', 'lon_index']] = MU_data['STATIONS'].apply(
                    lambda row: find_indices_of_the_grid(float(row["LATITUDE"]), float(row["LONGITUDE"]), 
                                                         satellite_coordinates['lat'], satellite_coordinates['lon'],
                                                         grid_size = self.grid_size), 
                                                    axis=1, result_type="expand")
                                    
                MU_databases_of_the_case = pool.starmap(do_the_match_up_for_one_date, [(satellite_date, satellite_files, 
                                                                                MU_data, MU_database_of_the_case, 
                                                                                MU_stations, info,
                                                                                self.grid_size, where_are_saved_satellite_data) 
                                                           for satellite_date in satellite_dates ])
                
                MU_databases.append( merge_dicts( MU_databases_of_the_case ) ) 
                        
                MU_data['STATIONS'] = MU_data['STATIONS'].drop(['lat_index', 'lon_index'], axis = 1)
            
        self.MU_database = merge_dicts(MU_databases)
        
        dump(self.MU_database, self.path_to_the_MU_database , compress=3)
        
        self.any_modification_has_been_done = True
        
    def Complete_the_MU_database(self, where_are_saved_satellite_data) : 
        
        paths_already_filled_in_the_database = set( get_non_empty_paths(self.MU_database) )
        MU_databases_to_add = []
        
        MU_stations = self.MU_database['STATIONS'].loc[:,list(self.MU_database['STATIONS'].columns)[:-1]]

        # Precompute filled paths to avoid redundant calculations
        filled_destination_paths = {
            i: fill_the_sat_paths(
                info, path_to_fill_to_where_to_save_satellite_files(where_are_saved_satellite_data),
                local_path=True,                 
                dates = get_the_values_from_a_list_comprehension( 
                                    [self.MU_database['DATES'][the_param] 
                                     for the_param in get_the_corresponding_insitu_parameter_name(info.Satellite_variable)],
                                    True)
            )
            for i, info in self.cases_to_process.iterrows()
        }
        
        self.any_modification_has_been_done = False
        
        # pool = multiprocessing.Pool(self.nb_of_cores_to_use)
        with multiprocess.Pool(self.nb_of_cores_to_use) as pool:

            for i, info in self.cases_to_process.iterrows() : 
                                              
                # info = self.cases_to_process.iloc[0]
                
                print(f'{i} over {self.cases_to_process.shape[0]-1} ({info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable})')
                
                satellite_files = find_sat_data_files(info, path_to_sat_data = filled_destination_paths[i])
                if not satellite_files : 
                    continue     
                
                paths_to_process = pool.starmap(find_the_path_to_satellite_MU_in_the_dict, 
                                                [(satellite_file, where_are_saved_satellite_data) for satellite_file in satellite_files ])
    
                index_to_process = np.where( ~ np.isin( paths_to_process, list(paths_already_filled_in_the_database) ) )[0]
                if not len(index_to_process) :
                    continue    
                
                paths_to_process = list( np.array(paths_to_process)[ index_to_process ] )
                satellite_files = list( np.array(satellite_files)[ index_to_process ] )
                
                MU_database_of_the_case = create_arborescence( paths_to_process )
     
                satellite_dates = np.unique( [match_obj.group(1) for file in satellite_files if (match_obj := re.search(r"(\d{4}/\d{2}/\d{2})", file))] )
                
                with xr.open_dataset(satellite_files[0]) as ds :
                    satellite_coordinates = {dim: arr.values for dim, arr in ds.coords.items()}
                        
                self.MU_database['STATIONS'][['lat_index', 'lon_index']] = self.MU_database['STATIONS'].apply(
                    lambda row: find_indices_of_the_grid(float(row["LATITUDE"]), float(row["LONGITUDE"]), 
                                                         satellite_coordinates['lat'], satellite_coordinates['lon'],
                                                         grid_size = self.grid_size), 
                                                    axis=1, result_type="expand")
                                    
                MU_databases_of_the_case = pool.starmap(do_the_match_up_for_one_date, [(satellite_date, satellite_files, 
                                                                                self.MU_database, MU_database_of_the_case, 
                                                                                MU_stations, info,
                                                                                self.grid_size, where_are_saved_satellite_data) 
                                                           for satellite_date in satellite_dates ])
    
                merged_MU_databases_of_the_case = merge_dicts( MU_databases_of_the_case )
                if "" in merged_MU_databases_of_the_case : 
                    del merged_MU_databases_of_the_case['']
                
                MU_databases_to_add.append( merged_MU_databases_of_the_case ) 
                        
                self.MU_database['STATIONS'] = self.MU_database['STATIONS'].drop(['lat_index', 'lon_index'], axis = 1)
                self.any_modification_has_been_done = True
            
        if self.any_modification_has_been_done : 
            
            MU_databases_to_add = merge_dicts( MU_databases_to_add )
            self.MU_database = merge_dicts( [MU_databases_to_add, self.MU_database] )
            
            dump(self.MU_database, self.path_to_the_MU_database, compress=3)
        
    # def Perform_matchups
    def Summarize_MU_with_statistics(self) : 

        if self.any_modification_has_been_done == False :          

            self.MU_summary = load( self.path_to_the_MU_database.replace('database', 'summary') )
            return
            
        paths_filled_in_the_database = set( [x for x in get_non_empty_paths(self.MU_database) if x.count('/') > 2] )
                
        with multiprocess.Pool(self.nb_of_cores_to_use) as pool:
            MU_summaries = pool.starmap(Summarize_the_matchup_database, [(path, 
                                                                          self.MU_database) 
                                                       for path in paths_filled_in_the_database ])  
            
        # For debugging
        # for path in paths_filled_in_the_database  : 
        #     print(path)
        #     a = Summarize_the_matchup_database(path, self.MU_database) 
            
        MU_summary = merge_dicts( [x for x in MU_summaries if x is not None ] ) 
        
        MU_summary = {**self.MU_database, **MU_summary}
        
        dump(MU_summary, self.path_to_the_MU_database.replace('database', 'summary'), compress=3)
        
        self.MU_summary = MU_summary
        
        
    def Compile_and_Save_MU_summary(self, Data_sources) :
        
        def process_data(data_dict):
            df = {k: data_dict[k] for k in Data_sources if k in data_dict}
            df = extract_dataframes_iterative(df)       
            df = list(df)
            return pd.concat(df)
        
        base_path_to_save_data = "/".join(self.path_to_the_MU_database.split('/')[:-1])
        
        if self.any_modification_has_been_done == False :          

            self.MU_database_df = pd.read_csv(f"{base_path_to_save_data}/database.csv", low_memory=False)
            self.MU_summary_df = pd.read_csv(f"{base_path_to_save_data}/summary.csv", low_memory=False)
            return

        MU_summary_df = process_data(self.MU_summary)
        MU_database_df = process_data(self.MU_database).reindex(MU_summary_df.index)
    
        mask = np.isfinite(MU_summary_df.index.get_level_values('Insitu_value'))
        MU_summary_df, MU_database_df = MU_summary_df[mask], MU_database_df[mask]
        
        self.MU_summary_df = MU_summary_df
        self.MU_database_df = MU_database_df
    
        MU_database_df.to_csv(f"{base_path_to_save_data}/database.csv")
        MU_summary_df.to_csv(f"{base_path_to_save_data}/summary.csv")
        
        
    def Make_scatterplot_and_save_statistics(self) : 
        
        # MU_summary_df = self.MU_summary_df.reset_index(drop = False)
        MU_summary_df = pd.read_csv("/".join(self.path_to_the_MU_database.split('/')[:-1]) + '/summary.csv')
        where_to_save_Match_Up_plots = self.where_to_save_Match_Up_data
        cases_to_process = self.cases_to_process
        zones = self.zones
        
        with multiprocess.Pool(self.nb_of_cores_to_use) as pool:
            
            pool.starmap(scatterplot_and_save_statistics, 
                         [(MU_summary_df, case_to_process, where_to_save_Match_Up_plots, zones) 
                          for _, case_to_process in cases_to_process.iterrows() ])
        
        # Get all statistics .csv files
        
        # Make the final table with them 


# =============================================================================
#### Main functions
# =============================================================================


def Match_up_with_insitu_measurements(core_arguments, zones, redo_the_MU_database, nb_of_cores_to_use,
                                      where_are_saved_satellite_data, where_to_save_Match_Up_data) : 
             
    core_arguments.update({'Temporal_resolution' : ['DAILY']})
    
    cases_to_process = get_all_cases_to_process(core_arguments)
               
    MU_database = MU_database_processing(where_to_save_Match_Up_data = where_to_save_Match_Up_data, 
                              cases_to_process = cases_to_process,
                              zones = zones,
                              redo_the_MU_database = redo_the_MU_database,
                              nb_of_cores_to_use = nb_of_cores_to_use)
    
    if len(MU_database.MU_database) == 0 :
        MU_database.Create_the_MU_database(where_are_saved_satellite_data)
    else : 
        MU_database.Complete_the_MU_database(where_are_saved_satellite_data)
      
    MU_database.Summarize_MU_with_statistics()
    
    MU_database.Compile_and_Save_MU_summary( core_arguments['Data_sources'] )
    # print("Not broken yet")
    MU_database.Make_scatterplot_and_save_statistics()
    
    # MU_database.Make_time_series()

