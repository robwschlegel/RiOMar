#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import utils and data downloading functions
# TODO: Get this to run without needing to call the exact root path
exec(open("/Users/rws/RiOMar/code/0_data_management.py").read())

# Additional packages
# import pickle # unused
from joblib import dump, load

# Static directory addresses
# TODO: Change these static pathways once this issue has been resolved
path_to_SOMLIT_insitu_data = '/Users/rws/RiOMar/data/IN_SITU/Somlit.csv'
path_to_REPHY_insitu_data = '/Users/rws/RiOMar/data/REPHY/Table1_REPHY_hydro_RIOMAR.csv'


# =============================================================================
#### Functions
# =============================================================================


def get_insitu_measurements(path_to_SOMLIT_insitu_data = None, path_to_REPHY_insitu_data = None) : 
    
    SOMLIT_data_filtered = None; REPHY_data_filtered = None
    
    if path_to_SOMLIT_insitu_data is not None : 
    
        QC_values_to_keep = [2,6,7]
        Depth_threshold = 10    
    
        SOMLIT_data = (pd.read_csv(path_to_SOMLIT_insitu_data, sep = ";", header = 2).iloc[1:]
                            .rename(columns = {'gpsLat*':'LATITUDE', 
                                               'gpsLong*':'LONGITUDE',
                                               'nomSite*':"Site"}))
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
                                .loc[:,['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE', 'TIME', 'TEMP', 'SAL', 'POC', 'SPM']])
        
        # SOMLIT_stations = SOMLIT_data_filtered.loc[:,["ID_SITE", "Site", "DATE", "HEURE", "Latitude", "Longitude"]].drop_duplicates().reset_index(drop = True)

    if path_to_REPHY_insitu_data is not None : 

        var_names = ['CHLOROA', 'TURB', 'TURB-FNU', 'TEMP', 'SALI']
        QC_values_to_keep = ['Bon']
        Depth_threshold = 10

        REPHY_data = pd.read_csv(path_to_REPHY_insitu_data, sep = ";", header = 0, encoding="ISO-8859-1")            
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
                                   .loc[:,['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE', 'TIME', 'VARIABLE', 'VALUE', 'QC']] 
                                   .pivot_table(index=['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE', 'TIME'], 
                                                columns="VARIABLE", values="VALUE", aggfunc="mean")
                                   .rename(columns = {'CHLOROA':'CHLA', 
                                                      'SALI':'SAL'})
                                   .reset_index(drop = False))
        
        REPHY_data_filtered['TURB_all'] = REPHY_data_filtered[["TURB-FNU", "TURB"]].mean(axis=1)
        REPHY_data_filtered['TURB'] = REPHY_data_filtered['TURB']
        REPHY_data_filtered = REPHY_data_filtered.drop(['TURB_all', 'TURB-FNU'], axis = 1)
        
    INSITU_measurements = pd.concat([SOMLIT_data_filtered, REPHY_data_filtered], ignore_index = True)
    INSITU_measurements['DATE'] = INSITU_measurements['TIME'].dt.date.astype(str)
    INSITU_measurements['TIME'] = INSITU_measurements['TIME'].dt.time.astype(str)
    
    INSITU_stations = (INSITU_measurements
                       .groupby(["SOURCE", "SITE", "LATITUDE", "LONGITUDE"], as_index = False)
                       .agg({"DATE": list}))
    # INSITU_stations = INSITU_measurements.loc[:,["SOURCE", "SITE", "LATITUDE", "LONGITUDE"]].drop_duplicates().reset_index(drop = True)

    return INSITU_measurements, INSITU_stations


def find_indices_of_the_grid(lat_station, lon_station, lat_map, lon_map, grid_size):
    
    half_grid_size = int( (grid_size - 1) / 2 )
    
    lat_idx = np.abs(lat_map - lat_station).argmin()
    lon_idx = np.abs(lon_map - lon_station).argmin()
    
    lat_idx_all = np.arange(lat_idx-half_grid_size, lat_idx+half_grid_size+1)
    lon_idx_all = np.arange(lon_idx-half_grid_size, lon_idx+half_grid_size+1)
    
    return {'lat_index' : lat_idx_all, 'lon_index' : lon_idx_all}


def do_the_match_up_for_one_date(satellite_date, satellite_files, MU_data, MU_database_of_the_case, 
                                 MU_stations, info, grid_size, where_are_saved_sat_files) : 

    date_pattern = re.compile(satellite_date)
    
    satellite_files_of_the_day = [file for file in satellite_files if date_pattern.search(file)]
    
    satellite_date_formatted = pd.to_datetime(satellite_date).strftime('%Y-%m-%d')
        
    for satellite_file in satellite_files_of_the_day : 
    
        try : 
            map_ini = xr.open_dataset(satellite_file, decode_times=True)
        except Exception as e:
            
            add_array_to_dict(dictionary = MU_database_of_the_case,
                              path = satellite_file.replace(f'{where_are_saved_sat_files}/', ''),
                              array = f"❌ Impossible to load the nc file {satellite_file} - Error : {e}")
            
            continue
        
        if ( (len(map_ini) == 0) or 
             (len(np.unique(map_ini.lat)) == 1) or 
             (len(np.unique(map_ini.lon)) == 1) ) : 
            
            add_array_to_dict(dictionary = MU_database_of_the_case,
                              path = satellite_file.replace(f'{where_are_saved_sat_files}/', ''),
                              array = f"❌ Wrong file format for {satellite_file}")
            # return 
            continue
        
        the_var_name = next(iter(map_ini.keys()))
        
        time_value = extract_the_time_from_the_satellite_file(map_ini)
        
        MU_values = [
            map_ini.isel(lat=station_info['lat_index'], lon=station_info['lon_index'])[the_var_name].values.flatten()
            if satellite_date_formatted in station_info['DATE'] else np.array([])
            for _, station_info in MU_data['STATIONS'].iterrows()
        ]
        
        MU_values = pd.DataFrame(MU_values)
        MU_values.columns = [f"{i}_{j}" for i in range(1, grid_size+1) for j in range(1, grid_size+1)]     
        MU_values.index = pd.MultiIndex.from_frame(
                                MU_stations.assign(**pd.concat([info,
                                                                pd.Series({'Satellite_algorithm' : return_the_parameter_name_based_on_file_name(satellite_file),
                                                                           'DATE' : satellite_date_formatted,
                                                                           'Satellite_time' : time_value})])))
                 
        map_ini.close()
        
        add_array_to_dict(dictionary = MU_database_of_the_case,
                          path = satellite_file.replace(f'{where_are_saved_sat_files}/', ''),
                          array = MU_values)
            
    return MU_database_of_the_case


def find_the_path_to_satellite_MU_in_the_dict(path_to_satellite_file, where_are_saved_sat_files) : 
    
    filename = path_to_satellite_file.split("/")[-1]
    
    # Extract parameter name (e.g., SPM-G)
    param = return_the_parameter_name_based_on_file_name(filename)
    path_in_the_dict = path_to_satellite_file.replace(filename, param).replace(where_are_saved_sat_files, '')
    
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
                         f'{grid_size}x{grid_size}_std': selected_data.std(axis=1, skipna = True),
                         f'{grid_size}x{grid_size}_n': selected_data.count(axis=1, numeric_only=True)})


def Summarize_the_matchup_database(path, MU_database) : 
    
    MU_values = access_item_in_a_dictionnary(MU_database, path)
    
    MU_summary_df = {size: compute_grid_stats(MU_values, size) for size in np.arange(1,9,2)}
    MU_summary_df = pd.DataFrame({key: value for res in MU_summary_df.values() for key, value in res.items()})
    
    var_name_map = {"SPM": "SPM", "CHL": "CHLA"}
    var_name_in_INSITU_dtb = next((v for k, v in var_name_map.items() if k in path), None)
    
    INSITU_values = ( MU_database['INSITU']
                     .xs( MU_summary_df.index.get_level_values('DATE')[0], level = 'DATE')
                     .rename(columns={var_name_in_INSITU_dtb: 'Insitu_value'})
                     .rename_axis(index={'TIME': 'Insitu_time'})[['Insitu_value']] )
    
    merged_df = (MU_summary_df.reset_index()
                 .merge(INSITU_values.reset_index(), 
                        on=list(set(MU_summary_df.index.names) & set(INSITU_values.index.names)), how='outer')
                 .set_index(list(set(MU_summary_df.index.names) | set(INSITU_values.index.names) | set(INSITU_values.columns)))
                 .reorder_levels( list(MU_summary_df.index.names) + ['Insitu_time', 'Insitu_value'] ) )
    
    MU_summary = create_arborescence( [path[1:]] )

    add_array_to_dict(MU_summary, path[1:], merged_df)
    
    return MU_summary

# =============================================================================
#### Classes 
# =============================================================================

class MU_database_processing : 
    
    def __init__(self, where_to_save_Match_Up_data, cases_to_process, nb_of_cores_to_use, redo_the_MU_database = False) :
        
        grid_size = 9 # 9x9
        
        os.makedirs( where_to_save_Match_Up_data, exist_ok=True)
                     
        path_to_the_MU_database = where_to_save_Match_Up_data + '/MU_database.joblib'
        
        if os.path.isfile(path_to_the_MU_database) and redo_the_MU_database == False : 
            MU_database = load( path_to_the_MU_database )        
        else :
            MU_database = []           
            
        self.MU_database = MU_database
        self.where_to_save_Match_Up_data = where_to_save_Match_Up_data
        self.path_to_the_MU_database = path_to_the_MU_database
        self.grid_size = grid_size
        self.cases_to_process = cases_to_process
        self.nb_of_cores_to_use = nb_of_cores_to_use
            
    def Create_the_MU_database(self, where_are_saved_sat_files) : 
        
        INSITU_data, INSITU_stations = get_insitu_measurements(path_to_SOMLIT_insitu_data, 
                                                               path_to_REPHY_insitu_data)

        MU_data = {'STATIONS' : INSITU_stations,
                   'DATES' : np.unique(INSITU_data['DATE']),
                   'INSITU' : INSITU_data.set_index(['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE', 'DATE', 'TIME'])}
                   
        MU_stations = MU_data['STATIONS'].loc[:,['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE']]
        
        MU_databases = [MU_data]
        filled_base_path = path_to_fill_to_where_to_save_satellite_files(where_are_saved_sat_files)
        
        # Parallel processing with context manager
        # pool = multiprocessing.Pool(self.nb_of_cores_to_use)
        with multiprocessing.Pool(self.nb_of_cores_to_use) as pool:
            
            for i, info in self.cases_to_process.iterrows() : 
                        
                # info = cases_to_process.iloc[i].copy()
                                          
                print(f'{i} over {self.cases_to_process.shape[0]-1} ({info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable})')
                        
                filled_destination_path = fill_the_sat_paths(info = info, 
                                                             path_to_fill = filled_base_path, 
                                                             local_path = True,
                                                             dates = MU_data['DATES'])
                
                MU_database_of_the_case = create_arborescence( [x.replace(f'{where_are_saved_sat_files}/', '') for x in filled_destination_path] )
                
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
                                                                                self.grid_size, where_are_saved_sat_files) 
                                                           for satellite_date in satellite_dates ])
                
                MU_databases.append( merge_dicts( MU_databases_of_the_case ) ) 
                        
                MU_data['STATIONS'] = MU_data['STATIONS'].drop(['lat_index', 'lon_index'], axis = 1)
            
        self.MU_database = merge_dicts(MU_databases)
        
        dump(self.MU_database, self.path_to_the_MU_database , compress=3)
        
    def Complete_the_MU_database(self, where_are_saved_sat_files) : 
        
        paths_already_filled_in_the_database = set( get_non_empty_paths(self.MU_database) )
        MU_databases_to_add = []
        
        MU_stations = self.MU_database['STATIONS'].loc[:,['SOURCE', 'SITE', 'LATITUDE', 'LONGITUDE']]

        # Precompute filled paths to avoid redundant calculations
        filled_destination_paths = {
            i: fill_the_sat_paths(
                info, path_to_fill_to_where_to_save_satellite_files(where_are_saved_sat_files),
                local_path=True, dates=self.MU_database['DATES']
            )
            for i, info in self.cases_to_process.iterrows()
        }
        
        # pool = multiprocessing.Pool(self.nb_of_cores_to_use)
        with multiprocessing.Pool(self.nb_of_cores_to_use) as pool:

            for i, info in self.cases_to_process.iterrows() : 
                                              
                # info = self.cases_to_process.iloc[0]
                
                print(f'{i} over {self.cases_to_process.shape[0]-1} ({info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Satellite_variable})')
                
                satellite_files = find_sat_data_files(info, path_to_sat_data = filled_destination_paths[i])
                
                if not satellite_files : 
                    continue     
                
                paths_to_process = pool.starmap(find_the_path_to_satellite_MU_in_the_dict, 
                                                [(satellite_file, where_are_saved_sat_files) for satellite_file in satellite_files ])
    
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
                                                                                self.grid_size, where_are_saved_sat_files) 
                                                           for satellite_date in satellite_dates ])
    
                merged_MU_databases_of_the_case = merge_dicts( MU_databases_of_the_case )
                if "" in merged_MU_databases_of_the_case : 
                    del merged_MU_databases_of_the_case['']
                
                MU_databases_to_add.append( merged_MU_databases_of_the_case ) 
                        
                self.MU_database['STATIONS'] = self.MU_database['STATIONS'].drop(['lat_index', 'lon_index'], axis = 1)
            
        MU_databases_to_add = merge_dicts( MU_databases_to_add )
        self.MU_database = merge_dicts( [MU_databases_to_add, self.MU_database] )
        
        dump(self.MU_database, self.path_to_the_MU_database, compress=3)
        
    # def Perform_matchups
    def Summarize_MU_with_statistics(self) : 

        paths_filled_in_the_database = set( [x for x in get_non_empty_paths(self.MU_database) if x.count('/') > 1] )
                
        with multiprocessing.Pool(self.nb_of_cores_to_use) as pool:
            MU_summaries = pool.starmap(Summarize_the_matchup_database, [(path, 
                                                                          self.MU_database) 
                                                       for path in paths_filled_in_the_database ])           
            
        MU_summary = merge_dicts( [x for x in MU_summaries if x is not None ] ) 
        
        MU_summary = {**self.MU_database, **MU_summary}
        
        dump(MU_summary, self.path_to_the_MU_database.replace('database', 'summary'), compress=3)
        
        self.MU_summary = MU_summary
        
        
    def Compile_and_Save_MU_summary(self, Data_sources) :
        
        MU_summary_df = {key: self.MU_summary[key] for key in Data_sources if key in self.MU_summary}
        
        MU_summary_df = list( extract_dataframes_iterative(MU_summary_df) )
        
        pd.concat( MU_summary_df ).to_csv(self.path_to_the_MU_database.replace('database.joblib', 'summary.csv'))


## Main functions

def Match_up_with_insitu_measurements(core_arguments, redo_the_MU_database, nb_of_cores_to_use,
                                      where_to_save_satellite_data, where_to_save_Match_Up_data) : 
             
    cases_to_process = get_all_cases_to_process(core_arguments)
               
    MU_database = MU_database_processing(where_to_save_Match_Up_data = where_to_save_Match_Up_data, 
                              cases_to_process = cases_to_process,
                              redo_the_MU_database = redo_the_MU_database,
                              nb_of_cores_to_use = nb_of_cores_to_use)
    
    if len(MU_database.MU_database) == 0 :

        # MU_database.Create_the_MU_database(path_to_SOMLIT_insitu_data, path_to_REPHY_insitu_data, where_to_save_satellite_data)
        MU_database.Create_the_MU_database(where_to_save_satellite_data)
    else : 
        MU_database.Complete_the_MU_database(where_to_save_satellite_data)
                
    MU_database.Summarize_MU_with_statistics()
    
    MU_database.Compile_and_Save_MU_summary(Data_sources)
    
    del MU_database
    gc.collect()