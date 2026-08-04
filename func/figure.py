#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import os, sys, pickle, re, glob, subprocess
import numpy as np
import xarray as xr
import pandas as pd
import rpy2.robjects as robjects
from functools import reduce
from PIL import Image

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (load_csv_files, path_to_fill_to_where_to_save_satellite_files,
                  align_bathymetry, define_parameters, order_zones)
# Create_the_plume_mask/delineate_plume_pipeline/create_polygon_mask
# used to be imported from this repo's own func/plume.py, which forked the
# plume-detection algorithm and had drifted out of sync with the version
# actually used to produce every other result in this project (e.g. it was
# missing the near-mouth-quantile threshold-bound estimation panache gained
# on 17 Jul). Now imported directly from the installed panache package so
# Figure 3/4's methodology panels reflect the real, current algorithm.
from panache.plume_algorithm import Create_the_plume_mask, delineate_plume_pipeline, create_polygon_mask
# panache.utils.define_parameters is a separately-maintained, cleaned-up
# fork of this repo's own util.py::define_parameters -- it's authoritative
# for anything the *algorithm* reads (thresholds, searching_strategies as
# preset names rather than raw grids, etc; confirmed via diff: RiOMar's
# copy has zero keys panache lacks, only 6 figure-display-only extras --
# lat/lon_new_resolution, lat/lon_range_of_the_area_to_check_for_clouds,
# lat/lon_range_of_the_map_to_plot). build_plume_parameters() below merges
# panache's authoritative values over RiOMar's, keeping only those 6 extras
# from the local copy.
from panache.utils import define_parameters as panache_define_parameters


def build_plume_parameters(Zone):
    """
    RiOMar's own figure-display parameters (map-crop extent, resolution
    reduction), overlaid with panache's authoritative algorithm parameters
    -- so Figure_3_panels()/Figure_3_zone_maps() detect plumes exactly as the real pipeline
    does, not as RiOMar's older, now-unmaintained local fork of the
    algorithm would.
    """
    parameters = define_parameters(Zone)
    parameters.update(panache_define_parameters(Zone))
    return parameters


def preprocess_annual_dataset_and_compute_land_mask(path_to_annual_ds, parameters) :

    """
    Load the annual map and use it to generate a land mask.

    This function processes an annual dataset by loading it and generating
    a land mask based on null values, at native resolution (plume detection
    no longer coarsens the data -- see panache.plume_algorithm.main_process,
    which sets ds_reduced = ds). RiOMar-specific figure-prep utility with no
    equivalent in panache.

    Parameters
    ----------
    path_to_annual_ds : str
        the path to the annual (yearly) averaged dataset
    parameters : dict
        Configuration parameters for plume detection. It contains the limits of the polygon defining the searching area.

    Returns
    -------
    tuple
        A tuple containing:
        - xarray.DataArray: Subset of the annual dataset for cloud coverage checks.
        - xarray.DataArray: Boolean mask where True represents land and False represents water.
    """

    # Load the annual dataset for plume detection
    with open(path_to_annual_ds, 'rb') as f:

        multi_annual_ds = pickle.load(f)['Basin_map']['map_data']
        multi_annual_ds.values[multi_annual_ds.values < 0] = np.nan

    # Select a subset of the annual dataset based on the area to check for cloud coverage
    multi_annual_ds_subset = multi_annual_ds.sel(lat=slice(parameters['lat_range_of_the_area_to_check_for_clouds'][0],
                                                           parameters['lat_range_of_the_area_to_check_for_clouds'][1]),
                                                lon=slice(parameters['lon_range_of_the_area_to_check_for_clouds'][0],
                                                          parameters['lon_range_of_the_area_to_check_for_clouds'][1]))

    # Create a land mask based on null values (land = True), at native resolution
    land_mask = multi_annual_ds.isnull()

    return multi_annual_ds_subset, land_mask


# =============================================================================
#### Utility functions
# =============================================================================


def do_R_plot(the_plume, where_to_save_the_plot, name_of_the_plot):
    """
    Convert a panache Create_the_plume_mask instance's SPM map and plume mask
    to an R dataframe and plot using ggplot2. Extracted out of the class
    itself (it used to be a bound method, the_plume.do_R_plot(...)) because
    it's a RiOMar-manuscript-specific plotting helper, not something the
    general-purpose panache library needs -- panache's own version of the
    class doesn't have it, only reads the same public attributes it sets.
    """

    folder_where_to_save_data = os.path.join(where_to_save_the_plot, 'DATA')
    os.makedirs( folder_where_to_save_data, exist_ok=True )

    # Create a Pandas DataFrame
    lat_values = the_plume.spm_map.lat.values
    lon_values = the_plume.spm_map.lon.values

    df = pd.DataFrame({
        'lat': np.repeat(lat_values, len(lon_values)),
        'lon': np.tile(lon_values, len(lat_values)),
        'analysed_spim': the_plume.spm_map.values.flatten()
    })

    if name_of_the_plot == 'B' :

        latitudes_used_for_finding_SPM_threshold = lat_values[the_plume.points_used_for_finding_SPM_threshold[:,1].astype(int)]
        longitudes_used_for_finding_SPM_threshold = lon_values[the_plume.points_used_for_finding_SPM_threshold[:,0].astype(int)]
        points_used_for_finding_SPM_threshold = pd.DataFrame({'longitude': longitudes_used_for_finding_SPM_threshold,
                                                              'latitude': latitudes_used_for_finding_SPM_threshold})
        points_used_for_finding_SPM_threshold.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}_points_used_for_finding_SPM_threshold.csv") )

        all_latitudes_used_for_finding_SPM_threshold = lat_values[the_plume.all_points_tested_for_finding_SPM_threshold[:,1].astype(int)]
        all_longitudes_used_for_finding_SPM_threshold = lon_values[the_plume.all_points_tested_for_finding_SPM_threshold[:,0].astype(int)]
        all_points_used_for_finding_SPM_threshold = pd.DataFrame({'longitude': all_longitudes_used_for_finding_SPM_threshold,
                                                              'latitude': all_latitudes_used_for_finding_SPM_threshold})
        all_points_used_for_finding_SPM_threshold.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}_all_points_used_for_finding_SPM_threshold.csv") )


    if 'plume_mask' in vars(the_plume) :
        df['plume'] = the_plume.plume_mask.values.flatten()

    index_to_keep = np.where((df.lat >= the_plume.parameters['lat_range_of_the_map_to_plot'][0]) &
                             (df.lat <= the_plume.parameters['lat_range_of_the_map_to_plot'][1]) &
                             (df.lon >= the_plume.parameters['lon_range_of_the_map_to_plot'][0]) &
                             (df.lon <= the_plume.parameters['lon_range_of_the_map_to_plot'][1]) )

    df = df.iloc[index_to_keep]

    df.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}.csv") )

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_3_panel']

    # Call the R function
    r_function(
        where_to_save_the_figure = robjects.StrVector([where_to_save_the_plot]),
        name_of_the_plot = robjects.StrVector([name_of_the_plot])
    )


def save_files_for_Figure_1(where_are_saved_satellite_data, where_to_save_the_figure,
                            mean_spm_path, coordinates_of_the_map) :
    # Main-panel background changed 2026-08-01 (per Robert) from a single
    # hardcoded snapshot date (2011/02/02) to the pixel-wise temporal mean
    # SPM field across the whole SEXTANT archive (func/compute_mean_spm.py,
    # output/STATS/mean_spm_national.nc) -- a proper climatological
    # background rather than one arbitrary day's cloud/turbidity pattern.

    folder_where_to_save_Figure_1_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURE_1', 'DATA')
    os.makedirs(folder_where_to_save_Figure_1_data, exist_ok = True)

    with xr.open_dataset(mean_spm_path) as ds :
        SPM_map = (ds['mean_analysed_spim']
                   .sel(lat=slice(coordinates_of_the_map['lat_min'], coordinates_of_the_map['lat_max']),
                        lon=slice(coordinates_of_the_map['lon_min'], coordinates_of_the_map['lon_max']))
                   .to_dataframe()
                   .reset_index()
                   .rename(columns={'mean_analysed_spim': 'analysed_spim'}))

        SPM_map.to_csv(folder_where_to_save_Figure_1_data + "/SPM_map.csv")

    save_insitu_stations_for_plot(folder_where_to_save_Figure_1_data)
    
    
def load_and_save_regional_maps_for_plot(where_are_saved_regional_maps, where_to_save_the_figure, dates_for_each_zone) :

    # NB: these per-zone regional maps are the insets combined into Figure 1
    # (with SOMLIT/REPHY stations) -- so they are written into FIGURE_1's own
    # DATA folder, not a separate FIGURE_2 (manuscript Figure 2 is the
    # satellite-vs-in-situ validation scatterplot from func/validate.py).
    folder_where_to_save_regional_zone_maps_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURE_1', 'DATA')
    os.makedirs(folder_where_to_save_regional_zone_maps_data, exist_ok = True)
    
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
            
            SPM_map.to_csv(folder_where_to_save_regional_zone_maps_data + f"/{key}.csv")
            
    save_insitu_stations_for_plot(folder_where_to_save_regional_zone_maps_data)


def save_insitu_stations_for_plot(folder_where_to_save_Figure_data) :     

    coordinates_of_the_RIOMARS = { zone_name : {'lat_min' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][0],
                                                'lat_max' : define_parameters(zone_name)['lat_range_of_the_map_to_plot'][1],
                                                'lon_min' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][0],
                                                'lon_max' : define_parameters(zone_name)['lon_range_of_the_map_to_plot'][1]}
                                  for zone_name in order_zones(['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY', 'BAY_OF_SEINE']) }
    
    pd.DataFrame.from_dict(coordinates_of_the_RIOMARS).to_csv(folder_where_to_save_Figure_data + "/RIOMAR_limits.csv")
    
    # Station coordinates -- written by func/validate.R's in-situ prep
    # section (metadata/in_situ_site_list.csv), which replaced the retired
    # func/validate.py::get_insitu_measurements() as the source of this data.
    insitu_stations = pd.read_csv(os.path.join(proj_dir, 'metadata', 'in_situ_site_list.csv'))
    insitu_stations = insitu_stations.rename(columns={'source': 'SOURCE', 'site': 'SITE',
                                                       'lat': 'LATITUDE', 'lon': 'LONGITUDE',
                                                       'zone': 'REGION'})
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
    # Dict order here is the panel order used by Figure 1's insets and
    # Figure_3_zone_maps() -- kept north to south (util.ZONE_ORDER) rather
    # than alphabetical or ad hoc, so multi-panel figures read consistently
    # top to bottom.
    #
    # Dates updated 2026-08-01 (per Robert): each zone's date of maximum
    # dynamic-threshold plume area (output/panache/dynamic/<zone>/Results.csv),
    # with func/util.R's plume_area_ceiling artefact screen applied first --
    # the raw per-zone maxima before that screen (Bay of Seine 2014-02-02 at
    # 9563 km^2; Gulf of Lion 2005-05-12 at 22359 km^2) are known flood-fill
    # spillover artefacts, both far above their zone's ceiling (2500 and
    # 10000 km^2 respectively), not real plume days.
    dates = {'GULF_OF_LION' : '2005-01-27',
            'BAY_OF_BISCAY' : '2018-01-17',
            'SOUTHERN_BRITTANY' : '2007-12-09',
            'BAY_OF_SEINE' : '2001-02-19'}
    return {zone : dates[zone] for zone in order_zones(dates.keys())}


# =============================================================================
#### Main functions
# =============================================================================


def Figure_1(where_are_saved_satellite_data, where_are_saved_regional_maps, where_to_save_the_figure):
    save_files_for_Figure_1(where_are_saved_satellite_data,
                            where_to_save_the_figure,
                            mean_spm_path=os.path.join(proj_dir, "output", "STATS", "mean_spm_national.nc"),
                            coordinates_of_the_map={"lat_min": 42, "lat_max": 51.5, "lon_min": -6, "lon_max": 8})

    # The four regional zoomed insets (the former standalone Figure_2(), with
    # SOMLIT/REPHY station overlays) are built inside Figure_1() in figure.R.
    # Changed 2026-08-01 (per Robert): each inset now uses its own zone's date
    # of maximum plume area (dates_for_each_zone()) instead of being cropped
    # from the national map's single date/mean -- the point of the change is
    # to actually show each zone's plume at its most visually striking, so a
    # climatological mean (now used for the national backdrop) or a single
    # shared arbitrary date would defeat that. Manuscript Figure 2 is the
    # satellite-vs-in-situ validation scatterplot instead, produced
    # separately by func/validate.py during 1_validate.py.
    load_and_save_regional_maps_for_plot(where_are_saved_regional_maps,
                                                      where_to_save_the_figure,
                                                      dates_for_each_zone())

    # Run as a standalone Rscript process rather than via in-process rpy2
    # (like every other figure function below): Figure_1() is the only
    # figure that loads the 'sf' package (for the zoomed insets' high-res
    # coastline), and 'sf' links against the system's GDAL/GEOS/PROJ, which
    # conflict with the conda geospatial stack (shapely/geopandas/pyproj)
    # already loaded in this process via panache -- two incompatible
    # GEOS/PROJ/HDF5 builds in one process crashes R. A separate process has
    # its own library address space, so the conflict can't occur.
    subprocess.run(
        ['Rscript', 'func/run_figure_1.R', where_to_save_the_figure],
        cwd=proj_dir,
        check=True,
    )


def regional_zone_maps(where_are_saved_regional_maps, where_to_save_the_figure, include_station_points=True):
    # Standalone regional-zone-maps figure, kept for Figure_5()'s
    # without-stations reference-map byproduct. The with-stations version is
    # no longer produced standalone -- it is embedded directly in Figure_1().

    the_dates_for_each_zone = dates_for_each_zone()

    load_and_save_regional_maps_for_plot(where_are_saved_regional_maps,
                                                      where_to_save_the_figure,
                                                      the_dates_for_each_zone)

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['regional_zone_maps']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               include_station_points=robjects.BoolVector([include_station_points]))


def Figure_2(where_to_save_the_figure):
    # Manuscript Figure 2: satellite-vs-in-situ validation scatterplots.
    # Combines the two already-rendered, 3x3 grid-size SEXTANT scatterplots
    # (spatially averaged over the 3x3 pixel window around each in situ
    # station -- not the 1x1 single-native-pixel version) from
    # func/validate.py's match-up pipeline (run during 1_validate.py). This
    # function does not regenerate the scatterplots themselves, so tweaking
    # them and re-running 1_validate.py is picked up automatically next time.
    # Panel b changed 2026-08-01 (per Robert) from the Chl-a scatterplot to
    # the SPM-vs-turbidity one: REPHY (the larger in situ network) measures
    # only turbidity, not SPM (only SOMLIT measures SPM directly), so this
    # panel is what actually makes use of the REPHY match-ups -- see
    # Section 2.3.1's discussion of Doxaran et al. (2009) for why comparing
    # satellite SPM against in situ turbidity is a defensible substitution.
    spm_scatterplot_path = os.path.join(where_to_save_the_figure, 'validation', 'scatterplot',
                                        'SEXTANT_SPM_SPM_Standard_OC5_3x3.png')
    turb_scatterplot_path = os.path.join(where_to_save_the_figure, 'validation', 'scatterplot',
                                         'SEXTANT_TUR_SPM_Standard_OC5_3x3.png')

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_2']

    # Call the R function
    r_function(spm_scatterplot_path=robjects.StrVector([spm_scatterplot_path]),
               turb_scatterplot_path=robjects.StrVector([turb_scatterplot_path]),
               where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def Figure_3_panels(where_are_saved_regional_maps, where_to_save_the_figure):
    # Renamed from Figure_4 (2026-08-01): produces the A-E methodology panels
    # for manuscript Figure 3, never manuscript Figure 4 (that name was left
    # over from before the manuscript was renumbered).

    # Static date for each zone to illustrate the plume detection steps
    Zone = 'GULF_OF_LION'
    plume_name = 'Grand Rhone'
    Date = '2014-01-04'

    parameters = build_plume_parameters(Zone)

    path_to_the_satellite_file_to_use = os.path.join(where_are_saved_regional_maps, 'REGIONAL_MAPS', Zone, 'SEXTANT', 'SPM',
                                                     'merged',
                                                     'Standard', 'MAPS', 'DAILY', Date[:4], f'{Date}.pkl')

    # Open and load the file (binary file assumed to contain data)
    with open(path_to_the_satellite_file_to_use, 'rb') as f:
        ds = pickle.load(f)['Basin_map']['map_data']

    # Plume detection now runs at native resolution (panache no longer
    # coarsens internally), so no resolution reduction happens here either.
    ds_reduced = ds

    bathymetry_data_aligned_to_reduced_map = align_bathymetry(ds_reduced,
                                                                            f'{where_are_saved_regional_maps}/REGIONAL_MAPS/{Zone}/Bathy_data.pkl')

    (_, land_mask) = preprocess_annual_dataset_and_compute_land_mask(
        (path_to_fill_to_where_to_save_satellite_files(where_are_saved_regional_maps + "/REGIONAL_MAPS/" + Zone)
         .replace('[DATA_SOURCE]/[PARAMETER]/[SENSOR]/[ATMOSPHERIC_CORRECTION]/[TIME_FREQUENCY]',
                  'SEXTANT/SPM/merged/Standard/MAPS/MULTIYEAR')
         .replace('[YEAR]/[MONTH]/[DAY]', 'Averaged_over_multi-years.pkl')
         ), parameters)

    inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)

    # Panels A-E feed the Figure 3 composite only (methodology figure); they
    # are not manuscript Figure 4 (that's Figure_4_timeseries()'s
    # Figure_4.png, saved into its own FIGURE_4/ folder). Writing straight
    # into FIGURE_3/ avoids the folder-name/manuscript-number mismatch this
    # used to have (2026-08-01) -- same fix already applied to
    # Figure_3_zone_maps()'s zone-maps panel.
    where_to_save_the_figure_4 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURE_3")

    the_plume = Create_the_plume_mask(ds_reduced,
                                      bathymetry_data_aligned_to_reduced_map,
                                      land_mask,
                                      parameters,
                                      plume_name)
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_4,
             name_of_the_plot='A')

    # Real scene-specific dynamic threshold (panache's current
    # determine_SPM_threshold, incl. its near-mouth-quantile bound
    # estimation when minimal/maximal_threshold aren't fixed) -- no longer
    # overridden with a hardcoded placeholder value.
    the_plume.determine_SPM_threshold(dynamic_determination_of_SPM_threshold=True)
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_4,
             name_of_the_plot='B')

    the_plume.do_a_raw_plume_detection()
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_4,
             name_of_the_plot='C')

    the_plume.include_cloudy_regions()

    the_plume.remove_resuspension_areas(
        maximal_bathymetry=parameters['maximal_bathymetric_for_zone_with_resuspension'][plume_name],
        minimal_distance_from_estuary=parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][
            plume_name])
    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
    #                    name_of_the_plot='before_shallow_water_removal')
    ##

    the_plume.remove_shallow_waters()
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_4,
             name_of_the_plot='D')

    the_plume.remove_close_river_mouth(the_plume.parameters['pixel_starting_points_close_river_mouth'])

    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
    #                    name_of_the_plot='before_')
    ##

    the_plume.dilate_and_merge_close_areas()

    the_plume.remove_undersized_shapes()

    the_plume.set_false_outside_search_area(inside_polygon_mask)

    the_plume.select_shape_by_core()

    the_plume.remove_shallow_waters()

    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_4,
    #                    name_of_the_plot='before_shrink_widen')
    ##

    if not np.isin(plume_name, ['Seine']):
        the_plume.remove_post_shrink_widening()

    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_4,
             name_of_the_plot='E')


def Figure_3_zone_maps(where_are_saved_regional_maps, where_to_save_the_figure):
    # Renamed from Figure_5 (2026-08-01): produces zone_maps_panel.png, an
    # input to manuscript Figure 3, never manuscript Figure 5 (that's
    # Figure_5_driver_comparison()'s Figure_5.png).
    the_dates_for_each_zone = dates_for_each_zone()

    # Folded into Figure 3 (methodology composite) rather than a standalone
    # manuscript figure -- this panel is now an input to Figure_3(), not
    # FIGURE_5 (which is manuscript Figure 5, the plume-vs-flow driver
    # comparison produced by Figure_5_driver_comparison()).
    where_to_save_the_figure_5 = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURE_3")
    os.makedirs(os.path.join(where_to_save_the_figure_5, "DATA"), exist_ok=True)

    for Zone, Date in the_dates_for_each_zone.items():

        parameters = build_plume_parameters(Zone)

        path_to_the_satellite_file_to_use = os.path.join(where_are_saved_regional_maps, 'REGIONAL_MAPS', Zone, 'SEXTANT',
                                                         'SPM', 'merged',
                                                         'Standard', 'MAPS', 'DAILY', Date[:4], f'{Date}.pkl')

        # Open and load the file (binary file assumed to contain data)
        with open(path_to_the_satellite_file_to_use, 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data']

            # Plume detection now runs at native resolution (panache no
            # longer coarsens internally), so no resolution reduction
            # happens here either.
        ds_reduced = ds

        bathymetry_data_aligned_to_reduced_map = align_bathymetry(ds_reduced,
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
        # Loop through each plume starting point to process plume detection.
        # Always dynamic threshold (matching Figure_3_panels() and every
        # zone's real metadata/zone_config_dynamic_*.json, all of which set
        # dynamic_threshold: True) -- this used to override several river
        # mouths (Seine, Grand Rhone, Petit Rhone, Loire, Sevre) with fixed
        # thresholds hand-tuned for a since-replaced set of dates, which
        # under-detected the plume on the zones' actual record dynamic-
        # threshold days (Robert, 2026-08-01: Bay of Seine and Gulf of Lion's
        # Figure 3 panels looked implausibly small for their record area).
        for plume_name, starting_point in parameters['starting_points'].items():

            the_plume = delineate_plume_pipeline(ds_reduced,
                                                        bathymetry_data_aligned_to_reduced_map,
                                                        land_mask,
                                                        parameters,
                                                        plume_name,
                                                        inside_polygon_mask,
                                                        dynamic_thresh=True)

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

        SPM_map.to_csv(where_to_save_the_figure_5 + f"/DATA/{Zone}.csv")

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

    r_function = robjects.r['Figure_3_zone_maps']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure_5]))

    regional_zone_maps(where_are_saved_regional_maps="output",
                       where_to_save_the_figure="figures",
                       include_station_points=False)


def Figure_3(where_to_save_the_figure):
    """
    Assembles manuscript Figure 3 (plume-detection methodology): panels A-D
    from Figure_3_panels() side by side on top, the per-zone plume maps panel
    from Figure_3_zone_maps() (zone_maps_panel.png) stacked below. Both write
    their intermediates straight into FIGURE_3/ (neither is a standalone
    manuscript figure itself), so this just reads them back out and
    composites -- must be called after both.

    Panels A-D are resized down to the zone-maps panel's native width
    (preserving aspect ratio) before the horizontal/vertical stack, rather
    than upscaling the zone-maps panel to match a wider top row -- avoids
    blurring the source PNGs.
    """
    figure_3_dir = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURE_3")
    os.makedirs(figure_3_dir, exist_ok=True)

    panel_paths = [os.path.join(figure_3_dir, f"{letter}.png") for letter in "ABCD"]
    zone_maps_path = os.path.join(figure_3_dir, "zone_maps_panel.png")
    missing = [p for p in panel_paths + [zone_maps_path] if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            "Figure_3() inputs missing: " + ", ".join(missing) +
            " -- run Figure_3_panels() and Figure_3_zone_maps() first."
        )

    zone_maps_panel = Image.open(zone_maps_path)
    target_width = zone_maps_panel.width
    tile_width = target_width // 4

    tiles = []
    for p in panel_paths:
        panel = Image.open(p)
        tile_height = round(tile_width * panel.height / panel.width)
        tiles.append(panel.resize((tile_width, tile_height), Image.LANCZOS))

    top_row = Image.new("RGB", (tile_width * 4, tiles[0].height), "white")
    for i, tile in enumerate(tiles):
        top_row.paste(tile, (tile_width * i, 0))

    composite = Image.new("RGB", (target_width, top_row.height + zone_maps_panel.height), "white")
    composite.paste(top_row, (0, 0))
    composite.paste(zone_maps_panel, (0, top_row.height))

    composite.save(os.path.join(figure_3_dir, "Figure_3.png"))


def Figure_5_driver_comparison(where_to_save_the_figure, max_lag_daily=14):
    """manuscript Figure 5: plume-area-vs-flow scatter + lagged correlation, one row per zone.
    Distinct from Figure_5() above (regional zone maps feeding the Figure 3
    composite) -- see func/figure.R::Figure_5_driver_comparison() for the
    naming-collision note. Only wireable here since the tidync/RNetCDF ->
    ncdf4 migration (2026-07-31) fixed an rpy2-embedded-R conflict that used
    to break this call.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_5_driver_comparison']
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               max_lag_daily=robjects.IntVector([max_lag_daily]))


def Figure_7_driver_rose(where_to_save_the_figure, n_sectors=16):
    """manuscript Figure 7: wind/wave direction-magnitude roses, coloured by
    the flow-controlled plume-area response, one row per zone. Replaces the
    original Figure 7/8 concept (X11-decomposed wind/wave magnitude time
    series) -- direction matters as much or more than magnitude for both.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_7_driver_rose']
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               n_sectors=robjects.IntVector([n_sectors]))


def Figure_8_driver_category(where_to_save_the_figure):
    """manuscript Figure 8: flow-controlled plume-area residual vs. wave
    height, coloured by on/off-shore wind category, one panel per zone.
    Generalises the Rhone-only rhone_wind_wave_effect() analysis.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_8_driver_category']
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def Figure_9_gam_partial(where_to_save_the_figure, stats_dir="output/STATS"):
    """manuscript Figure 9: GAM partial-dependence curves for flow, wind,
    wave, and current (tide intentionally excluded from the plot -- it stays
    in the underlying GAM/Table 6 stats, just not visualised), one row per
    zone. Refits func/driver_interactions.R::fit_gam() from the already-saved
    daily_driver_matrix_<zone>.csv (Stage 4 output) rather than a separate
    model or a full pipeline rerun.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_9_gam_partial']
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               stats_dir=robjects.StrVector([stats_dir]))


def Figure_4_S1_timeseries(where_are_saved_plume_results_with_dynamic_threshold,
                           where_are_saved_plume_results_with_fixed_threshold,
                           where_to_save_the_figure):
    # Renamed from Figure_6_7 (2026-08-01): that name matched neither of the
    # two manuscript figures it actually produces (4 and S1) -- left over
    # from before the manuscript was renumbered. Data prep here is shared
    # (both dynamic- and fixed-threshold Results.csv, all zones); the R side
    # is now two separate functions, Figure_4_timeseries() and
    # Figure_S1_thresholds(), each reading the same ts_data.csv this writes.
    figure_4_dir = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURE_4")
    os.makedirs(figure_4_dir + '/DATA', exist_ok=True)

    # Source util.R (via figure.R) to grab the same zone-specific
    # plume_area_ceiling used by load_plume_ts() -- this function reads
    # Results.csv directly with pandas rather than going through
    # load_plume_ts(), so it was never getting the ceiling's spike removal.
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    plume_area_ceiling = dict(zip(robjects.r['plume_area_ceiling'].names,
                                  list(robjects.r['plume_area_ceiling'])))

    ts_files_with_dynamic_threshold = glob.glob(os.path.join(
        where_are_saved_plume_results_with_dynamic_threshold, '*', 'Results.csv'))

    def _load_results(ts_file):
        df = pd.read_csv(ts_file)
        zone = os.path.basename(os.path.dirname(ts_file))
        df['Zone'] = zone
        df['date'] = pd.to_datetime(df['date']).dt.date
        df['Years'] = pd.to_datetime(df['date']).dt.year
        df['Satellite_sensor'] = 'merged'
        df.loc[df['area_of_the_plume_mask_in_km2'] > plume_area_ceiling[zone],
              'area_of_the_plume_mask_in_km2'] = np.nan
        return df

    ts_data_with_dynamic_threshold = pd.concat(
        [_load_results(f) for f in ts_files_with_dynamic_threshold])
    ts_data_with_dynamic_threshold['Dynamic_threshold'] = True

    ts_files_with_fixed_threshold = glob.glob(os.path.join(
        where_are_saved_plume_results_with_fixed_threshold, '*', 'Results.csv'))

    ts_data_with_fixed_threshold = pd.concat(
        [_load_results(f) for f in ts_files_with_fixed_threshold])
    ts_data_with_fixed_threshold['Dynamic_threshold'] = False

    pd.concat([ts_data_with_dynamic_threshold, ts_data_with_fixed_threshold]).to_csv(
        os.path.join(figure_4_dir, 'DATA', 'ts_data.csv'))

    # Both R functions take the top-level "figures" path and build their own
    # ARTICLE/FIGURE_* subfolder internally (matching every other figure
    # function's convention) -- Figure_S1_thresholds() reads the same
    # ts_data.csv from FIGURE_4/DATA/ that Figure_4_timeseries() just wrote.
    robjects.r['Figure_4_timeseries'](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))
    robjects.r['Figure_S1_thresholds'](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def _prep_x11_weekly_data(where_are_saved_X11_results, data_dir):
    """
    Reads the weekly X11-decomposed plume-area/river-flow time series for
    all 4 zones and writes the combined ts_plume_data.csv/ts_river_data.csv
    that figure.R's compute_x11_zone_plots() reads. Shared by the dynamic
    and static calls in Figure_X11_weekly_results() below.
    """
    os.makedirs(os.path.join(data_dir, 'DATA'), exist_ok=True)

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

    pd.concat(ts_plume_data).to_csv(os.path.join(data_dir, 'DATA', 'ts_plume_data.csv'))
    pd.concat(ts_river_data).to_csv(os.path.join(data_dir, 'DATA', 'ts_river_data.csv'))


def Figure_X11_weekly_results(where_are_saved_X11_results_dynamic, where_are_saved_X11_results_static,
                              where_to_save_the_figure):
    """
    Replaces the deprecated Figure_8_9_10()/Figures_8_9_10(), which produced
    3 plots (Seasonal/Interannual/Residual) per threshold but only wired the
    dynamic-threshold ones into the manuscript (Figure 6, Figure S2) --
    its static-threshold counterpart was computed and saved but never cited
    anywhere (2026-08-01 audit). This version wires all 4 real manuscript
    figures: Figure 6 + Figure S2 (dynamic), Figure S5 + Figure S6 (static,
    supplementary robustness check mirroring 6/S2). Renamed 2026-07-31 (S3->S2,
    S6->S5, S7->S6) when the daily-vs-weekly Figure S2 was removed from the
    manuscript and every later supplementary figure renumbered down by one.
    """
    figure_6_dir = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURE_6")
    figure_s5_dir = os.path.join(where_to_save_the_figure, "ARTICLE", "FIGURE_S5")
    _prep_x11_weekly_data(where_are_saved_X11_results_dynamic, figure_6_dir)
    _prep_x11_weekly_data(where_are_saved_X11_results_static, figure_s5_dir)

    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    for r_func_name in ['Figure_6_x11_interannual', 'Figure_S2_x11_components',
                       'Figure_S5_x11_interannual_static', 'Figure_S6_x11_components_static']:
        try:
            robjects.r[r_func_name](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))
        except Exception as e:
            print(f"Warning: {r_func_name} R plot failed: {e}. Skipping.")


def Figure_S3_seasonal_boxplots(where_to_save_the_figure):
    """
    Migrated 2026-08-01 from manuscript/make_figures_tables.R's
    generate_figure_s4_seasonal_thresholds() into the real pipeline, so it
    writes straight to figures/ARTICLE/FIGURE_S3/ instead of via the
    manuscript/figures/ copy step. No Python-side data prep needed -- the R
    function reads output/panache/{dynamic,static}/{zone}/Results.csv directly.
    Renamed 2026-07-31 from Figure_S4_seasonal_boxplots(), see
    Figure_X11_weekly_results()'s note above.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    robjects.r['Figure_S3_seasonal_boxplots'](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))

