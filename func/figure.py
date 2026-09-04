#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import os, sys, re, glob, subprocess, json
import numpy as np
import xarray as xr
import pandas as pd
import rpy2.robjects as robjects
from functools import reduce
from PIL import Image

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (load_csv_files, order_zones, get_registry_row, registry_filename)
from panache.plume_algorithm import (Create_the_plume_mask, delineate_plume_pipeline, create_polygon_mask,
                                     derive_masks_from_bathymetry, estimate_near_mouth_bounds)
from panache.utils import define_parameters, align_bathymetry
from panache.io import load_map_data


# =============================================================================
#### panache bridging functions
# =============================================================================


def build_plume_parameters(Zone):
    """
    panache's authoritative algorithm parameters, overlaid with the
    zone-tuned near-mouth-quantile settings from
    metadata/zone_config_dynamic_<Zone>.json -- so Figure_3_panels()/
    Figure_3_zone_maps() detect plumes exactly as the real pipeline does.

    near_mouth_lower_quantile/near_mouth_upper_quantile/gradient_steepness_fraction
    are set only in that JSON, read here and overlaid last -- panache's
    define_parameters() does not set them at all, so without this they
    silently fell back to panache's hardcoded defaults (0.25/0.75/0.9)
    instead of each zone's tuned values, which only the operational panache
    CLI path (panache.runner) was reading.
    """
    parameters = define_parameters(Zone)

    zone_config_path = os.path.join(proj_dir, "metadata", f"zone_config_dynamic_{Zone}.json")
    with open(zone_config_path) as f:
        zone_config = json.load(f)
    for key in ("near_mouth_lower_quantile", "near_mouth_upper_quantile", "gradient_steepness_fraction"):
        if key in zone_config:
            parameters[key] = zone_config[key]

    return parameters


# =============================================================================
#### Utility functions
# =============================================================================


def do_R_plot(the_plume, where_to_save_the_plot, name_of_the_plot):
    """
    Convert a panache Create_the_plume_mask instance's SPM map and plume mask
    to a CSV for plot_methodology_worked_example_panel() (func/figure.R) to plot later. The R call
    itself is deferred to a single batched Rscript subprocess in
    Figure_3_panels() (run_figure_3_panels.R) rather than done here
    in-process, since plot_methodology_worked_example_panel() now renders a high-res coastline via
    sf, which conflicts with the conda geospatial stack already loaded in
    this process via panache -- same workaround as Figure_1().
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

    # Same shallow-water criterion remove_shallow_waters() applies to the
    # plume mask: lets plot_methodology_worked_example_panel() draw
    # the bathymetric exclusion boundary actually used to produce this
    # panel. bathymetric_threshold is 0 at the three zones that don't use
    # this general exclusion (only the Gulf of Lion does), so this column is
    # all-False and draws nothing there.
    bathymetric_threshold = the_plume.parameters.get('bathymetric_threshold', 0)
    df['shallow'] = the_plume.bathymetry_map.values.flatten() > -bathymetric_threshold

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

        # Export the SPM value found at every point tested along the
        # search-direction transects, plus the near-mouth min/max bounds and
        # the final derived SPM_threshold, so plot_methodology_worked_example_panel() can draw a
        # panel showing how the gradient-cutoff filtering (kept vs. rejected
        # transect points) and the near-mouth quantile bounds combine to
        # pick the plume-edge threshold.
        # minimal_threshold/maximal_threshold are not stored on `the_plume`
        # by determine_SPM_threshold(), so they're recomputed here the same
        # way panache computes them internally (deterministic given the same
        # scene/parameters) rather than duplicated inside panache itself.
        plume_name = the_plume.plume_name
        lower_q = the_plume.parameters.get('near_mouth_lower_quantile', 0.25)
        upper_q = the_plume.parameters.get('near_mouth_upper_quantile', 0.75)
        minimal_threshold, maximal_threshold = estimate_near_mouth_bounds(
            the_plume.spm_map.values, the_plume.land_mask.values,
            the_plume.parameters['pixel_starting_points'][plume_name],
            lower_quantile=lower_q, upper_quantile=upper_q,
            radius_pixels=the_plume.parameters.get('near_mouth_radius_pixels', 15),
            window_shape=the_plume.parameters.get('near_mouth_window_shape', 'circle'),
        )

        start_y, start_x = the_plume.parameters['pixel_starting_points'][plume_name]
        all_pts = the_plume.all_points_tested_for_finding_SPM_threshold
        x_idx = all_pts[:, 0].astype(int)
        y_idx = all_pts[:, 1].astype(int)

        # Pixel-index Euclidean distance from the river mouth, converted to
        # km via the scene's own lat/lon grid spacing.
        lat_res_km = abs(lat_values[1] - lat_values[0]) * 111.0
        lon_res_km = abs(lon_values[1] - lon_values[0]) * 111.0 * np.cos(np.deg2rad(np.mean(lat_values)))
        distance_km = np.hypot((x_idx - start_x) * lon_res_km, (y_idx - start_y) * lat_res_km)

        kept_set = set(map(tuple, the_plume.points_used_for_finding_SPM_threshold.astype(int)))
        kept = np.array([(x, y) in kept_set for x, y in zip(x_idx, y_idx)])

        transect_values = pd.DataFrame({
            'longitude': lon_values[x_idx],
            'latitude': lat_values[y_idx],
            'distance_km': distance_km,
            'analysed_spim': the_plume.spm_map.values[y_idx, x_idx],
            'kept': kept,
        })
        transect_values.to_csv(os.path.join(folder_where_to_save_data, f"{name_of_the_plot}_transect_values.csv"))

        threshold_values = pd.DataFrame({
            'minimal_threshold': [minimal_threshold],
            'maximal_threshold': [maximal_threshold],
            'SPM_threshold': [the_plume.SPM_threshold],
            'gradient_steepness_fraction': [the_plume.parameters.get('gradient_steepness_fraction', 0.9)],
            'quantile_to_use': [the_plume.parameters['quantile_to_use'][plume_name]],
        })
        threshold_values.to_csv(os.path.join(folder_where_to_save_data, f"{name_of_the_plot}_threshold_values.csv"))


    if 'plume_mask' in vars(the_plume) :
        df['plume'] = the_plume.plume_mask.values.flatten()

    index_to_keep = np.where((df.lat >= the_plume.parameters['lat_range_of_plume_area'][0]) &
                             (df.lat <= the_plume.parameters['lat_range_of_plume_area'][1]) &
                             (df.lon >= the_plume.parameters['lon_range_of_plume_area'][0]) &
                             (df.lon <= the_plume.parameters['lon_range_of_plume_area'][1]) )

    df = df.iloc[index_to_keep]

    df.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}.csv") )


def save_files_for_Figure_1(where_are_saved_satellite_data, where_to_save_the_figure,
                            mean_spm_path, coordinates_of_the_map) :

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
    
    
def load_and_save_regional_maps_for_plot(where_to_save_the_figure, dates_for_each_zone) :

    folder_where_to_save_regional_zone_maps_data = os.path.join(where_to_save_the_figure, 'ARTICLE', 'FIGURE_1', 'DATA')
    os.makedirs(folder_where_to_save_regional_zone_maps_data, exist_ok = True)

    for key, date in dates_for_each_zone.items() :

        coordinates_of_the_map = define_parameters(key)

        zone_config_path = os.path.join(proj_dir, "metadata", f"zone_config_dynamic_{key}.json")
        with open(zone_config_path) as f:
            zone_config = json.load(f)

        date_ts = pd.Timestamp(date)
        nc_path = glob.glob(os.path.join(zone_config['input_path'], f'{date_ts.year}', f'{date_ts.month:02d}', f'{date_ts.day:02d}', '*.nc'))[0]

        SPM_map_da = load_map_data(nc_path,
                                   lon_range=tuple(coordinates_of_the_map['lon_range_of_plume_area']),
                                   lat_range=tuple(coordinates_of_the_map['lat_range_of_plume_area']),
                                   variable_name=zone_config['variable_name'])

        SPM_map = SPM_map_da.to_dataframe().reset_index()

        SPM_map.to_csv(folder_where_to_save_regional_zone_maps_data + f"/{key}.csv")

    save_insitu_stations_for_plot(folder_where_to_save_regional_zone_maps_data)


def save_insitu_stations_for_plot(folder_where_to_save_Figure_data) :     

    coordinates_of_the_RIOMARS = { zone_name : {'lat_min' : define_parameters(zone_name)['lat_range_of_plume_area'][0],
                                                'lat_max' : define_parameters(zone_name)['lat_range_of_plume_area'][1],
                                                'lon_min' : define_parameters(zone_name)['lon_range_of_plume_area'][0],
                                                'lon_max' : define_parameters(zone_name)['lon_range_of_plume_area'][1]}
                                  for zone_name in order_zones(['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY', 'BAY_OF_SEINE']) }
    
    pd.DataFrame.from_dict(coordinates_of_the_RIOMARS).to_csv(folder_where_to_save_Figure_data + "/RIOMAR_limits.csv")
    
    # Station coordinates -- written by func/validate.R
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

    dates = {'GULF_OF_LION' : '2025-01-11',
            'BAY_OF_BISCAY' : '2025-01-11',
            'SOUTHERN_BRITTANY' : '2025-01-11',
            'BAY_OF_SEINE' : '2025-01-11'}
    return {zone : dates[zone] for zone in order_zones(dates.keys())}


# =============================================================================
#### Main functions
# =============================================================================


def Figure_1(where_are_saved_satellite_data, where_to_save_the_figure):
    save_files_for_Figure_1(where_are_saved_satellite_data,
                            where_to_save_the_figure,
                            mean_spm_path=os.path.join(proj_dir, "output", "STATS", "mean_spm_national.nc"),
                            coordinates_of_the_map={"lat_min": 42, "lat_max": 51.5, "lon_min": -6, "lon_max": 8})

    # The four regional zoomed insets are built inside Figure_1() in figure.R.
    load_and_save_regional_maps_for_plot(where_to_save_the_figure,
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


def regional_zone_maps(where_to_save_the_figure, include_station_points=True):

    the_dates_for_each_zone = dates_for_each_zone()

    load_and_save_regional_maps_for_plot(where_to_save_the_figure,
                                         the_dates_for_each_zone)

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['regional_zone_maps']

    # Call the R function
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               include_station_points=robjects.BoolVector([include_station_points]))


def Figure_2(where_to_save_the_figure):
    # Manuscript slot "validation_scatterplot_panel" -- see
    # manuscript/figure_table_registry.csv (util.get_registry_row()) for its
    # current figure number/output folder and the R function that renders it.
    # This Python entry point's own name is an unchanging handle and is not
    # kept in sync with the manuscript number.

    spm_scatterplot_path = os.path.join(where_to_save_the_figure, 'validation', 'scatterplot',
                                        'SEXTANT_SPM_SPM_Standard_OC5_3x3.png')
    turb_scatterplot_path = os.path.join(where_to_save_the_figure, 'validation', 'scatterplot',
                                         'SEXTANT_TUR_SPM_Standard_OC5_3x3.png')

    registry_row = get_registry_row('validation_scatterplot_panel')

    # Source the R script
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r[registry_row['r_function']]

    # Call the R function
    r_function(spm_scatterplot_path=robjects.StrVector([spm_scatterplot_path]),
               turb_scatterplot_path=robjects.StrVector([turb_scatterplot_path]),
               where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def Figure_3_panels(where_are_saved_regional_maps, where_to_save_the_figure):
    # Produces the A-E methodology panels

    # Static date for each zone to illustrate the plume detection steps
    Zone = 'GULF_OF_LION'
    plume_name = 'Grand Rhone'
    Date = '2025-01-11'

    parameters = build_plume_parameters(Zone)

    zone_config_path = os.path.join(proj_dir, "metadata", f"zone_config_dynamic_{Zone}.json")
    with open(zone_config_path) as f:
        zone_config = json.load(f)

    date_ts = pd.Timestamp(Date)
    nc_path = glob.glob(os.path.join(zone_config['input_path'], f'{date_ts.year}', f'{date_ts.month:02d}', f'{date_ts.day:02d}', '*.nc'))[0]

    # Crop to panache's own plume-area bbox (not regmap.py's wider
    # Basin_limits) -- every bbox used anywhere in this pipeline comes from
    # panache, no exceptions.
    ds = load_map_data(nc_path,
                       lon_range=tuple(parameters['lon_range_of_plume_area']),
                       lat_range=tuple(parameters['lat_range_of_plume_area']),
                       variable_name=zone_config['variable_name'])

    # Plume detection now runs at native resolution (panache no longer
    # coarsens internally), so no resolution reduction happens here either.
    ds_reduced = ds

    bathymetry_data_aligned_to_reduced_map = align_bathymetry(ds_reduced, f'{where_are_saved_regional_maps}/REGIONAL_MAPS/{Zone}/Bathy_data.pkl')

    (_, land_mask) = derive_masks_from_bathymetry(bathymetry_data_aligned_to_reduced_map, parameters)

    inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)

    # Panels A-E feed the plume_methodology_panel composite -- see
    # manuscript/figure_table_registry.csv for its current figure number.
    where_to_save_the_figure_3 = os.path.join(
        where_to_save_the_figure, "ARTICLE", get_registry_row("plume_methodology_panel")['output_subdir'])

    the_plume = Create_the_plume_mask(ds_reduced,
                                      bathymetry_data_aligned_to_reduced_map,
                                      land_mask,
                                      parameters,
                                      plume_name)
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_3,
              name_of_the_plot='A')

    the_plume.determine_SPM_threshold(dynamic_determination_of_SPM_threshold=True)
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_3,
             name_of_the_plot='B')

    the_plume.do_a_raw_plume_detection()
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_3,
             name_of_the_plot='C')

    the_plume.include_cloudy_regions()

    the_plume.remove_resuspension_areas(
        maximal_bathymetry=parameters['maximal_bathymetric_for_zone_with_resuspension'][plume_name],
        minimal_distance_from_estuary=parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][
            plume_name])
    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_3,
    #                    name_of_the_plot='before_shallow_water_removal')
    ##

    the_plume.remove_shallow_waters()
    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_3,
             name_of_the_plot='D')

    the_plume.remove_close_river_mouth(the_plume.parameters['pixel_starting_points_close_river_mouth'])

    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_3,
    #                    name_of_the_plot='before_')
    ##

    the_plume.dilate_and_merge_close_areas()

    the_plume.remove_undersized_shapes()

    the_plume.set_false_outside_search_area(inside_polygon_mask)

    the_plume.select_shape_by_core()

    the_plume.remove_shallow_waters()

    ##
    # the_plume.do_R_plot(where_to_save_the_plot=where_to_save_the_figure_3,
    #                    name_of_the_plot='before_shrink_widen')
    ##

    if not np.isin(plume_name, ['Seine']):
        the_plume.remove_post_shrink_widening()

    do_R_plot(the_plume, where_to_save_the_plot=where_to_save_the_figure_3,
              name_of_the_plot='E')

    # Run as a standalone Rscript process rather than via in-process rpy2:
    # plot_methodology_worked_example_panel() (func/figure.R) now renders a high-res coastline via
    # sf, which conflicts with the conda geospatial stack already loaded in
    # this process via panache -- same workaround as Figure_1().
    subprocess.run(
        ['Rscript', 'func/run_figure_3_panels.R', where_to_save_the_figure_3],
        cwd=proj_dir,
        check=True,
    )


def Figure_3_zone_maps(where_are_saved_regional_maps, where_to_save_the_figure):

    the_dates_for_each_zone = dates_for_each_zone()

    # Folded into the plume_methodology_panel composite

    where_to_save_the_figure_3 = os.path.join(
        where_to_save_the_figure, "ARTICLE", get_registry_row("plume_methodology_panel")['output_subdir'])
    os.makedirs(os.path.join(where_to_save_the_figure_3, "DATA"), exist_ok=True)

    for Zone, Date in the_dates_for_each_zone.items():

        parameters = define_parameters(Zone)

        zone_config_path = os.path.join(proj_dir, "metadata", f"zone_config_dynamic_{Zone}.json")
        with open(zone_config_path) as f:
            zone_config = json.load(f)

        date = pd.Timestamp(Date)
        nc_path = glob.glob(os.path.join(zone_config['input_path'], f'{date.year}', f'{date.month:02d}', f'{date.day:02d}', '*.nc'))[0]

        SPM_map_da = load_map_data(nc_path,
                                   lon_range=tuple(parameters['lon_range_of_plume_area']),
                                   lat_range=tuple(parameters['lat_range_of_plume_area']),
                                   variable_name=zone_config['variable_name'])

        plume_masks_path = os.path.join(where_are_saved_regional_maps, 'panache', 'dynamic', Zone, 'PlumeMasks.nc')
        with xr.open_dataset(plume_masks_path) as mask_ds:
            plume_mask = (mask_ds['plume_mask']
                          .sel(time=Date, river='ALL')
                          .reindex(lat=SPM_map_da.lat, lon=SPM_map_da.lon, method='nearest')
                          .load())

        SPM_map = SPM_map_da.to_dataframe().reset_index()
        SPM_map['plume'] = plume_mask.values.astype(bool).flatten()

        # Lets plot_methodology_zone_maps_panel() draw the bathymetric
        # exclusion boundary that actually matches panache's algorithm at
        # each zone. `bathymetric_threshold` (remove_shallow_waters()'s
        # general exclusion) is 0 at three of the four zones, so it would
        # draw nothing there -- instead this uses
        # maximal_bathymetric_for_zone_with_resuspension, the depth
        # remove_resuspension_pixels() actually applies at every zone
        # (uniform across all rivers within a given zone: Seine 20 m,
        # Gironde/Charente/Sevre 10 m, Loire/Vilaine 10 m, Grand/Petit
        # Rhone 20 m), so any one river's value represents its zone.
        bathymetry_map = align_bathymetry(SPM_map_da, f'{where_are_saved_regional_maps}/REGIONAL_MAPS/{Zone}/Bathy_data.pkl')
        resuspension_threshold = next(iter(parameters['maximal_bathymetric_for_zone_with_resuspension'].values()))
        SPM_map['shallow'] = (bathymetry_map.values.flatten() > -resuspension_threshold)

        SPM_map.to_csv(where_to_save_the_figure_3 + f"/DATA/{Zone}.csv")

    # Run as a standalone Rscript process rather than via in-process rpy2:
    # plot_methodology_zone_maps_panel() (func/figure.R) now renders a high-res coastline
    # via sf, which conflicts with the conda geospatial stack already loaded
    # in this process via panache -- same workaround as Figure_1().
    subprocess.run(
        ['Rscript', 'func/run_figure_3_zone_maps.R', where_to_save_the_figure_3],
        cwd=proj_dir,
        check=True,
    )

    regional_zone_maps(where_to_save_the_figure="figures",
                       include_station_points=False)


def Figure_3(where_to_save_the_figure):
    """
    Assembles the plume_methodology_panel figure (plume-detection
    methodology): panels A-D from Figure_3_panels() side by side on top,
    panel e) (transect_panel.png,
    also from Figure_3_panels()) as a full-width row in the middle, and the
    per-zone plume maps panel from Figure_3_zone_maps() (zone_maps_panel.png,
    panels f-i) stacked below. All three write their intermediates straight
    into the plume_methodology_panel slot's output folder (none is a
    standalone manuscript figure itself), so this just reads them back out
    and composites -- must be called after both Figure_3_panels() and
    Figure_3_zone_maps(). See manuscript/figure_table_registry.csv for the
    slot's current figure number.

    Panels A-D and e) are resized down to the zone-maps panel's native width
    (preserving aspect ratio) before the horizontal/vertical stack, rather
    than upscaling the zone-maps panel to match a wider top row -- avoids
    blurring the source PNGs.
    """
    output_subdir = get_registry_row("plume_methodology_panel")['output_subdir']
    figure_3_dir = os.path.join(where_to_save_the_figure, "ARTICLE", output_subdir)
    os.makedirs(figure_3_dir, exist_ok=True)

    panel_paths = [os.path.join(figure_3_dir, f"{letter}.png") for letter in "ABCD"]
    transect_panel_path = os.path.join(figure_3_dir, "transect_panel.png")
    zone_maps_path = os.path.join(figure_3_dir, "zone_maps_panel.png")
    missing = [p for p in panel_paths + [transect_panel_path, zone_maps_path] if not os.path.exists(p)]
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

    transect_panel = Image.open(transect_panel_path)
    transect_row_height = round(target_width * transect_panel.height / transect_panel.width)
    transect_row = transect_panel.resize((target_width, transect_row_height), Image.LANCZOS)

    composite = Image.new("RGB", (target_width, top_row.height + transect_row.height + zone_maps_panel.height), "white")
    composite.paste(top_row, (0, 0))
    composite.paste(transect_row, (0, top_row.height))
    composite.paste(zone_maps_panel, (0, top_row.height + transect_row.height))

    composite.save(os.path.join(figure_3_dir, registry_filename(output_subdir)))


def Figure_4_S1_timeseries(where_are_saved_plume_results_with_dynamic_threshold,
                           where_are_saved_plume_results_with_fixed_threshold,
                           where_to_save_the_figure):

    plume_area_timeseries_row = get_registry_row("plume_area_timeseries")
    figure_4_dir = os.path.join(where_to_save_the_figure, "ARTICLE", plume_area_timeseries_row['output_subdir'])
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
        # panache v5.0.0+ writes one row per individual river mouth plus a
        # combined 'ALL' row per date -- keep only the zone-level union.
        df = df[df['river'] == 'ALL']
        zone = os.path.basename(os.path.dirname(ts_file))
        df['date'] = pd.to_datetime(df['date']).dt.date
        # Null every measurement column (not just area) on a screened day --
        # mirrors util.R::load_plume_ts()'s plume_area_ceiling_exceeded fix,
        # extended here now that the plume_area_timeseries figure also plots
        # SPM mass, which was previously left unscreened on the same anomalous days.
        ceiling_exceeded = df['area_of_the_plume_mask_in_km2'] > plume_area_ceiling[zone]
        df.loc[ceiling_exceeded, df.columns.difference(['date'])] = np.nan

        # Fill gaps in the date range with explicit NaN rows
        # so the plume_area_timeseries figure plots real data gaps (cloud cover, etc.)
        # as breaks rather than geom_path/geom_point silently drawing a
        # straight line across them. This reads Results.csv directly rather
        # than going through util.R::load_plume_ts() (which already does
        # this gap-fill), specifically to keep the ceiling-based spike
        # removal above -- so the gap-fill needs redoing here too.
        df = df.set_index(pd.DatetimeIndex(pd.to_datetime(df['date'])))
        full_range = pd.date_range(df.index.min(), df.index.max(), freq='D')
        df = df.reindex(full_range)
        df['date'] = df.index.date

        df['Zone'] = zone
        df['Years'] = pd.to_datetime(df['date']).dt.year
        df['Satellite_sensor'] = 'merged'
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
    # function's convention, looked up from the registry) -- the
    # thresholds_comparison slot's R function reads the same ts_data.csv from
    # this slot's DATA/ that plot_plume_area_timeseries() just wrote.
    robjects.r[plume_area_timeseries_row['r_function']](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))
    robjects.r[get_registry_row("thresholds_comparison")['r_function']](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def Figure_5_seasonal_analysis(where_are_saved_plume_results_with_dynamic_threshold,
                               where_are_saved_plume_results_with_static_threshold,
                               where_to_save_the_figure):
    """The seasonal_boxplot_heatmap figure (sec:results_seasonal, see
    manuscript/figure_table_registry.csv for its current figure number):
    heatmap of monthly median values (% of each zone's own observed range)
    for all four plume properties and five drivers, per zone, dynamic
    threshold -- one small zone x month heatmap per variable, 9 in total in a
    single 3x3 grid.
    All data loading (incl. the along-coast PCA projection, which needs
    func/multi.R) and plotting happens in R.
    See func/figure.R::plot_seasonal_boxplot_heatmap()
    following the same no-Python-prep pattern already used by
    Figure_S3_seasonal_boxplots() below. Writes a shared
    <output_subdir>/DATA/monthly_boxplot_data.csv (both thresholds,
    full daily-level detail, not just the medians plotted here) that the
    updated Figure_S3_seasonal_boxplots() reads back in, so the dynamic/
    static computation is only done once. Writes its PNG directly (a
    single figure now, no Python-side compositing needed).
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    robjects.r[get_registry_row("seasonal_boxplot_heatmap")['r_function']](
        where_are_saved_plume_results_with_dynamic_threshold=robjects.StrVector([where_are_saved_plume_results_with_dynamic_threshold]),
        where_are_saved_plume_results_with_static_threshold=robjects.StrVector([where_are_saved_plume_results_with_static_threshold]),
        where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


def Figure_S_daily_flow(where_to_save_the_figure, max_lag_daily=14):
    """Supplementary "Sx. Lagged daily correlations" figure (fig:daily_flow):
    daily plume area vs. river flow scatter + lagged correlation, per zone.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r[get_registry_row("daily_flow_lagged_correlation")['r_function']]
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               max_lag_daily=robjects.IntVector([max_lag_daily]))


def Figure_7_driver_rose(where_to_save_the_figure, n_sectors=8):
    """manuscript figure: wind/wave direction-magnitude roses, coloured by
    the flow-controlled plume-area response, one row per zone. See
    manuscript/figure_table_registry.csv (slot "driver_rose_diagram") for
    its current figure number.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r[get_registry_row("driver_rose_diagram")['r_function']]
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               n_sectors=robjects.IntVector([n_sectors]))


def Figure_8_gam_partial(where_to_save_the_figure, stats_dir="output/STATS"):
    """GAM partial-dependence curves for flow, wind, wave, and current (tide
    intentionally excluded from the plot -- it stays in the underlying
    GAM/driver_stats_table stats, just not visualised), one row per zone. Refits
    func/driver_interactions.R::fit_gam() from the already-saved
    daily_driver_matrix_<zone>.csv (Stage 4 output) rather than a separate
    model or a full pipeline rerun. See manuscript/figure_table_registry.csv
    (slot "gam_partial_effects") for its current figure number.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r[get_registry_row("gam_partial_effects")['r_function']]
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]),
               stats_dir=robjects.StrVector([stats_dir]))


def _prep_x11_weekly_data(where_are_saved_X11_results, data_dir):
    """
    Reads the weekly X11-decomposed plume-area/river-flow time series for
    all 4 zones and writes the combined ts_plume_data.csv/ts_river_data.csv
    that figure.R's compute_x11_zone_plots() reads. Shared by the dynamic
    and static calls in Figure_X11_weekly_results() below.
    """
    os.makedirs(os.path.join(data_dir, 'DATA'), exist_ok=True)

    regions = ["BAY_OF_BISCAY", "GULF_OF_LION", "BAY_OF_SEINE", "SOUTHERN_BRITTANY"]

    # Globbed per exact zone directory (not '*') so the per-river X11 output
    # directories (e.g. BAY_OF_SEINE_Seine) added alongside the zone-level
    # ones don't get pulled in here too -- a '*' wildcard + substring region
    # match previously matched both, duplicating dates per zone.
    ts_plume_data = []
    for region in regions:
        for ts_file in glob.glob(
                os.path.join(where_are_saved_X11_results, region, 'X11_ANALYSIS', 'area_of_the_plume_mask_in_km2',
                             'SEXTANT_merged_Standard_WEEKLY.csv')):
            ts_data = pd.read_csv(ts_file)
            ts_data['Zone'] = region
            ts_plume_data.append(ts_data)

    ts_river_data = []
    for region in regions:
        for ts_file in glob.glob(
                os.path.join(where_are_saved_X11_results, region, 'X11_ANALYSIS', 'river_flow',
                             'River_flow___WEEKLY.csv')):
            ts_data = pd.read_csv(ts_file)
            ts_data['Zone'] = region
            ts_river_data.append(ts_data)

    pd.concat(ts_plume_data).to_csv(os.path.join(data_dir, 'DATA', 'ts_plume_data.csv'))
    pd.concat(ts_river_data).to_csv(os.path.join(data_dir, 'DATA', 'ts_river_data.csv'))


def _prep_x11_dynamic_vs_static_data(where_are_saved_X11_results_dynamic, where_are_saved_X11_results_static, data_dir):
    """
    Reads the weekly X11-decomposed plume-area time series for all 4 zones
    under both detection thresholds and writes the combined
    ts_plume_dynamic_vs_static.csv that figure.R's
    compute_x11_dynamic_vs_static_plots() reads -- pairs the SAME variable
    (plume area) under the two thresholds, rather than plume area against
    river flow (which _prep_x11_weekly_data() above does for the
    x11_interannual_river_flow/x11_components_dynamic figures).
    """
    os.makedirs(os.path.join(data_dir, 'DATA'), exist_ok=True)
    regions = ["BAY_OF_BISCAY", "GULF_OF_LION", "BAY_OF_SEINE", "SOUTHERN_BRITTANY"]

    def _load_plume(where_are_saved_X11_results, threshold_label):
        dfs = []
        for region in regions:
            for ts_file in glob.glob(
                    os.path.join(where_are_saved_X11_results, region, 'X11_ANALYSIS',
                                 'area_of_the_plume_mask_in_km2', 'SEXTANT_merged_Standard_WEEKLY.csv')):
                ts_data = pd.read_csv(ts_file)
                ts_data['Zone'] = region
                ts_data['threshold'] = threshold_label
                dfs.append(ts_data)
        return pd.concat(dfs)

    dynamic_df = _load_plume(where_are_saved_X11_results_dynamic, 'dynamic')
    static_df = _load_plume(where_are_saved_X11_results_static, 'static')

    pd.concat([dynamic_df, static_df]).to_csv(os.path.join(data_dir, 'DATA', 'ts_plume_dynamic_vs_static.csv'))


def Figure_X11_weekly_results(where_are_saved_X11_results_dynamic, where_are_saved_X11_results_static,
                              where_to_save_the_figure):
    """
    Wires all 6 real manuscript figures (slot_key -> current number, see
    manuscript/figure_table_registry.csv): x11_interannual_river_flow +
    x11_seasonal_river_flow + x11_residual_river_flow (dynamic threshold,
    plume area vs. river flow, split into interannual/seasonal/residual --
    each its own standalone 4-zone-panel figure, never combined into one
    image), x11_interannual_dynamic_vs_static + x11_seasonal_dynamic_vs_static
    + x11_residual_dynamic_vs_static (interannual/seasonal/residual, dynamic
    vs. static threshold comparison of plume area itself).
    """
    x11_river_flow_subdir = get_registry_row("x11_interannual_river_flow")['output_subdir']
    x11_dyn_vs_static_subdir = get_registry_row("x11_interannual_dynamic_vs_static")['output_subdir']
    x11_river_flow_dir = os.path.join(where_to_save_the_figure, "ARTICLE", x11_river_flow_subdir)
    figure_s5_dir = os.path.join(where_to_save_the_figure, "ARTICLE", x11_dyn_vs_static_subdir)
    _prep_x11_weekly_data(where_are_saved_X11_results_dynamic, x11_river_flow_dir)
    _prep_x11_dynamic_vs_static_data(where_are_saved_X11_results_dynamic, where_are_saved_X11_results_static, figure_s5_dir)

    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    for slot_key in ['x11_interannual_river_flow', 'x11_seasonal_river_flow',
                     'x11_residual_river_flow',
                     'x11_interannual_dynamic_vs_static',
                     'x11_seasonal_dynamic_vs_static',
                     'x11_residual_dynamic_vs_static']:
        r_func_name = get_registry_row(slot_key)['r_function']
        try:
            robjects.r[r_func_name](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))
        except Exception as e:
            print(f"Warning: {r_func_name} ({slot_key}) R plot failed: {e}. Skipping.")


def Figure_S3_seasonal_boxplots(where_to_save_the_figure):
    """
    Migrated from manuscript/make_figures_tables.R's
    generate_figure_s4_seasonal_thresholds() into the real pipeline, so it
    writes straight to the seasonal_boxplots_dynamic_vs_static slot's output
    folder (see manuscript/figure_table_registry.csv) instead of via the
    manuscript/figures/ copy step. No Python-side data prep needed -- the R
    function reads output/panache/{dynamic,static}/{zone}/Results.csv directly.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)
    robjects.r[get_registry_row("seasonal_boxplots_dynamic_vs_static")['r_function']](where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))


# =============================================================================
#### Deprecated functions
# =============================================================================

def Figure_8_driver_category(where_to_save_the_figure):
    """manuscript Figure 8: flow-controlled plume-area residual vs. wave
    height, coloured by on/off-shore wind category, one panel per zone.
    Generalises the Rhone-only rhone_wind_wave_effect() analysis.
    """
    figure_R_path = os.path.join(func_dir, 'figure.R')
    robjects.r['source'](figure_R_path)

    r_function = robjects.r['Figure_8_driver_category']
    r_function(where_to_save_the_figure=robjects.StrVector([where_to_save_the_figure]))

