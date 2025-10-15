#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================


import os, sys, pickle, glob, multiprocess, imageio, re, gc
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import rpy2.robjects as robjects
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.ndimage import label, binary_dilation, center_of_mass, distance_transform_edt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from collections import deque
from shapely.geometry import Polygon
from shapely.vectorized import contains
from skimage import morphology
from functools import reduce

multiprocess.set_start_method('spawn', force = True) # MacOS friendly multiprocessing

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (load_file, coordinates_of_pixels_to_inspect, define_parameters, 
                  align_bathymetry_to_resolution, unique_years_between_two_dates, 
                  path_to_fill_to_where_to_save_satellite_files, load_shapefile_data,
                  fill_the_sat_paths, get_all_cases_to_process_for_regional_maps_or_plumes_or_X11)


# =============================================================================
#### Utility function
# =============================================================================


def reduce_resolution(ds, lat_bin_size_in_degree, lon_bin_size_in_degree):
        
    """
    Reduce the spatial resolution of a dataset by aggregating it into larger bins.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset with latitude (`lat`) and longitude (`lon`) coordinates.
    lat_bin_size_in_degree : float
        The desired bin size in degrees for the latitude dimension.
    lon_bin_size_in_degree : float
        The desired bin size in degrees for the longitude dimension.

    Returns
    -------
    xarray.Dataset
        A new dataset with reduced resolution, aggregated over the specified bin sizes.

    Notes
    -----
    The function calculates the average spacing of latitude and longitude points in
    the dataset and determines the aggregation factor required to achieve the specified
    bin sizes. The `xarray.Dataset.coarsen` method is used to perform the aggregation.

    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> # Create a sample dataset
    >>> lat = np.arange(-90, 91, 1)
    >>> lon = np.arange(-180, 181, 1)
    >>> data = np.random.rand(len(lat), len(lon))
    >>> ds = xr.Dataset({"values": (["lat", "lon"], data)}, coords={"lat": lat, "lon": lon})
    >>> # Reduce the resolution
    >>> ds_reduced = reduce_resolution(ds, lat_bin_size_in_degree=5, lon_bin_size_in_degree=5)
    >>> ds_reduced
    <xarray.Dataset>
    Dimensions:  (lat: 36, lon: 72)
    Coordinates:
      * lat      (lat) float64 -87.5 -82.5 -77.5 ... 77.5 82.5 87.5
      * lon      (lon) float64 -177.5 -172.5 -167.5 ... 172.5 177.5
    Data variables:
        values   (lat, lon) float64 ...
    """
    
    diff_lat = np.diff(np.unique(ds.lat)).mean()
    diff_lon = np.diff(np.unique(ds.lon)).mean()
    
    lat_factor = round( lat_bin_size_in_degree / diff_lat )
    lon_factor = round( lon_bin_size_in_degree / diff_lon )
    
    ds_reduced = ds.coarsen(lat=lat_factor, lon=lon_factor, boundary='trim').mean()
    
    return ds_reduced


def flood_fill(data, start, SPM_threshold, directions):
    
    """
    Perform a flood fill on a 2D array, starting from a given point and marking cells
    based on a specified threshold and direction vectors.

    Parameters
    ----------
    data : numpy.ndarray
        A 2D array representing the input data (e.g., satellite or environmental data).
    start : tuple of int
        The starting position for the flood fill as (row, column).
    SPM_threshold : float
        Threshold value. Cells with values greater than or equal to this threshold
        will be marked as part of the mask.
    directions : list of tuple of int
        List of direction vectors (e.g., [(0, 1), (1, 0), ...]) that determine how
        neighboring cells are inspected during the flood fill process.

    Returns
    -------
    mask : numpy.ndarray
        A boolean 2D array with `True` indicating cells that meet the threshold and
        are part of the flood-filled region.
    done_pixels : numpy.ndarray
        A boolean 2D array with `True` indicating all cells that were visited during
        the flood fill process.

    Notes
    -----
    - The flood fill checks neighboring cells based on the specified directions.
    - Handles special cases where neighboring cells are NaN by searching for the
      closest finite values.
    - Ensures continuity in marking when there are gaps exceeding the threshold
      within the same direction.

    Examples
    --------
    >>> import numpy as np
    >>> from collections import deque
    >>> # Create sample data
    >>> data = np.array([
    ...     [0, 0, 3, 4, 5],
    ...     [0, 2, 4, 5, 6],
    ...     [1, 3, 5, 7, 8],
    ...     [0, 1, 2, 3, 4]
    ... ])
    >>> start = (2, 2)
    >>> SPM_threshold = 4
    >>> directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    >>> mask, done_pixels = flood_fill(data, start, SPM_threshold, directions)
    >>> print(mask)
    [[False False False  True  True]
     [False False  True  True  True]
     [False False  True  True  True]
     [False False False False False]]
    >>> print(done_pixels)
    [[False False  True  True  True]
     [False  True  True  True  True]
     [ True  True  True  True  True]
     [False  True  True  True  True]]
    """
        
    # Initialize a stack with the starting point
    stack = deque([start])
    
    # Initialize the mask array indicating which pixels are part of the plume (True = plume).
    mask = np.zeros(data.shape, dtype=bool)
    
    # Initialize an empty set to keep track of already visited cells
    done = set()
    
    # Extract the number of rows and columns from the data array
    n_rows, n_cols = data.shape

    # Begin the flood fill process
    while stack:

        # Pop a coordinate (lat_i, lon_i) from the stack
        lat_i, lon_i = stack.pop()

        # Skip if the current coordinate has already been visited
        if (lat_i, lon_i) in done:
            continue

        # Mark this cell as visited by adding it to the 'done' set
        done.add((lat_i, lon_i))

        # If the value at the current cell is below the threshold, continue to the next cell
        if data[lat_i, lon_i] < SPM_threshold:
            
            # If this is the first cell visited ...
            if len(done) <= 1 :
                
                # ... Add the neighboring cells (in the given directions) to the stack for further exploration
                [ stack.append((lat_i + x[0], lon_i + x[1])) for x in directions]                
                
            continue
    
        # If the value at the current cell meets or exceeds the threshold, mark it in the mask as True
        if data[lat_i, lon_i] >= SPM_threshold:
            mask[lat_i, lon_i] = True
                    
        # Generate the coordinates of neighboring cells based on the current location and directions
        coordinates_to_inspect = [ (lat_i + d_lat, lon_i + d_lon) for d_lat, d_lon in directions]
        
        # If any neighboring cell is out of bounds (either row or column), skip to the next iteration
        if ( any( [ (x[0] > data.shape[0]-1) or (x[0] < 0)  for x in coordinates_to_inspect] ) or 
             any( [ (x[1] > data.shape[1]-1) or (x[1] < 0) for x in coordinates_to_inspect] ) ) :
            continue
        
        # If all neighboring cells are NaN, handle the special case
        if all( np.isnan( [ data[ coordinates ] for coordinates in coordinates_to_inspect ] ) ) :
            
            # If this is the first cell, search for the closest finite values
            if len(done) == 1 : 
                closest_finite_values = []
                closest_coordinates = coordinates_to_inspect
                
                # Continue searching for finite values within the data array until at least 10 are found
                while (len(closest_finite_values) < 10) and (len(closest_coordinates) < 40000000) :
                    # print(len(closest_coordinates))
                    closest_coordinates = [ (lat_coord + d_lat, lon_coord + d_lon) for d_lat, d_lon in directions for lat_coord, lon_coord in closest_coordinates]
                    closest_finite_values = [ data[ coordinates ] for coordinates in closest_coordinates if np.isfinite( data[ coordinates ] ) ]
                    
                # If more than 80% of the closest finite values exceed the threshold, mark the mask as True
                if np.sum(np.array(closest_finite_values) > SPM_threshold) > len(closest_finite_values)*0.8 : 
                    closest_rows, closest_columns = zip(*closest_coordinates)
                    mask[closest_rows, closest_columns] = True
                    
                    # Add the closest coordinates to the stack for further exploration
                    [ stack.append(x) for x in closest_coordinates ]
            
            continue

        # For each neighboring cell, check if it is valid and continue the flood fill
        for d_lat, d_lon in directions:
            n_lat, n_lon = lat_i + d_lat, lon_i + d_lon

            # Skip if the neighboring cell is out of bounds
            if n_lat < 0 or n_lat >= n_rows or n_lon < 0 or n_lon >= n_cols:
                continue

            # If the neighboring cell contains a NaN value, add it to the stack for further processing
            if np.isnan(data[n_lat, n_lon]):
                stack.append((n_lat, n_lon))
                continue

            # If the neighboring cell's value exceeds the threshold, add it to the stack
            if data[n_lat, n_lon] > SPM_threshold:
                stack.append((n_lat, n_lon))

                # If the latitude difference is greater than 1, mark all cells between as True in the mask
                if abs(d_lat) > 1:
                    if d_lat > 1:
                        lat_index_to_True = np.arange(lat_i, n_lat + 1)
                    else:
                        lat_index_to_True = np.arange(n_lat, lat_i + 1)
                    mask[lat_index_to_True, n_lon] = True

                # If the longitude difference is greater than 1, mark all cells between as True in the mask
                if abs(d_lon) > 1:
                    if d_lon > 1:
                        lon_index_to_True = np.arange(lon_i, n_lon + 1)
                    else:
                        lon_index_to_True = np.arange(n_lon, lon_i + 1)
                    mask[n_lat, lon_index_to_True] = True

    # Create a boolean array indicating the pixels that have been visited (done)
    done_pixels = np.zeros(mask.shape, dtype=bool)  
    x_coords_of_done_pixels, y_coords_of_done_pixels = zip(*done)
    done_pixels[np.array(x_coords_of_done_pixels), np.array(y_coords_of_done_pixels)] = True
    
    return mask, done_pixels


def find_high_value_pixels(data, center_lat, center_lon, radius_km, SPM_threshold):
        
    """
    Identify pixels in a 2D dataset that are within a specified radius from a given center point 
    and have values exceeding a threshold.

    Parameters
    ----------
    data : numpy.ndarray
        A 2D array representing the input data (e.g., satellite data).
    center_lat : float
        Latitude of the center point in degrees.
    center_lon : float
        Longitude of the center point in degrees.
    radius_km : float
        Radius within which pixels are considered, in kilometers.
    SPM_threshold : float
        Threshold value. Pixels with values greater than this threshold will be included in the mask.

    Returns
    -------
    numpy.ndarray
        A boolean 2D array where `True` indicates pixels that are within the specified radius 
        and exceed the threshold.

    Notes
    -----
    The function uses the Haversine formula to calculate the great-circle distance between
    each pixel and the specified center point. Only pixels within the radius and above the
    threshold are included in the returned mask.

    Examples
    --------
    >>> import numpy as np
    >>> # Create a sample 2D dataset with latitudes and longitudes
    >>> lat = np.linspace(-10, 10, 5)  # 5 latitudes from -10 to 10 degrees
    >>> lon = np.linspace(-20, 20, 5)  # 5 longitudes from -20 to 20 degrees
    >>> data = np.random.rand(len(lat), len(lon)) * 100  # Random data values
    >>> # Define center and parameters
    >>> center_lat = 0
    >>> center_lon = 0
    >>> radius_km = 1000
    >>> SPM_threshold = 50
    >>> mask = find_high_value_pixels(data, center_lat, center_lon, radius_km, SPM_threshold)
    >>> print(mask)
    [[False  True False ...  True]
     [ True False False ...]]
    """
    
    # Earth's radius in kilometers
    earth_radius_km = 6371.0

    # Convert latitude and longitude to radians
    lat_radians = np.deg2rad(data.lat)
    lon_radians = np.deg2rad(data.lon)
    
    # Center point in radians
    center_lat_rad = np.deg2rad(center_lat)
    center_lon_rad = np.deg2rad(center_lon)
    
    # Calculate distances using Haversine formula
    dlat = lat_radians - center_lat_rad
    dlon = lon_radians - center_lon_rad
    a = np.sin(dlat / 2.0)**2 + np.cos(center_lat_rad) * np.cos(lat_radians) * np.sin(dlon / 2.0)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distances = earth_radius_km * c
    
    # Create a mask for pixels within the specified radius and above the SPM_threshold
    mask_radius = distances <= radius_km
    mask_SPM_threshold = data > SPM_threshold
    mask = mask_radius & mask_SPM_threshold
    
    return mask


def make_the_plot(path_to_the_figure_file_to_save, ds, ds_reduced, france_shapefile,
                  lon_range_of_the_map_to_plot, lat_range_of_the_map_to_plot, 
                  bathymetric_threshold, bathymetry_data_aligned_to_reduced_map, thresholds,
                  plot_the_plume_area,
                  mask_area = None, close_river_mouth_area = None,
                  pixel_done = None,
                  show_bathymetric_mask = False) : 
    
    """
   Generate and save a visualization comparing the original and processed datasets, including optional plume and bathymetric masks.

   Parameters
   ----------
   path_to_the_figure_file_to_save : str
       Path to save the generated figure as a PNG file.
   ds : xarray.DataArray
       Original dataset to be plotted on the first subplot.
   ds_reduced : xarray.DataArray
       Processed or reduced dataset used to determine color bar limits and optional masking.
   france_shapefile : geopandas.GeoDataFrame
       Shapefile of France's boundary used for overlaying on the plots.
   lon_range_of_the_map_to_plot : tuple of float
       Longitude range (min, max) for both subplots.
   lat_range_of_the_map_to_plot : tuple of float
       Latitude range (min, max) for both subplots.
   bathymetric_threshold : float
       Depth threshold (in meters) used to create a bathymetric mask.
   bathymetry_data_aligned_to_reduced_map : xarray.DataArray
       Bathymetry dataset aligned with the reduced dataset grid.
   thresholds : dict
       Dictionary containing SPM thresholds to describe the plume area in the title (e.g., {"threshold1": value1, "threshold2": value2}).
   plot_the_plume_area : bool
       If `True`, plots the plume area on the second subplot; otherwise, the second subplot contains a placeholder message.
   mask_area : xarray.DataArray, optional
       Mask indicating the plume area to highlight, with `True` values for plume pixels. Default is `None`.
   pixel_done : xarray.DataArray, optional
       Boolean mask indicating pixels that were processed, with `True` for processed pixels. Default is `None`.
   show_bathymetric_mask : bool, optional
       If `True`, overlays a bathymetric mask on the second subplot. Default is `False`.

   Returns
   -------
   None
       The function saves the figure to the specified path and does not return any value.

   Notes
   -----
   - The function uses logarithmic normalization for the color scale of the plotted datasets.
   - The bathymetric mask is created based on the specified depth threshold.
   - The title of the second subplot dynamically includes SPM thresholds and bathymetric information.
   - If `plot_the_plume_area` is `False`, the second subplot shows a placeholder message indicating no plume detection.

   Examples
   --------
   >>> import xarray as xr
   >>> import geopandas as gpd
   >>> # Example inputs
   >>> ds = xr.DataArray(...)
   >>> ds_reduced = xr.DataArray(...)
   >>> france_shapefile = gpd.read_file("path_to_shapefile.shp")
   >>> lon_range = (-10, 10)
   >>> lat_range = (40, 50)
   >>> bathymetry = xr.DataArray(...)
   >>> thresholds = {"threshold1": 0.5, "threshold2": 1.0}
   >>> make_the_plot(
   ...     "output_figure",
   ...     ds,
   ...     ds_reduced,
   ...     france_shapefile,
   ...     lon_range,
   ...     lat_range,
   ...     bathymetric_threshold=50,
   ...     bathymetry_data_aligned_to_reduced_map=bathymetry,
   ...     thresholds=thresholds,
   ...     plot_the_plume_area=True
   ... )
   """
    
    # Determine color bar limits based on the dataset and optional mask
    if mask_area is not None : 
        # Use the mask area to calculate the color bar range
        
        # max_color_bar = np.max( [np.nanquantile(ds.values, 0.99), 1] )
        max_color_bar = np.nanmax( [np.nanquantile(ds_reduced.values[np.where(mask_area.values)], 0.85), 1] )
        # max_color_bar = np.nanquantile(ds.values, 0.99)
        min_color_bar = np.nanmax( [np.nanmin(ds_reduced.values[np.where(mask_area.values)]), 0.1] )
        
    else : 
        # Calculate the color bar range without a mask

        max_color_bar = np.max( [np.nanquantile(ds.values, 0.95), 1] )
        # max_color_bar = np.nanquantile(ds_reduced.values, 0.99)
        # max_color_bar = np.nanquantile(ds.values, 0.99)
        min_color_bar = np.max([np.nanmin(ds_reduced.values), 0.1]) 
    
    # max_color_bar = 1.5
    # min_color_bar = 0.6
    
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # Plot the original dataset on the first subplot
    plt.subplot(1, 2, 1)
    ds.plot(vmin = min_color_bar, vmax = max_color_bar, norm=colors.LogNorm())
    france_shapefile.boundary.plot(ax=ax1, linewidth=1, edgecolor='black')

    # Set plot titles and axis labels for the first subplot
    ax1.set_title('Original map')
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.set_xlim([lon_range_of_the_map_to_plot[0], lon_range_of_the_map_to_plot[1]])
    ax1.set_ylim([lat_range_of_the_map_to_plot[0], lat_range_of_the_map_to_plot[1]])
    
    # Plot the reduced dataset and additional elements on the second subplot
    plt.subplot(1, 2, 2)
    
    if plot_the_plume_area : 
        
        # Plot the reduced dataset with the same color bar limits
        ds.plot(vmin = min_color_bar, vmax = max_color_bar, norm=colors.LogNorm())
       
        if show_bathymetric_mask : 
           
            # Create a bathymetric mask based on the threshold
            bathymetric_mask = xr.DataArray(np.zeros(ds_reduced.shape, dtype=bool), coords=ds_reduced.coords, dims=ds_reduced.dims)
            
            # Identify pixels above the bathymetric threshold
            bathymetric_mask.values[ (bathymetry_data_aligned_to_reduced_map.values > -bathymetric_threshold)] = True
            
            # Overlay the bathymetric mask on the plot
            ax2.contourf(bathymetric_mask.lon, bathymetric_mask.lat, bathymetric_mask.values, levels=[0.5, 1], colors=['lightgray']) 

        if pixel_done is not None : 
            
            # Overlay processed pixels on the plot
            ax2.contourf(pixel_done.lon, pixel_done.lat, pixel_done.values, levels=[0.5, 1], colors=['black'], alpha=1)
        
        # Overlay the plume mask on the plot
        ax2.contourf(mask_area.lon, mask_area.lat, mask_area.values, levels=[0.5, 1], colors=['red'])
       
        if close_river_mouth_area is not None : 
            # Overlay the mask of close river mouth 
            # index_True_values = np.where(close_river_mouth_area.values)
            ax2.contourf(close_river_mouth_area.lon, close_river_mouth_area.lat, close_river_mouth_area.values, 
                         levels=[0.5, 1], colors=['rosybrown']) 
       
        # Overlay France's boundary on the second subplot
        france_shapefile.boundary.plot(ax=ax2, linewidth=1, edgecolor='black')
        
        # Set axis limits for the second subplot
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        ax2.set_xlim([lon_range_of_the_map_to_plot[0], lon_range_of_the_map_to_plot[1]])
        ax2.set_ylim([lat_range_of_the_map_to_plot[0], lat_range_of_the_map_to_plot[1]])
        
        # Add a title summarizing the plume detection
        fig.suptitle(f'Plume detection ({os.path.basename(path_to_the_figure_file_to_save)})', fontsize = 20)
       
    else : 
        
        # Add a placeholder title when no plume is detected
        fig.suptitle(f'No plume detection because too much nan in the searching area ({os.path.basename(path_to_the_figure_file_to_save)})', fontsize = 20)
        
    # Add a title to the second subplot with bathymetric and SPM threshold details
    ax2.set_title(f'Area with bathy > {bathymetric_threshold}m and SPM {"; ".join([f"> {round(value, 1)} ({key})" for key, value in thresholds.items() if value is not None])} g m-3')
       
    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the figure to the specified file
    plt.savefig(f'{path_to_the_figure_file_to_save}.png')
            
    # Close the figure to free up memory
    plt.close(fig)
    
    
def identify_the_shape_label_corresponding_to_the_plume(mask_area, core_of_the_plume) :

    """
    Assigns a label to the main plume area, which is the one covering the coordinates of the plume core.
    
    Parameters
    ----------
    mask_area : xarray.DataArray
        A boolean array indicating the pixels inside (True) or outside (False) the plume area.
    core_of_the_plume : tuple of float
        A tuple containing the (latitude, longitude) coordinates of the plume's core.
        
    Returns
    -------
    label_of_the_plume_shape : int
        The label corresponding to the main plume area.
    labeled_array : numpy.ndarray
        A 2D array where each connected shape in the mask is assigned a unique label.
    num_features : int
        The total number of unique connected shapes (features) identified in the mask.
    
    Notes
    -----
    If the core pixel does not belong to any labeled shape (label is 0), the function identifies the closest shape
    to the plume core based on centroid distance.
    """  
    
    # Find the index of the latitude closest to the core's latitude
    lat_core_idx = np.abs(mask_area.lat - core_of_the_plume[0]).argmin()
    # Find the index of the longitude closest to the core's longitude
    lon_core_idx = np.abs(mask_area.lon - core_of_the_plume[1]).argmin()
    # Combine the indices into a tuple representing the core pixel
    pixel_core_of_the_plume = (int(lat_core_idx), int(lon_core_idx)) 
   
    # Label connected regions in the mask and count the total number of distinct shapes
    labeled_array, num_features = label(mask_area.values)
    
    # Get the label of the shape at the core pixel
    label_of_the_plume_shape = labeled_array[pixel_core_of_the_plume]
    
    # If the core pixel does not belong to any shape (label is 0)
    if label_of_the_plume_shape == 0 : 
        
        # Calculate centroids of all labeled shapes in the mask
        centroids = center_of_mass(mask_area.values, labeled_array, range(1, num_features + 1))
        
        # Calculate distances from the core pixel to each shape's centroid
        distances = cdist([pixel_core_of_the_plume], centroids)
        
        # Identify the closest shape's label (adjusted for 1-based labeling)
        label_of_the_plume_shape = np.argmin(distances) + 1  # labels are 1-indexed    
    
    # Return the label of the main plume shape, the labeled array, and the total number of features
    return label_of_the_plume_shape, labeled_array, num_features


def merge_plume_shape_with_close_shapes(mask_area, core_of_the_plume, land_mask, structure_of_the_dilation) :

    """
    Merges the main plume shape in a mask with nearby shapes if they are close enough. This process
    adjusts the mask to include neighboring shapes that, when connected, form a single larger shape.
    
    Parameters
    ----------
    mask_area : xarray.DataArray
        A boolean array (mask) indicating the pixels inside (True) or outside (False) the plume area.
    core_of_the_plume : tuple of float
        Coordinates (latitude, longitude) identifying the core of the plume shape to be preserved 
        (i.e., the main plume shape).
    land_mask : xarray.DataArray
        A boolean array where True indicates land areas that must not be part of the plume shape.
    structure_of_the_dilation : numpy.ndarray
        Structuring element used for the morphological dilation of the main plume shape to connect 
        it to neighboring plume areas if they are close enough.
        
    Returns
    -------
    mask_area : xarray.DataArray
        The updated mask where the plume shape has potentially been merged with nearby shapes.
    """
    
    # Identify and label the main plume area, returning its label and the labeled array
    label_of_the_shape_to_keep, labeled_array, num_features = identify_the_shape_label_corresponding_to_the_plume(mask_area, core_of_the_plume)
    
    # Create a new mask with all False values to store the main plume area
    mask_zone_to_keep = xr.zeros_like(mask_area).astype(bool)
    
    # Set True for pixels belonging to the main plume shape
    mask_zone_to_keep.values[labeled_array == label_of_the_shape_to_keep] = True
    
    # Dilate the main plume shape to extend its boundary and include neighboring pixels
    dilated_mask_zone_to_keep = binary_dilation(mask_zone_to_keep.values, structure=structure_of_the_dilation)
    
    # Identify new candidate pixels that are in the dilated area but not yet in the main plume shape
    new_candidates = np.where( dilated_mask_zone_to_keep & (mask_zone_to_keep.values == False) )
    new_candidates_coordinates = list(zip(new_candidates[0], new_candidates[1]))
     
    # Define a function to check if adding a candidate pixel merges separate shapes
    def check_if_joins_shapes(labeled_array, x, y):
        
        """
        Checks if setting a specific pixel to True will merge two separate shapes.
        
        Parameters
        ----------
        labeled_array : numpy.ndarray
            Array where each connected shape is labeled with a unique integer.
        x, y : int
            Indices of the pixel to test.
        
        Returns
        -------
        bool
            True if setting the pixel to True reduces the total number of shapes.
        """
        
        # Create a temporary copy of the labeled array
        temp_array = labeled_array.copy()
        # Set the pixel at (x, y) to 1 (mark it as part of a shape)
        temp_array[x, y] = 1  
        
        # Re-label the shapes in the modified array
        new_labeled_array, new_num_features = label(temp_array)
        
        # If the number of shapes decreases, adding this pixel merges two shapes
        return new_num_features < num_features
    
    # Check if any candidate pixel, when set to True, will join separate shapes
    if any( [check_if_joins_shapes(labeled_array, x, y) for x, y in new_candidates_coordinates] ) :
        
        # If so, update the mask_area by setting all candidate pixels to True
        rows_idd, columns_idd = zip(*new_candidates_coordinates)
        mask_area.values[rows_idd, columns_idd] = True
        
        # Ensure land pixels remain False in the updated mask
        mask_area.values[land_mask.values == True] = False
         
    # Return the updated mask
    return mask_area


def Set_cloudy_regions_to_True(ds_reduced, mask_area, land_mask, SPM_threshold) :
            
    """
    Refines the plume shape by marking cloudy areas as part of the plume 
    if they are enclosed by SPM concentrations higher than the SPM_threshold.
    This function helps address the issue of cloudy regions disrupting plume detection.

    Parameters
    ----------
    ds_reduced : xarray.DataArray
        The Suspended Particulate Matter (SPM) map containing concentration values.
    mask_area : xarray.DataArray
        A boolean mask indicating the plume area (True for plume, False otherwise).
    land_mask : xarray.DataArray
        A mask indicating land areas (True for land, False for ocean).
    SPM_threshold : float
        The SPM threshold value used to detect and confirm the plume area.

    Returns
    -------
    mask_area : xarray.DataArray
        Updated boolean mask where additional cloudy areas have been marked as True.
    """
        
    # Create a copy of the original mask to prevent unintended modifications
    mask_area_to_use = mask_area.copy()
    
    # Mark all land regions as True since they are irrelevant for plume detection
    mask_area_to_use.values[ land_mask.values | np.isfinite(ds_reduced.values) ] = True
        
    # Label connected regions of False values in the mask
    labeled_false_array, num_false_features = label(~mask_area_to_use)
    
    # Define a dilation strategy that expands regions by one cell in all directions
    dilation_strategy = np.array([[True, True, True],
                                [True, True, True],
                                [True, True, True]])
    
    # Iterate through each labeled region of False values (non-plume areas)
    for false_area_idd in range(num_false_features+1) : 
        
        # Create a boolean mask for the current false area
        false_area = (labeled_false_array == false_area_idd)
        
        # Skip areas that touch the edge of the map as they are not enclosed
        if (any( false_area[0,:] ) or any( false_area[:,0] ) or any( false_area[-1,:] ) or any( false_area[:,-1] )) : 
           continue      
        
        # Dilate the false area to expand its boundary
        false_area_diluted = binary_dilation(false_area, structure=dilation_strategy)
        
        # Identify new candidate pixels added by the dilation process
        new_candidates = np.where( false_area_diluted & (false_area == False) )
        
        # Convert new candidate coordinates into a list of tuples
        new_candidates_coordinates = list(zip(new_candidates[0], new_candidates[1]))
        
        # Check if each candidate satisfies one of the conditions:
        # 1. The SPM concentration is above the threshold
        # 2. The pixel is over land
        test = [bool((ds_reduced[(x_i,x_j)] > SPM_threshold) or (land_mask[(x_i,x_j)] == True)) for x_i, x_j in new_candidates_coordinates]
        
        # Calculate the percentage of candidates that meet the conditions
        percentage_of_test_ok = 100 * len([x for x in test if x == True]) / len(test)
        
        # If more than 90% of the candidates meet the conditions, mark the false area as plume
        if percentage_of_test_ok > 90 : 
            mask_area.values[labeled_false_array == false_area_idd] = True
            
    # Return the updated mask with the cloudy regions marked as part of the plume
    return mask_area


def load_and_resize_files(file_name_config) : 
    
    """
    Loads a file and resizes its map data based on the provided configuration.

    Parameters
    ----------
    file_name_config : list
        A list where:
        - file_name_config[0] (str): The file path to load.
        - file_name_config[1] (int): The desired width for resizing.
        - file_name_config[2] (int): The desired height for resizing.

    Returns
    -------
    numpy.ndarray
        A resized version of the map data.
    
    Notes
    -----
    The function assumes the input file contains a dictionary-like structure
    with a "Basin_map" key or directly contains "map_data". If "Basin_map" is present,
    its "map_data" is used.
    """
    
    # Load the file specified in the configuration.
    data = load_file(file_name_config[0])
    
    # Check if the file contains a "Basin_map" key. If so, extract its value.
    if 'Basin_map' in data : 
        data = data['Basin_map']
        
    # Extract the "map_data" from the file contents.
    data = data['map_data']
    
    # Resize the map data to the specified width and height.
    resized_data = reduce_resolution(data, file_name_config[1], file_name_config[2])
    
    # Return the resized map data.
    return resized_data


def compute_gradient_with_directions_vectorized(spm_map, start_point, directions, max_steps, 
                                                lower_high_values_to = None,
                                                create_X_intermediates_between_each_direction = 1):
    
    """
    Computes the gradient in multiple directions from a starting point using vectorized operations.

    Parameters
    ----------
    spm_map : ndarray
        A 2D array of Suspended Particulate Matter (SPM) concentrations.
    start_point : tuple
        Coordinates (y, x) of the starting point in the map.
    directions : list of tuple
        List of direction vectors (dy, dx) for gradient computation.
    max_steps : int, optional
        Maximum number of steps to compute in each direction (default is 100).
    create_X_intermediates_between_each_direction : int, optional
        Number of intermediate direction vectors to generate between existing ones (default is 1).

    Returns
    -------
    relative_gradient_values : ndarray
        Array of normalized gradient values for each direction and step.
    gradient_points : ndarray
        Array of coordinates (y, x) corresponding to the gradient values.
    spm_values : ndarray
        SPM values at the computed gradient points.
    """
        
    # Get the dimensions of the map and unpack the starting point coordinates
    rows, cols = spm_map.shape
    start_y, start_x = start_point
    
    # Convert the list of directions into a NumPy array for easier manipulation
    directions = np.array(directions)
    
    # Generate intermediate direction vectors
    all_directions = []
    for rep in range(create_X_intermediates_between_each_direction) : 
        
        if rep == 1 :
        
            # Compute new directions by averaging adjacent vectors
            new_directions = (all_directions[:-1] + all_directions[1:]) / 2.0

            # Define the order of elements and assemble all directions
            order_of_elements = np.append( np.tile([1,2], int( (len(all_directions) + len(new_directions) -1) / 2) ), [1] )   
            all_directions = assemble_and_order_directions(all_directions, new_directions, order_of_elements) # [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1])
                                   
        else : 
            
            # Compute intermediate directions from initial ones
            intermediate_directions = (directions[:-1] + directions[1:]) / 2.0
            order_of_elements = np.append( np.tile([1,2], int( (len(directions) + len(intermediate_directions) -1) / 2) ), [1] )   
            all_directions = assemble_and_order_directions(directions, intermediate_directions, order_of_elements)
        
    # Create steps (range of step numbers)
    max_steps = 50 if max_steps is None else max_steps
    if len(spm_map.shape) == 2 :
        steps = np.arange(0, max_steps).reshape(-1, 1) 
    elif len(spm_map.shape) == 3 : 
        steps = np.arange(0, max_steps).reshape(-1, 1, 1) 
    
    # Compute new positions by adding scaled direction vectors to the start point
    new_y = (start_y + steps * all_directions[:, 0]).astype(int).transpose()  # Cast to int
    new_x = (start_x + steps * all_directions[:, 1]).astype(int).transpose()  # Cast to int
    
    # Ensure coordinates are within map bounds
    valid_mask = (new_y >= 0) & (new_y < rows) & (new_x >= 0) & (new_x < cols)
    new_y = np.where(valid_mask, new_y, -1)  # Mark invalid coordinates as -1 (invalid index)
    new_x = np.where(valid_mask, new_x, -1)
    index_invalid_coordinates = np.where((new_y == -1) | (new_x == -1))
    
    # Extract SPM values at the computed positions
    spm_values = np.take(spm_map, new_y * cols + new_x, mode='clip')  # Flatten the 2D array into 1D
    spm_values[ index_invalid_coordinates ] = np.nan # Assign nan to the coordinates -1
    
    # Apply upper limit to avoid bias from extreme values 
    if lower_high_values_to is not None : 
        spm_values = xr.where(spm_values > lower_high_values_to, lower_high_values_to, spm_values)
    
    # Compute gradient values as differences between consecutive steps
    gradient_values = np.diff(spm_values, axis=-1, prepend=np.full((spm_values.shape[0], 1), np.nan))
    
    if lower_high_values_to is not None : 
        gradient_values[np.where(spm_values == lower_high_values_to)] = float('inf') # np.nanmax(gradient_values) 
    
    # Normalize gradient values using a robust range
    relative_gradient_values = (gradient_values / (np.nanmax(spm_values) - np.nanmin(spm_values))
                                # (gradient_values / (np.nanquantile(spm_values, 0.9) - np.nanquantile(spm_values, 0.1))
                                # if (np.nanquantile(spm_values, 0.9) - np.nanquantile(spm_values, 0.1)) > 0 
                                if (np.nanmax(spm_values) > np.nanmin(spm_values)) 
                                else None)
    
    # Create an array of gradient points, ignoring invalid coordinates
    gradient_points = np.stack((new_x, new_y), axis=-1).astype(float)
    gradient_points[ index_invalid_coordinates ] = [np.nan, np.nan] # Assign nan to the coordinates -1
    
    return relative_gradient_values, gradient_points, spm_values


def filter_gradient_points_vectorized(gradient_values, gradient_points, absolute_values, land_mask, threshold_on_gradient_values):
    
    """
    Filters gradient points based on the gradient values and proximity to land.

    Parameters
    ----------
    gradient_values : ndarray
        2D array of gradient values for each direction and step.
    gradient_points : ndarray
        Array of coordinates (y, x) corresponding to the gradient values.
    absolute_values : ndarray
        Absolute values associated with the gradient points.
    land_mask : ndarray
        Binary mask indicating land (1) and water (0) areas.
    threshold_on_gradient_values : float
        Threshold for gradient value filtering.

    Returns
    -------
    gradient_points_to_keep : ndarray
        Filtered coordinates of gradient points that meet the criteria.
    gradient_values_to_keep : ndarray
        Filtered gradient values corresponding to the kept points.
    absolute_values_to_keep : ndarray
        Filtered absolute values corresponding to the kept points.
    """
            
    # Identify points that are far from land
    are_pixels_far_from_land = pixels_far_from_land(land_mask, pixel_positions = gradient_points, distance_threshold = 20)
    
    # # Plotting the points where the gradient was computed
    # plt.imshow(spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = 8)
    # for x, y in gradient_points[are_pixels_far_from_land].reshape(-1, gradient_points.shape[-1]) :
    #     plt.scatter(x, y, s=1)
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
        
    # # Get the index of the first finite values for each axis = -1
    # first_finite_mask = np.isfinite(gradient_values)
    # first_finite_indices = np.argmax( first_finite_mask , axis=-1)
    # no_finite_elements = first_finite_mask.any(axis=-1) == False
    # first_finite_indices[no_finite_elements] = 999 
    # mask_values_after_the_first_finite = np.arange(gradient_values.shape[1])[None, :] > first_finite_indices[:, None]  
    
    # Step 1: Create a mask for valid gradient values
    # mask_higher_than_threshold = (abs(gradient_values) > threshold_on_gradient_values) & are_pixels_far_from_land & mask_values_after_the_first_finite
    mask_higher_than_threshold = (abs(gradient_values) > threshold_on_gradient_values) & are_pixels_far_from_land
        
    # Step 2: Compute the maximum step index where the gradient exceeds the threshold
    # max_steps_higher_than_threshold = mask_higher_than_threshold.shape[1] - np.argmax(mask_higher_than_threshold, axis=-1) 
    max_steps_higher_than_threshold = np.argmax( np.cumsum(mask_higher_than_threshold, axis = -1) , axis = -1 )
    # max_steps_higher_than_threshold[no_finite_elements] = -1
        
    # Step 3: Create a mask for valid steps based on the max step index
    step_mask = (max_steps_higher_than_threshold)[:, None] > np.arange(gradient_values.shape[1])
    
    # # Plotting the points where the gradient was computed
    # plt.imshow(spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = 8)
    # for x, y in gradient_points[step_mask].reshape(-1, gradient_points.shape[-1]) :
    #     plt.scatter(x, y, s=1)
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
    
    # Compute thresholds for absolute values
    superior_limit = np.append( np.nanquantile(absolute_values[np.where(are_pixels_far_from_land)], 0.95), [7] ).min()
    inferior_limit = np.nanquantile( absolute_values[np.where(are_pixels_far_from_land)], 0.05 )
    threshold_on_absolute_values = 10 ** ( np.log10( inferior_limit ) + (( np.log10( superior_limit ) - np.log10( inferior_limit ) ) / 1.5) )
        
    # Combine all masks
    mask_on_absolute_values = absolute_values > threshold_on_absolute_values # 0.75
    combined_mask = step_mask & mask_on_absolute_values & are_pixels_far_from_land

    # # Plotting the points where the gradient was computed
    # plt.imshow(spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = 8)
    # for x, y in gradient_points[combined_mask].reshape(-1, gradient_points.shape[-1]) :
    #     plt.scatter(x, y, s=1)
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
            
    # Filter values based on the combined mask
    gradient_values_to_keep = gradient_values[np.where(combined_mask)]
    gradient_points_to_keep = gradient_points[np.where(combined_mask)] 
    absolute_values_to_keep = absolute_values[np.where(combined_mask)] 
    
    return gradient_points_to_keep, gradient_values_to_keep, absolute_values_to_keep


def find_SPM_threshold(spm_map, land_mask, start_point, directions, max_steps, maximal_threshold, minimal_threshold, quantile_to_use) : 

    """
    Finds the Suspended Particulate Matter (SPM) threshold by analyzing gradients 
    in specified directions and filtering gradient points.

    Parameters
    ----------
    spm_map : ndarray
        2D array representing SPM concentrations.
    land_mask : ndarray
        Boolean mask where True indicates land and False indicates water.
    start_point : tuple
        Starting point (y, x) for gradient computation.
    directions : list
        List of direction vectors for gradient computation.

    Returns
    -------
    SPM_threshold : float
        The computed SPM threshold value.
    """
    
    # nb_of_ocean_pixels = len(np.where(land_mask.values == False)[0])

    # Compute gradients and gradient points in specified directions
    gradient_values, gradient_points, absolute_values = compute_gradient_with_directions_vectorized(spm_map = spm_map, 
                                                                                   start_point = start_point, 
                                                                                   directions = directions,
                                                                                   max_steps = max_steps, # Max steps for gradient expansion
                                                                                   lower_high_values_to=maximal_threshold, # 7
                                                                                   create_X_intermediates_between_each_direction = 2) # Number of intermediate directions

    # threshold_on_gradient_values = 0.35 # 0.2 # 0.15
    threshold_on_gradient_values = np.nanmax(gradient_values[np.isfinite(gradient_values)]) * 0.9

    # are_pixels_far_from_land = pixels_far_from_land(land_mask, pixel_positions = gradient_points, distance_threshold = 20)

    # # Plotting the points where the gradient was computed
    # plt.imshow(spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = maximal_threshold)
    # for x, y in gradient_points.reshape(-1, gradient_points.shape[-1]) :
    #     plt.scatter(x, y, s=1)
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
        
    # # Plot the gradient values over distance
    # for i in range(gradient_values.shape[0]):  # Loop over the 3rd dimension (directions)
    #     plt.plot(gradient_values[i, :], label=f'Direction {i+1}')
    # plt.hlines([-threshold_on_gradient_values, threshold_on_gradient_values], xmin = 0, xmax = max_steps)
    # plt.title("Gradient Computation in Specified Directions")
    # plt.xlabel("Step")
    # plt.ylabel("Gradient of SPM")
    # plt.show()   

    # plt.hist(gradient_values.flatten(), bins=10, edgecolor='black', alpha=0.7)
    # plt.title('Distribution of Gradient Values')
    # plt.xlabel('Gradient Value')
    # plt.ylabel('Frequency')
    # plt.show()

    # Filter gradient points based on specified criteria
    filtered_points, _, filtered_values = filter_gradient_points_vectorized(gradient_values, gradient_points, absolute_values, land_mask,
                                                                            threshold_on_gradient_values = threshold_on_gradient_values) 

    # plt.imshow(spm_map[:, :], cmap='viridis', origin = "lower", vmin = 0.1, vmax = maximal_threshold)
    # for coords in [x for x in filtered_points]:
    #     plt.scatter(coords[0], coords[1], color='red', s=1)  # Mark the points where gradient was computed
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()

    # Compute the SPM threshold as the 10th percentile of filtered values
    # SPM_threshold = np.nanmin( filtered_values )
    SPM_threshold = np.max( [np.nanquantile(filtered_values, quantile_to_use), minimal_threshold] ) # 0.25
    
    return SPM_threshold, filtered_points, gradient_points


def find_first_nan_after_finite(arr):
    
    """
    Detects the first NaN value that occurs after finite values in a 3D array.

    Parameters
    ----------
    arr : ndarray
        3D array of shape (n_steps, 1, n_directions), where the second dimension is 1.

    Returns
    -------
    first_nan_indices : ndarray
        Array containing the index of the first NaN value after finite values for each direction.
        If no such NaN exists, returns -1 for that direction.
    """
    
    # Check if the input array has shape (100, 1, 5)
    # assert arr.shape[1] == 1, "The second dimension must be 1."
    
    # Flatten the array along the second dimension
    flattened = arr[:, :]
    
    # Masks for NaN and finite values
    nan_mask = np.isnan(flattened)
    finite_mask = np.isfinite(flattened)

    # Find the index of the first finite value for each direction
    first_finite_idx = np.argmax(finite_mask, axis=0)

    # Mask to exclude leading NaNs and consider only after the first finite value
    mask_after_finite = np.arange(flattened.shape[0])[:, None] > first_finite_idx[None, :]
    
    # Mask for NaNs that occur after finite values
    nan_after_finite = nan_mask & mask_after_finite

    # Find the first NaN index after finite values for each direction
    first_nan_idx = np.where(np.any(nan_after_finite, axis=0), np.argmax(nan_after_finite, axis=0), -1)
    
    return first_nan_idx


def pixels_far_from_land(land_mask, pixel_positions, distance_threshold):
    
    """
    Filters pixel positions based on their proximity to land.

    Parameters
    ----------
    land_mask : ndarray
        Boolean array where True indicates land and False indicates water.
    pixel_positions : ndarray
        Array of shape (n_positions, 2) representing pixel coordinates [(y, x), ...].
    distance_threshold : float
        Minimum allowable distance from land in pixels.

    Returns
    -------
    filtered_positions : ndarray
        Boolean array of shape (n_positions,) indicating whether each pixel is sufficiently far from land.
    """   
    
    # Compute Euclidean distance transform from land
    distance_map = distance_transform_edt(~land_mask)  # EDT treats False as "objects" for distance computation
    
    # Calculate distances for the given pixel positions
    distance_from_land = np.array( [ distance_map[int(y), int(x)]
                                      if ( (np.isfinite(x)) and (np.isfinite(y)) ) else -1
                                      for index, pixel_position in enumerate(pixel_positions) 
                                      for (x,y) in pixel_position ] ).reshape(pixel_positions.shape[:-1])
            
    # Compute differences in distances for validation
    diff_distance_from_land = np.hstack((np.full((distance_from_land.shape[0], 1), 999), 
                                         np.diff(distance_from_land)))
    
    # Filter positions based on distance threshold and distance differences
    filtered_positions = (distance_from_land > distance_threshold) | (diff_distance_from_land > 0)
    filtered_positions[:int(filtered_positions.shape[0]/2)] = distance_from_land[:int(filtered_positions.shape[0]/2)] > 0
        
    return filtered_positions


def first_true_block(arr):
        
    """
    Finds the position of the first contiguous block of True values in a boolean array.

    Parameters
    ----------
    arr : np.ndarray
        A 1D boolean array.

    Returns
    -------
    tuple
        A tuple (start, end) indicating the start and end indices of the first block of True values.
        If no True values are found, returns (-1, -1).
    """
    
    # Find the indices where the array is True
    true_indices = np.where(arr)[0]
    
    if len(true_indices) == 0:
        return -1, -1  # Return (-1, -1) if no True values are found
    
    # Calculate differences between consecutive True indices
    diffs = np.diff(true_indices)
    
    # Find boundaries where consecutive indices are interrupted
    block_boundaries = np.where(diffs != 1)[0]
    
    if len(block_boundaries) == 0:  # Single contiguous block case
        return true_indices[0], true_indices[-1]  # Return the first block
    
    # Handle multiple blocks, returning the first one
    block_start = true_indices[0]
    block_end = true_indices[block_boundaries[0]]  # End of the first block
    
    return block_start, block_end


def last_true_block(arr):
    
    """
    Finds the position of the last contiguous block of True values in a boolean array.

    Parameters
    ----------
    arr : np.ndarray
        A 1D boolean array.

    Returns
    -------
    tuple
        A tuple (start, end) indicating the start and end indices of the last block of True values.
        If no True values are found, returns (-1, -1).
    """
    
    if not np.any(arr):  # Check if the array contains any True values
        return -1, -1

    # Pad the array with False to detect transitions more easily
    padded = np.r_[arr, False] 
    
    # Compute differences to identify boundaries of True blocks
    diff = np.diff(padded.astype(int))

    # Identify the end of the last block of True values
    last_end = np.where(diff == -1)[0]
    last_end = last_end[-1] 

    # Find the start of the last block of True values
    last_start_candidates = np.where(diff[:last_end + 1] == 1)[0]
    if len(last_start_candidates) > 0:
        last_start = last_start_candidates[-1]
    else:
        last_start = 0  

    return last_start, last_end


def load_and_filter_arrays(file_name):
    
    """
    Loads a file and filters its arrays based on specific criteria.

    Parameters
    ----------
    file_name : str
        Path to the file to load.

    Returns
    -------
    dict
        A dictionary containing filtered arrays and other values from the file.
    """
    
    # Load the file into a dictionary
    d = load_file(file_name)
    
    # Initialize an empty dictionary for filtered results
    filtered_dict = {}
    
    # Iterate through the dictionary items
    for key, value in d.items():
        
        # Check if the value is a NumPy array or a Pandas DatetimeIndex
        if isinstance(value, np.ndarray) or isinstance(value, pd.DatetimeIndex) :
            if len(value) <= 10:  # Include arrays with length <= 10
                filtered_dict[key] = value
        else:
            filtered_dict[key] = value # Include non-NumPy values
            
    return filtered_dict


def assemble_and_order_directions(directions_1, directions_2, order_of_elements) : 
        
    """
    Assembles and orders direction vectors based on a specified order.

    Parameters
    ----------
    directions_1 : np.ndarray
        Array of direction vectors labeled as group 1.
    directions_2 : np.ndarray
        Array of direction vectors labeled as group 2.
    order_of_elements : list
        A list specifying the order in which directions should be assembled, where 
        elements are labeled as 1 (group 1) or 2 (group 2).

    Returns
    -------
    np.ndarray
        A 2D array containing assembled and ordered direction vectors.
    """
    
    # Initialize an empty array to store ordered directions
    all_directions = np.empty((len(order_of_elements), 2))
        
    # Loop through the order array and assign directions accordingly
    for i, idx in enumerate(order_of_elements):
        index = np.where( np.where(np.array(order_of_elements )== idx)[0] == i )[0]
        if idx == 1:
            all_directions[i] = directions_1[index]  # Assign from directions_1 for label 1
        elif idx == 2:
            all_directions[i] = directions_2[index]  # Assign from directions_2 for label 2
            
    return(all_directions)


def set_mask_area_values_to_False_based_on_an_index_object(mask_area, index_object, gradient_points) :

    """
    Modify a mask array by setting values to False for points outside a polygon defined by gradient points.

    Parameters
    ----------
    mask_area : numpy.ndarray
        A 2D array representing the mask area where True values indicate active regions.
    index_object : numpy.ndarray
        Array of indices representing the points used to define the polygon.
    gradient_points : numpy.ndarray
        Array of points representing gradient vectors.

    Returns
    -------
    numpy.ndarray
        Updated mask area with False values outside the polygon.
    """

    # # Plotting the points where the gradient was computed
    # plt.imshow(mask_area, cmap='viridis', origin = "lower")
    # for x, y in gradient_points[index_to_keep].reshape(-1, gradient_points[index_to_keep].shape[-1]):
    #     # plt.scatter(x, y, color='red', s=1)  # Mark the points where gradient was computed
    #     plt.scatter(x, y, s=1)  # Mark the points where gradient was computed
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()

    # Combine gradient points and ensure no NaNs exist
    points_in_the_polygon = np.vstack( (gradient_points[0], gradient_points[index_object][::-1], gradient_points[0][0]) )
    # polygon_boundaries = np.vstack( (gradient_points[0][0], gradient_points[index_object][::-1], gradient_points[0][0]) )
    points_in_the_polygon = points_in_the_polygon[ ~np.isnan(points_in_the_polygon).any(axis=1) ]
    
    if (np.unique( points_in_the_polygon[:,0] ).size == 1) or (np.unique( points_in_the_polygon[:,1] ).size == 1) : 
        return mask_area
    
    # Create a convex hull polygon from the points
    polygon = ConvexHull( points_in_the_polygon )
    polygon_boundaries = points_in_the_polygon[ np.append(polygon.vertices, polygon.vertices[0]) ]
    
    ### Convave version
    # polygon_boundaries = concave_hull(points_in_the_polygon)
    
    # Create a path object for the polygon
    polygon_path = Path(polygon_boundaries, closed=True)
    
    
    # # Plotting the points where the gradient was computed
    # plt.imshow(mask_area, cmap='viridis', origin = "lower")
    # plt.title("Gradient Points Following Direction Vectors")
    # ax = plt.gca()  # Get the current Axes instance
    # ax.add_patch(PathPatch(polygon_path, facecolor='none', edgecolor='red', lw=2))  # Add the polygon path to the plot
    # plt.show()
    
    # Identify True points in the mask and their coordinates
    true_indices = np.where(mask_area.values)
    points = np.column_stack((true_indices[1], true_indices[0]))  # shape (N, 2), where N is number of True points
   
    # Check which points are inside the polygon
    inside_polygon = polygon_path.contains_points(points)   
    
    # # Plotting the points where the gradient was computed
    # plt.imshow(mask_area, cmap='viridis', origin = "lower")
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.scatter( [x[0] for x in points[inside_polygon]] ,
    #               [x[1] for x in points[inside_polygon]], 
    #               color='red')
    # plt.show()
   
    # Update the mask by setting points outside the polygon to False
    for i, (row, col) in enumerate(points):
        if not inside_polygon[i]:
            mask_area.values[col, row] = False  # Set to False if the point is outside the polygon
            
    return mask_area


def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points on the Earth's surface.

    Parameters
    ----------
    lat1 : float
        Latitude of the first point in degrees.
    lon1 : float
        Longitude of the first point in degrees.
    lat2 : float
        Latitude of the second point in degrees.
    lon2 : float
        Longitude of the second point in degrees.

    Returns
    -------
    float
        Distance between the two points in kilometers.
    """
        
    EARTH_RADIUS = 6371.0  # Earth's radius in kilometers
    
    # Convert degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return EARTH_RADIUS * c
        

def remove_coastal_areas_with_sediment_resuspension(ds_reduced, mask_area, land_mask, bathymetry_data_aligned_to_reduced_map,
                                                    SPM_threshold, 
                                                    core_of_the_plume,
                                                    bathymetric_threshold_for_zone_with_resuspension,
                                                    distance_from_estuary_threshold) : 
        
    """
    Remove coastal areas influenced by sediment resuspension based on given criteria.

    Parameters
    ----------
    ds_reduced : xr.DataArray
        Reduced spatial dataset representing suspended particulate matter (SPM).
    mask_area : xr.DataArray
        Boolean mask indicating areas of interest.
    land_mask : xr.DataArray
        Boolean mask indicating land areas.
    bathymetry_data_aligned_to_reduced_map : xr.DataArray
        Aligned bathymetry data for the reduced dataset.
    SPM_threshold : float
        Threshold for SPM above which areas are considered for resuspension.
    core_of_the_plume : tuple
        Latitude and longitude of the plume core (estuary location).
    bathymetric_threshold_for_zone_with_resuspension : float
        Depth threshold (negative values) for zones with sediment resuspension.
    distance_from_estuary_threshold : float
        Maximum allowable distance (in kilometers) from the estuary for resuspension zones.

    Returns
    -------
    xr.DataArray
        Updated mask with resuspension areas removed.
    """
    
    # Identify areas with SPM greater than the specified threshold
    SPM_criterion = ds_reduced > SPM_threshold
    
    # Identify areas with bathymetry shallower than the threshold
    bathymetry_criterion = bathymetry_data_aligned_to_reduced_map.values > -bathymetric_threshold_for_zone_with_resuspension
    
    # Calculate distance from the estuary core
    lat_grid, lon_grid = np.meshgrid(mask_area.lat, mask_area.lon, indexing='ij')
    distance_from_estuary = haversine(lat_grid, lon_grid, core_of_the_plume[0], core_of_the_plume[1])
    distance_from_estuary_criterion = distance_from_estuary > distance_from_estuary_threshold
        
    # Combine criteria into a single mask
    combined_criterion = distance_from_estuary_criterion & SPM_criterion & bathymetry_criterion

    # Find connected shapes in the combined mask that overlap with the land mask
    connected_shapes = find_connected_shapes(xr.DataArray(combined_criterion, coords=ds_reduced.coords, dims=ds_reduced.dims), 
                                             land_mask)
    
    # Update the mask to exclude the identified connected shapes
    mask_area.values[connected_shapes.values] = False
    
    return mask_area
    

def find_connected_shapes(shape_1, shape_2):
     
    """
    Identify connected shapes in one mask that overlap with another mask.

    Parameters
    ----------
    shape_1 : xr.DataArray
        First boolean mask for connected shape identification.
    shape_2 : xr.DataArray
        Second boolean mask indicating overlap regions.

    Returns
    -------
    xr.DataArray
        Boolean mask of connected shapes in shape_1 overlapping with shape_2.
    """
    
    # Label connected regions in both masks
    shape_1_labeled, num_shape_1_labels = label(shape_1)
    shape_2_labeled, num_shape_2_labels = label(shape_2)
    
    # Dilate the second mask to check adjacency
    shape_2_dilated = binary_dilation(shape_2)
    
    # Identify labels in the first mask that overlap with the dilated second mask
    connected_labels = np.unique(shape_1_labeled[shape_2_dilated])
    
    # Remove the background label (0)
    connected_labels = connected_labels[connected_labels > 0]
    
    # Create a mask for connected shapes
    connected_mask = np.isin(shape_1_labeled, connected_labels)
    
    # Convert the result to an xarray.DataArray
    connected_shapes = xr.DataArray(
        connected_mask,
        coords=shape_1.coords,
        dims=shape_1.dims,
        name="connected_shapes"
    )
    
    return connected_shapes


def find_index_and_values_of_multiple_directions_in_the_plume_area(mask_area, ds_reduced,
                                                                   pixel_starting_point, searching_strategy_direction,
                                                                   max_steps) : 
        
    """
    Identify and retrieve index and values along multiple directions in the plume area.

    Parameters
    ----------
    mask_area : xr.DataArray
        Boolean mask indicating the plume area.
    ds_reduced : xr.DataArray
        Reduced spatial dataset for analysis.
    pixel_starting_point : tuple
        Starting pixel coordinates (row, column) for the search.
    searching_strategy_direction : np.ndarray
        Array of direction vectors for searching in the plume area.

    Returns
    -------
    tuple
        - boolean_values_in_the_area_of_the_plume : list
          Boolean values of pixels in the identified plume area.
        - values_in_the_area_of_the_plume : list
          Data values of the identified plume area.
        - direction_points : np.ndarray
          Coordinates of the points in the identified directions.
    """
    
    # Compute gradients along specified directions starting from the pixel
    _, direction_points, direction_boolean = compute_gradient_with_directions_vectorized(spm_map = mask_area.astype(float), 
                                                                                    start_point = pixel_starting_point, 
                                                                                    directions = searching_strategy_direction,
                                                                                    max_steps = max_steps,
                                                                                    create_X_intermediates_between_each_direction = 2) # 35 for the Gulf of Lion
    
    # # Plotting the points where the gradient was computed
    # plt.imshow(mask_area, cmap='viridis', origin = "lower")
    # for x, y in direction_points.reshape(-1, direction_points.shape[-1]):
    #     # plt.scatter(x, y, color='red', s=1)  # Mark the points where gradient was computed
    #     plt.scatter(x, y, s=1)  # Mark the points where gradient was computed
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
    
    # Recompute gradients for the reduced dataset
    _, direction_points, direction_values = compute_gradient_with_directions_vectorized(spm_map = ds_reduced, 
                                                                                    start_point = pixel_starting_point, 
                                                                                    directions = searching_strategy_direction,
                                                                                    max_steps = max_steps,
                                                                                    create_X_intermediates_between_each_direction = 2) # 35 for the Gulf of Lion
    
    # plt.imshow(ds_reduced, cmap='viridis', origin = "lower", vmin = 0.1, vmax = 10)
    # for x, y in direction_points.reshape(-1, direction_points.shape[-1]):
    #     # plt.scatter(x, y, color='red', s=1)  # Mark the points where gradient was computed
    #     plt.scatter(x, y, s=1)  # Mark the points where gradient was computed
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
    
    # Check if all values are zero, indicating no plume area
    if (direction_values == 0).all() or (direction_boolean == 0).all() :                         
        
        boolean_values_in_the_area_of_the_plume = []
        values_in_the_area_of_the_plume = []
    
    else : 
        
        # Extract valid data based on the boolean array
        boolean_values_in_the_area_of_the_plume = direction_boolean.values[ :np.max(np.where(direction_boolean.values.any(axis=0)))+1 ]
        values_in_the_area_of_the_plume = direction_values.values[ :np.max(np.where(direction_boolean.values.any(axis=0)))+1 ]
        
    return boolean_values_in_the_area_of_the_plume, values_in_the_area_of_the_plume, direction_points


def find_the_index_of_the_plume_starting_point(ds_reduced, starting_point):
        
    """
    Locate the index of the plume's starting point in the reduced dataset.
    
    Parameters
    ----------
    ds_reduced : xr.DataArray
        Reduced spatial dataset.
    starting_point : tuple
        Latitude and longitude of the plume's starting point.
    
    Returns
    -------
    tuple
        Index of the starting point as (row, column).
    """
   
    # Find the index of the closest latitude and longitude to the starting point
    lat_idx = np.abs(ds_reduced.lat - starting_point[0]).argmin()
    lon_idx = np.abs(ds_reduced.lon - starting_point[1]).argmin()
    pixel_starting_point = (int(lat_idx), int(lon_idx))
    
    return pixel_starting_point


def return_stats_dictionnary(final_mask_area, spm_reduced_map, spm_map, parameters, thresholds, return_empty_dict = False) :
    
    """
    Compute or return an empty dictionary with statistics about the plume area.

    Parameters
    ----------
    final_mask_area : xr.DataArray
        Boolean mask of the final plume area.
    spm_reduced_map : xr.DataArray
        SPM dataset at reduced resolution.
    spm_map : xr.DataArray
        Original SPM dataset.
    parameters : dict
        Dictionary of map parameters (e.g., resolution, ranges).
    thresholds : dict
        Dictionary of thresholds for SPM.
    return_empty_dict : bool, optional
        If True, return an empty dictionary with NaN values.

    Returns
    -------
    dict
        Dictionary containing plume statistics or an empty dictionary with NaN values.
    """
        
    def make_an_empty_dict() : 
        
        # Initialize an empty dictionary with NaN values
        data_to_return = {'date' : np.array(spm_reduced_map.date_for_plot)}
        data_to_return.update({f'{stat_name}': np.nan for stat_name in ['n_pixel_in_the_plume_area', 'area_of_the_plume_mask_in_km2', 
                                                                        'mean_SPM_in_the_plume_area', 'sd_SPM_in_the_plume_area', 
                                                                        'mass_SPM_in_the_plume_area_in_g_m', 'lat_centroid_of_the_plume_area',
                                                                        'lon_centroid_of_the_plume_area', 'lat_weighted_centroid_of_the_plume_area',
                                                                        'lon_weighted_centroid_of_the_plume_area', 'confidence_index_in_perc']})
        data_to_return.update({f'SPM_threshold_{plume_name}': np.nan for plume_name, threshold in thresholds.items()})
        
        return data_to_return

    
    if return_empty_dict : 
        return make_an_empty_dict()
     
    # Calculate the number of pixels in the plume area
    n_pixel_in_the_plume_area = np.sum(final_mask_area) # Number of pixels in the plume area
    
    # Calculate the area of one pixel in km²
    area_of_one_pixel = (parameters['lon_new_resolution'] * 111.32 * np.cos(np.radians(np.mean(parameters['lat_range_of_the_map_to_plot'])))) * (parameters['lat_new_resolution'] * 111.32) # Area of one pixel in km²
    area_of_the_plume_mask = np.sum(final_mask_area) * area_of_one_pixel # Total plume area in km²
    
    # Compute mean and standard deviation of SPM in the plume area
    mean_SPM_in_the_plume_area = np.nanmean(spm_reduced_map.values[final_mask_area]) # Mean SPM in the plume area
    sd_SPM_in_the_plume_area = np.nanstd(spm_reduced_map.values[final_mask_area]) # Standard deviation of SPM
    
    # Calculate the mass of SPM in the plume area
    mass_SPM_in_the_plume_area = mean_SPM_in_the_plume_area * area_of_the_plume_mask * 1000 # Calculate mass of SPM
    
    # Identify pixels in the plume area
    plume_area_pixels = np.array(np.where(final_mask_area & np.isfinite(spm_reduced_map)))
    if plume_area_pixels.size == 0 : 
        return make_an_empty_dict()
    lat_plume_area_pixels, lon_plume_area_pixels = final_mask_area.lat[plume_area_pixels[0]], final_mask_area.lon[plume_area_pixels[1]]
        
    # Calculate centroids of the plume area
    centroid_of_the_plume_area = [ float(lat_plume_area_pixels.mean()) , float(lon_plume_area_pixels.mean()) ]
    weighted_centroid_of_the_plume_mask = [ np.average(lat_plume_area_pixels, weights = spm_reduced_map.values[plume_area_pixels[0], plume_area_pixels[1]]) , 
                                            np.average(lon_plume_area_pixels, weights = spm_reduced_map.values[plume_area_pixels[0], plume_area_pixels[1]]) ]
        
    # Align the final mask to the original dataset and calculate confidence index
    mask_aligned_to_ds = final_mask_area.astype(int).interp_like(spm_map, method='nearest', kwargs={"fill_value": 0}).astype(bool)
    spm_data_in_the_plume_area = spm_map.values[mask_aligned_to_ds.values]

    # Calculate the confidence index as a percentage of non-NaN values in the plume area
    nb_of_nan_in_the_plume_area = np.isnan( spm_data_in_the_plume_area ).sum()
    confidence_index = 100 * (1 - nb_of_nan_in_the_plume_area / spm_data_in_the_plume_area.size)
        
    # Prepare the statistics dictionary   
    data_to_return = {'date' : np.array(spm_map.date_for_plot),
                      'n_pixel_in_the_plume_area' : float(n_pixel_in_the_plume_area),
                      'area_of_the_plume_mask_in_km2' : float(area_of_the_plume_mask),
                      'mean_SPM_in_the_plume_area' : float(mean_SPM_in_the_plume_area),
                      'sd_SPM_in_the_plume_area' : float(sd_SPM_in_the_plume_area),
                      'mass_SPM_in_the_plume_area_in_g_m' : float(mass_SPM_in_the_plume_area),
                      'lat_centroid_of_the_plume_area' : centroid_of_the_plume_area[0],
                      'lon_centroid_of_the_plume_area' : centroid_of_the_plume_area[1],
                      'lat_weighted_centroid_of_the_plume_area' : weighted_centroid_of_the_plume_mask[0],
                      'lon_weighted_centroid_of_the_plume_area' : weighted_centroid_of_the_plume_mask[1],
                      'confidence_index_in_perc' : confidence_index}
    
    data_to_return.update({f'SPM_threshold_{plume_name}': threshold for plume_name, threshold in thresholds.items()})
    
    return data_to_return


def create_polygon_mask(dataset, parameters) : 
    
    """
    Create a polygon mask that covers the searching area.

    Parameters
    ----------
    dataset : xarray.DataArray
        The input dataset to be used to extract polygon coordinates.
    parameters : dict
        Configuration parameters for plume detection. It contains the limits of the polygon defining the searching area.

    Returns
    -------
    xarray.DataArray
        Boolean mask where True represents areas inside the polygon.
    """
    
    # If the searching zone of the plume area (lat_range_to_search_plume_area & lon_range_to_search_plume_area) has more than two coordinates, create a polygon mask
    if len(parameters['lat_range_to_search_plume_area']) > 2 : 

        # Create polygon coordinates from the lat/lon ranges provided
        polygon_coordinates = [(lon, lat) for lat, lon in zip( [*parameters['lat_range_to_search_plume_area'], parameters['lat_range_to_search_plume_area'][0]],
                                                               [*parameters['lon_range_to_search_plume_area'], parameters['lon_range_to_search_plume_area'][0]])]
        
        # Create the polygon object
        polygon = Polygon(polygon_coordinates)
    
        # Create a grid of lon/lat and check if points are inside the polygon
        lon, lat = np.meshgrid(dataset.coords['lon'], dataset.coords['lat'])
        inside_polygon_mask = contains(polygon, lon, lat)
        
    else : 
        
         # If it's a simple range, create a rectangular mask based on the lat/lon boundaries
         inside_polygon_mask = xr.zeros_like(dataset).astype(bool)
         inside_polygon_mask.values = ( (inside_polygon_mask.lat >= parameters['lat_range_to_search_plume_area'][0]) & 
                                        (inside_polygon_mask.lat <= parameters['lat_range_to_search_plume_area'][1]) &
                                        (inside_polygon_mask.lon >= parameters['lon_range_to_search_plume_area'][0]) & 
                                        (inside_polygon_mask.lon <= parameters['lon_range_to_search_plume_area'][1]) ) 
         
    return inside_polygon_mask


def preprocess_annual_dataset_and_compute_land_mask(path_to_annual_ds, parameters) :

    """
    Load the annual map and use it to generate a land mask.

    This function processes an annual dataset by loading it, reducing its resolution, 
    and generating a land mask based on null values.

    Parameters
    ----------
    path_to_annual_ds : str
        the path to the annual (yearly) averaged dataset
    parameters : dict
        Configuration parameters for plume detection. It contains the limits of the polygon defining the searching area.
    is_SEXTANT_file : bool
        Boolean value indicating if the data source is from SEXTANT databases.

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
    
    # Reduce the resolution of the annual dataset
    if parameters['lat_new_resolution'] is not None : 
        multi_annual_ds_reduced = reduce_resolution(multi_annual_ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])

    # Create a land mask based on null values (land = True)
    land_mask = multi_annual_ds_reduced.isnull()
    
    return multi_annual_ds_subset, land_mask


def Check_if_the_area_is_too_cloudy(dataset, map_wo_clouds, parameters) :

    """
    Check if a specific area within the dataset is too cloudy based on a threshold.

    This function evaluates cloud coverage within a defined sub-area of the dataset
    by comparing cloudy pixel counts to a user-defined threshold.

    Parameters
    ----------
    dataset : xarray.DataArray
        The dataset containing the data to evaluate for cloud coverage.
    map_wo_clouds : xarray.DataArray
        A reference dataset that indicates the valid data points for comparison.
    parameters : dict
        Dictionary of parameters, including:
        - 'lat_range_of_the_area_to_check_for_clouds': tuple of float
            Latitude range of the sub-area to evaluate.
        - 'lon_range_of_the_area_to_check_for_clouds': tuple of float
            Longitude range of the sub-area to evaluate.
        - 'threshold_of_cloud_coverage_in_percentage': float
            Percentage threshold to determine if the area is too cloudy.

    Returns
    -------
    bool
        True if the area is too cloudy (exceeds the threshold), False otherwise.
    """    

    # Select the sub-area of the dataset for cloud cover evaluation based on latitude and longitude ranges.
    sub_area_to_check_for_zeros = dataset.sel(lat=slice(parameters['lat_range_of_the_area_to_check_for_clouds'][0],
                                                        parameters['lat_range_of_the_area_to_check_for_clouds'][1]), 
                                              lon=slice(parameters['lon_range_of_the_area_to_check_for_clouds'][0], 
                                                        parameters['lon_range_of_the_area_to_check_for_clouds'][1]))
        
    # Count the number of cloudy pixels within the sub-area.
    # Cloudy pixels are null in the current dataset but valid in the annual dataset.    
    n_cloudy_pixels = (sub_area_to_check_for_zeros.isnull() & map_wo_clouds.notnull()).sum().item()
    
    # Count the total number of valid pixels in the reference dataset.
    n_total_pixel = map_wo_clouds.notnull().sum().item()
    
    # Calculate the cloud coverage percentage and compare it to the threshold.
    # If the percentage of cloudy pixels exceeds the threshold, return True.
    # if n_total_pixel == 0 :
    #     test = True
    # else :
    #     test = (100 * n_cloudy_pixels / n_total_pixel) > parameters['threshold_of_cloud_coverage_in_percentage']
    cloud_test = (100 * n_cloudy_pixels / n_total_pixel) > parameters['threshold_of_cloud_coverage_in_percentage']
    
    # Return the result indicating whether the area is too cloudy.
    return cloud_test


def fast_delimitation_of_a_river_plume_area(spm_map, land_mask, start_point, SPM_threshold, maximal_threshold, max_steps = 20) : 
    
    """
    Delimits a river plume area based on a Suspended Particulate Matter (SPM) map 
    by computing gradients in specified directions and creating a polygon boundary.

    Parameters
    ----------
    spm_map : xarray.DataArray
        A 2D array representing the SPM concentration map.
    land_mask : xarray.DataArray
        A boolean mask where True indicates land areas.
    start_point : tuple
        Coordinates of the river mouth (row, column) to begin plume delimitation.
    max_steps : int, optional
        Maximum number of steps to compute gradients in each direction, by default 20.
    SPM_threshold : float, optional
        Threshold for SPM concentration to distinguish plume areas, by default 0.5.

    Returns
    -------
    matplotlib.path.Path
        A polygon path object representing the delimited river plume area.

    Notes
    -----
    - The function relies on computing gradients in multiple directions from the 
      start point to determine the boundaries of the plume.
    - Gradients are analyzed to determine areas where SPM concentration decreases sharply.
    - A convex hull is constructed to define the plume boundary.

    Examples
    --------
    >>> plume_path = fast_delimitation_of_a_river_plume_area(
    ...     spm_map=spm_map,
    ...     land_mask=land_mask,
    ...     start_point=(50, 50),
    ...     max_steps=30,
    ...     SPM_threshold=0.5
    ... )
    >>> plt.imshow(spm_map, cmap='viridis')
    >>> plt.plot(plume_path.vertices[:, 0], plume_path.vertices[:, 1], 'r-')
    >>> plt.show()
    """
    
    # Define directions for inspecting neighboring pixels
    directions = coordinates_of_pixels_to_inspect( {"directions" : {'grid' : np.array([ [True,  True,  True],
                                                                                        [True,  True,  True],
                                                                                        [True,  True,  True],
                                                                                      ]),
                                                                    'coordinates_of_center' : (1,1)}} )
    directions['directions'].append(directions['directions'][0]) # Ensure loop closure for directions

    # Compute gradients and points along specified directions
    gradient_values, direction_points, absolute_values = compute_gradient_with_directions_vectorized(spm_map = spm_map, 
                                                                                   start_point = start_point, 
                                                                                   directions = directions['directions'],
                                                                                   max_steps = max_steps, # Max steps for gradient expansion
                                                                                   lower_high_values_to = np.min( [SPM_threshold*5, maximal_threshold] ),
                                                                                   create_X_intermediates_between_each_direction = 2) # Number of intermediate directions

    if gradient_values is None : 
        # mask_values_equal_to_SPM_threshold = absolute_values.values == SPM_threshold*1.5
        # max_steps_for_each_direction = np.array([ first_true_block(arr)[1] for arr in mask_values_equal_to_SPM_threshold ])
        return None  

    # Define a threshold for detecting significant changes in SPM gradient
    # threshold_on_gradient_values = -0.15 # 0.2
    threshold_on_gradient_values = np.nanmax(gradient_values[np.isfinite(gradient_values)]) * 0.9

    # # Plot the direction points on the SPM map
    # plt.imshow(spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = 8)
    # for x, y in direction_points.reshape(-1, direction_points.shape[-1]) :
    #     plt.scatter(x, y, s=1)
    # plt.title("Gradient Points Following Direction Vectors")
    # plt.show()
        
    # # Plot gradient values over distance for each direction
    # for i in range(gradient_values.shape[0]):  # Loop over the 3rd dimension (directions)
    #     plt.plot(gradient_values[i, :], label=f'Direction {i+1}')
    # plt.hlines(threshold_on_gradient_values, xmin = 0, xmax = max_steps)
    # plt.title("Gradient Computation in Specified Directions")
    # plt.xlabel("Step")
    # plt.ylabel("Gradient of SPM")
    # plt.show()   
    
    # Mask 1: Detect areas below the gradient threshold                
    mask_lower_than_threshold = ( gradient_values < threshold_on_gradient_values)
    test_to_select_the_max_step = np.cumsum(mask_lower_than_threshold, axis = -1)
    max_steps_for_each_direction = np.argmax( test_to_select_the_max_step == 1 , axis = -1 ) 
       
    # test_to_select_the_max_step[test_to_select_the_max_step > 2] = 2
    

    # Mask 2: Adjust max steps for directions touching land or exceeding SPM threshold
    for index in np.where(max_steps_for_each_direction == 0)[0] :
        
        valid_points = direction_points[index][(direction_points[index] >= 0).all(axis=1)]                    
        y_coords, x_coords = valid_points[:, 1].astype(int), valid_points[:, 0].astype(int)
    
        if any(land_mask.values[y_coords, x_coords]) : 
    
            combined_condition = (land_mask.values[y_coords, x_coords] == False) & \
                                 (spm_map.values[y_coords, x_coords] > SPM_threshold)
        
            max_steps_for_each_direction[index] = (np.min( [last_true_block(combined_condition)[1] + 1, len(combined_condition) -1] ) 
                                                   if combined_condition.any() else 0)
    
    # Compute polygon boundaries based on max steps
    polygon_boundaries = np.array([ list(direction_points[index][value].astype(int)) for index, value in enumerate(max_steps_for_each_direction)])
    
    if len(np.unique(polygon_boundaries)) <= 4 : 
        return None
        
    # Create a convex hull from the boundary points
    polygon = ConvexHull( polygon_boundaries )
    polygon_boundaries = polygon_boundaries[ np.append(polygon.vertices, polygon.vertices[0]) ]
    
    # Create a polygon path object
    polygon_path = Path(polygon_boundaries, closed=True)
    
    # # Plot the final plume boundary
    # plt.imshow(spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = 8)
    # plt.title("Gradient Points Following Direction Vectors")
    # ax = plt.gca()  # Get the current Axes instance
    # ax.add_patch(PathPatch(polygon_path, facecolor='none', edgecolor='red', lw=2))  # Add the polygon path to the plot
    # plt.show()
                        
    return polygon_path


# =============================================================================
#### Main function
# =============================================================================

def main_process(file_name, file_names_pattern,
                 parameters,
                 bathymetry_data_aligned_to_reduced_map,
                 france_shapefile,
                 map_wo_clouds, land_mask,
                 inside_polygon_mask,
                 plume_dir,
                 regional_map_dir,
                 dynamic_thresh):
    """
    Process a single satellite data file for plume detection.

    This function detects and analyzes plumes in a satellite data file by
    applying thresholds, masks, and other criteria.

    Parameters
    ----------
    file_name : str
        Path to the satellite data file.
    file_names_pattern : str
        Pattern for matching files in directories.
    parameters : dict
        Configuration parameters for plume detection.
    bathymetry_data_aligned_to_reduced_map : xr.DataArray
        Bathymetric data aligned to the reduced resolution.
    is_SEXTANT_file : bool
        Indicates if the file comes from the SEXTANT database.
    france_shapefile : gpd.GeoDataFrame
        Shapefile of France for plotting.
    map_wo_clouds : xr.DataArray
        Dataset for cloud coverage checks.
    land_mask : xr.DataArray
        Mask identifying land pixels.
    inside_polygon_mask : xr.DataArray
        Mask to limit processing area inside a specific polygon.

    Returns
    -------
    dict
        Processed results, including plume statistics and confidence index.
    """
    # file_name = "/home/calanus/RiOMar/output/REGIONAL_MAPS/GULF_OF_LION/SEXTANT/SPM/merged/Standard/MAPS/WEEKLY/1998/Averaged_over_01_01.pkl"
    # file_name = "/home/calanus/RiOMar/output/REGIONAL_MAPS/GULF_OF_LION/SEXTANT/SPM/merged/Standard/MAPS/WEEKLY/2013/Averaged_over_01_01.pkl"
    # print(file_name)

    # Define path to save the figure generated during processing
    path_to_the_figure_file_to_save = f'{os.path.dirname(file_name).replace("MAPS", "PLUME_DETECTION").replace(
        regional_map_dir, plume_dir)}/MAPS/{os.path.basename(file_name)}'.replace('.pkl', '')

    # Ensure the directory for saving figures exists
    os.makedirs(f'{os.path.dirname(path_to_the_figure_file_to_save)}', exist_ok=True)

    # Open and load the file (binary file assumed to contain data)
    with open(file_name, 'rb') as f:
        ds = pickle.load(f)['Basin_map']['map_data']

    # Reduce the resolution of the dataset to the specified latitude and longitude resolutions
    ds_reduced = (reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])
                  if parameters['lat_new_resolution'] is not None
                  else ds)

    # Initialize a list to store masks for all plume areas
    all_mask_area = []
    all_river_mouth_to_remove = []
    thresholds = {key: None for key in parameters['starting_points']}

    # If the percentage of cloud coverage exceeds the specified threshold, return default values without processing the plume area
    if (Check_if_the_area_is_too_cloudy(ds, map_wo_clouds, parameters)):
        # Plot the map with no plume area (due to clouds)
        make_the_plot(path_to_the_figure_file_to_save, ds, ds_reduced, france_shapefile,
                      parameters['lon_range_of_the_map_to_plot'],
                      parameters['lat_range_of_the_map_to_plot'],
                      parameters['bathymetric_threshold'],
                      bathymetry_data_aligned_to_reduced_map,
                      thresholds,
                      plot_the_plume_area=False)

        data_to_return = return_stats_dictionnary(None, ds_reduced, ds, parameters,
                                                  thresholds, return_empty_dict=True)

        return data_to_return

    # Loop through each plume starting point to process plume detection
    for plume_name, starting_point in parameters['starting_points'].items():

        the_plume = Pipeline_to_delineate_the_plume(ds_reduced,
                                                    bathymetry_data_aligned_to_reduced_map,
                                                    land_mask,
                                                    parameters,
                                                    plume_name,
                                                    inside_polygon_mask,
                                                    dynamic_thresh)

        if the_plume is None:
            continue

        thresholds[plume_name] = the_plume.SPM_threshold
        all_mask_area.append(the_plume.plume_mask)
        if "close_river_mouth_mask" in vars(the_plume):
            all_river_mouth_to_remove.append(the_plume.close_river_mouth_mask)

    # If no valid plume area is detected, plot and return default values
    if (len(all_mask_area) == 0) or (any([x.values.any() for x in all_mask_area]) == False):
        make_the_plot(path_to_the_figure_file_to_save, ds, ds_reduced, france_shapefile,
                      lon_range_of_the_map_to_plot=parameters['lon_range_of_the_map_to_plot'],
                      lat_range_of_the_map_to_plot=parameters['lat_range_of_the_map_to_plot'],
                      bathymetric_threshold=parameters['bathymetric_threshold'],
                      bathymetry_data_aligned_to_reduced_map=bathymetry_data_aligned_to_reduced_map,
                      thresholds=thresholds,
                      plot_the_plume_area=False)

        data_to_return = return_stats_dictionnary(None, ds_reduced, ds, parameters,
                                                  thresholds, return_empty_dict=True)

        del all_mask_area
        gc.collect()

        return data_to_return

    # Combine all detected plume areas using logical OR
    final_mask_area = reduce(np.logical_or, all_mask_area)
    final_close_river_mouth_area = reduce(np.logical_or, all_river_mouth_to_remove)

    data_to_return = return_stats_dictionnary(final_mask_area, ds_reduced, ds, parameters, thresholds)

    # Plot the final map with the plume area
    make_the_plot(path_to_the_figure_file_to_save, ds, ds_reduced, france_shapefile,
                  lon_range_of_the_map_to_plot=parameters['lon_range_of_the_map_to_plot'],
                  lat_range_of_the_map_to_plot=parameters['lat_range_of_the_map_to_plot'],
                  bathymetric_threshold=parameters['bathymetric_threshold'],
                  bathymetry_data_aligned_to_reduced_map=bathymetry_data_aligned_to_reduced_map,
                  thresholds=thresholds,
                  plot_the_plume_area=True,
                  mask_area=final_mask_area,
                  close_river_mouth_area=final_close_river_mouth_area,
                  pixel_done=None,
                  show_bathymetric_mask=True)

    # Cleanup to free memory
    del final_mask_area, all_mask_area
    gc.collect()

    return data_to_return


# =============================================================================
#### Classes 
# =============================================================================

class Create_the_plume_mask : 
    
    """
    A class to detect and refine plume areas based on SPM, bathymetry, and land masks.

    Attributes
    ----------
    spm_map : xr.DataArray
        Map of [SPM] concentrations.
    bathymetry_map : xr.DataArray
        Aligned map of bathymetry values.
    land_mask : xr.DataArray
        Boolean map indicating land regions (True = land).
    protocol : list of str
        List of executed methods in order.
    plume_name : str
        Name of the plume being processed.
    parameters : dict
        Configuration and parameters for plume detection.

    Methods
    -------
    determine_SPM_threshold(manual_determination_of_SPM_threshold=None)
        Determines the SPM threshold for plume detection.
    do_a_raw_plume_detection()
        Performs initial plume area detection.
    include_cloudy_regions_to_plume_area()
        Expands the plume area to include cloudy regions.
    remove_the_areas_with_sediment_resuspension(maximal_bathymetry=None, minimal_distance_from_estuary=None)
        Removes areas with sediment resuspension.
    remove_shallow_waters(bathymetric_threshold=None)
        Removes shallow water areas from the plume.
    dilate_the_main_plume_area_to_merge_close_plume_areas()
        Merges closely located plume areas.
    remove_small_shapes_that_do_not_meet_a_minimum_size_criterion(minimum_size_threshold=3)
        Removes small shapes that do not meet a size threshold.
    set_pixels_to_False_if_outside_of_the_searching_area(searching_area)
        Masks plume pixels outside the specified searching area.
    identify_the_main_plume_shape_based_on_the_plume_core_location(plume_core_location=None)
        Identifies the main plume shape based on core location.
    remove_parts_of_the_plume_area_identified_only_on_the_edge_of_the_searching_area()
        Removes parts of the plume area identified on the search area edge.
    remove_parts_of_the_plume_area_with_very_high_SPM_on_the_edge_of_the_searching_zone()
        Removes plume areas with high SPM values near the edge of the search zone.
    remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase()
        Removes plume areas that widen after the shrinking phase.
    """
    
    def __init__(self, spm_reduced_map, bathymetry_map_aligned_to_spm_reduced_map, land_mask, parameters, plume_name) :
        
        """
        Initialize the plume detection class.

        Parameters
        ----------
        spm_reduced_map : xr.DataArray
            Map of [SPM] concentrations.
        bathymetry_map_aligned_to_spm_reduced_map : xr.DataArray
            Map of bathymetry values aligned to the SPM map.
        land_mask : xr.DataArray
            Boolean map where True indicates land.
        parameters : dict
            Configuration and parameters for plume detection.
        plume_name : str
            Name of the plume being processed.
        """
        
        self.spm_map = spm_reduced_map
        self.bathymetry_map = bathymetry_map_aligned_to_spm_reduced_map
        self.land_mask = land_mask
        self.protocol = []
        self.plume_name = plume_name
        self.parameters = parameters
        self.parameters["pixel_starting_points"] = {plume_name : find_the_index_of_the_plume_starting_point(self.spm_map, 
                                                                                                           self.parameters['starting_points'][plume_name])}
        self.parameters["pixel_starting_points_close_river_mouth"] = {river_mouth : find_the_index_of_the_plume_starting_point(self.spm_map, coordinates) 
                                                                          for river_mouth, coordinates in self.parameters['river_mouth_to_exclude'].items()}

    def determine_SPM_threshold(self, dynamic_determination_of_SPM_threshold) : 
        
        """
        Determine the SPM threshold for plume detection.

        Parameters
        ----------
        manual_determination_of_SPM_threshold : float, optional
            User-specified SPM threshold. If False, the threshold is determined automatically.
        """
        
        # spm_map = self.spm_map
        # land_mask = self.land_mask
        # start_point = self.parameters['pixel_starting_points'][self.plume_name]
        # directions = self.parameters['searching_strategy_directions'][self.plume_name]
        # max_steps = self.parameters['max_steps_for_the_directions'][self.plume_name]
        # maximal_threshold = self.parameters['maximal_threshold'][self.plume_name]
        # minimal_threshold = self.parameters['minimal_threshold'][self.plume_name]
        # quantile_to_use = self.parameters['quantile_to_use'][self.plume_name]
        
        if dynamic_determination_of_SPM_threshold : 
            SPM_threshold, points_used_for_finding_SPM_threshold, all_points_tested = find_SPM_threshold(spm_map = self.spm_map, 
                                                land_mask = self.land_mask,
                                                start_point = self.parameters['pixel_starting_points'][self.plume_name], 
                                                directions = self.parameters['searching_strategy_directions'][self.plume_name],
                                                max_steps = self.parameters['max_steps_for_the_directions'][self.plume_name],
                                                maximal_threshold = self.parameters['maximal_threshold'][self.plume_name],
                                                minimal_threshold = self.parameters['minimal_threshold'][self.plume_name],
                                                quantile_to_use = self.parameters['quantile_to_use'][self.plume_name])
            self.points_used_for_finding_SPM_threshold = points_used_for_finding_SPM_threshold
            
            all_points_tested = all_points_tested.reshape(-1,2)
            all_points_tested = all_points_tested[~np.isnan(all_points_tested).any(axis=1)]
            self.all_points_tested_for_finding_SPM_threshold = all_points_tested
            
        else : 
            SPM_threshold = self.parameters['fixed_threshold'][self.plume_name]
                
        self.SPM_threshold = SPM_threshold
        self.protocol.append(f'{len(self.protocol)} : determine_SPM_threshold')
                
        
    def do_a_raw_plume_detection(self) : 
        
        """
        Perform the initial detection of the plume area.
        """
        
        # data = self.spm_map.values
        # start = self.parameters['pixel_starting_points'][self.plume_name]
        # SPM_threshold = self.SPM_threshold
        # directions = self.parameters['searching_strategy_directions'][self.plume_name]
        
        mask, pixel_done = flood_fill(data = self.spm_map.values, 
                                      start = self.parameters['pixel_starting_points'][self.plume_name], 
                                      SPM_threshold = self.SPM_threshold,
                                      directions = self.parameters['searching_strategy_directions'][self.plume_name])  
        
        self.plume_mask = xr.DataArray(mask, coords=self.spm_map.coords, dims=self.spm_map.dims)
        self.protocol.append(f'{len(self.protocol)} : do_a_raw_plume_detection')
        
    def include_cloudy_regions_to_plume_area(self) : 
        
        """
        Expand the detected plume area to include cloudy regions.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        self.plume_mask = Set_cloudy_regions_to_True(self.spm_map, self.plume_mask, self.land_mask, self.SPM_threshold) 
        
        self.protocol.append(f'{len(self.protocol)} : include_cloudy_regions_to_plume_area')
        
    def remove_the_areas_with_sediment_resuspension(self, maximal_bathymetry = None, minimal_distance_from_estuary = None) : 
        
        """
        Remove areas with sediment resuspension from the plume mask.

        Parameters
        ----------
        maximal_bathymetry : float, optional
            Maximum bathymetry value to consider as part of sediment resuspension zones.
            If None, it is retrieved from the parameters dictionary.
        minimal_distance_from_estuary : float, optional
            Minimum distance from the estuary to consider as part of sediment resuspension zones.
            If None, it is retrieved from the parameters dictionary.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        if maximal_bathymetry is None : 
            maximal_bathymetry = self.parameters['maximal_bathymetric_for_zone_with_resuspension'][self.plume_name]
        
        if minimal_distance_from_estuary is None : 
            minimal_distance_from_estuary = self.parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][self.plume_name]
        
        self.plume_mask = remove_coastal_areas_with_sediment_resuspension(self.spm_map, self.plume_mask, 
                                                                    self.land_mask, self.bathymetry_map,
                                                                    self.SPM_threshold, 
                                                                    self.parameters['core_of_the_plumes'][self.plume_name],
                                                                    maximal_bathymetry,
                                                                    minimal_distance_from_estuary)
    
        self.protocol.append(f'{len(self.protocol)} : remove_the_areas_with_sediment_resuspension')
    
    def remove_shallow_waters(self, bathymetric_threshold = None) :
        
        """
        Remove shallow water regions from the plume mask.

        Parameters
        ----------
        bathymetric_threshold : float, optional
            Bathymetric depth threshold. Pixels with bathymetry values greater than
            the negative of this threshold are excluded from the plume mask.
            If None, the value is retrieved from the parameters dictionary.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        if bathymetric_threshold is None :
            bathymetric_threshold = self.parameters['bathymetric_threshold']
        
        self.plume_mask.values[self.bathymetry_map.values > -bathymetric_threshold] = False
        self.protocol.append(f'{len(self.protocol)} : remove_shallow_waters')
    
    
    def dilate_the_main_plume_area_to_merge_close_plume_areas(self) :
    
        """
        Dilate the plume mask to merge closely located plume areas.
        """    
    
        if self.plume_mask.any() == False : 
            return None
        
        self.plume_mask = merge_plume_shape_with_close_shapes(self.plume_mask, 
                                                        self.parameters['core_of_the_plumes'][self.plume_name],
                                                        self.land_mask,
                                                        structure_of_the_dilation = np.array([[True, True, True, True, True],
                                                                                              [True, True, True, True, True],
                                                                                              [True, True, True, True, True],
                                                                                              [True, True, True, True, True],
                                                                                              [True, True, True, True, True]]))
        
        self.protocol.append(f'{len(self.protocol)} : Dilate_the_main_plume_area_to_merge_close_plume_areas')
        
    
    def remove_small_shapes_that_do_not_meet_a_minimum_size_criterion(self, minimum_size_threshold = 3) : 
        
        """
        Remove small shapes from the plume mask based on a size criterion.

        Parameters
        ----------
        minimum_size_threshold : int, optional
            Minimum size (number of pixels) for a shape to be retained in the plume mask.
            Default is 3.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        self.plume_mask.values = morphology.remove_small_objects(self.plume_mask.values, minimum_size_threshold)     
        self.protocol.append(f'{len(self.protocol)} : Remove_small_shapes_that_do_not_meet_a_minimum_size_criterion')
        
    def set_pixels_to_False_if_outside_of_the_searching_area(self, searching_area) :
    
        """
        Exclude pixels outside the specified searching area from the plume mask.

        Parameters
        ----------
        searching_area : xr.DataArray
            Boolean mask indicating the valid searching area (True = valid).
        """    
    
        if self.plume_mask.any() == False : 
            return None
        
        self.plume_mask = self.plume_mask.where(searching_area, other=False)
        self.protocol.append(f'{len(self.protocol)} : Set_pixels_to_False_if_outside_of_the_searching_area')
        
    def identify_the_main_plume_shape_based_on_the_plume_core_location(self, plume_core_location = None) :
        
        """
        Identify the main plume shape using the plume core location.

        Parameters
        ----------
        plume_core_location : tuple of int, optional
            Coordinates of the plume core. If None, the value is retrieved
            from the parameters dictionary.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        if plume_core_location is None : 
            plume_core_location = self.parameters['core_of_the_plumes'][self.plume_name]
        
        label_of_the_shape_to_keep, labeled_array, num_features = identify_the_shape_label_corresponding_to_the_plume(self.plume_mask, plume_core_location)
        self.plume_mask.values = (labeled_array == label_of_the_shape_to_keep)
        self.protocol.append(f'{len(self.protocol)} : identify_the_main_plume_shape_based_on_the_plume_core_location')

    def remove_parts_of_the_plume_area_identified_only_on_the_edge_of_the_searching_area(self) :
        
        """
        Remove parts of the plume area that are identified only on the edge of the searching area.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        boolean_values_in_the_area_of_the_plume, \
            values_in_the_area_of_the_plume, \
            direction_points = find_index_and_values_of_multiple_directions_in_the_plume_area(mask_area = self.plume_mask, 
                                                                                              ds_reduced = self.spm_map,
                                                                                              pixel_starting_point = self.parameters['pixel_starting_points'][self.plume_name], 
                                                                                              searching_strategy_direction = self.parameters['searching_strategy_directions'][self.plume_name],
                                                                                              max_steps = self.parameters['max_steps_for_the_directions'][self.plume_name])
        
        if len(boolean_values_in_the_area_of_the_plume) == 0 : 
            return None
        
        edge_idd = int( round( boolean_values_in_the_area_of_the_plume.shape[0] / 3 ) / 2 )
        non_zero_in_middle = np.any( boolean_values_in_the_area_of_the_plume[edge_idd:-edge_idd, :] != 0, axis=0)
        non_zero_in_middle = first_true_block(non_zero_in_middle)
        boolean_values_in_the_area_of_the_plume[ :, (non_zero_in_middle[1] + 1):] = 0
        
        refine_removal_of_edge_effect = np.any(boolean_values_in_the_area_of_the_plume[[2,-2], :] != 0, axis=0)
        refine_removal_of_edge_effect = first_true_block(refine_removal_of_edge_effect)
        boolean_values_in_the_area_of_the_plume[ [0,1,-1,-2], (refine_removal_of_edge_effect[1] + 1):] = 0
       
        index_to_keep = np.where(boolean_values_in_the_area_of_the_plume == 1)
        
        if len(index_to_keep[0]) == 0 :
            return None

        self.plume_mask = set_mask_area_values_to_False_based_on_an_index_object(self.plume_mask, 
                                                                           index_object = index_to_keep, 
                                                                           gradient_points = direction_points)
        
        self.protocol.append(f'{len(self.protocol)} : remove_parts_of_the_plume_area_which_have_very_high_SPM_on_the_edge_of_the_searching_area')
        
    def remove_parts_of_the_plume_area_with_very_high_SPM_on_the_edge_of_the_searching_zone(self) :
    
        """
        Remove parts of the plume area with very high SPM values near the edge of the searching zone.
        """    
    
        if self.plume_mask.any() == False : 
            return None
    
        boolean_values_in_the_area_of_the_plume, \
            values_in_the_area_of_the_plume, \
            direction_points = find_index_and_values_of_multiple_directions_in_the_plume_area(self.plume_mask, 
                                                                                                self.spm_map,
                                                                                                self.parameters['pixel_starting_points'][self.plume_name], 
                                                                                                self.parameters['searching_strategy_directions'][self.plume_name],
                                                                                                self.parameters['max_steps_for_the_directions'][self.plume_name])
    
        if len(boolean_values_in_the_area_of_the_plume) == 0 : 
            return None
        
        SPM_threshold_on_the_edge = self.SPM_threshold * 1.5
        test_values = np.array([ (values_in_the_area_of_the_plume[0][i] > SPM_threshold_on_the_edge) or 
                                (values_in_the_area_of_the_plume[-1][i] > SPM_threshold_on_the_edge) for  
                                i in np.arange( values_in_the_area_of_the_plume.shape[1] ) ])
        test_values[ :first_true_block( test_values[5:] )[1]+1 +5 ] = False
        
        # Process each selected row
        if test_values.any() : 
            
            for index in np.arange( np.where(test_values)[0][0], len(test_values) ) :
                
                sequence_values = values_in_the_area_of_the_plume[:, index]
                sequence_plume_idd = boolean_values_in_the_area_of_the_plume[:, index]
                
                # Find the indices where the values exceed the threshold
                # above_threshold = sequence_values > SPM_threshold
                above_threshold = sequence_plume_idd == 1
                     
                first_block_index = first_true_block( above_threshold )
                last_block_index = last_true_block( above_threshold )
                
                if first_block_index[0] == 0 : 
                    boolean_values_in_the_area_of_the_plume[:first_block_index[1]+1, index] = 0
                
                if last_block_index[1] == len(sequence_values)-1 : 
                    boolean_values_in_the_area_of_the_plume[last_block_index[0]:, index] = 0
            
            index_to_keep = np.where(boolean_values_in_the_area_of_the_plume == 1)
    
            self.plume_mask = set_mask_area_values_to_False_based_on_an_index_object(self.plume_mask, 
                                                                                     index_object = index_to_keep, 
                                                                                     gradient_points = direction_points)
            
        self.protocol.append(f'{len(self.protocol)} : remove_parts_of_the_plume_area_with_very_high_SPM_on_the_edge_of_the_searching_zone')

    def remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase(self) : 
         
        """
        Remove parts of the plume area that widen after the shrinking phase.
        """
        
        if self.plume_mask.any() == False : 
            return None
        
        boolean_values_in_the_area_of_the_plume, \
            values_in_the_area_of_the_plume, \
            direction_points = find_index_and_values_of_multiple_directions_in_the_plume_area(self.plume_mask, 
                                                                                                self.spm_map,
                                                                                                self.parameters['pixel_starting_points'][self.plume_name], 
                                                                                                self.parameters['searching_strategy_directions'][self.plume_name],
                                                                                                self.parameters['max_steps_for_the_directions'][self.plume_name])
        
        if len(boolean_values_in_the_area_of_the_plume) == 0 : 
            return None
        
        row_sums = np.sum(boolean_values_in_the_area_of_the_plume == 1, axis=0)
        increase_indices = np.where(np.diff(row_sums) >= 1)[0] + 1
        indices_to_keep = np.where(increase_indices > (boolean_values_in_the_area_of_the_plume.shape[1] * 2/3))[0]
        
        if len(indices_to_keep) > 0 : 

            increase_indices_to_keep = np.where( np.full((direction_points.shape[0], 1), increase_indices[ np.min( indices_to_keep ) ] - 1) 
                                                    > np.arange(direction_points.shape[1]) )
            
            self.plume_mask = set_mask_area_values_to_False_based_on_an_index_object(self.plume_mask, increase_indices_to_keep, direction_points)
            
        self.protocol.append(f'{len(self.protocol)} : remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase')
        
    def remove_close_river_mouth(self, positions_of_river_mouths = {}) :
        
        """
        Removes the area of a close river plume from the current plume mask.
    
        This function identifies and removes regions of overlapping plumes associated 
        with nearby river mouths. It uses the plume boundary determined by 
        `fast_delimitation_of_a_river_plume_area` and updates the plume mask accordingly.
    
        Parameters
        ----------
        positions_of_river_mouths : dict, optional
            A dictionary of positions of river mouths. Each key corresponds to a river mouth's
            identifier, and the value is a tuple of coordinates (row, column).
    
        Returns
        -------
        None
            The method modifies `self.plume_mask` in place, removing the areas of close river plumes.
            It also updates the `self.protocol` attribute to log the operation.
    
        Notes
        -----
        - This function assumes the presence of `self.spm_map`, `self.land_mask`, 
          `self.plume_mask`, and `self.SPM_threshold` attributes.
        - The plume boundaries are determined using a convex hull polygon.
        """
           
        if self.plume_mask.any() == False : 
            return None
        
        self.close_river_mouth_mask = xr.zeros_like(self.plume_mask)
        
        # Loop through each provided river mouth position
        for i, position_of_a_river_mouth in positions_of_river_mouths.items() :  
            
            # spm_map = self.spm_map
            # land_mask = self.land_mask
            # start_point = position_of_a_river_mouth
            # SPM_threshold = self.SPM_threshold
            # maximal_threshold = self.parameters['maximal_threshold'][self.plume_name]
            # max_steps = 20
            
            # Delimit the plume area for the current river mouth
            delimitation_of_the_plume =  fast_delimitation_of_a_river_plume_area(spm_map = self.spm_map, 
                                                                                land_mask = self.land_mask,
                                                                                start_point = position_of_a_river_mouth,
                                                                                SPM_threshold = self.SPM_threshold,
                                                                                maximal_threshold = self.parameters['maximal_threshold'][self.plume_name],
                                                                                max_steps = 20)
            
            if delimitation_of_the_plume is None : 
                continue
                                
            # # Optional: Uncomment the lines below to visualize the plume boundary
            # plt.imshow(self.spm_map, cmap='viridis', origin = "lower", vmin = 0.5, vmax = 8)
            # plt.title("Gradient Points Following Direction Vectors")
            # ax = plt.gca()  # Get the current Axes instance
            # ax.add_patch(PathPatch(delimitation_of_the_plume, facecolor='none', edgecolor='red', lw=2))  # Add the polygon path to the plot
            # plt.show()
            
            # Identify the coordinates of plume pixels currently marked as True
            true_indices = np.where(self.plume_mask.values)
            points = np.column_stack((true_indices[1], true_indices[0])) # Convert to (x, y) coordinates
            
            # Check if these points are inside the current plume polygon
            inside_polygon = delimitation_of_the_plume.contains_points(points) 
                
            # Update the plume mask to exclude points inside the polygon
            index_corresponding_to_the_close_river_mouth = points[ np.where(inside_polygon)[0] ][:,1], points[ np.where(inside_polygon)[0] ][:,0]
            self.plume_mask.values[ index_corresponding_to_the_close_river_mouth ] = False
            self.close_river_mouth_mask.values[ index_corresponding_to_the_close_river_mouth ] = True

        # Log the operation in the protocol
        self.protocol.append(f'{len(self.protocol)} : remove_the_river_plume_from_the_mouth_of_the_neighboring_river')


    def do_R_plot(self, where_to_save_the_plot, name_of_the_plot):
        """
        Convert the SPM map and plume mask to an R dataframe and plot using ggplot2.
        """
          
        folder_where_to_save_data = os.path.join(where_to_save_the_plot, 'DATA')
        os.makedirs( folder_where_to_save_data, exist_ok=True )
          
        # Create a Pandas DataFrame
        lat_values = self.spm_map.lat.values
        lon_values = self.spm_map.lon.values
        
        df = pd.DataFrame({
            'lat': np.repeat(lat_values, len(lon_values)),
            'lon': np.tile(lon_values, len(lat_values)),
            'analysed_spim': self.spm_map.values.flatten()
        })
        
        if name_of_the_plot == 'B' :
            
            latitudes_used_for_finding_SPM_threshold = lat_values[self.points_used_for_finding_SPM_threshold[:,1].astype(int)]
            longitudes_used_for_finding_SPM_threshold = lon_values[self.points_used_for_finding_SPM_threshold[:,0].astype(int)]
            points_used_for_finding_SPM_threshold = pd.DataFrame({'longitude': longitudes_used_for_finding_SPM_threshold, 
                                                                  'latitude': latitudes_used_for_finding_SPM_threshold})
            points_used_for_finding_SPM_threshold.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}_points_used_for_finding_SPM_threshold.csv") )

            all_latitudes_used_for_finding_SPM_threshold = lat_values[self.all_points_tested_for_finding_SPM_threshold[:,1].astype(int)]
            all_longitudes_used_for_finding_SPM_threshold = lon_values[self.all_points_tested_for_finding_SPM_threshold[:,0].astype(int)]
            all_points_used_for_finding_SPM_threshold = pd.DataFrame({'longitude': all_longitudes_used_for_finding_SPM_threshold, 
                                                                  'latitude': all_latitudes_used_for_finding_SPM_threshold})
            all_points_used_for_finding_SPM_threshold.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}_all_points_used_for_finding_SPM_threshold.csv") )
                    
        
        if 'plume_mask' in vars(self) :
            df['plume'] = self.plume_mask.values.flatten()
        
        index_to_keep = np.where((df.lat >= self.parameters['lat_range_of_the_map_to_plot'][0]) & 
                                 (df.lat <= self.parameters['lat_range_of_the_map_to_plot'][1]) & 
                                 (df.lon >= self.parameters['lon_range_of_the_map_to_plot'][0]) & 
                                 (df.lon <= self.parameters['lon_range_of_the_map_to_plot'][1]) )
        
        df = df.iloc[index_to_keep]

        df.to_csv( os.path.join(folder_where_to_save_data, f"{name_of_the_plot}.csv") )
        
        # Source the R script
        figure_R_path = os.path.join(func_dir, 'figure.R')
        robjects.r['source'](figure_R_path)
        # robjects.r['source']("myRIOMAR_dev/_5_Figures_for_article/utils.R")

        r_function = robjects.r['Figure_4']

        # Call the R function
        r_function(
            where_to_save_the_figure = robjects.StrVector([where_to_save_the_plot]),
            name_of_the_plot = robjects.StrVector([name_of_the_plot])
        )


# =============================================================================
#### Pipeline function
# =============================================================================


def Pipeline_to_delineate_the_plume(ds_reduced,
                                    bathymetry_data_aligned_to_reduced_map,
                                    land_mask,
                                    parameters,
                                    plume_name,
                                    inside_polygon_mask,
                                    dynamic_thresh):
    the_plume = Create_the_plume_mask(ds_reduced,
                                      bathymetry_data_aligned_to_reduced_map,
                                      land_mask,
                                      parameters,
                                      plume_name)

    the_plume.determine_SPM_threshold(dynamic_thresh)
    # the_plume.SPM_threshold = 10

    the_plume.do_a_raw_plume_detection()
    # the_plume.plume_mask.plot()

    the_plume.include_cloudy_regions_to_plume_area()

    the_plume.remove_the_areas_with_sediment_resuspension(
        maximal_bathymetry=parameters['maximal_bathymetric_for_zone_with_resuspension'][plume_name],
        minimal_distance_from_estuary=parameters['minimal_distance_from_estuary_for_zone_with_resuspension'][
            plume_name])

    the_plume.remove_shallow_waters()

    the_plume.remove_close_river_mouth(the_plume.parameters['pixel_starting_points_close_river_mouth'])

    the_plume.dilate_the_main_plume_area_to_merge_close_plume_areas()

    the_plume.remove_small_shapes_that_do_not_meet_a_minimum_size_criterion()

    the_plume.set_pixels_to_False_if_outside_of_the_searching_area(inside_polygon_mask)

    the_plume.identify_the_main_plume_shape_based_on_the_plume_core_location()

    the_plume.remove_shallow_waters()

    # the_plume.remove_parts_of_the_plume_area_identified_only_on_the_edge_of_the_searching_area()

    # the_plume.remove_parts_of_the_plume_area_with_very_high_SPM_on_the_edge_of_the_searching_zone()

    if not np.isin(plume_name, ['Seine']):
        the_plume.remove_parts_of_the_plume_area_that_widden_after_the_shrinking_phase()

    # the_plume.plume_mask.plot()
    # the_plume.close_river_mouth_mask.plot()
    # the_plume.protocol
    return the_plume


# =============================================================================
#### Workflow functions
# =============================================================================


def apply_plume_mask(core_arguments, Zones, time_step, nb_cores, 
                     dynamic_thresh, regional_map_dir, plume_dir):
      
    """
    Apply plume detection and analysis to satellite data.

    This function processes satellite images to detect and analyze river plumes. It handles data loading,  
    plume detection, statistical analysis, and visualization of results.

    Parameters
    ----------
    core_arguments : dict
        Core configuration parameters including 'start_day' and 'end_day'
    Zones : list of str
        List of geographic regions for plume detection
    time_step : str or list of str
        Temporal resolution(s) of the data (e.g., 'DAILY', 'MONTHLY')
    nb_cores : int
        Number of CPU cores to use for parallel processing
    dynamic_thresh : bool
        Whether to use dynamic threshold determination for plume detection
    regional_map_dir : str
        Base directory containing the regional satellite maps
    plume_dir : str
        Directory where plume detection results will be saved

    Returns
    -------
    None
        Results are saved to disk as:
        - CSV files containing plume statistics
        - PNG images showing detected plumes
        - Animated GIFs combining the detection results

    Notes
    -----
    The function performs the following steps:
    1. Loads and preprocesses satellite data
    2. Aligns bathymetry data to the satellite resolution
    3. Creates masks for land and search areas
    4. Detects plumes using parallel processing
    5. Saves results as statistics and visualizations

    The results include:
    - Plume area and statistics
    - SPM concentration thresholds
    - Confidence indices
    - Visualization maps
    """

    core_arguments.update(
        {'Years': unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
         'Zones': Zones,
         'Satellite_variables': ['SPM'],
         'Temporal_resolution': ([time_step]
                                 if isinstance(time_step, str)
                                 else time_step)})

    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)

    france_shapefile = load_shapefile_data()

    for i, info in cases_to_process.iterrows():

        # info = cases_to_process.iloc[i]

        print(
            f'{i} over {cases_to_process.shape[0] - 1} ({info.Zone} / {info.Data_source} / {info.sensor_name} / {info.atmospheric_correction} / {info.Year} / {info.Temporal_resolution})')

        # Retrieve specific parameters based on the selected zone.
        parameters = define_parameters(info.Zone)

        # Build the file pattern to locate the satellite data files.
        file_names_pattern = fill_the_sat_paths(info,
                                                path_to_fill_to_where_to_save_satellite_files(
                                                    regional_map_dir + "/" + info.Zone).replace(
                                                    '[TIME_FREQUENCY]', ''),
                                                local_path=True).replace('/*/*/*',
                                                                         f'MAPS/{info.Temporal_resolution}/{info.Year}/*.pkl')

        # Check if the directory for the year's data exists; if not, skip processing.
        if not os.path.exists(os.path.dirname(file_names_pattern)):
            print(f"Missing satellite data here : {file_names_pattern}")
            print("Skip to the next iterate")
            continue

        # Find all files matching the specified pattern.
        file_names = glob.glob(file_names_pattern)

        # Open the first file to extract the dataset structure.
        with open(file_names[0], 'rb') as f:
            ds = pickle.load(f)['Basin_map']['map_data']

        # Reduce the spatial resolution of the dataset to match the new resolution
        ds_reduced = (reduce_resolution(ds, parameters['lat_new_resolution'], parameters['lon_new_resolution'])
                      if parameters['lat_new_resolution'] is not None
                      else ds)

        # Align bathymetry data to the dataset resolution.
        bathymetry_data_aligned_to_reduced_map = align_bathymetry_to_resolution(ds_reduced,
                                                                                f'{regional_map_dir}/{info.Zone}/Bathy_data.pkl')

        # Create a mask that delineate the searching area.
        inside_polygon_mask = create_polygon_mask(ds_reduced, parameters)

        # Create an annual dataset used for validating cloud-free areas, and generate a land mask
        (map_wo_clouds, land_mask) = preprocess_annual_dataset_and_compute_land_mask(
            (re.sub(r'/[0-9]{2,}', '', file_names_pattern)
             .replace(info.Temporal_resolution, 'MULTIYEAR')
             .replace('*.pkl', 'Averaged_over_multi-years.pkl')),
            parameters)

        # Process each file in parallel
        with multiprocess.Pool(nb_cores) as pool:

            results = pool.starmap(main_process,
                                   [(file_name,
                                     file_names_pattern,
                                     parameters,
                                     bathymetry_data_aligned_to_reduced_map,
                                     france_shapefile,
                                     map_wo_clouds,
                                     land_mask,
                                     inside_polygon_mask,
                                     plume_dir,
                                     regional_map_dir,
                                     dynamic_thresh) for file_name in file_names])
        # results = []
        # for file_name in file_names:
        #     result = main_process(file_name,
        #             file_names_pattern, 
        #             parameters,
        #             bathymetry_data_aligned_to_reduced_map,
        #             france_shapefile, 
        #             map_wo_clouds,
        #             land_mask,
        #             inside_polygon_mask,
        #             plume_dir,
        #             regional_map_dir,
        #             dynamic_thresh)
        #     results.append(result)

        # Create a DataFrame from the results, sort by date, and save it to a CSV
        statistics = pd.DataFrame([x for x in results if x is not None]).sort_values('date').reset_index(drop=True)
        folder_name = os.path.dirname(file_names[0]).replace(regional_map_dir,
                                                             plume_dir).replace("MAPS", "PLUME_DETECTION")
        os.makedirs(folder_name, exist_ok=True)
        statistics.to_csv(f'{folder_name}/Results.csv', index=False)

        # Create a GIF from the saved maps by combining all PNG images
        saved_maps = sorted(glob.glob(f'{folder_name}/MAPS/*.png'))
        with imageio.get_writer(f'{folder_name}/GIF.gif', mode='I', fps=1) as writer:
            for figure_file in saved_maps:
                image = imageio.imread(figure_file)
                writer.append_data(image)

        del results, ds_reduced, bathymetry_data_aligned_to_reduced_map, inside_polygon_mask, map_wo_clouds, land_mask
        gc.collect()

        # For debugging
        for file_name in file_names[0:10] :
            print(file_name)
            main_process(file_name,
                  file_names_pattern,
                  parameters,
                  bathymetry_data_aligned_to_reduced_map,
                  france_shapefile,
                  map_wo_clouds,
                  land_mask,
                  inside_polygon_mask,
                  plume_dir,
                  regional_map_dir, 
                  dynamic_thresh)

    # global_cases_to_process = cases_to_process.drop(['Year'], axis = 1).drop_duplicates().reset_index(drop = True)

    # for i in range(global_cases_to_process.shape[0]) :

    #     info = global_cases_to_process.iloc[i].copy()

    #     plot_time_series_of_plume_areas(working_directory = working_directory,
    #                                     Zone = info.Zone,
    #                                     Data_source = info.Data_source,
    #                                     Satellite_sensor = info.Satellite_sensor,
    #                                     atmospheric_correction = info.atmospheric_correction,
    #                                     Temporal_resolution = info.Temporal_resolution,
    #                                     Years = np.arange(1998,2025))


def make_and_plot_time_series_of_plume_areas(core_arguments, Zones, 
                                             nb_cores, time_step,
                                             plume_dir_in, plume_dir_out):
    """
    Calls the R function `plot_time_series_of_plume_area_and_thresholds` from Python.

    Args:
        working_directory (str): Working directory path.
        Zone (str): Zone name.
        Data_source (str): Data source name.
        Satellite_sensor (list): List of satellite sensors.
        atmospheric_correction (str): Atmospheric correction type.
        Temporal_resolution (str): Time resolution.
        Years (list): List of years.
    """

    core_arguments.update(
        {'Years': unique_years_between_two_dates(core_arguments['start_day'], core_arguments['end_day']),
         'Zones': Zones,
         'Satellite_variables': ['SPM'],
         'Temporal_resolution': ([time_step]
                                 if isinstance(time_step, str)
                                 else time_step)})

    Plumes_per_zone = {Zone: list(define_parameters(Zone)['core_of_the_plumes'].keys()) for Zone in Zones}

    # Source the R script
    plume_R_path = os.path.join(func_dir, 'plume.R')
    robjects.r['source'](plume_R_path)
    # robjects.r['source']("myRIOMAR_dev/_3_plume_detection/utils.R")

    r_function = robjects.r['plot_time_series_of_plume_area_and_thresholds']

    # Call the R function
    r_function(
        where_are_saved_plume_results=robjects.StrVector([plume_dir_in]),
        where_to_save_plume_time_series=robjects.StrVector([plume_dir_out]),
        Zone=robjects.StrVector(core_arguments['Zones']),
        Data_source=robjects.StrVector(core_arguments['Data_sources']),
        Satellite_sensor=robjects.StrVector(core_arguments['Sensor_names']),
        atmospheric_correction=robjects.StrVector(core_arguments['Atmospheric_corrections']),
        Temporal_resolution=robjects.StrVector(core_arguments['Temporal_resolution']),
        Years=robjects.IntVector(core_arguments['Years']),
        Plumes=robjects.ListVector(Plumes_per_zone),
        nb_cores=robjects.IntVector([nb_cores])
    )

