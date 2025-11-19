#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import urllib.request
import bz2
# import netCDF4 as nc
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# import numpy as np
from datetime import datetime, timedelta

def download_and_plot(dl_var, dl_date, bbox, output_dir, overwrite=False):
    """
    Downloads a NetCDF file from a specified URL and plots a variable as a map.

    Parameters:
    - dl_var: Variable to download (e.g., "SPM", "SPIM", "CHLA").
    - dl_date: Date in YYYY-MM-DD format.
    - bbox: Bounding box as [lon_min, lon_max, lat_min, lat_max].
    - output_dir: Directory to save the downloaded file and plot.
    - overwrite: Whether to overwrite existing files.
    """

    # Download code
    print("Downloading file...")

    # Prep date strings
    dl_date_flat = dl_date.replace("-", "")
    url_year_doy = f"{dl_date[:4]}/{datetime.strptime(dl_date, '%Y-%m-%d').strftime('%j')}"

    # Get general URL based on desired variables
    if dl_var.upper() in ["SPM", "SPIM", "CHLA"]:
        url_base = "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
    else:
        raise ValueError("Variable not yet available")

    # Get product specifics
    if dl_var.upper() in ["SPM", "SPIM"]:
        file_name = f"{dl_date_flat}-EUR-L4-SPIM-ATL-v01-fv01-OI.nc.bz2"
        url_product = "EUR-L4-SPIM-ATL-v01"
        nc_var_name = "analysed_spim"
        var_label = "SPM [g m-3]"
        nc_file = os.path.join(output_dir, f"{dl_date_flat}-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")
    elif dl_var.upper() == "CHLA":
        file_name = f"{dl_date_flat}-EUR-L4-CHL-ATL-v01-fv01-OI.nc.bz2"
        url_product = "EUR-L4-CHL-ATL-v01"
        nc_var_name = "analysed_chl_a"
        var_label = "chl a [mg m-3]"
        nc_file = os.path.join(output_dir, f"{dl_date_flat}-EUR-L4-CHL-ATL-v01-fv01-OI.nc")
    else:
        raise ValueError("Variable not yet available")

    # Assemble final URL
    url_final = f"{url_base}/{url_product}/{url_year_doy}/{file_name}"
    file_name_full = os.path.join(output_dir, file_name)

    # Fetch file
    if os.path.exists(file_name_full) and not overwrite:
        print(f"{file_name_full} already exists. Set --overwrite True to force the download.")
    else:
        urllib.request.urlretrieve(url_final, file_name_full)
        print(f"File downloaded at: {file_name_full}")

        # Extract the .bz2 file
        with open(file_name_full, 'rb') as source, open(file_name_full[:-4], 'wb') as dest:
            dest.write(bz2.decompress(source.read()))
        print(f"File extracted at: {file_name_full[:-4]}")

    # Plotting code
    print("Plotting...")

    # Open the NetCDF file
    # nc_data = nc.Dataset(nc_file)

    # Extract longitude, latitude, and the specified variable
    # lon = nc_data.variables['lon'][:]
    # lat = nc_data.variables['lat'][:]
    # time = nc_data.variables['time'][:]
    # var = nc_data.variables[nc_var_name][:]

    # Close the NetCDF file
    # nc_data.close()

    # Create a grid of longitude and latitude
    # lon_grid, lat_grid = np.meshgrid(lon, lat)

    # Filter data to bounding box
    # lon_mask = (lon_grid >= bbox[0]) & (lon_grid <= bbox[1])
    # lat_mask = (lat_grid >= bbox[2]) & (lat_grid <= bbox[3])
    # var_filtered = np.ma.masked_where(~(lon_mask & lat_mask), var)

    # Open the NetCDF file using xarray
    ds = xr.open_dataset(nc_file)
    print(ds)

    # Extract the specified variable and coordinates
    # var = ds[nc_var_name]
    # lon = ds.lon
    # lat = ds.lat
    time = ds.time
    date_value = ds.time.values[0]
    date_datetime = pd.to_datetime(date_value).to_pydatetime()
    print(date_value)
    # lon_grid, lat_grid = np.meshgrid(lon, lat)
    
    # Select data within the bounding box
    # var_filtered = var.where(
    #     (lon >= bbox[0]) & (lon <= bbox[1]) &
    #     (lat >= bbox[2]) & (lat <= bbox[3]),
    #     drop=True
    # )
    var_subset = ds[nc_var_name].where(
        (ds.lon >= bbox[0]) & (ds.lon <= bbox[1]) &
        (ds.lat >= bbox[2]) & (ds.lat <= bbox[3]),
        drop=True
    )

    # Create a figure and axis with a Plate Carree projection
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Plot the data
    raster = var_subset.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='viridis',
        add_colorbar=True,
        cbar_kwargs={'label': 'Variable Units'}
    )

    # Add geospatial features
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.gridlines(draw_labels=True)

    # Add a title
    plt.title(f'Map of {nc_var_name} on {time[0]}')

    # Show the plot
    plt.show()

    # %%
    # var_filtered.plot()
    # var_filtered.plot.pcolormesh(x="lon", y="lat");

    # Close the dataset (optional, as xarray handles file closing when the object is garbage collected)
    # ds.close()

    # Get date for plot label
    # plot_date = datetime(1998, 1, 1) + timedelta(seconds=1234567890)

    # Get the filtered lon and lat coordinates from var_filtered
    # filtered_lon = var_filtered.lon
    # filtered_lat = var_filtered.lat

    # %%
    # Create a meshgrid using the filtered coordinates
    # lon_grid, lat_grid = np.meshgrid(filtered_lon, filtered_lat)
    # lon_grid, lat_grid = np.meshgrid(var_filtered.lon, var_filtered.lat)

    # Transpose var_filtered if necessary to match the orientation of lon_grid and lat_grid
    # var_filtered_values = var_filtered.values
    # if var_filtered_values.shape != lon_grid.shape:
    #     var_filtered_values = var_filtered_values.T

    # Plot using matplotlib with pcolormesh
    # plt.figure(figsize=(10, 6))
    # raster = plt.pcolormesh(lon_grid, lat_grid, var_filtered, cmap='viridis', shading='auto')

    # Add colorbar and labels
    # plt.colorbar(raster, label=var_label)
    # plt.title(f"Map of {nc_var_name} on {time}")
    # plt.xlabel("Longitude")
    # plt.ylabel("Latitude")  

    # Save the plot
    # plot_name = os.path.join(output_dir, f"{nc_var_name}_{dl_date}.png")
    # plt.savefig(plot_name, bbox_inches='tight', dpi=300)
    # print(f"Image saved at: {plot_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download a NetCDF file and plot a variable as a map.")
    parser.add_argument("--variable", type=str, required=True, help="Variable to download (e.g., 'SPM', 'SPIM', 'CHLA')")
    parser.add_argument("--date", type=str, required=True, help="Date in YYYY-MM-DD format")
    parser.add_argument("--boundingbox", type=float, nargs=4, required=True, help="Bounding box as lon_min lon_max lat_min lat_max")
    parser.add_argument("--outputdir", type=str, required=True, help="Directory to save the downloaded file and plot")
    parser.add_argument("--overwrite", type=bool, default=False, help="Overwrite existing files. Default = False")

    args = parser.parse_args()
    download_and_plot(args.variable, args.date, args.boundingbox, args.outputdir, args.overwrite)

