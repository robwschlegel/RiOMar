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
import numpy as np
import pandas as pd
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

    # Open the NetCDF file using xarray
    ds = xr.open_dataset(nc_file)

    # Get date value as a label for plotting
    date_value = ds.time.values[0]  # Extract the numpy.datetime64 value
    date_time_obj = pd.to_datetime(date_value)
    date_label = date_time_obj.strftime('%Y-%m-%d')

    # Subset NetCDF file to desired bounding box
    var_subset = ds[nc_var_name].where(
        (ds.lon >= bbox[0]) & (ds.lon <= bbox[1]) &
        (ds.lat >= bbox[2]) & (ds.lat <= bbox[3]),
        drop=True
    )

    # Create a figure and axes with a Plate Carree projection
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
    plt.title(f'Map of {nc_var_name} on {date_label}')

    # Show the plot
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download a NetCDF file and plot a variable as a map.")
    parser.add_argument("--variable", type=str, required=True, help="Variable to download (e.g., 'SPM', 'SPIM', 'CHLA')")
    parser.add_argument("--date", type=str, required=True, help="Date in YYYY-MM-DD format")
    parser.add_argument("--boundingbox", type=float, nargs=4, required=True, help="Bounding box as lon_min lon_max lat_min lat_max")
    parser.add_argument("--outputdir", type=str, required=True, help="Directory to save the downloaded file and plot")
    parser.add_argument("--overwrite", type=bool, default=False, help="Overwrite existing files. Default = False")

    args = parser.parse_args()
    download_and_plot(args.variable, args.date, args.boundingbox, args.outputdir, args.overwrite)

