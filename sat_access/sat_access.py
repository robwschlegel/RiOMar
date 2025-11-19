#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ftplib
import netCDF4 as nc
import matplotlib.pyplot as plt
import argparse

def download_and_plot(ftp_server, ftp_path, username, password, local_filename, variable_name):
    
    # Download the NetCDF file from FTP
    ftp = ftplib.FTP(ftp_server)
    ftp.login(user=username, passwd=password)
    with open(local_filename, 'wb') as f:
        ftp.retrbinary('RETR ' + ftp_path, f.write)
    ftp.quit()

    # Load the NetCDF file
    data = nc.Dataset(local_filename)

    # Extract variables
    lons = data.variables['longitude'][:]
    lats = data.variables['latitude'][:]
    var = data.variables[variable_name][:]

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.contourf(lons, lats, var, levels=20, cmap='viridis')
    plt.colorbar(label=variable_name)
    plt.title(f'Map of {variable_name}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download a NetCDF file from FTP and plot a variable as a map.')
    parser.add_argument('--server', required=True, help='FTP server address')
    parser.add_argument('--path', required=True, help='Path to the NetCDF file on the FTP server')
    parser.add_argument('--username', required=True, help='FTP username')
    parser.add_argument('--password', required=True, help='FTP password')
    parser.add_argument('--output', default='downloaded_file.nc', help='Local filename to save the NetCDF file')
    parser.add_argument('--variable', required=True, help='Variable name in the NetCDF file to plot')

    args = parser.parse_args()

    download_and_plot(args.server, args.path, args.username, args.password, args.output, args.variable)

