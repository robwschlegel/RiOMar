# Overview

This Python script allows you to download a NetCDF file from an FTP server and visualize a specified variable as a geographical map. 
It is designed specifically to work for the suite of ODATIS products put into production for/during the RiOMar project.

# features

- Downloads a NetCDF file from an FTP server using provided credentials.
- Extracts and visualizes a specified variable as a contour map.
- Supports custom output filenames for the downloaded file.
- Uses argparse for easy command-line argument handling.

# Requirements

## Python

- Python 3.x
- Packages:
-- netCDF4
-- matplotlib

Install the required packages using:

```
pip install netCDF4 matplotlib
```

## R

- R 4.x
- Packages
-- argparse
-- ncdf4
-- httr
-- dplyr

# Usage

## Command-Line Arguments

| Argument | Description | Required | Example Value |
|----|----|----|----|
| `--server` | FTP server address | Yes | ftp.example.com |
| `--path` | Path to the NetCDF file on the FTP server | Yes | /data/file.nc |
| `--username` | FTP username | Yes | user123 |
| `--password` | FTP password | Yes | pass123 |
| `--variable` | The data layer to visualise | Yes |temperature |
| `--output` | Local filename to save the NetCDF file | No | downloaded_data.nc |

## Example Usage

### Basic Example

```
python sat_access.py \
    --server ftp.example.com \
    --path /data/ocean_data.nc \
    --username myusername \
    --password mypassword \
    --variable temperature
```

### Custom Output Filename

```
python sat_access.py \
    --server ftp.example.com \
    --path /data/ocean_data.nc \
    --username myusername \
    --password mypassword \
    --output my_ocean_data.nc \
    --variable salinity
```

# How It Works

1. FTP Download:
  The script connects to the specified FTP server using the provided credentials and downloads the NetCDF file to your local machine.
2. Data Extraction:
  It loads the NetCDF file and extracts the longitude, latitude, and the specified variable (e.g., temperature, salinity).
3. Visualization:
  The script generates a contour map of the specified variable using matplotlib, with longitude and latitude as the axes.

# Assumptions

- The NetCDF file contains variables named longitude and latitude. If your file uses different names, you will need to modify the script accordingly.
- The specified variable (e.g., temperature) exists in the NetCDF file.

# Troubleshooting

- FTP Connection Issues: Ensure the FTP server address, path, username, and password are correct. Check your internet connection and firewall settings.
- Missing Variables: If the script fails to find longitude, latitude, or the specified variable, verify the variable names in your NetCDF file using a tool like Panoply.
- Module Errors: Ensure all required Python packages are installed. Use `pip install netCDF4 matplotlib` if you encounter import errors.

# License

This script is provided as-is and is free to use. No warranty is provided.

# Support

For questions or issues, please contact Robert at: robert.schlegel@imev-mer.fr

