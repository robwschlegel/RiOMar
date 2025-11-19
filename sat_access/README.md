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
-- xarray
-- matplotlib
-- cartopy

Install the required packages from a terminal with:

```
pip install xarray matplotlib cartopy
```

## R

- R 4.x
- Packages
-- argparse
-- ncdf4
-- curl
-- ggplot2
-- reshape2

Install the required packages from an R terminal with:

```
install.packages(c("argparse", "ncdf4", "curl", "ggplot2", "reshape2"))
```

### Windows

Make the script executable by write clicking on it and checking the box allowing it to be run as a program.

### Linux

Make the script executable (run within the same location):

```
chmod +x sat_access.R
```

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

The following examples assume the user is in a terminal at the same location as the script.

### Python

```
python sat_access.py --variable chla --date 2025-10-15 --boundingbox 4 6 42 44 --outputdir downloads --overwrite TRUE
```

### R

```
./sat_access.R --variable chla --date 2025-09-15 --boundingbox 4 6 42 44 --outputdir downloads --overwrite TRUE
```

__NB:__ The `./` before `sat_access.R` is necessary for bash to understand that this R script is meant to be run as an executable.

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

