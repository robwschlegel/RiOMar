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

To check if you have Python installed, open a terminal and type:

```
python --version
```

If your machine cannot find an active version of python, but you should have one, see the troubleshooting below.

Otherwise, once you have an active Python environment available in your terminal, install the required packages from a terminal with:

```
pip install xarray matplotlib cartopy
```

Note that if you have not used `cartopy` before, the first time you run the `sat_access.py` script it will download a few shape files. 
This may take a few minutes on a slow internet connection. During which time the map menu will appear to be hanging, but wait it out for a while.
You will see activity in the console as it downloads the necessary shapefiles.

## R

- R 4.x
- Packages
-- argparse
-- ncdf4
-- curl
-- ggplot2
-- reshape2

To check if you have R installed, open a terminal and type:

```
R --version
```

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
| `--variable` | The data layer to download | Yes | SPM |
| `--daterange` | Date range for desired data in YYYY-MM-DD format | Yes | 2025-10-15 |
| `--outputdir` | Local folder to save the NetCDF file | Yes | downloads |
| `--overwrite` | Overwrite existing files; Default = False | No | True |
| `--plot` | Whether or not to plot the downloaded data. Default = False | No | True |
| `--boundingbox` | Bounding box as lon_min lon_max lat_min lat_max | No | 4 6 42 44 |

## Example Usage

The following examples assume the user is in a terminal at the same location as the script.

### Python

Download a single chl a file

```
python sat_access.py --variable chla --date 2025-10-15 --outputdir .
```

 Download multiple SPM files
 
 ```
python sat_access.py --variable SPM --date 2025-10-15 2025-10-17 --outputdir .
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

## Windows

If this returns nothing:

```
python --version
```

But you are certain python is installed (e.g. via anaconda), you may need to add python to your PATH.

In windows:

1. Open System Environment Variables:

  - Press Win + S, type "Environment Variables", and select "Edit the system environment variables".
  - Click "Environment Variables".

2. Edit the PATH Variable:

  - Under "User variables" or "System variables", find the Path variable and click "Edit".
  - Add the following paths (adjust for your Anaconda installation):
  - `C:\Users\YourUsername\anaconda3`
  - `C:\Users\YourUsername\anaconda3\Scripts`
  - `C:\Users\YourUsername\anaconda3\Library\bin`
  - Click OK to save.

3. Restart PowerShell:

  - Close and reopen PowerShell for the changes to take effect.

Alternatively, if you already have a virtual environment installed on your system it would be preferable to activate it before running the script.

To check the virtual environments installed on your system (with anaconda/miniconda):

```
conda env list

```

To activate the environment of choice:


```
conda activate your_env_name

```

# License

This script is provided as-is and is free to use. No warranty is provided.

# Support

For questions or issues, please contact Robert at: robert.schlegel@imev-mer.fr

