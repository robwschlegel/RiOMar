# Overview

This Python script allows you to download a series NetCDF files for chl a and/or SPM from an FTP server and visualize a specified variable as a geographical map. 
It is designed specifically to work for the suite of ODATIS products put into production for/during the RiOMar project.

# features

- Uses argparse for easy command-line argument handling
- Downloads a range of NetCDF file from an FTP server
- Extracts and visualizes a specified date and variable as a map
- Available for both `Python` and `R`

# Requirements

## Python

- Python: 3.x
- Packages:
  - xarray
  - matplotlib
  - cartopy

To check if you have Python installed, open a terminal and type:

```
python --version
```

Or if you would like to use an existing virtual environment (recommended):

```
conda activate your_env_name
```

If your machine cannot find an active version of python, but you should have one, see the troubleshooting below.

Otherwise, once you have an active Python environment available in your terminal, install the required packages from the terminal with:

```
pip install xarray matplotlib cartopy
```

Note that if you have not used `cartopy` before, the first time you run the `sat_access.py` script it will download a few shape files. 
This may take a few minutes on a slow internet connection. During which time the map menu will appear to be hanging.
You will know it is working if you see activity in the console as it downloads the necessary shapefiles.

## R

- R: 4.x
- Packages:
  - argparse
  - ncdf4
  - curl
  - ggplot2
  - reshape2

To check if you have R installed, open a terminal and type:

```
R --version
```

Install the required packages from an the terminal with:

```
Rscript -e "install.packages(c('argparse', 'ncdf4', 'curl', 'ggplot2', 'reshape2'), repos='https://cloud.r-project.org/')"
```

### Windows

Make the script executable by write clicking on it and checking the box allowing it to be run as a program.

### Linux

Make the script executable (run within the same location as the script):

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
| `--plot` | Whether or not to plot the downloaded data; Default = False | No | True |
| `--boundingbox` | Bounding box as: lon_min lon_max lat_min lat_max | No | 4 6 42 44 |

## Example Usage

The following examples assume the user is in a terminal at the same location as the script.

### Python

Download a single chl a file

```
python sat_access.py --variable chla --date 2025-10-15 --outputdir downloads
```

Download multiple SPM files and plot one
 
```
python sat_access.py --variable SPM --date 2025-10-15 2025-10-17 --outputdir downloads --plot = True --boundingbox 4 6 42 44
```

### R

Download a single chl a file and plot it:

```
./sat_access.R --variable chla --date 2025-09-15 --outputdir downloads --plot = TRUE --boundingbox 4 6 42 44
```

Download multiple SPM files but plot none:

```
./sat_access.R --variable SPM --date 2025-11-15 2025-11-17 --outputdir downloads --plot = FALSE
```

__NB:__ The `./` before `sat_access.R` is necessary for bash to understand that this R script is meant to be run as an executable.

# How It Works

1. FTP Download:
  - The script connects to the specified FTP server based on the requested variable and date range and downloads the NetCDF file to your local machine.
2. Visualisation:
  - It loads the NetCDF file and clips the data to the specified bounding box the longitude, latitude, and the specified variable before plotting it as a map 
  using __`matplotlib`__ (`Python`) or __`ggplot2`__ (`R`).
3. Saving plots:
  - The `R` script creates and saves the plot directly to the specified `outputdir` folder, while the `Python` script displays the plot in an interactive window.

# Troubleshooting

- FTP Connection Issues: Check your internet connection and firewall settings. Or wait a few seconds and try again.
- Missing Variables: If not all required variables have been imputed in the terminal, an automated error message should be generated to help you along.
- Module Errors: Ensure all required `Python` or `R` packages are installed.
- `Python` or `R` Not Found: See additional steps below.

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

