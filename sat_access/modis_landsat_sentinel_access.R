# modis_landsat_sentinel_access.R

# This script demonstrates a workflow that can be adapted to access data from:
# MODIS, Landsat, and Sentinel platforms

# Meaning that it provides access to all products from the catalog's for:
# NASA, USGS, Copernicus


# Libraries ---------------------------------------------------------------

# After a lot of research it was decided that the package 'rsat' provides
# the best all-round access to marine satellite products.
# See here for their explanation why: https://docs.ropensci.org/rsat/

earth_up <- read_csv("~/pCloudDrive/Documents/info/earthdata_pswd.csv")
