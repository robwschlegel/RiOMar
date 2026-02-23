# modis_landsat_sentinel_access.R

# This script demonstrates a workflow that can be adapted to access data from:
# MODIS, Landsat, and Sentinel platforms

# Meaning that it provides access to all products from the catalog's for:
# NASA, USGS, Copernicus

# After a lot of research it was decided that the package 'rsat' provides
# the best all-round access to marine satellite products in R.
# See here for their explanation why: https://docs.ropensci.org/rsat/

# NB: On Linux it is necessary to ensure the following libraries are installed.
## Copy and paste the following two lines into a terminal (without the #) and run them:
# sudo apt update
# sudo apt install r-cran-rcpp gdal-bin libgdal-dev libproj-dev openssl libssl-dev xml2 libxml2-dev libmagick++-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev libfribidi-dev


# Libraries ---------------------------------------------------------------

# Start with the development version of the package:
# remotes::install_github("spatialstatisticsupna/rsat", build_vignettes=TRUE)

library(tidyverse)
library(rsat)
library(sf)


# Custom functions --------------------------------------------------------

# NB: There is a bug in the rsat_download funciton...
rsat_download2 <- function(x, db_path, verbose = FALSE, parallel = FALSE, ...) {
  require(rsat)
  
  args <- list(...)
  
  if(missing(db_path)){
    db_path <- get_database(x)
    if(db_path == ""){
      stop("db_path or global environment database needed for image downloading.")
    }
  }
  
  # filter records
    x <- records(x)
  dataspace <- x[get_api_name(x) %in% "dataspace"]
  usgs <- x[get_api_name(x) %in% "usgs"]
  lpdaac <- x[get_api_name(x) %in% "lpdaac"]
  
  if(parallel){
    functions_list <- list(
      list(func = connection$getApi("lpdaac")$download_lpdaac_records,
           args = list(lpdaac_records = lpdaac, db_path = db_path,verbose = verbose,...)),
      list(func = rsat:::connection$getApi("dataspace")$dataspace_download_records,
           args = list(records = dataspace, db_path = db_path, verbose = verbose,...)),
      list(func = connection$getApi("usgs")$espa_order_and_download,
           args = list(usgs = usgs,db_path = db_path, verbose = verbose,...))
    )
    null.list <- mclapply(functions_list, function(entry) {
      do.call(entry$func, entry$args)
    }, mc.cores = 3)
  } else {
    functions_list <- list(
      list(func = rsat:::connection$getApi("usgs")$order_usgs_records,
           args = list(espa_orders = usgs,db_path = db_path,verbose = verbose,...)),
      list(func = rsat:::connection$getApi("lpdaac")$download_lpdaac_records,
           args = list(lpdaac_records = lpdaac,db_path = db_path,verbose = verbose,...)),
      list(func = rsat:::connection$getApi("dataspace")$dataspace_download_records,
           args = list(records = dataspace,db_path = db_path,verbose = verbose,...)),
      list(func = rsat:::connection$getApi("usgs")$download_espa_orders,
           args = list(espa.orders = usgs, db_path = db_path, verbose = verbose,...))
    )
    null.list <- lapply(functions_list, function(entry) {
      do.call(entry$func, entry$args)
    })
  }
}


# Credentials -------------------------------------------------------------

# NB: The following credentials are stored in .csv files that are not shared on GitHub.
# To create these files go to the websites listed below and create an account. 
# Then, create a .csv file with two columns and one row. 
# The headers of the columns should be:
# usrname , psswrd
# with the row underneath containing the username and password for each API.

# NASA Earth Data
# https://www.earthdata.nasa.gov/
earth_data <- read_csv("~/pCloudDrive/Documents/info/earthdata_pswd.csv")
# Copernicus Open Access Hub
# https://browser.stac.dataspace.copernicus.eu/
copernicus <- read_csv("~/pCloudDrive/Documents/info/copernicus_pswd.csv")
# USGS Earth Explorer
# https://earthexplorer.usgs.gov/
usgs_ee <- read_csv("~/pCloudDrive/Documents/info/usgs_ee_pswd.csv")

# Set credentials
### NB: The documentatioon does not match the source code / expeted behaviour
## NASA
set_credentials(earth_data$usrname, earth_data$psswrd, credential = "earthdata")
## Copernicus
set_credentials(copernicus$usrname, copernicus$psswrd, credential = "dataspace")
## USGS
## NB: It appears that the rsat package is attempting to connect to USGS via the Earth Data API
# So ideally one can set the same username and password for both NASA and USGS
# set_credentials(usgs_ee$usrname, usgs_ee$psswrd, credential = "earthdata")

# Check that all credentials have been entered correctly:
print_credentials()


# Product IDs -------------------------------------------------------------

# See here for a deeper explanation:
# https://docs.ropensci.org/rsat/articles/rsat1_search.html

## For the Landsat program:
# Landsat 4-5 mission: "LANDSAT_TM_C1"
# Landsat-7 mission: "LANDSAT_ETM_C1"
# Landsat-8 mission: "LANDSAT_8_C1"

## For the MODIS program:
# Terra satellite: "MOD09GA"
# Aqua satellite: "MYD09GA"

## For the Sentinel program:
# Sentinel-2 mission: "S2MSI2A"
# Sentinel-3 mission: "SY_2_SYN___"

# All listed products
rsat::rsat_products()


#  Searching --------------------------------------------------------------

# First establish a bounding box of interest
study_bbox <- st_sf(st_as_sfc(st_bbox(c(
  xmin = 6.9,
  xmax = 7.7,
  ymin = 43.4,
  ymax = 43.8 
), crs = 4326)))

# Then the desired date range
study_time <- seq(as.Date("2020-09-26"),as.Date("2020-10-10"))

# Set the desired download and processing folder paths
# NB: Change according to your desired folder setup
db.base <- "~/data/rsat"
db.path <- file.path(db.base,"database")
ds.path <- file.path(db.base,"datasets")

# Create the folders if they don't yet exist
if(!file.exists(db.base)) dir.create(db.base)
if(!file.exists(db.path)) dir.create(db.path)
if(!file.exists(ds.path)) dir.create(ds.path)

# Create the rtoi object used by rsat to store necessary info
# NB: This is a sort of database and does not behave as one might expect
# One should only create it once, and not overwrite it
# Generally the package will prevent this from happening
# study_rtoi <- new_rtoi(name = "study_site",
#                        region = study_bbox,
#                        db_path = db.path,
#                        rtoi_path = ds.path)

# Load it if it has already been created
study_rtoi <- read_rtoi(file.path(ds.path, "study_site")) 

# View
print(study_rtoi)

# Search for some MODIS Aqua, Sentinel-2, and Sentinel-3 data
rsat_search(region = study_rtoi,
            product = c("MYD09GA", "S2MSI2A", "SY_2_SYN___"),
            dates = study_time)

# Look at the files it found
print(study_rtoi)

# Look at the records more closely
## NB: By default this only shows the first six, but there are many
study_records <- records(study_rtoi)
study_records


# Previewing --------------------------------------------------------------

# The plot function adapts automagically to the rtoi

# Plot the dates
plot(study_rtoi, "dates")

# Here we are previewing the Sentinel-3 data
# NB: Doesn't seem to work as expected...
plot(study_rtoi,
     "preview",
     product = "SY_2_SYN___",
     dates = as.Date("2020-10-10"))


# Downloading -------------------------------------------------------------

# Authenticate credentials
set_credentials(copernicus$usrname, copernicus$psswrd)
print_credentials()

# All downloads are then managed directly with one function
rsat_download2(study_rtoi)

