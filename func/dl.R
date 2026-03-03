# func/dl.R

# R code for running downloads

# TODO: Export as much of this to 0_download_data.py as possible


# Setup -------------------------------------------------------------------

# Get internal functions, meta-data, etc.
source("func/util.R")

# Repo available here:
# https://github.com/RiOMar-projet/sat_access
source("~/sat_access/sat_access_script.R")

# AVISO+ credentials
aviso_plus_cred <- read.csv("~/pCloudDrive/Documents/info/aviso_plus_pswd.csv")

# Time steps to walk through
vec_time_steps <- c("daily", "weekly")#, "monthly")

# The variables we want in the two different ODATIS-MR correction types
vars_nirswir <- c("CDOM", "CHL", "SPM", "TUR", "RRS", "SST")
vars_polymer <- c("CDOM", "CHL", "SPM", "TUR", "RRS")


# Download all SEXTANT ----------------------------------------------------

# Unnecessary as this is all run via the main python script
# Though I would like to remove the file structure complexity in the project generally...


# Download ODATIS-MR ------------------------------------------------------

## MODIS ------------------------------------------------------------------

# TODO: Could convert this to it's own wrapper function for convenience and to reduce lines of code / complexity
# OR just the whole thing into a massive ridiculous nested for loop...
for(i in nrow(zones_bbox)){
  
  # Get the target zone
  zone <- zones_bbox[i,]
  
  for(j in length(vec_time_steps)){
    for(k in length(vars_nirswir)){
      
      # Download the data
      download_nc(
        dl_var = vars_nirswir[k],
        dl_dates = c("2002-07-04"), # TODO: THis will need to be handled differently for monthly data
        # dl_dates = c("2002-07-04", "2024-12-31"),
        dl_product = "ODATIS-MR",
        dl_sensor = "MODIS",
        dl_correction = "nirswir", 
        dl_time_step = vec_time_steps[j],
        dl_bbox = c(zone$lon_min, zone$lon_max, 
                    zone$lat_min, zone$lat_max),
        username = aviso_plus_cred$usrname,
        password = aviso_plus_cred$psswrd,
        output_dir = file.path("~/data/ODATIS-MR/MODIS", zone$zone, vec_time_steps[j]),
        overwrite = FALSE
      )
      
    }
  }
}

