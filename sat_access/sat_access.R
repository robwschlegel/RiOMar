#!/usr/bin/env Rscript


# TODO:
# Check for download directory and either create one or return error that it doesn't exist


# Libraries ---------------------------------------------------------------

library(argparse) # For parsing arguments from the command line
library(ncdf4)    # For reading NetCDF files
suppressPackageStartupMessages(library(curl))     # For FTP download
library(ggplot2)  # For visualization
library(reshape2) # For data reshaping

# message("All libraries loaded.")


# Parse arguments ---------------------------------------------------------

# Create the parser
parser <- ArgumentParser(description = "Download a NetCDF file from FTP and plot a variable as a map.")

# Add arguments
parser$add_argument("-v", "--variable", type = "character", required = TRUE, help = "Surface variable to fetch and plot")
parser$add_argument("-d", "--date", type = "character", required = TRUE, help = "Date of the desired variable in YYYY-MM-DD")
parser$add_argument("-bbox", "--boundingbox", nargs = 4, type = "double", 
                    help = "The bounding box for plotting. Must be given as: lonmin, lonmax, latmin, latmax")
parser$add_argument("-od", "--outputdir", type = "character", required = TRUE, help = "Location to save the NetCDF file and output image")
parser$add_argument("-ov", "--overwrite", type = "logical", default = FALSE, 
                    help = "Whether to overwrite an existing file or not. Default = FALSE")
# parser$add_argument("--server", type = character, required = TRUE, help = "FTP server address")
# parser$add_argument("--path", type = character, required = TRUE, help = "Path to the NetCDF file on the FTP server")
# parser$add_argument("--username", type = character, required = TRUE, help = "FTP username")
# parser$add_argument("--password", type = character, required = TRUE, help = "FTP password")

# Create the function
args <- parser$parse_args()

# message("All arguments parsed.")

# testers...
# dl_var = "SPM"
# dl_var = "CHLA"
# dl_var = "SST"
# dl_date = "2025-11-16"
# bbox = c(4, 6, 42, 44)
# output_dir = "sat_access/downloads"
# overwrite = FALSE


# The function to call ----------------------------------------------------

# Main function to download and plot
download_and_plot <- function(dl_var, dl_date, bbox, output_dir, overwrite) {

  
  ## Download code -----------------------------------------------------------

  message("Downloading file...")
  
  # Prep date strings
  dl_date_flat <- gsub("-", "", dl_date)
  url_year_doy <- paste0(substr(dl_date, start = 1, stop = 4),"/",strftime(as.Date(dl_date), format = "%j"))
  
  # Get general URL based on desired variables
  if(toupper(dl_var) %in% c("SPM", "SPIM", "CHLA")){
    url_base <- "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
  } else {
    stop("Variable not yet available")
  }
  
  # Get product specifics
  if(toupper(dl_var) %in% c("SPM", "SPIM")){
    file_name <- paste0(dl_date_flat,"-EUR-L4-SPIM-ATL-v01-fv01-OI.nc.bz2")
    url_product <- "EUR-L4-SPIM-ATL-v01"
    nc_file <- file.path(output_dir, paste0(dl_date_flat,"-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"))
    nc_var_name <- "analysed_spim"
    var_label <- "SPM [g m-3]"
  } else if(toupper(dl_var) == "CHLA"){
    file_name <- paste0(dl_date_flat,"-EUR-L4-CHL-ATL-v01-fv01-OI.nc.bz2")
    url_product <- "EUR-L4-CHL-ATL-v01"
    nc_file <- file.path(output_dir, paste0(dl_date_flat,"-EUR-L4-CHL-ATL-v01-fv01-OI.nc"))
    nc_var_name <- "analysed_chl_a"
    var_label <- "chl a [mg m-3]"
  } else {
    stop("Variable not yet available")
  }
  
  # Assemble final URL
  url_final <- paste(url_base, url_product, url_year_doy, file_name, sep = "/")
  file_name_full <- file.path(output_dir, file_name)
  
  # Fetch file
  # GET(url_product, write_disk(output_file_ext, overwrite = TRUE), progress(), content = "raw")
  if(file.exists(file_name_full) & !overwrite){
    message(paste0(file_name_full," already exists. Set --overwrite TRUE to force the download."))
    # return()
  } else {
    curl::curl_download(url_final, destfile = file_name_full)
    
    # On Linux + Mac
    system(paste("bunzip2 -k", file_name_full))
    
    # On windows
    # Example: Unzip a .bz2 file using 7z
    # bz2_file <- "path/to/your/file.bz2"
    # output_file <- "path/to/your/output_file"  # Remove the .bz2 extension
    
    # Run the 7z command- On Windows, you can use 7z (from 7-Zip) if it is installed:
    # system(paste0("7z e ", bz2_file, " -o", file.path(dirname(bz2_file), "output_dir")))
    
    message(paste0("File downloaded at: ",file_name_full))
  }
  

  ## Plotting code -----------------------------------------------------------

  message("Plotting...")
  
  # Set plot name
  plot_name <- file.path(output_dir, paste0(nc_var_name,"_",dl_date,".png"))
  
  # Open the NetCDF file
  nc_data <- nc_open(nc_file)
  
  # Extract longitude, latitude, and the specified variable
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  time <- ncvar_get(nc_data, "time")
  var <- ncvar_get(nc_data, nc_var_name)
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Reshape data for ggplot (if necessary)
  var_df <- melt(var)
  names(var_df) <- c("lon_idx", "lat_idx", "value")
  var_df$lon <- lon[var_df$lon_idx]
  var_df$lat <- lat[var_df$lat_idx]
  
  # Filter df to bounding box
  var_df_sub <- var_df[var_df$lon >= bbox[1] & var_df$lon <= bbox[2],]
  var_df_sub <- var_df_sub[var_df_sub$lat > bbox[3] & var_df_sub$lat <= bbox[4], ]
  
  # Get date for plot label
  plot_date <- as.Date(as.POSIXct(time, origin = "1998-01-01"))
  
  # Plot using ggplot2
  p <- ggplot(var_df_sub, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(var_label) +
    labs(title = paste("Map of", nc_var_name, "on", plot_date),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(filename = plot_name, plot = p)
  message(paste0("Image saved at: ",plot_name))
}


# Parse the args ----------------------------------------------------------

# Call the main function with parsed arguments
download_and_plot(
  dl_var = args$variable,
  dl_date = args$date,
  bbox = args$boundingbox,
  output_dir = args$outputdir,
  overwrite = args$overwrite
)

