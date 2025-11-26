#!/usr/bin/env Rscript


# Libraries ---------------------------------------------------------------

# Check for missing libraries and install them if necessary
if (!all(c("argparse", "ncdf4", "curl", "reshape2", "ggplot2") %in% installed.packages())) {
  install.packages(c("argparse", "ncdf4", "curl", "reshape2", "ggplot2"), repos = "https://cloud.r-project.org/")
}

# Uncomment this line if R cannot find python on Windows:
# options(python_cmd = "C:/Path/To/Your/Python/python.exe")

# Activate libraries
library(argparse) # For parsing arguments from the command line
library(ncdf4)    # For reading NetCDF files
suppressPackageStartupMessages(library(curl)) # For FTP download
library(reshape2) # For data reshaping
library(ggplot2)  # For visualization


# Parse arguments ---------------------------------------------------------

# NB: If you have opened this script in R/RStudio and intend to run it manually, this section can be ignored.

# Create the parser
parser <- ArgumentParser(description = "Download a NetCDF file from FTP and plot a variable as a map.")

# Add arguments
parser$add_argument("-v", "--variable", type = "character", required = TRUE, help = "Surface variable to fetch and plot")
parser$add_argument("-d", "--daterange", type = "character", nargs = '+', required = TRUE, 
                    help = "Date range of the desired variable in YYYY-MM-DD. Provide only one or two values. 
                    Note that if two dates are provided, only the first will be used for plotting.")
parser$add_argument("-od", "--outputdir", type = "character", required = TRUE, help = "Location to save the NetCDF file and output image")
parser$add_argument("-ov", "--overwrite", type = "logical", default = FALSE, 
                    help = "Whether to overwrite an existing file or not. Default = FALSE")
parser$add_argument("-p", "--plot", type = "logical", default = FALSE, 
                    help = "Whether or not to plot the downloaded data. Default = False")
parser$add_argument("-bbox", "--boundingbox", nargs = 4, type = "double", required = FALSE,
                    help = "The bounding box for plotting. Must be given as: lonmin, lonmax, latmin, latmax")

# Create the function
args <- parser$parse_args()


# The function to call ----------------------------------------------------

# Main function to download and plot
download_and_plot <- function(dl_var, dl_dates, bbox, output_dir, plot_var, overwrite) {
  
  ## Download code -----------------------------------------------------------

  # Check date range
  if(length(dl_dates) > 2){
    stop("Please provide only one or two dates for the date range.")
  }
  
  message("Downloading files...")
  
  # Set start and end dates for download
  start_date <- as.Date(dl_dates[1])
  if(length(dl_dates) == 2){
    end_date <- as.Date(dl_dates[2])
  } else{
    end_date <- start_date
  }
  
  # Iterate over each date in the range
  current_date <- start_date
  while(current_date <= end_date){
    
    # Prep date strings
    dl_date <- current_date
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
    if(file.exists(file_name_full) & !overwrite){
      
      message(paste0(file_name_full," already exists. Set --overwrite TRUE to force the download."))
      
    } else {
      
      # Ensure output directory exists
      if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
      }
      
      # Download
      curl::curl_download(url_final, destfile = file_name_full)
      message(paste0("File downloaded at: ",file_name_full))
      
      # Unzip
      system(paste("bunzip2 -k -f", file_name_full))
      message("File unzipped at: ", gsub(".bz2","",file_name_full))
    }
    
    # Move to the next day
    current_date <- current_date + 1
  }


  ## Plotting code -----------------------------------------------------------

  if(!plot_var){
    message("Plotting skipped as per user request.")
  } else {
    message("Plotting...")
    
    # User warning about date range
    if(length(dl_dates) == 2){
      message("Two dates provided; the last date will be used for plotting.")
    }
    
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
    
    # Reshape data for ggplot
    var_df <- melt(var)
    names(var_df) <- c("lon_idx", "lat_idx", "value")
    var_df$lon <- lon[var_df$lon_idx]
    var_df$lat <- lat[var_df$lat_idx]
    
    if(is.null(bbox)){
      bbox <- c(
        min(lon, na.rm = TRUE),
        max(lon, na.rm = TRUE),
        min(lat, na.rm = TRUE),
        max(lat, na.rm = TRUE)
      )
      message(paste("No bounding box provided, using full extent:", paste(round(bbox, 4), collapse = ", ")))
    }
    
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
}


# Parse the args ----------------------------------------------------------

# NB: If you have opened this script in R/RStudio and intend to run it manually, 
# uncomment the following chunk of code and add your desired values directly:
# args <- list(
#   variable = "SPM",
#   daterange = c("2023-01-01", "2023-01-05"),
#   outputdir = "path/to/output/dir",
#   overwrite = TRUE,
#   plot = TRUE,
#   boundingbox = c(-10, 10, 35, 45)
# )

# Call the main function with parsed arguments
download_and_plot(
  dl_var = args$variable,
  dl_dates = args$daterange,
  output_dir = args$outputdir,
  overwrite = args$overwrite,
  plot_var = args$plot,
  bbox = args$boundingbox
)

