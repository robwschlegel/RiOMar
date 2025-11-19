#!/usr/bin/env Rscript

# Load required libraries
library(argparse)
library(ncdf4)    # For reading NetCDF files
library(httr)     # For FTP download
library(ggplot2)  # For visualization
library(reshape2) # For data reshaping (if needed)

# Define the argument parser
parser <- ArgumentParser(description = "Download a NetCDF file from FTP and plot a variable as a map.")

parser$add_argument("--server", type = character, required = TRUE, help = "FTP server address")
parser$add_argument("--path", type = character, required = TRUE, help = "Path to the NetCDF file on the FTP server")
parser$add_argument("--username", type = character, required = TRUE, help = "FTP username")
parser$add_argument("--password", type = character, required = TRUE, help = "FTP password")
parser$add_argument("--output", type = character, default = "downloaded_data.nc", help = "Local filename to save the NetCDF file")
parser$add_argument("--variable", type = character, required = TRUE, help = "Variable name in the NetCDF file to plot")

# Parse arguments
args <- parser$parse_args()

# Function to download a file from FTP
download_ftp_file <- function(server, path, username, password, output_file) {
  url <- paste0("ftp://", username, ":", password, "@", server, path)
  GET(url, write_disk(output_file, overwrite = TRUE), progress())
}

# Function to plot a NetCDF variable as a map
plot_ncdf_map <- function(nc_file, variable_name) {
  # Open the NetCDF file
  nc_data <- nc_open(nc_file)
  
  # Extract longitude, latitude, and the specified variable
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude")
  var <- ncvar_get(nc_data, variable_name)
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Reshape data for ggplot (if necessary)
  var_df <- melt(var)
  names(var_df) <- c("lon_idx", "lat_idx", "value")
  var_df$lon <- lon[var_df$lon_idx]
  var_df$lat <- lat[var_df$lat_idx]
  
  # Plot using ggplot2
  p <- ggplot(var_df, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = paste("Map of", variable_name),
         x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  print(p)
}

# Main function to download and plot
download_and_plot <- function(server, path, username, password, output_file, variable_name) {
  message("Downloading file from FTP...")
  download_ftp_file(server, path, username, password, output_file)
  
  message("Plotting data...")
  plot_ncdf_map(output_file, variable_name)
}

# Call the main function with parsed arguments
download_and_plot(
  server = args$server,
  path = args$path,
  username = args$username,
  password = args$password,
  output_file = args$output,
  variable_name = args$variable
)

