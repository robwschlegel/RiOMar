# func/analyse_spatiotemporal.R
# 2026-02-23

# This script will download satellite data
# Then load it in bite sized pieces
# Then perform an time series analysis


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# Get satellite download function
source("~/sat_access/sat_access_script.R")

# lon lat ranges
lon_range <- c(7.03, 7.42)
lat_range <- c(43.44, 43.73)


# Functions ---------------------------------------------------------------

# Load and prep one day of data
# testers...
# file_name <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
# lon_range <- c(1, 4)
load_sextant <- function(file_name, lon_range, lat_range){

  # Find the date
  sextant_one_date <- as.Date(tidync(file_name)[["attribute"]][["value"]][["start_date"]])
  
  # The necessary code
  sextant_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = sextant_one_date) |> 
    dplyr::select(lon, lat, date, mask, analysed_spim)
  
  # Exit
  return(sextant_one)
}


# Download data -----------------------------------------------------------

# NB: Uncomment and run the following lines to download data

# A few days of SEXTANT data
# download_nc(
#   dl_var = "SPM",
#   dl_dates = c("2024-09-01", "2024-09-05"),
#   output_dir = "~/Downloads/SEXTANT", # Change as desired/required
#   overwrite = FALSE # Change to TRUE to force downloads
# )

# A few days of MODIS data, cut to a bounding box
# download_nc(
#   dl_var = "SPM",
#   dl_dates = c("2008-12-12", "2008-12-31"),
#   dl_product = "ODATIS-MR",
#   dl_sensor = "MODIS",
#   dl_bbox = c(3, 4, 42.5, 44),
#   output_dir = "~/Downloads/MODIS", # Change as desired/required
#   overwrite = TRUE # Change to TRUE to force downloads
# )


# Load data ---------------------------------------------------------------

# All of the SEXTANT 1998 files
# sextant_1998_dir <- dir("~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# sextant_1999_dir <- dir("~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1999", pattern = ".nc", recursive = TRUE, full.names = TRUE)

# Load and combine
# system.time(
# sextant_1998 <- plyr::ldply(sextant_1998_dir, load_sextant, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# )
# sextant_1999 <- plyr::ldply(sextant_1999_dir, load_sextant, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# Combine and save
# sextant_1998_1999 <- rbind(sextant_1998, sextant_1999)
# save(sextant_1998, file = "data/SEXTANT/sextant_1998.RData")

# Load the data
load("data/SEXTANT/sextant_1998.RData")


# Spatial analysis --------------------------------------------------------

# Look at the monthly average SPM 
sextant_1998_monthly <- sextant_1998 |> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) |> 
  summarise(mean_spm = mean(analysed_spim, na.rm = TRUE), .by = c("lon", "lat", "year", "month"))

# Map of the 12 months
ggplot(data = sextant_1998_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  # annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  facet_wrap(~month)
  # facet_grid(year~month)


# Temporal analysis -------------------------------------------------------

# Time series analysis
sextant_1998_ts <- sextant_1998 |> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) |> 
  summarise(mean_spm = mean(analysed_spim, na.rm = TRUE), .by = c("year", "month", "date"))

# Time series plot of the monthly average SPM
ggplot(data = sextant_1998_ts, aes(x = date, y = mean_spm)) +
  geom_line() +
  geom_point() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = "Date", y = "Mean SPM", title = "Monthly Average SPM in 1998") +
  theme_minimal()

# Boxplot of the SPM per month
ggplot(data = sextant_1998_ts, aes(x = factor(month), y = mean_spm)) +
  geom_boxplot() +
  labs(x = "Month", y = "Mean SPM", title = "Monthly Average SPM in 1998") +
  theme_minimal()

# Boxplots per year
ggplot(data = sextant_1998_ts, aes(x = factor(year), y = mean_spm)) +
  geom_boxplot() +
  labs(x = "Year", y = "Mean SPM", title = "Monthly Average SPM in 1998") +
  theme_minimal()


# Animation ---------------------------------------------------------------

# Once a brick of data is loaded, it is possible to walk through one day at a time

# Animate using gganimate to show the months as the time steps
p_map <- ggplot(data = sextant_1998_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  # facet_wrap(~month) +
  labs(title = "SEXTANT SPM data for 1998", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  scale_fill_viridis_c(option = "D")

# Add animation
animated_plot <- p_map +
  transition_states(
    month,
    transition_length = 1,
    state_length = 1
  ) +
  enter_fade() +
  exit_fade()

# Render the animation
animate(animated_plot, fps = 10, duration = 10, 
        renderer = gifski_renderer(file = "animations/sextant_1998.gif",
                                   width = 1200,    # Increase width in pixels
                                   height = 1000))   # Increase height in pixels

