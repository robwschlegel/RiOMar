# func/validate.R

# This script contains the code necessary to run validation on a number of
# satellite products and variables against multiple sources of in situ data


# Setup -------------------------------------------------------------------

# The packages needed for this script
library(tidyverse)
library(raster) # Used here to select specific pixels in a .nc file
library(ncdf4)
library(sf) # Used for complex shape files
library(scales) # For better plot labels
library(ggExtra) # For histogram border plots
library(gt) # For fancy tables
library(geosphere) # For pixel distances
library(doParallel); registerDoParallel(cores = detectCores()-2)

# The shared functions
source("func/util.R")

## NB: File pathways no called within pipeline function
# Get all SEXTANT file names
# files_SEXTANT_SPM <- dir("~/pCloudDrive/data/SEXTANT/SPM", pattern = ".nc", full.names = TRUE, recursive = TRUE)
# files_SEXTANT_CHL <- dir("~/pCloudDrive/data/SEXTANT/CHLA", pattern = ".nc", full.names = TRUE, recursive = TRUE)

# Get all MODIS file names
# files_MODIS_path <- "/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/"
# files_MODIS_SPM <- dir(files_MODIS_path, recursive = TRUE, full.names = TRUE, pattern = "SPM")
# files_MODIS_CHL <- dir(files_MODIS_path, recursive = TRUE, full.names = TRUE, pattern = "CHL")
# files_MODIS_TUR <- dir(files_MODIS_path, recursive = TRUE, full.names = TRUE, pattern = "TUR")
# files_MODIS_CDOM <- dir(files_MODIS_path, recursive = TRUE, full.names = TRUE, pattern = "CDOM")
# files_MODIS_RRS <- dir(files_MODIS_path, recursive = TRUE, full.names = TRUE, pattern = "RRS")
# files_MODIS_SST <- dir(files_MODIS_path, recursive = TRUE, full.names = TRUE, pattern = "SST")


# Functions ---------------------------------------------------------------

## Stats ------------------------------------------------------------------

# The full suite of stats to calculate
compute_stats <- function(x_vec, y_vec){
  
  if(!is.numeric(x_vec)) stop("x_vec is not numeric")
  if(!is.numeric(y_vec)) stop("y_vec is not numeric")
  
  if(length(x_vec) < 3){
    return(data.frame(row.names = NULL,
                      n = length(x_vec),
                      Slope = NA,
                      Slope_log = NA,
                      RMSE = NA,
                      MSA = NA,
                      MAPE = NA,
                      Bias = NA,
                      Error = NA))
  }
  
  # Calculate RMSE (Root Mean Square Error)
  rmse <- sqrt(mean((y_vec - x_vec)^2, na.rm = TRUE))
  
  # Calculate MAPE (Mean Absolute Percentage Error)
  mape <- mean(abs((y_vec - x_vec) / x_vec), na.rm = TRUE) * 100
  
  # Calculate MSA (Mean Squared Adjustment)
  msa <- mean(abs(y_vec - x_vec), na.rm = TRUE)
  
  # Calculate linear slope
  lin_fit <- lm(y_vec ~ x_vec)
  slope <- coef(lin_fit)[2]
  
  # Calculate log-log linear slope
  log_lin_fit <- lm(log10(y_vec) ~ log10(x_vec))
  log_slope <- coef(log_lin_fit)[2]
  
  # Calculate Bias
  log_ratio <- log10(y_vec / x_vec)
  log_ratio_median <- median(log_ratio, na.rm = TRUE)
  bias_perc <- 100 * (sign(log_ratio_median) * (10^abs(log_ratio_median) - 1))
  
  # Calculate error
  log_ratio_median_abs <- median(abs(log_ratio), na.rm = TRUE)
  error_perc <- 100 * (10^log_ratio_median_abs - 1)
  
  # Combine int data.frame and exit
  return(data.frame(row.names = NULL,
                    n = length(x_vec),
                    Slope = round(slope, 2),
                    Slope_log = round(log_slope, 2),
                    RMSE = round(rmse, 6),
                    MSA = round(msa, 6),
                    MAPE = round(mape, 2),
                    Bias = round(bias_perc, 2),
                    Error = round(error_perc, 2)))
}


## Plots -----------------------------------------------------------------

# It does what it says on the tin
var_labels <- function(var_name) {
  
  if(grepl("CHL", var_name)){
    unit_x <- expression(mg~m^-3)
    unit_y <- expression(mg~m^-3)
    axis_limits <- c(10^-2, 10^2)
    var_colour <- "green4"
  } else if(grepl("SPM|TUR", var_name)){
    axis_limits <- c(10^-2, 10^3)
    var_colour<- "red4"
    if(grepl("TUR", var_name)){
      unit_x <- expression(NTU)
      unit_y <- expression(g~m^-3)
    } else {
      unit_x <- expression(g~m^-3)
      unit_y <- expression(g~m^-3)
    }
  } else if(grepl("TEMP", variable)) {
    unit_x <- expression('"°C"')
    unit_y <- expression('"°C"')
    axis_limits <- c(3, 30)
    var_colour <- "orange4"
  } else {
    stop(paste("Could not find ", variable))
  }
  
  # List and exit
  to_return <- list("unit_x" = unit_x,
                    "unit_y" = unit_y,
                    "var_colour" = var_colour, 
                    "axis_limits" = axis_limits)
  return(to_return)
}

# Convenience colour wrapper
colours_of_stations <- function() {
  
  # colour_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
  #                  "Point C" = mako(n = 1,begin = 0.85,end = 0.85), 
  #                  "Luc-sur-Mer" = mako(n = 1,begin = 0.90,end = 0.9), 
  #                  "Smile" = mako(n = 1,begin = 0.95,end = 0.95),
  #                  
  #                  "Bizeux" = viridis(n = 1,begin = 1,end = 1), # Bretagne
  #                  "Le Buron" = viridis(n = 1,begin = 0.925,end = 0.925), 
  #                  "Cézembre" = viridis(n = 1,begin = 0.85,end = 0.85), 
  #                  "Estacade" = viridis(n = 1,begin = 0.775,end = 0.775), 
  #                  "Astan" = viridis(n = 1,begin = 0.7,end = 0.7), 
  #                  "Portzic" = viridis(n = 1,begin = 0.625,end = 0.625), 
  #                  
  #                  "Antioche" = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
  #                  "pk 86" = plasma(n = 1,begin = 0.1,end = 0.1), 
  #                  "pk 52" = plasma(n = 1,begin = 0.15,end = 0.15),
  #                  "pk 30" = plasma(n = 1,begin = 0.2,end = 0.2),
  #                  "Comprian" = plasma(n = 1,begin = 0.25,end = 0.25), 
  #                  "Eyrac" = plasma(n = 1,begin = 0.3,end = 0.3), 
  #                  "Bouee 13" = plasma(n = 1,begin = 0.35,end = 0.35), 
  #                  
  #                  "Sola" = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
  #                  "Sete" = rocket(n = 1,begin = 0.85,end = 0.85), 
  #                  "Frioul" = rocket(n = 1,begin = 0.90,end = 0.90),
  #                  "Point B" = rocket(n = 1,begin = 0.95,end = 0.95)
  # ) 
  
  # colour_values = c('Manche orientale - Mer du Nord' = mako(n = 1,begin = 0.8,end = 0.8), # Manche
  #                  'Baie de Seine' = mako(n = 1,begin = 0.95,end = 0.95),
  # 
  #                  'Manche occidentale' = viridis(n = 1,begin = 1,end = 1), # Bretagne
  #                  'Bretagne Sud' = viridis(n = 1,begin = 0.625,end = 0.625),
  # 
  #                  'Pays de la Loire - Pertuis' = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
  #                  'Sud Golfe de Gascogne' = plasma(n = 1,begin = 0.35,end = 0.35),
  # 
  #                  'Golfe du Lion' = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
  #                  'Mer ligurienne - Corse' = rocket(n = 1,begin = 0.95,end = 0.95)
  # )
  
  # colour_values = c('BAY OF SEINE' = viridis::mako(n = 1, begin = 0.8,end = 0.8),
  #                  'SOUTHERN BRITTANY' = viridis::viridis(n = 1, begin = 0.9,end = 0.9),
  #                  'GULF OF BISCAY' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
  #                  'GULF OF LION' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  colour_values = c('Bay of Seine' = viridis::mako(n = 1, begin = 0.8,end = 0.8),
                   'S. Brittany' = viridis::viridis(n = 1, begin = 0.9,end = 0.9),
                   'Bay of Biscay' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
                   'Gulf of Lion' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  
  return(colour_values)
}

# The figure code wrapper
# var_name = "SPM"; match_up_df = zone_in_situ_SEXTANT; match_up_stats = stats_SEXTANT; sat_name = "SEXTANT"
validation_plots <- function(var_name, sat_name, match_up_df, match_up_stats) {
  
  # Get axis labels
  plot_meta_data <- var_labels(var_name)
  
  # Get 1:1 line limits
  identity_line <- data.frame(x = plot_meta_data$axis_limits, 
                              y = plot_meta_data$axis_limits)
  
  # Subset datasets for chosen variable
  match_up_df_var <- match_up_df |> filter(variable == var_name) |> 
    mutate(zone_pretty = factor(zone, 
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")))
  match_up_stats_var <- match_up_stats |> filter(variable == var_name)
  
  # Get colour values
  colour_values <- colours_of_stations()
  
  # Get stats to plot
  if(str_detect(var_name, 'SST')){
    # Error_value <- match_up_stats_var$error_Pahlevan_linear
    # Bias_value <- match_up_stats_var$bias_Pahlevan_linear
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope
  } else {
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope_log
  }
  
  # Correct var_name for satellite
  if(sat_name == "SEXTANT"){
    if(var_name %in% c("TUR", "SPM")){
      var_name_sat <- "SPIM"
    } else if(var_name == "CHLA"){
      var_name_sat <- "CHL"
    } else {
      stop("Impossible satellite variable matchup")
    }
    # TODO: Add other satellites to this
  } else {
    var_name_sat <- var_name
  }
  
  
  # Create title and subtitle
  plot_title = paste("In situ", var_name, "Vs.", sat_name, var_name_sat)
  
  # Create the plot
  scatterplot <- ggplot(data = match_up_df_var, aes(x = value_in_situ, y = value_satellite)) + 
    
    geom_point(aes(colour = zone_pretty, shape = source), size = 6, show.legend = TRUE) + 
    geom_line(data = identity_line, aes(x = x, y = y), linetype = "dashed", show.legend = FALSE) +
    
    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('In~situ~measurements~(', plot_meta_data$unit_x, ')'))) +
    
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Satellite~estimates~(', plot_meta_data$unit_y, ')'))) + 
    
    coord_equal(xlim = plot_meta_data$axis_limits, 
                ylim = plot_meta_data$axis_limits) +
    
    annotate(geom = 'text', x = plot_meta_data$axis_limits[1], y = plot_meta_data$axis_limits[2], 
             hjust = 0, vjust = 1, color = "black", size = 12.5,
             label = paste('Error = ', round(ifelse(Error_value %>% is.numeric(),Error_value, NA), 1), "%\n",
                           'Bias = ', round(ifelse(Bias_value %>% is.numeric(),Bias_value, NA), 1), " %\n",
                           # 'r² (linearity) = ', round(statistics_values$r2_log, 2),"\n",
                           # 'Slope = ', round(ifelse(Slope_value %>% is.numeric(),Slope_value, NA), 2),"\n",
                           'n = ', nrow(match_up_df_var), sep = "")) + 
    
    labs(title = plot_title) +
    
    scale_color_manual(name = "zone", values = colour_values, drop = FALSE) +
    
    # scale_linetype_manual(values = c("Identity line" = "dashed",
    #                                  "Linear regression" = "solid"), name = "") +
    
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 10), order = 1),
           shape = guide_legend(ncol = 1, override.aes = list(size = 10), order = 2),
           # linetype = guide_legend(override.aes = list(color = c("black"), shape = c(NA), linetype = c("dashed")), ncol = 2, order = 2)
           ) +
    
    ggplot_theme() + 
    
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 30),
          legend.margin = margin(5, 10, 5, 5),
          plot.subtitle = element_text(hjust = 0.5, color = "black", face = "bold.italic"),
          plot.title = element_text(color = plot_meta_data$var_colour, face = "bold", size = 35))
  # scatterplot
  
  if(var_name %>% str_detect("SST")){ 
    scatterplot <- scatterplot +  
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      scale_x_continuous(name = parse(text = paste0('In~situ~measurements~(', plot_meta_data$unit_x, ')'))) +
      scale_y_continuous(name = parse(text = paste0('Satellite~estimates~(', plot_meta_data$unit_y, ')')))
  } else {
    scatterplot <- scatterplot + 
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      annotation_logticks()
  }
  # scatterplot
  
  # Add histograms to x and y axes
  scatterplot_with_side_hist <- ggMarginal(scatterplot, type = "histogram", groupFill = TRUE, alpha = 1)
  # scatterplot_with_side_hist
  ggsave(paste0("figures/validation/scatterplot_",sat_name,"_",var_name,".png"), 
         plot = scatterplot_with_side_hist, height = 16, width = 16, bg = "white")
  
  # Barplot of frequency per year
  bar_plot_freq_per_year <- match_up_df_var |> 
    mutate(Year = year(date)) |> 
    dplyr::count(Year) |> 
    ggplot(aes(x = Year)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_data$var_colour, linewidth = 1.5) + 
    scale_x_continuous(breaks = seq(1998, 2025, by = 3), name = "", labels = function(x) substring(x, 3, 4)) +
    scale_y_continuous(expand = c(0,0), name = "n per year") +
    labs(title = plot_title) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_data$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_year
  ggsave(paste0("figures/validation/barplot_annual_",sat_name,"_",var_name,".png"), 
         plot = bar_plot_freq_per_year, height = 10, width = 14, bg = "white")
  
  # Barplot of monthly counts
  bar_plot_freq_per_month <- match_up_df_var |> 
    mutate(Month = month(date)) |> 
    dplyr::count(Month) |> 
    ggplot(aes(x = Month)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_data$var_colour, linewidth = 2) + 
    scale_x_continuous(breaks = 1:12, name = "", labels = function(x) month.abb[x]) +
    scale_y_continuous(expand = c(0,0), name = "n per month") +
    coord_cartesian(xlim = c(1,12)) +
    labs(title = plot_title) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_data$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_month
  ggsave(paste0("figures/validation/barplot_monthly_",sat_name,"_",var_name,".png"), 
         plot = bar_plot_freq_per_month, height = 10, width = 14, bg = "white")
  return()
}


## Tables ---------------------------------------------------------------

# Create pretty tables
# file_path = "output/MATCH_UP_DATA/FRANCE/STATISTICS/"; sat_name = "SEXTANT"
# TODO: This will need to be modified to work with other satellite products
validation_tables <- function(file_path, sat_name) {
  
  # Load all output stats
  # TODO: Clean up zone names and order by latitude
  stat_files <- map_dfr(dir(file_path, pattern = sat_name, full.names = TRUE), read_csv) |> 
    dplyr::select(zone, Bias, Error, n, variable) |> 
    mutate_at(vars(-zone, -variable), ~ round(., 1)) |> 
    mutate(var_label = paste("When compared with in situ", ifelse(variable == 'TUR', 'TURB', variable))) |>
    dplyr::select(-variable)
    # dplyr::rename(`Bias (%)` = Bias, `Error (%)` = Error)

  desired_colnames <- names(stat_files) |> 
    str_replace_all("Bias", "Bias (%)") |> 
    str_replace_all("Error", "Error (%)") |> 
    str_remove_all("var_label")
  names(desired_colnames) <- names(stat_files)
  
  # Create the SPM/TUR table
  # TODO: Wrap this into a function to be called per variable
  table_SPM_TUR <- stat_files |>
    filter(!grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite SPM', subtitle = 'Compared with SPM and Turbidity in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    # Pretty tabs
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3), # A bit more padding for summaries
                row_group.padding = px(4)) |> # And even more for our groups
    # More styling
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Create the CHLA table
  table_CHLA <- stat_files |>
    filter(grepl("CHL", var_label)) |> 
    mutate(zone = paste0("**", zone, "**")) |> 
    gt(rowname_col = 'zone', groupname_col = 'var_label', process_md = TRUE) |> 
    cols_label(.list = desired_colnames) |> 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) |> 
    tab_header(title = 'Performances of satellite CHL', subtitle = 'Compared with CHLA in situ measurements.') |> 
    sub_missing(missing_text = "-") |> 
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3),
                row_group.padding = px(4)) |> 
    opt_stylize(style = 6, color = 'gray') |> 
    tab_style(style = cell_text(align = "center"),
              locations = cells_column_labels())
  
  # Save and exit
  ## SPM/TUR
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.html"), inline_css = TRUE)
  gtsave(table_SPM_TUR, filename = paste0("figures/validation/table_",sat_name,"_SPM_TUR.png"), expand = 10)
  ## CHLA
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.html"), inline_css = TRUE)
  gtsave(table_CHLA, filename = paste0("figures/validation/table_",sat_name,"_CHLA.png"), expand = 10)
}


# Load in situ ------------------------------------------------------------

# Look at times when MODIS, MERIS, OLCI-A, and OLCI-B are overhead in the ODATIS-MR L3 products
get_start_end_time("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_BISCAY/daily/L3m_20020704__FRANCE_03_MOD_CDOM-NS_DAY_00.nc")
get_start_end_time("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/GULF_OF_LION/daily/L3m_20020704__FRANCE_03_MOD_NRRS555-NS_DAY_00.nc")
get_start_end_time("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/SOUTHERN_BRITTANY/daily/L3m_20030912__FRANCE_03_MER_CDOM-PO_DAY_00.nc")

# Get all L3 satellite daily start and end times when overhead France
## MODIS
times_MODIS <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/BAY_OF_SEINE/daily", pattern = "SPM", full.names = TRUE),
                           get_start_end_time, .parallel = TRUE)
# Double check that times are the same for any region
times_MODIS_check <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MODIS/GULF_OF_LION/daily", pattern = "SPM", full.names = TRUE),
                                 get_start_end_time, .parallel = TRUE)
times_MODIS_test <- bind_cols(times_MODIS, times_MODIS_check) |> 
  mutate(start_diff = start_time.x - start_time.y, end_diff = end_time.x - end_time.y)

## MERIS
times_MERIS <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/MERIS/BAY_OF_SEINE/daily", pattern = "SPM", full.names = TRUE),
                           get_start_end_time, .parallel = TRUE)

## OLCI-A
times_OLCI_A <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-A/BAY_OF_SEINE/daily", pattern = "SPM", full.names = TRUE),
                           get_start_end_time, .parallel = TRUE)

## OLCI-B
times_OLCI_B <- plyr::ldply(dir("/media/calanus/HDD2TB/home/calanus/data/ODATIS-MR/OLCI-B/BAY_OF_SEINE/daily", pattern = "SPM", full.names = TRUE),
                           get_start_end_time, .parallel = TRUE)

# REPHY
## NB: The script used to process the data, and necessary files, were provided by Victor Pochic
# "data/INSITU_data/REPHY/REPHY_Dataset_20250408.R"
REPHY <- read.csv2("data/INSITU_data/REPHY/Table1_REPHY_hydro_20250408.csv", fileEncoding = "ISO-8859-1") |> 
  mutate(lon = as.numeric(lon), lat = as.numeric(lat), source = "REPHY") |> 
  dplyr::rename(site = Code_point_Libelle, date = Date, variable = Code.parametre, value = Valeur_mesure)

# SOMLIT
# "data/INSITU_data/SOMLIT/SOMLIT_prep.R"
SOMLIT <- read_csv("data/INSITU_data/SOMLIT/Somlit_clean.csv") |> 
  mutate(source = "SOMLIT")

# Filter out sites that are within the RiOMar regions
# in_situ_site_list <- bind_rows(dplyr::select(REPHY, source, lon, lat, site),
#                                dplyr::select(SOMLIT, source, lon, lat, site)) |>
#   distinct() |>
#   summarise(lon = mean(lon), lat = mean(lat), .by = c("source", "site")) |>
#   mutate(zone = case_when(lon >= zones_bbox$lon_min[1] & lon <= zones_bbox$lon_max[1] &
#                             lat >= zones_bbox$lat_min[1] & lat <= zones_bbox$lat_max[1] ~ "GULF_OF_LION",
#                           lon >= zones_bbox$lon_min[2] & lon <= zones_bbox$lon_max[2] &
#                             lat >= zones_bbox$lat_min[2] & lat <= zones_bbox$lat_max[2] ~ "BAY_OF_SEINE",
#                           lon >= zones_bbox$lon_min[3] & lon <= zones_bbox$lon_max[3] &
#                             lat >= zones_bbox$lat_min[3] & lat <= zones_bbox$lat_max[3] ~ "BAY_OF_BISCAY",
#                           lon >= zones_bbox$lon_min[4] & lon <= zones_bbox$lon_max[4] &
#                             lat >= zones_bbox$lat_min[4] & lat <= zones_bbox$lat_max[4] ~ "SOUTHERN_BRITTANY")) |> 
#   mutate(zone_pretty = factor(zone, 
#                               levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
#                               labels = c("Bay of Seine", "S. Brittany", "Bay of Biscay", "Gulf of Lion")), .after = "zone") |> 
#   mutate(source = factor(source, levels = c("SOMLIT", "REPHY")))
# write_csv(in_situ_site_list, "metadata/in_situ_site_list.csv")
in_situ_site_list <- read_csv("metadata/in_situ_site_list.csv")

# Filter in situ stations to just those within a zone
zone_sites <- in_situ_site_list |> filter(!is.na(zone))

# Clean up and combine all in situ data into one dataframe
## REPHY
clean_REPHY <- right_join(REPHY, zone_sites, by = c("source", "site", "lon", "lat")) |> 
  filter(Qualite.resultat == "Bon") |> # Only keep 'good' measurements
  filter(as.numeric(Profondeur.metre) <= 10) |>  # Only keep values within 10 meters of surface
  dplyr::select(source, site, lon, lat, date, variable, value) |> 
  mutate(variable = case_when(variable == "SALI" ~ "SAL",
                              variable == "TURB" ~ "TUR",
                              variable == "CHLOROA" ~ "CHLA", TRUE ~ variable),
         date = as.Date(date)) |> 
  filter(value >= 0)

## SOMLIT
clean_SOMLIT <- right_join(SOMLIT, zone_sites, by = c("source", "site", "lon", "lat")) |> 
  filter(prof_num <= 10) |> # Filter data deeper than 10 meters
  # Filter bad flags
  mutate(TEMP = case_when(temp_QC %in% c(2, 6, 7) ~ temp),
         SAL = case_when(sal_QC %in% c(2, 6, 7) ~ sal),
         POC = case_when(POC_QC %in% c(2, 6, 7) ~ POC),
         SPM = case_when(SPM_QC %in% c(2, 6, 7) ~ SPM),
         CHLA = case_when(CHLA_QC %in% c(2, 6, 7) ~ CHLA)) |> 
  dplyr::select(source, site, lon, lat, date, TEMP, SAL, POC, SPM, CHLA) |> 
  pivot_longer(TEMP:CHLA, values_to = "value", names_to = "variable") |>
  filter(value >= 0) #|> 
  # summarise(value = mean(value, na.rm = TRUE), .by = c("source", "site", "lon", "lat", "date", "variable"))

## Combine
zone_data_in_situ <- bind_rows(clean_REPHY, clean_SOMLIT) |> 
  # Get only variables of interest
  filter(variable %in% c('TEMP', 'SAL', 'POC', 'SPM', 'CHLA', 'TUR')) |> 
  # Create daily means
  summarise(value = mean(value, na.rm = TRUE), .by = c("source", "site", "lon", "lat", "date", "variable"))


# Map In situ -------------------------------------------------------------

# Load high-res shape files
# Tuto here : https://www.etiennebacher.com/posts/2021-12-27-mapping-french-rivers-network/
if(!exists("borders_FR")) borders_FR <- read_sf("data/FRANCE_shapefile/gadm41_FRA_0.shp")
if(!exists("rivers_FR")) rivers_FR <- st_intersection(read_sf("data/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"), borders_FR)

# Map all in situ stations, highlighting the zones and the stations used
# TODO: Add text labels showing the number of REPHY and SOMLIT stations per zone
in_situ_station_map <- ggplot() +
  geom_sf(data = borders_FR, color = "black", fill = "sienna4", inherit.aes = FALSE) +
  geom_sf(data = rivers_FR, color = "lightblue", inherit.aes = FALSE, linewidth = 0.2) +
  geom_rect(data = zones_bbox, fill = NA, linewidth = 2, show.legend = FALSE,
            aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max, colour = zone_pretty)) +
  # coord_map("moll") +
  coord_sf(xlim = c(-5, 10), ylim = c(41.5, 51)) +
  geom_point(data = in_situ_site_list,
             aes(x = lon, y = lat, shape = source), color = "black", size = 4) +
  geom_point(data = filter(in_situ_site_list, is.na(zone)),
             aes(x = lon, y = lat, shape = source), color = "red", size = 3) +
  geom_point(data = filter(in_situ_site_list, !is.na(zone)),
             aes(x = lon, y = lat, shape = source, color = zone_pretty), size = 3) +
  # ggrepel::geom_text_repel(data =  filter(in_situ_site_list, source == "SOMLIT"), aes(x = lon, y = lat, label = site),
  #                 color = "red", size = 5, max.overlaps = 20) +
  # labs(title = paste(nrow(SOMLIT_stations), "SOMLIT stations"), x = NULL, y = NULL) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(labels = scales::unit_format(unit = "°E")) +
  scale_y_continuous(labels = scales::unit_format(unit = "°N")) +
  scale_color_manual(name = "zone", values = colours_of_stations(), drop = FALSE) +
  cowplot::theme_map() +
  ggplot_theme() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.73, 0.8),
        legend.box.background = element_rect(colour = "black", fill = "white"),
        legend.box.margin = margin(5, 5, 5, 5),
        axis.text = element_text(size = 20))
# in_situ_station_map
ggsave("figures/map_in_situ_stations.png", height = 14, width = 15.5, bg = "white")


# Prep satellite pixels ---------------------------------------------------

# Create the data.frames of pixel matchups
## SEXTANT
write_pixels(zone_sites, "SEXTANT")

## MODIS
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "MODIS")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "MODIS")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "MODIS")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "MODIS")

## MERIS
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "MERIS")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "MERIS")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "MERIS")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "MERIS")

## OLCI-A
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "OLCI-A")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "OLCI-A")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "OLCI-A")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "OLCI-A")

## OLCI-B
write_pixels(filter(zone_sites, zone == "BAY_OF_SEINE"), "OLCI-B")
write_pixels(filter(zone_sites, zone == "SOUTHERN_BRITTANY"), "OLCI-B")
write_pixels(filter(zone_sites, zone == "BAY_OF_BISCAY"), "OLCI-B")
write_pixels(filter(zone_sites, zone == "GULF_OF_LION"), "OLCI-B")


# Extract satellite data --------------------------------------------------

# SEXTANT
extract_pixels_all("SEXTANT")

# MODIS
# TODO: This can be optimized with a map or walk function
# But it isn't bad to have more control over the process...
# Though, one could walk through this entire section, all sensors and zones, in one function call...
# Could also use m_ply()
extract_pixels_all("MODIS", "BAY_OF_SEINE"); gc()
extract_pixels_all("MODIS", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("MODIS", "BAY_OF_BISCAY"); gc()
extract_pixels_all("MODIS", "GULF_OF_LION"); gc()

# MERIS
extract_pixels_all("MERIS", "BAY_OF_SEINE"); gc()
extract_pixels_all("MERIS", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("MERIS", "BAY_OF_BISCAY"); gc()
extract_pixels_all("MERIS", "GULF_OF_LION"); gc()

# OLCI-A
extract_pixels_all("OLCI-A", "BAY_OF_SEINE"); gc()
extract_pixels_all("OLCI-A", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("OLCI-A", "BAY_OF_BISCAY"); gc()
extract_pixels_all("OLCI-A", "GULF_OF_LION"); gc()

# OLCI-B
extract_pixels_all("OLCI-B", "BAY_OF_SEINE"); gc()
extract_pixels_all("OLCI-B", "SOUTHERN_BRITTANY"); gc()
extract_pixels_all("OLCI-B", "BAY_OF_BISCAY"); gc()
extract_pixels_all("OLCI-B", "GULF_OF_LION"); gc()


# Validation stats --------------------------------------------------------

# TODO: To turn this section into a pipeline function

# Load prepped SEXTANT data
zone_median_SEXTANT <- data.table::fread("output/MATCH_UP_DATA/FRANCE/zone_median_SEXTANT.csv") |> 
  mutate(date = as.Date(date))

# Combine extracted sat data with in situ
zone_in_situ_SEXTANT <- zone_data_in_situ |> 
  left_join(zone_median_SEXTANT, by = c("source", "site", "date", "variable")) |> 
  filter(value > 0, median > 0) |> 
  dplyr::rename(value_in_situ = value, value_satellite = median) |> 
  mutate(season = case_when(
    month(date) %in% c(12, 1, 2) ~ "Winter", month(date) %in% 3:5  ~ "Spring",
    month(date) %in% 6:8  ~ "Summer", month(date) %in% 9:11 ~ "Autumn"), .after = "date") |> 
  dplyr::select(zone, dplyr::everything())

# Create big grid of sites to look for obvious outliers
# ggplot(data = zone_in_situ_SEXTANT, aes(x = value_in_situ, y = value_satellite)) +
#   geom_point(aes(colour = source, shape = variable)) +
#   facet_wrap(~site, scales = "free")

# TODO: Think of a way to optimise this
# Stats for all groups together by variable
zone_in_situ_SEXTANT_stats_01 <- zone_in_situ_SEXTANT |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL", season = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all groups together by variable and season
zone_in_situ_SEXTANT_stats_02 <- zone_in_situ_SEXTANT |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all zones by variable
zone_in_situ_SEXTANT_stats_03 <- zone_in_situ_SEXTANT |> mutate(source = "ALL", site = "ALL", season = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all zones by variable and season
zone_in_situ_SEXTANT_stats_04 <- zone_in_situ_SEXTANT |> mutate(source = "ALL", site = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all sources by variable
zone_in_situ_SEXTANT_stats_05 <- zone_in_situ_SEXTANT |> mutate(zone = "GLOBAL", site = "ALL", season = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all sources by variable and season
zone_in_situ_SEXTANT_stats_06 <- zone_in_situ_SEXTANT |> mutate(zone = "GLOBAL", site = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all zones and sources by variable
zone_in_situ_SEXTANT_stats_07 <- zone_in_situ_SEXTANT |> mutate(site = "ALL", season = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all zones and sources by variable and season
zone_in_situ_SEXTANT_stats_08 <- zone_in_situ_SEXTANT |> mutate(site = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all sites by variable
zone_in_situ_SEXTANT_stats_09 <- zone_in_situ_SEXTANT |> mutate(site = "ALL", season = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for all sites by variable and season
zone_in_situ_SEXTANT_stats_10 <- zone_in_situ_SEXTANT |> mutate(site = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for each site by variable
zone_in_situ_SEXTANT_stats_11 <- zone_in_situ_SEXTANT |> mutate(season = "ALL") |> 
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Stats for each site by variable and season
zone_in_situ_SEXTANT_stats_12 <- zone_in_situ_SEXTANT |>
  summarise(compute_stats(value_in_situ, value_satellite), .by = c("zone", "source", "site", "season", "variable"))

# Bind all together
zone_in_situ_SEXTANT_stats <- bind_rows(zone_in_situ_SEXTANT_stats_01, zone_in_situ_SEXTANT_stats_02, zone_in_situ_SEXTANT_stats_03,
                                        zone_in_situ_SEXTANT_stats_04, zone_in_situ_SEXTANT_stats_05, zone_in_situ_SEXTANT_stats_06,
                                        zone_in_situ_SEXTANT_stats_07, zone_in_situ_SEXTANT_stats_08, zone_in_situ_SEXTANT_stats_09,
                                        zone_in_situ_SEXTANT_stats_10, zone_in_situ_SEXTANT_stats_11, zone_in_situ_SEXTANT_stats_12) |> 
  mutate(sensor = "SEXTANT", .before = "zone")

# Save results
write_csv(zone_in_situ_SEXTANT_stats, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/","SEXTANT","_stats.csv"))

# Plot to look for odd results
# TODO: Consider running ANOVAs on the base values per the different groupings
# Or rather run them on the resulting stats visualised here
## Histogram
# ggplot(data = zone_in_situ_SEXTANT_stats, aes(x = Bias)) +
#   geom_histogram(aes(fill = variable))#, position = "dodge")
# ggplot(data = zone_in_situ_SEXTANT_stats, aes(x = Bias, y = Error)) +
#   geom_point(aes(colour = zone, shape = variable, size = source)) +
#   scale_x_continuous(limits = c(-300, 300)) +
#   scale_y_continuous(limits = c(0, 300))
# ggplot(data = zone_in_situ_SEXTANT_stats, aes(y = Error)) +
#   geom_boxplot(aes(fill = zone)) +
#   facet_grid(variable~season, scales = "free_y") +
#   coord_cartesian(ylim = c(0, 200))
#   # coord_cartesian(ylim = c(-200, 200))


#  Validation figures -----------------------------------------------------

# TODO: To turn this section into a pipeline function

# Reload base data matchup
zone_median_SEXTANT <- data.table::fread("output/MATCH_UP_DATA/FRANCE/zone_median_SEXTANT.csv")

# Reload stats
stats_SEXTANT <- map_dfr(dir("output/MATCH_UP_DATA/FRANCE/STATISTICS/", 
                             pattern = "SEXTANT", full.names = TRUE), read_csv)

# Create figures
plyr::l_ply(unique(stats_SEXTANT$variable), validation_plots, .parallel = TRUE,
            sat_name = "SEXTANT", match_up_df = zone_in_situ_SEXTANT, match_up_stats = stats_SEXTANT)


# Validation tables -------------------------------------------------------

validation_tables("output/MATCH_UP_DATA/FRANCE/STATISTICS/", "SEXTANT")

