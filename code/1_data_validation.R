library(plyr)
library(tidyverse)
library(scales) 
library(ggExtra)
library(ggpubr)
library(patchwork)
library(viridis)
library(gridExtra)
library(ggpmisc)
library(maps)
library(ggthemes)
library(mapproj)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(MASS)
library(Metrics)
library(doParallel); registerDoParallel(cores = detectCores()/2)

#### Define path ####

work_dir <- "/home/terrats/Desktop/RIOMAR/"
path_to_SOMLIT_data_with_MU_results <- file.path(work_dir, "DATA/IN_SITU/")
path_to_river_data <- file.path(work_dir, "DATA/RIVERS/FOR_MAPS/FILES/")
where_to_save_MU_results <- file.path(work_dir, "RESULTS/MATCH_UP/")

#### Configuration ####

max_CV <- 30 # 30
min_n <- NA # 5
Maximum_time_difference = NA # 3 # in hour
vars_to_plot <- c("SPM-G", "CHLA") # c("SST-DAY", 'SST-NIGHT') 
# stations_to_remove <- c('Sola', 'Bouee 13', 'Eyrac', 'Comprian')
stations_to_remove <- c()
Version_best_spatio_temporal_grid <- 1 # V1 show good results
SOMLIT_QC_codes_to_keep <- c(2,6,7)

#### Load data ####

df <- load_SOMLIT_data(path_to_SOMLIT_data_with_MU_results, stations_to_remove)

color_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
                 "Point C" = mako(n = 1,begin = 0.85,end = 0.85), 
                 "Luc-sur-Mer" = mako(n = 1,begin = 0.90,end = 0.9), 
                 "Smile" = mako(n = 1,begin = 0.95,end = 0.95),
                 
                 "Bizeux" = viridis(n = 1,begin = 1,end = 1), # Bretagne
                 "Le Buron" = viridis(n = 1,begin = 0.925,end = 0.925), 
                 "Cézembre" = viridis(n = 1,begin = 0.85,end = 0.85), 
                 "Estacade" = viridis(n = 1,begin = 0.775,end = 0.775), 
                 "Astan" = viridis(n = 1,begin = 0.7,end = 0.7), 
                 "Portzic" = viridis(n = 1,begin = 0.625,end = 0.625), 
                 
                 "Antioche" = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
                 "pk 86" = plasma(n = 1,begin = 0.1,end = 0.1), 
                 "pk 52" = plasma(n = 1,begin = 0.15,end = 0.15),
                 "pk 30" = plasma(n = 1,begin = 0.2,end = 0.2),
                 "Comprian" = plasma(n = 1,begin = 0.25,end = 0.25), 
                 "Eyrac" = plasma(n = 1,begin = 0.3,end = 0.3), 
                 "Bouee 13" = plasma(n = 1,begin = 0.35,end = 0.35), 
                 
                 "Sola" = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
                 "Sete" = rocket(n = 1,begin = 0.85,end = 0.85), 
                 "Frioul" = rocket(n = 1,begin = 0.90,end = 0.90),
                 "Point B" = rocket(n = 1,begin = 0.95,end = 0.95)
) 

SOMLIT_stations <- df %>% distinct(Latitude, Longitude, Site)

# world_map = map_data("world") %>% filter(! long > 180)
# # France_boundaries <- ne_countries(country = "France", scale = 'medium', type = 'map_units', returnclass = 'sf')  
# France_boundaries <- read_sf(file.path(work_dir, "DATA", "FRANCE_shapefile", "gadm41_FRA_0.shp"))  
# riversData <- read_sf(file.path(path_to_river_data, "HydroRIVERS_v10_eu.shp")) 
# riversData <- riversData %>% st_intersection(France_boundaries)
# # Tuto here : https://www.etiennebacher.com/posts/2021-12-27-mapping-french-rivers-network/
# # 
# SOMLIT_station_map <- ggplot() +
#   # geom_polygon(data = world_map, aes(x=long, y = lat, group = group), fill = "black") +
#   geom_sf(data=France_boundaries, color="black", inherit.aes = FALSE) +
#   geom_sf(data=riversData, color="white", inherit.aes = FALSE) +
#   # geom_sf(data=riversData, color="white", inherit.aes = FALSE) +
#   # expand_limits(x = c(-8,11), y = c(52, 40)) +
#   # coord_map("moll") +
#   coord_sf(xlim = c(-8,11), ylim = c(40,52)) +
#   geom_point(data = SOMLIT_stations, aes(x = Longitude, y = Latitude), color = "red", size = 3) +
#   geom_text_repel(data = SOMLIT_stations, aes(x = Longitude, y = Latitude, label = Site),
#                   color = "red", size = 5, max.overlaps = 20) +
#   labs(title = paste(nrow(SOMLIT_stations), "SOMLIT stations"), x = "", y = "") +
#   scale_x_continuous(label = function(x) {paste(x, "°E", sep = "")}) +
#   scale_y_continuous(label = function(x) {paste(x, "°N", sep = "")}) +
#   theme_map() +
#   ggplot_theme() +
#   theme(plot.title = element_text(color = "red", size = 30),
#         axis.text = element_text(size = 20))
# # 
# # save_plot_as_png(SOMLIT_station_map, name = "SOMLIT_station_map", width = 16, height = 15, path = where_to_save_MU_results)

sat_products_to_process <- names(df) %>% 
  str_subset("_value$") %>% 
  str_replace("_value", "") %>% 
  str_replace("CHLA_|SPM-G_|SPM-R_|SST-DAY_|SST-NIGHT_", "") %>% 
  str_replace("_[0-9]x[0-9]", "") %>% 
  unique() %>% 
  sort()

valid_MU_criteria_for_title = paste("Match-up criteria : Coefficient of Variation ≤ ", max_CV, "% / n ≥ ", min_n, ' / time difference ≤ ', Maximum_time_difference, "h", sep = "")

statistics_and_scatterplots <- sat_products_to_process %>% llply(function(sat_product_to_process) {
  
  min_n_of_the_product <- ifelse(str_detect(sat_product_to_process, pattern = "SEXTANT_merged|polymer") & is.finite(min_n), 1, min_n)
  
  SPM_CHLA_plots <- vars_to_plot %>% llply(function(variable) {
    
    insitu_variable <- case_when(variable %in% c("SPM-G", "SPM-R") ~ "SPM", .default = variable) %>% str_replace("-DAY|-NIGHT", "")
    sat_product_full_name <- paste(variable, sat_product_to_process, sep = "_")
    
    if ( any(str_detect(names(df), sat_product_full_name)) == FALSE ) { return() }
    
    df_formatted <- SOMLIT_stations$Site %>% llply(function(Site_name) {
      
      df_Site <- get_sat_data_usign_the_optimal_grid(df = df, Site_name = Site_name, 
                                                     Version_optimized_grid = Version_best_spatio_temporal_grid,
                                                     sat_product_full_name = sat_product_full_name)
      
      return(df_Site)
      
    }, .inform = TRUE) %>% rbind.fill()
    
    df_of_the_sat_product <- df_formatted %>% 
      
      dplyr::rename(insitu_value := {{insitu_variable}},
                    insitu_qc := !!sym(paste("q", insitu_variable, sep = ""))) %>% 
      
      # mutate(sat_value = case_when( grepl("polymer", sat_product_to_process) ~ sat_value, 
      #                               .default = sat_median)) %>% # 1x1 grid only for polymer atm correction
      mutate(sat_value = sat_median) %>% # 1x1 grid only for polymer atm correction
      
      filter_at(vars(all_of(c("insitu_value", "sat_value"))), ~ is.finite(.)) %>% 
      filter(insitu_qc %in% SOMLIT_QC_codes_to_keep) %>% 
      
      mutate(sat_CV = 100 * sat_sd / sat_value,
             Year = lubridate::year(DATE),
             Month = lubridate::month(DATE))  
    
    if (max_CV %>% is.finite()) {df_of_the_sat_product <- df_of_the_sat_product %>% filter(sat_CV <= max_CV) }  
    if (min_n_of_the_product %>% is.finite()) {df_of_the_sat_product <- df_of_the_sat_product %>% filter(sat_n >= min_n_of_the_product) }
    
    if (nrow(df_of_the_sat_product) == 0) { return() }
    
    df_with_only_positive_values <- df_of_the_sat_product %>% filter_at(vars(ends_with("value")), ~ . > 0)
    
    statistics_df <- compute_stats(sat_values = df_with_only_positive_values$sat_value, 
                                   insitu_values = df_with_only_positive_values$insitu_value)
    
    if (str_detect(variable, 'SST')) {
      statistics_df['bias_Pahlevan'] = statistics_df['bias_Pahlevan_linear']
      statistics_df['error_Pahlevan'] = statistics_df['error_Pahlevan_linear']
      statistics_df['slope_log'] = statistics_df['slope']
    }
    
    Figures <- make_the_figures(df_with_only_positive_values, var_to_assess = variable, statistics_outputs = statistics_df)
    
    # if (variable %>% str_detect("SST")) { 
    save_plot_as_png(Figures$scatterplot_with_side_histograms, name = paste(variable, "V", Version_best_spatio_temporal_grid, sep = "_"), 
                     path = file.path(where_to_save_MU_results, "SCATTERPLOTS", "Per_sat_product"), 
                     width = 18, height = 14)
    # }
    
    return( list("scatterplot" = Figures$scatterplot, 
                 "scatterplot_with_side_hist" = Figures$scatterplot_with_side_hist, 
                 "freq_per_year" = Figures$bar_plot_freq_per_year,
                 "freq_per_month" = Figures$bar_plot_freq_per_month,
                 "statistics" = tibble(Sat_product = sat_product_to_process,
                                       Variable = variable,
                                       statistics_df %>% purrr::discard_at(vars(starts_with("lm_line"))) %>% data.frame())) )
    
  }, .inform = TRUE) %>% setNames(vars_to_plot)
  
  if (SPM_CHLA_plots %>% llply(is.null) %>% unlist() %>% all()) {
    return()
  }
  
  return( list("statistics" = ldply(SPM_CHLA_plots, function(x) x$statistics, .id = NULL),
               "scatterplots" = SPM_CHLA_plots %>% llply(function(x) {x$scatterplot}) %>% setNames(names(SPM_CHLA_plots))) )
  
}, .inform = TRUE, .parallel = TRUE) %>% 
  set_names(sat_products_to_process)



vars_to_plot %>% l_ply(function(var) {
  
  SEXTANT_plots <- c("SEXTANT_meris", "SEXTANT_modis", "SEXTANT_seawifs", "SEXTANT_viirs", "SEXTANT_merged") %>% paste("_Standard", sep = "") %>% 
    llply(function(x) {
      if (var %in% names(statistics_and_scatterplots[[x]]$scatterplots) == FALSE) {
        return()
      }
      statistics_and_scatterplots[[x]]$scatterplots[[var]] + 
        labs(title = str_replace_all(x, "SEXTANT|_|Standard", replacement = "")) + 
        guides(color=guide_legend(ncol=1), 
               linetype = guide_legend(ncol = 1,
                                       override.aes = list(color = c("black"),
                                                           shape = c(NA),
                                                           linetype = c("dashed")))) + 
        theme(legend.box.background = element_rect(color = "white"), plot.title = element_text(size = 50))}) 
  
  SEXTANT_plots <- ggarrange(plotlist = SEXTANT_plots, ncol = 4, nrow = 2, align = 'hv', common.legend = TRUE, legend = "right")
  
  SEXTANT_plots <- annotate_figure(
    annotate_figure(SEXTANT_plots, 
                    top=text_grob(valid_MU_criteria_for_title %>% paste(" -  For the merged product : n ≥ 1 / time difference ≤ 1d"), 
                                  face = "italic", size = 35 * 1.5, color = "black")),
    top=text_grob("SEXTANT", face = "bold", size = 45*1.5, color = "black"),
  )
  
  ODATIS_plots <- c("acolite", "polymer", "nirswir") %>% 
    
    llply(function(atm_correction) {
      
      sensor_plots <- c("meris", "modis", "olcia", "olcib") %>% llply(function(sensor_name) {
        
        product_name <- paste("ODATIS", sensor_name, atm_correction, sep = "_")
        
        if (product_name %in% names(statistics_and_scatterplots) == FALSE) {
          return(NULL)
        }
        
        if (var %in% names(statistics_and_scatterplots[[product_name]]$scatterplots) == FALSE) {
          return()
        }
        
        statistics_and_scatterplots[[product_name]]$scatterplots[[var]] + 
          labs(title = str_replace_all(product_name, "ODATIS|_", replacement = " ")) + 
          guides(color=guide_legend(ncol=1), 
                 linetype = guide_legend(ncol = 1,
                                         override.aes = list(color = c("black"),
                                                             shape = c(NA),
                                                             linetype = c("dashed")))) + 
          theme(legend.box.background = element_rect(color = "white"), plot.title = element_text(size = 40))  
        
      })
      
      return(sensor_plots)
      
    }) 
  
  ODATIS_plots <- ggarrange(plotlist = ODATIS_plots %>% flatten(), 
                            ncol = 4, nrow = 3, 
                            align = 'hv', common.legend = TRUE, legend = "right")
  
  ODATIS_plots <- annotate_figure(
    annotate_figure(ODATIS_plots, 
                    top=text_grob(valid_MU_criteria_for_title, face = "italic", size = 50, color = "black")),
    top=text_grob("ODATIS", face = "bold", size = 70, color = "black"),
  )
  
  final_plots = ggarrange(ODATIS_plots, SEXTANT_plots, ncol = 1, nrow = 2, align = "hv", common.legend = TRUE, heights = c(0.15, 0.1))
  
  save_plot_as_png(final_plots, name = paste("Product_comparison", var, "V", Version_best_spatio_temporal_grid, sep = "_"), 
                   path = file.path(where_to_save_MU_results, "SCATTERPLOTS"),
                   width = 60, height = 60, res = 100)
  
}, .inform = TRUE, .parallel = TRUE)


## Various functions to investigate


ggplot_theme <-   function() {
  theme(text = element_text(size=35, colour = "black"), #25
        plot.title = element_text(hjust = 0.5, size = 55),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x=element_text(angle=0),
        axis.ticks.length=unit(.25, "cm"))}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

save_plot_as_png <- function(plot, name = c(), width = 14, height = 8.27, path, res = 150) {
  
  graphics.off()
  
  if (name %>% length() == 1) {
    if (dir.exists(path) == FALSE) {dir.create(path, recursive = TRUE)}
    path <- file.path(path, paste(name, ".png", sep = ""))
  } else {
    path <- paste(path, ".png", sep = "")
  }
  
  if (grepl(pattern = ".png.png", path)) {path <- path %>% gsub(pattern = ".png.png", replacement = ".png", x = .)}
  
  png(path, width = width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
  
}

compute_Jamet_score <- function(stat_df, statistics) {
  
  all_score_values <- statistics %>% llply(function(statistic) {
    
    stat_values <- stat_df[,statistic]
    
    if (statistic %in% c("Slope")) {
      invert_abs_stat_values <- abs(1-stat_values)
      score_values <- (invert_abs_stat_values - max(invert_abs_stat_values)) / 
        (min(invert_abs_stat_values) - max(invert_abs_stat_values))
    }
    
    if (statistic %in% c("Intercept", "Bias", "Bias_in_perc")) {
      abs_stat_values <- abs(stat_values)
      score_values <- (abs_stat_values - max(abs_stat_values)) / 
        (min(abs_stat_values) - max(abs_stat_values))
    }
    
    if (statistic %in% c("R2")) {
      score_values <- (stat_values - min(stat_values)) / (max(stat_values) - min(stat_values))
    }
    
    if (statistic %in% c("RMSE", "MAPE", "nRMSE", "MedAE", "MedAPE", "MAE")) {
      score_values <- (stat_values - max(stat_values)) / (min(stat_values) - max(stat_values))
    }
    
    return( tibble(!!paste(statistic, "score", sep = "_") := score_values ) )
    
  }) %>% 
    as.data.frame() %>% 
    mutate(Total_score = rowSums(.))
  
  return( all_score_values )
  
} 


compute_stats <- function(sat_values, insitu_values, axis_limits = c(NA, NA)) {
  
  # outlier_threshold <- 1.25 * quantile(abs(sat_values - insitu_values), probs = 0.95)
  # index_to_keep <- which( abs(sat_values - insitu_values) < outlier_threshold )
  # insitu_values <- insitu_values[index_to_keep]
  # sat_values <- sat_values[index_to_keep]
  
  MedAE_multiplicative <- 10 ^ mean( abs( log10(sat_values) - log10(insitu_values) ) )
  MedAPE <- (MedAE_multiplicative - 1) * 100
  MedAE <- (MedAPE/100) * mean(insitu_values)
  
  Bias_multiplicative <- 10 ^ mean( log10(sat_values) - log10(insitu_values) )
  Bias_multiplicative_in_perc <- (Bias_multiplicative - 1)  * 100
  Bias_multiplicative <- (Bias_multiplicative_in_perc/100) * mean(insitu_values)  
  
  RMSE <- rmse(actual = insitu_values, predicted = sat_values)
  RMSE_in_perc = 100 * RMSE / mean(insitu_values)
  MAPE <- 100 * mape(actual = insitu_values, predicted = sat_values)
  
  log_sat_values <- log(sat_values)
  log_insitu_values <- log(insitu_values)
    
  weights_log <- rlm(log_sat_values ~ log_insitu_values)$weights
  index_to_keep_log <- which(weights_log == 1)
  
  lm_log <- lm(log_sat_values[index_to_keep_log] ~ log_insitu_values[index_to_keep_log]) %>% summary()
  slope_log <- lm_log$coefficients[2,1]
  intercept_log <- lm_log$coefficients[1,1]
  r2_log <- lm_log$r.squared
  
  lm_line_log <- data.frame(x = exp( log_insitu_values[index_to_keep_log] ) ,
                        y_pred = exp (log_insitu_values[index_to_keep_log]*slope_log) * exp( intercept_log ),
                        true_y = exp( log_sat_values[index_to_keep_log] ) ) %>% 
    add_row(data.frame(x = axis_limits, y_pred = exp( (log(axis_limits)*slope_log) ) * exp(intercept_log)))
  
  
  weights <- rlm(sat_values ~ insitu_values)$weights
  index_to_keep <- which(weights == 1)
  
  lm <- lm(sat_values[index_to_keep] ~ insitu_values[index_to_keep]) %>% summary()
  slope <- lm$coefficients[2,1]
  intercept <- lm$coefficients[1,1]
  r2 <- lm$r.squared
  
  lm_line <- data.frame(x = insitu_values[index_to_keep] ,
                        y_pred = insitu_values[index_to_keep]*slope + intercept,
                        true_y = sat_values[index_to_keep] ) %>% 
              add_row(data.frame(x = axis_limits, y_pred = axis_limits*slope + intercept))
  
  
  bias_Pahlevan = (sat_values / insitu_values) %>% log10() %>% median()
  bias_Pahlevan = 100 * sign(bias_Pahlevan) * (10^abs(bias_Pahlevan) - 1)
    
  error_Pahlevan = (sat_values / insitu_values) %>% log10() %>% abs() %>% median()
  error_Pahlevan = 100 * (10^error_Pahlevan - 1)
  
  bias_Pahlevan_linear = ((sat_values - insitu_values) / insitu_values) %>% median()
  bias_Pahlevan_linear = 100 * (bias_Pahlevan_linear)
  
  error_Pahlevan_linear = ((sat_values - insitu_values) / insitu_values) %>% abs() %>% median()
  error_Pahlevan_linear = 100 * (error_Pahlevan_linear)
  
  return(list(MedAE_multiplicative = MedAE_multiplicative,
              MedAPE = MedAPE,
              MedAE = MedAE,
              Bias_multiplicative = Bias_multiplicative,
              Bias_multiplicative_in_perc = Bias_multiplicative_in_perc,
              RMSE = RMSE,
              RMSE_in_perc = RMSE_in_perc,
              MAPE = MAPE,
              slope_log = slope_log,
              intercept_log = intercept_log,
              r2_log = r2_log,
              lm_line_log = lm_line_log,
              slope = slope,
              intercept = intercept,
              r2 = r2,
              lm_line = lm_line,
              bias_Pahlevan = bias_Pahlevan,
              error_Pahlevan = error_Pahlevan,
              bias_Pahlevan_linear = bias_Pahlevan_linear,
              error_Pahlevan_linear = error_Pahlevan_linear,
              n = length(insitu_values)))
  
}



get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


find_satellite_date_range <- function(satellite_name) {
  
  case_when(satellite_name == 'seawifs' ~ c("1997-08-01", "2010-31-12"),
            satellite_name == 'meris' ~ c("2002-04-28", "2012-04-09"),
            satellite_name == 'modis' ~ c("2002-05-04", Sys.Date() %>% as.character()),
            satellite_name == 'viirs' ~ c("2011-10-28", Sys.Date() %>% as.character()),
            satellite_name == 'olcia' ~ c("2016-02-16", Sys.Date() %>% as.character()),
            satellite_name == 'olcib' ~ c("2018-04-25", Sys.Date() %>% as.character()),
            satellite_name == 'merged' ~ c("1998-01-01", Sys.Date() %>% as.character()))
  
}


best_spatio_temporal_grid <- function(Site_name, Version) {
  
  if (Version == 0) {
  
    best_strategies <- list(
      "Antioche" = list("CHLA_ODATIS" = c(3, "3x3"), 
                        "CHLA_SEXTANT" = c(1, "3x3"),
                        "SPM-G_ODATIS" = c(9, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "3x3")),
      "Astan" = list( "CHLA_ODATIS" = c(9, "5x5"),
                      "CHLA_SEXTANT" = c(1, "3x3"),
                      "SPM-G_ODATIS" = c(1, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "3x3")),
      "Bizeux" = list( "CHLA_ODATIS" = c(1, "3x3"),
                       "CHLA_SEXTANT" = c(1, "3x3"),
                       "SPM-G_ODATIS" = c(9, "5x5"),
                       "SPM-G_SEXTANT" = c(1, "None")),
      "Bouee 13" = list( "CHLA_ODATIS" = c(9, "1x1"),
                         "CHLA_SEXTANT" = c(1, "None"),
                         "SPM-G_ODATIS" = c(1, "1x1"),
                         "SPM-G_SEXTANT" = c(1, "None")),
      "Cézembre" = list( "CHLA_ODATIS" = c(1, "5x5"),
                         "CHLA_SEXTANT" = c(1, "1x1"),
                         "SPM-G_ODATIS" = c(1, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "3x3")),
      "Comprian" = list( "CHLA_ODATIS" = c(3, "1x1"), 
                         "CHLA_SEXTANT" = c(1, "1x1"),
                         "SPM-G_ODATIS" = c(3, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "1x1")),
      "Estacade" = list( "CHLA_ODATIS" = c(9, "9x9"), 
                         "CHLA_SEXTANT" = c(1, "1x1"),
                         "SPM-G_ODATIS" = c(3, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "5x5")),
      "Eyrac" = list( "CHLA_ODATIS" = c(3, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(1, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "None")),
      "Frioul" = list( "CHLA_ODATIS" = c(1, "3x3"), 
                       "CHLA_SEXTANT" = c(1, "1x1"),
                       "SPM-G_ODATIS" = c(1, "3x3"),
                       "SPM-G_SEXTANT" = c(1, "3x3")),
      "Le Buron" = list( "CHLA_ODATIS" = c(3, "5x5"), 
                         "CHLA_SEXTANT" = c(1, "None"),
                         "SPM-G_ODATIS" = c(9, "1x1"),
                         "SPM-G_SEXTANT" = c(1, "None")),
      "Luc-sur-Mer" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                            "CHLA_SEXTANT" = c(1, "3x3"),
                            "SPM-G_ODATIS" = c(9, "1x1"),
                            "SPM-G_SEXTANT" = c(1, "3x3")),
      "pk 30" = list( "CHLA_ODATIS" = c(1, "3x3"), 
                      "CHLA_SEXTANT" = c(1, "3x3"),
                      "SPM-G_ODATIS" = c(1, "5x5"),
                      "SPM-G_SEXTANT" = c(1, "3x3")),
      "pk 52" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "3x3"),
                      "SPM-G_ODATIS" = c(1, "3x3"),
                      "SPM-G_SEXTANT" = c(1, "None")),
      "pk 86" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "3x3"),
                      "SPM-G_ODATIS" = c(1, "5x5"),
                      "SPM-G_SEXTANT" = c(1, "1x1")),
      "Point B" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "3x3"),
                        "SPM-G_ODATIS" = c(1, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "1x1")),
      "Point C" = list( "CHLA_ODATIS" = c(9, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "3x3"),
                        "SPM-G_ODATIS" = c(1, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "3x3")),
      "Point L" = list( "CHLA_ODATIS" = c(5, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "None"),
                        "SPM-G_ODATIS" = c(9, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "None")),
      "Portzic" = list( "CHLA_ODATIS" = c(9, "5x5"), 
                        "CHLA_SEXTANT" = c(1, "3x3"), 
                        "SPM-G_ODATIS" = c(5, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "None")),
      "Sete" = list( "CHLA_ODATIS" = c(1, "7x7"), 
                     "CHLA_SEXTANT" = c(1, "3x3"), 
                     "SPM-G_ODATIS" = c(1, "1x1"),
                     "SPM-G_SEXTANT" = c(1, "None")),
      "Smile" = list( "CHLA_ODATIS" = c(5, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "3x3"), 
                      "SPM-G_ODATIS" = c(9, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "1x1")),
      "Sola" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                     "CHLA_SEXTANT" = c(1, "1x1"), 
                     "SPM-G_ODATIS" = c(1, "3x3"),
                     "SPM-G_SEXTANT" = c(1, "1x1"))
    )  
    
  }
  
  if (Version == 1) {
    
    best_strategies <- list(
      "Antioche" = list("CHLA_ODATIS" = c(3, "3x3"), 
                        "CHLA_SEXTANT" = c(1, "9x9"),
                        "SPM-G_ODATIS" = c(9, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "3x3"),
                        "SST-DAY_ODATIS" = c(9, "1x1"),
                        "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Astan" = list( "CHLA_ODATIS" = c(9, "5x5"),
                      "CHLA_SEXTANT" = c(1, "9x9"),
                      "SPM-G_ODATIS" = c(1, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "9x9"),
                      "SST-DAY_ODATIS" = c(1, "1x1"),
                      "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Bizeux" = list( "CHLA_ODATIS" = c(1, "3x3"),
                       "CHLA_SEXTANT" = c(1, "3x3"),
                       "SPM-G_ODATIS" = c(9, "5x5"),
                       "SPM-G_SEXTANT" = c(1, "9x9"),
                       "SST-DAY_ODATIS" = c(9, "1x1"),
                       "SST-NIGHT_ODATIS" = c(20, "3x3")),
      "Bouee 13" = list( "CHLA_ODATIS" = c(9, "1x1"),
                         "CHLA_SEXTANT" = c(1, "9x9"),
                         "SPM-G_ODATIS" = c(1, "1x1"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST-DAY_ODATIS" = c(1, "3x3"),
                         "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Cézembre" = list( "CHLA_ODATIS" = c(1, "5x5"),
                         "CHLA_SEXTANT" = c(1, "9x9"),
                         "SPM-G_ODATIS" = c(1, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "3x3"),
                         "SST-DAY_ODATIS" = c(3, "1x1"),
                         "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Comprian" = list( "CHLA_ODATIS" = c(3, "1x1"), 
                         "CHLA_SEXTANT" = c(1, "7x7"),
                         "SPM-G_ODATIS" = c(3, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST-DAY_ODATIS" = c(1, "5x5"),
                         "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Estacade" = list( "CHLA_ODATIS" = c(9, "9x9"), 
                         "CHLA_SEXTANT" = c(1, "3x3"),
                         "SPM-G_ODATIS" = c(3, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "5x5"),
                         "SST-DAY_ODATIS" = c(1, "7x7"),
                         "SST-NIGHT_ODATIS" = c(20, "1x1")),
      "Eyrac" = list( "CHLA_ODATIS" = c(3, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(1, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "None"),
                      "SST-DAY_ODATIS" = c(3, "1x1"),
                      "SST-NIGHT_ODATIS" = c(10, "3x3")),
      "Frioul" = list( "CHLA_ODATIS" = c(1, "3x3"), 
                       "CHLA_SEXTANT" = c(1, "5x5"),
                       "SPM-G_ODATIS" = c(1, "3x3"),
                       "SPM-G_SEXTANT" = c(1, "5x5"),
                       "SST-DAY_ODATIS" = c(1, "3x3"),
                       "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Le Buron" = list( "CHLA_ODATIS" = c(3, "5x5"), 
                         "CHLA_SEXTANT" = c(1, "9x9"),
                         "SPM-G_ODATIS" = c(9, "1x1"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST-DAY_ODATIS" = c(3, "7x7"),
                         "SST-NIGHT_ODATIS" = c(10, "7x7")),
      "Luc-sur-Mer" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                            "CHLA_SEXTANT" = c(1, "9x9"),
                            "SPM-G_ODATIS" = c(9, "1x1"),
                            "SPM-G_SEXTANT" = c(1, "9x9"),
                            "SST-DAY_ODATIS" = c(1, "1x1"),
                            "SST-NIGHT_ODATIS" = c(20, "9x9")),
      "pk 30" = list( "CHLA_ODATIS" = c(1, "3x3"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(1, "5x5"),
                      "SPM-G_SEXTANT" = c(1, "None"),
                      "SST-DAY_ODATIS" = c(5, "1x1"),
                      "SST-NIGHT_ODATIS" = c(20, "5x5")),
      "pk 52" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(3, "9x9"),
                      "SPM-G_SEXTANT" = c(1, "None"),
                      "SST-DAY_ODATIS" = c(1, "1x1"),
                      "SST-NIGHT_ODATIS" = c(10, "3x3")),
      "pk 86" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "7x7"),
                      "SPM-G_ODATIS" = c(3, "5x5"),
                      "SPM-G_SEXTANT" = c(1, "9x9"),
                      "SST-DAY_ODATIS" = c(1, "9x9"),
                      "SST-NIGHT_ODATIS" = c(15, "1x1")),
      "Point B" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "3x3"),
                        "SPM-G_ODATIS" = c(1, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "1x1"),
                        "SST-DAY_ODATIS" = c(1, "1x1"),
                        "SST-NIGHT_ODATIS" = c(10, "7x7")),
      "Point C" = list( "CHLA_ODATIS" = c(9, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "9x9"),
                        "SPM-G_ODATIS" = c(1, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "5x5"),
                        "SST-DAY_ODATIS" = c(3, "1x1"),
                        "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Point L" = list( "CHLA_ODATIS" = c(5, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "9x9"),
                        "SPM-G_ODATIS" = c(9, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "9x9"),
                        "SST-DAY_ODATIS" = c(1, "1x1"),
                        "SST-NIGHT_ODATIS" = c(10, "5x5")),
      "Portzic" = list( "CHLA_ODATIS" = c(9, "5x5"), 
                        "CHLA_SEXTANT" = c(1, "3x3"), 
                        "SPM-G_ODATIS" = c(5, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "9x9"),
                        "SST-DAY_ODATIS" = c(3, "9x9"),
                        "SST-NIGHT_ODATIS" = c(10, "7x7")),
      "Sete" = list( "CHLA_ODATIS" = c(1, "7x7"), 
                     "CHLA_SEXTANT" = c(1, "7x7"), 
                     "SPM-G_ODATIS" = c(1, "1x1"),
                     "SPM-G_SEXTANT" = c(1, "9x9"),
                     "SST-DAY_ODATIS" = c(1, "1x1"),
                     "SST-NIGHT_ODATIS" = c(20, "1x1")),
      "Smile" = list( "CHLA_ODATIS" = c(5, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "7x7"), 
                      "SPM-G_ODATIS" = c(9, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "1x1"),
                      "SST-DAY_ODATIS" = c(9, "1x1"),
                      "SST-NIGHT_ODATIS" = c(10, "1x1")),
      "Sola" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                     "CHLA_SEXTANT" = c(1, "1x1"), 
                     "SPM-G_ODATIS" = c(1, "3x3"),
                     "SPM-G_SEXTANT" = c(1, "1x1"),
                     "SST-DAY_ODATIS" = c(3, "9x9"),
                     "SST-NIGHT_ODATIS" = c(20, "3x3"))
    )  
    
  }
  
  
  if (Version == 2) { # Done with median
    
    best_strategies <- list(
      "Antioche" = list("CHLA_ODATIS" = c(3, "7x7"), 
                        "CHLA_SEXTANT" = c(1, "9x9"),
                        "SPM-G_ODATIS" = c(9, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "5x5"),
                        "SST_ODATIS" = c(9, "1x1")),
      "Astan" = list( "CHLA_ODATIS" = c(9, "5x5"),
                      "CHLA_SEXTANT" = c(1, "9x9"),
                      "SPM-G_ODATIS" = c(9, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "9x9"),
                      "SST_ODATIS" = c(1, "1x1")),
      "Bizeux" = list( "CHLA_ODATIS" = c(1, "3x3"),
                       "CHLA_SEXTANT" = c(1, "3x3"),
                       "SPM-G_ODATIS" = c(9, "5x5"),
                       "SPM-G_SEXTANT" = c(1, "9x9"),
                       "SST_ODATIS" = c(9, "1x1")),
      "Bouee 13" = list( "CHLA_ODATIS" = c(9, "1x1"),
                         "CHLA_SEXTANT" = c(1, "9x9"),
                         "SPM-G_ODATIS" = c(1, "1x1"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST_ODATIS" = c(1, "3x3")),
      "Cézembre" = list( "CHLA_ODATIS" = c(1, "9x9"),
                         "CHLA_SEXTANT" = c(1, "1x1"),
                         "SPM-G_ODATIS" = c(1, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "3x3"),
                         "SST_ODATIS" = c(3, "1x1")),
      "Comprian" = list( "CHLA_ODATIS" = c(3, "1x1"), 
                         "CHLA_SEXTANT" = c(1, "7x7"),
                         "SPM-G_ODATIS" = c(3, "3x3"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST_ODATIS" = c(1, "5x5")),
      "Estacade" = list( "CHLA_ODATIS" = c(9, "9x9"), 
                         "CHLA_SEXTANT" = c(1, "3x3"),
                         "SPM-G_ODATIS" = c(3, "5x5"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST_ODATIS" = c(1, "7x7")),
      "Eyrac" = list( "CHLA_ODATIS" = c(3, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(9, "1x1"),
                      "SPM-G_SEXTANT" = c(1, "None"),
                      "SST_ODATIS" = c(3, "1x1")),
      "Frioul" = list( "CHLA_ODATIS" = c(1, "3x3"), 
                       "CHLA_SEXTANT" = c(1, "1x1"),
                       "SPM-G_ODATIS" = c(1, "3x3"),
                       "SPM-G_SEXTANT" = c(1, "5x5"),
                       "SST_ODATIS" = c(1, "3x3")),
      "Le Buron" = list( "CHLA_ODATIS" = c(3, "9x9"), 
                         "CHLA_SEXTANT" = c(1, "9x9"),
                         "SPM-G_ODATIS" = c(9, "9x9"),
                         "SPM-G_SEXTANT" = c(1, "9x9"),
                         "SST_ODATIS" = c(3, "7x7")),
      "Luc-sur-Mer" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                            "CHLA_SEXTANT" = c(1, "9x9"),
                            "SPM-G_ODATIS" = c(3, "1x1"),
                            "SPM-G_SEXTANT" = c(1, "9x9"),
                            "SST_ODATIS" = c(1, "1x1")),
      "pk 30" = list( "CHLA_ODATIS" = c(1, "9x9"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(1, "5x5"),
                      "SPM-G_SEXTANT" = c(1, "None"),
                      "SST_ODATIS" = c(5, "1x1")),
      "pk 52" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "None"),
                      "SPM-G_ODATIS" = c(1, "3x3"),
                      "SPM-G_SEXTANT" = c(1, "None"),
                      "SST_ODATIS" = c(1, "1x1")),
      "pk 86" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "9x9"),
                      "SPM-G_ODATIS" = c(1, "5x5"),
                      "SPM-G_SEXTANT" = c(1, "1x1"),
                      "SST_ODATIS" = c(1, "9x9")),
      "Point B" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "3x3"),
                        "SPM-G_ODATIS" = c(3, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "1x1"),
                        "SST_ODATIS" = c(1, "1x1")),
      "Point C" = list( "CHLA_ODATIS" = c(9, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "3x3"),
                        "SPM-G_ODATIS" = c(1, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "5x5"),
                        "SST_ODATIS" = c(3, "1x1")),
      "Point L" = list( "CHLA_ODATIS" = c(5, "1x1"), 
                        "CHLA_SEXTANT" = c(1, "9x9"),
                        "SPM-G_ODATIS" = c(9, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "9x9"),
                        "SST_ODATIS" = c(1, "1x1")),
      "Portzic" = list( "CHLA_ODATIS" = c(1, "9x9"), 
                        "CHLA_SEXTANT" = c(1, "3x3"), 
                        "SPM-G_ODATIS" = c(5, "1x1"),
                        "SPM-G_SEXTANT" = c(1, "9x9"),
                        "SST_ODATIS" = c(3, "9x9")),
      "Sete" = list( "CHLA_ODATIS" = c(1, "7x7"), 
                     "CHLA_SEXTANT" = c(1, "9x9"), 
                     "SPM-G_ODATIS" = c(1, "9x9"),
                     "SPM-G_SEXTANT" = c(1, "9x9"),
                     "SST_ODATIS" = c(1, "1x1")),
      "Smile" = list( "CHLA_ODATIS" = c(5, "1x1"), 
                      "CHLA_SEXTANT" = c(1, "7x7"), 
                      "SPM-G_ODATIS" = c(9, "3x3"),
                      "SPM-G_SEXTANT" = c(1, "1x1"),
                      "SST_ODATIS" = c(9, "1x1")),
      "Sola" = list( "CHLA_ODATIS" = c(1, "1x1"), 
                     "CHLA_SEXTANT" = c(1, "1x1"), 
                     "SPM-G_ODATIS" = c(1, "3x3"),
                     "SPM-G_SEXTANT" = c(1, "1x1"),
                     "SST_ODATIS" = c(3, "9x9"))
    )  
    
  }

  return(best_strategies[[Site_name]])
  
} 
  
  

get_sat_data_usign_the_optimal_grid <- function(df, Site_name, Version_optimized_grid, sat_product_full_name) {
  
  best_grid <- best_spatio_temporal_grid(Site_name, Version = Version_optimized_grid)[[paste( sat_product_full_name %>% str_split("_", simplify = TRUE) %>% .[1,1:2], collapse = "_" )]]
  
  time_diff_of_the_site <- best_grid[1] %>% as.numeric()
  grid_size_of_the_site <- best_grid[2]  
  
  if ((grid_size_of_the_site == "None") || (grid_size_of_the_site %>% is.null())) { return() }
  
  name_of_the_sat_product_of_the_site <- paste(sat_product_full_name, grid_size_of_the_site, sep = "_")
  
  if (any(str_detect(names(df), name_of_the_sat_product_of_the_site)) == FALSE) { return() }
  
  df_Site <- df %>%
    filter(Site == Site_name) %>% 
    dplyr::mutate(across(contains(name_of_the_sat_product_of_the_site), .names = "sat_{.col}")) %>% 
    rename_with(~ gsub(paste("sat", name_of_the_sat_product_of_the_site, sep = "_"), "sat", .), contains("sat_")) 
  
  if ( nrow(df_Site) == 0 ) { return() }
  
  if (str_detect(sat_product_full_name, pattern = "SEXTANT_merged") == FALSE) {
    
    index_wrong_time_format = grepl("^([01]?[0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9]", df_Site$sat_UTC_time) == FALSE
    df_Site$sat_UTC_time[ which(index_wrong_time_format) ] <- NA
    
    df_Site <- df_Site %>% 
      mutate( sat_UTC_time = lubridate::hms(sat_UTC_time) ) %>% 
      mutate(Time_difference = abs( lubridate::hour(HEURE) - lubridate::hour(sat_UTC_time) )) %>% 
      mutate( sat_UTC_time = as.character(sat_UTC_time) ) 
    
    index_to_remove <- which( df_Site$Time_difference > time_diff_of_the_site )
    df_Site$sat_mean[index_to_remove] <- NA; df_Site$sat_median[index_to_remove] <- NA; df_Site$sat_value[index_to_remove] <- NA
    
  }
  
  df_to_return <- df_Site %>% select_at(vars(DATE:qCHLA, starts_with("sat_"), contains('Site')))
  
  return(df_to_return)
  
}


unit_color_and_axis_limits_for_the_scatterplot <- function(variable) {
  
  if (variable == "CHLA") {
    unit <- expression(mg~m^-3)
    unit_text = "mg m-3"
    color = "green3"
    axis_limits <- c(10^-2, 10^2)
  } else if (variable %in% c("SPM-G", "SPM-R")) {
    unit <- expression(g~m^-3)
    unit_text = "g m-3"
    axis_limits <- c(10^-2, 10^3)
    if (variable == "SPM-G") {
      color = "red3"
    } else if (variable == "SPM-R") {
      color = "brown"
    }
  } else if (variable %>% str_detect("SST")) {
    unit <- expression('"°C"')
    unit_text = "°C"
    color = "orange"
    axis_limits <- c(3, 30)
  }
  
  to_return <- list("unit" = unit, 
                    "unit_text" = unit_text, 
                    "color" = color, 
                    "axis_limits" = axis_limits)
  
  return(to_return)
  
}


make_the_figures <- function(data, var_to_assess, statistics_outputs) {
  
  color_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
                   "Point C" = mako(n = 1,begin = 0.85,end = 0.85), 
                   "Luc-sur-Mer" = mako(n = 1,begin = 0.90,end = 0.9), 
                   "Smile" = mako(n = 1,begin = 0.95,end = 0.95),
                   
                   "Bizeux" = viridis(n = 1,begin = 1,end = 1), # Bretagne
                   "Le Buron" = viridis(n = 1,begin = 0.925,end = 0.925), 
                   "Cézembre" = viridis(n = 1,begin = 0.85,end = 0.85), 
                   "Estacade" = viridis(n = 1,begin = 0.775,end = 0.775), 
                   "Astan" = viridis(n = 1,begin = 0.7,end = 0.7), 
                   "Portzic" = viridis(n = 1,begin = 0.625,end = 0.625), 
                   
                   "Antioche" = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
                   "pk 86" = plasma(n = 1,begin = 0.1,end = 0.1), 
                   "pk 52" = plasma(n = 1,begin = 0.15,end = 0.15),
                   "pk 30" = plasma(n = 1,begin = 0.2,end = 0.2),
                   "Comprian" = plasma(n = 1,begin = 0.25,end = 0.25), 
                   "Eyrac" = plasma(n = 1,begin = 0.3,end = 0.3), 
                   "Bouee 13" = plasma(n = 1,begin = 0.35,end = 0.35), 
                   
                   "Sola" = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
                   "Sete" = rocket(n = 1,begin = 0.85,end = 0.85), 
                   "Frioul" = rocket(n = 1,begin = 0.90,end = 0.90),
                   "Point B" = rocket(n = 1,begin = 0.95,end = 0.95)
  ) 
  
  unit_color_and_axis_limits <- unit_color_and_axis_limits_for_the_scatterplot(variable = var_to_assess)
  
  identity_line <- data.frame(x = unit_color_and_axis_limits$axis_limits, 
                              y = unit_color_and_axis_limits$axis_limits)
  
  scatterplot <- ggplot(data = data, aes(x = insitu_value, y = sat_value)) + 
    
    geom_point(aes(color = Site_factor), size = 6, show.legend = TRUE) + 
    geom_line(data = identity_line, aes(x = x, y = y, linetype = "Identity line")) +

    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Insitu~measurements~(', unit_color_and_axis_limits$unit, ')'))) +
    
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Satellite~estimates~(', unit_color_and_axis_limits$unit, ')'))) + 
    
    coord_cartesian(xlim = unit_color_and_axis_limits$axis_limits, 
                    ylim = unit_color_and_axis_limits$axis_limits) +
    
    annotate(geom = 'text', x = unit_color_and_axis_limits$axis_limits[1], y = unit_color_and_axis_limits$axis_limits[2], 
             hjust = 0, vjust = 1, color = "black", size = 12.5,
             label = paste('Error = ', round(statistics_outputs$error_Pahlevan, 1), "%\n",
                           'Bias = ', round(statistics_outputs$bias_Pahlevan, 1), " %\n",
                           # 'r² (linearity) = ', round(statistics_outputs$r2_log, 2),"\n",
                           'Slope = ', round(statistics_outputs$slope_log, 2),"\n",
                           'n = ', nrow(data), sep = "")) + 
    
    labs(title = var_to_assess) +
    
    scale_color_manual(name = "", values = color_values, drop = FALSE) +
    
    scale_linetype_manual(values = c("Identity line" = "dashed",
                                     "Linear regression" = "solid"), name = "") +
    
    guides(color=guide_legend(ncol=2, 
                              override.aes = list(size = 5),
                              label.theme = element_text(size = 20))) +
    
    guides(linetype = guide_legend(override.aes = list(color = c("black"),
                                                       shape = c(NA),
                                                       linetype = c("dashed")),
                                   label.theme = element_text(size = 20),
                                   ncol = 2)) +
    
    ggplot_theme() + 
    
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.position = c(.8,.2),
          legend.margin=margin(c(-30,5,5,5)),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(color = unit_color_and_axis_limits$color, face = "bold", size = 35))
  
  if (var_to_assess %>% str_detect("SST")) { 
    scatterplot <- scatterplot +  
      geom_line(data = statistics_outputs$lm_line, aes(x = x, y = y_pred, linetype = "Linear regression"),
                color = "black", linewidth = 1.5) + 
      scale_x_continuous(name = parse(text = paste0('Insitu~measurements~(', unit_color_and_axis_limits$unit, ')'))) +
      scale_y_continuous(name = parse(text = paste0('Satellite~estimates~(', unit_color_and_axis_limits$unit, ')')))
  } else {
    scatterplot <- scatterplot + 
      geom_line(data = statistics_outputs$lm_line_log, aes(x = x, y = y_pred, linetype = "Linear regression"),
                color = "black", linewidth = 1.5) +
      annotation_logticks()
  }
  
  scatterplot_with_side_hist <- ggMarginal(scatterplot, type = "histogram", fill = "white")
  
  bar_plot_freq_per_year <- data %>% 
    dplyr::count(Year) %>% 
    ggplot(aes(x = Year)) + 
    geom_col(aes(y = n), fill = "white", color = unit_color_and_axis_limits$color, linewidth = 1.5) + 
    scale_x_continuous(breaks = 1998:2030, name = "", labels = function(x) substring(x, 3, 4)) +
    scale_y_continuous(expand = c(0,0), name = "n per year") +
    ggplot_theme() 
  
  bar_plot_freq_per_month <- data %>% 
    dplyr::count(Month) %>% 
    ggplot(aes(x = Month)) + 
    geom_col(aes(y = n), fill = "white", color = unit_color_and_axis_limits$color, linewidth = 2) + 
    scale_x_continuous(breaks = 1:12, name = "", labels = function(x) month.abb[x]) +
    scale_y_continuous(expand = c(0,0), name = "n per month") +
    coord_cartesian(xlim = c(1,12)) +
    ggplot_theme()
  
  return(list("scatterplot" = scatterplot, "scatterplot_with_side_histograms" = scatterplot_with_side_hist,
              "bar_plot_freq_per_year" = bar_plot_freq_per_year,
              "bar_plot_freq_per_month" = bar_plot_freq_per_month))
  
}



color_by_value_gt <- function(df, group_columns, columns_to_color) {
  
  # Normalize the values within each group defined by multiple columns
  df_with_color_values <- df %>%
    group_by(across(all_of(group_columns))) %>%
    mutate(across(all_of(columns_to_color), ~ scales::rescale(., to = c(0, 1), na.rm = TRUE))) %>%
    ungroup()
  
  # Create the gt table
  gt_table <- gt(df_with_color_values)
  
  # Apply continuous color scale to both "Bagdad" and "Liban" columns
  for (col in columns_to_color) {
    gt_table <- gt_table %>%
      data_color(
        columns = col,
        colors = scales::col_numeric(
          palette = viridis(100),
          domain = NULL  # Set to NULL as values are already rescaled
        )
      )
  }
  
  gt_table$`_data`[,columns_to_color] <- df[,columns_to_color]
  
  return(gt_table)
}


load_SOMLIT_data <- function(path_to_SOMLIT_data_with_MU_results, stations_to_remove = c()) {

  path_to_SOMLIT_data_with_MU_results %>% 
    list.files(pattern = "Somlit_with_MU", full.names = TRUE) %>% 
    llply(function(x) {read_csv(x, show_col_types = FALSE)[,-1]}) %>% 
    purrr::reduce(left_join, by = c("DATE", "HEURE_x", "COEF_MAREE", "MAREE", "PROF_TEXT", "PROF_NUM", "Site","Latitude", "Longitude", 
                                    "T",   "qT",  "S", "qS",  "COP", "qCOP", "MES", "qMES","CHLA", "qCHLA", "HEURE_y")) %>% 
    filter(PROF_TEXT == "S", Site %in% stations_to_remove == FALSE) %>% 
    mutate(Site_factor = factor(Site, levels = c("Point L", "Point C", # Manche 
                                                 "Luc-sur-Mer", "Smile",
                                                 "Bizeux", "Le Buron", "Cézembre", # Bretagne
                                                 "Estacade", "Astan",
                                                 "Portzic",
                                                 "Antioche", # Golfe de Gascogne
                                                 "pk 86", "pk 52", "pk 30",
                                                 "Comprian", "Eyrac", "Bouee 13", 
                                                 "Sola", "Sete", "Frioul", "Point B" # Golfe du Lion
    ))) %>% 
    dplyr::rename(HEURE = HEURE_x, SPM = MES, qSPM = qMES, SST = 'T', qSST = qT) %>% 
    rename_at(vars(starts_with('SST-NS_DAY')), function(x) {str_replace(x, 'SST-NS_DAY', 'SST-DAY')}) %>% 
    dplyr::select(- HEURE_y) %>% 
    distinct()  
  
}
