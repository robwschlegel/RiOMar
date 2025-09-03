list_of_packages <- c("plyr", "tidyverse", "scales", "viridis", "ggExtra", "Metrics", "MASS", "gt")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, require, character.only = TRUE)

#### Main function ####

# satellite_median = c(4,0.549999987706542,3,2,0.5999999865889549,1.1099999751895666,1,0.9699999783188104,3,5)
# satellite_n = c(0, 1, 0, 0, 1, 1, 0, 1, 0, 0)
# satellite_sd = c(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
# satellite_times = c(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
# insitu_variable = c('SPM', 'SPM', 'SPM', 'SPM', 'SPM', 'SPM', 'TURB', 'TURB', 'TURB', 'TURB')
# insitu_value = c(0.034, 0.331, 0.787, 0.148, 0.297, 0.217, 0.342, 0.342, 0.182, 0.354)
# insitu_time = c('11:15:00','10:10:00','11:40:00','14:00:00','10:30:00','11:05:00','11:35:00','10:45:00','10:10:00','11:20:00')
# insitu_Data_source = c('REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY')
# min_n = 1
# max_CV = NaN
# max_hour_diff = NaN
# grid_size = 1
# date = c('2018-01-02','2018-01-02','2018-01-08','2018-01-16','2018-01-16','2018-01-15','2018-01-15','2018-01-29','2018-01-29','2018-01-29')
# satellite_source = 'SEXTANT'
# satellite_sensor = 'modis'
# satellite_atm_corr = 'Standard'
# satellite_algorithm = 'suspended_matters'
# site_name = c('Marseillan (a)','SÃ¨te mer','Parc Leucate 2','Marseillan (a)','SÃ¨te mer','Barcares','Parc Leucate 2','Barcares','Marseillan (a)','Parc Leucate 2')
# region_name = c('Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion',
#                 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion')
# zones = c('GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY')
# where_to_save_MU_results <- '/home/terrats/Desktop/RIOMAR/TEST/MATCH_UP_DATA/'

Save_validation_scatterplots_and_stats <- function(satellite_median, satellite_n, satellite_sd, satellite_times, 
                                                   insitu_variable, insitu_value, insitu_time, insitu_Data_source, site_name, region_name, zones,
                          min_n, max_CV, max_hour_diff, grid_size, date,
                          satellite_source, satellite_sensor, satellite_atm_corr, satellite_algorithm,
                          where_to_save_MU_results) {
  
  args <- convert_nan(satellite_median, satellite_n, satellite_sd, satellite_times, 
                      insitu_variable, insitu_value, insitu_time, insitu_Data_source, site_name, region_name, zones,
                      min_n, max_CV, max_hour_diff, grid_size, date,
                      satellite_source, satellite_sensor, satellite_atm_corr, satellite_algorithm,
                      where_to_save_MU_results) %>% 
          setNames(c("satellite_median", "satellite_n", "satellite_sd", "satellite_times", 
                     "insitu_variable", "insitu_value", "insitu_time", "insitu_Data_source", "site_name", "region_name", "zones",
                     "min_n", "max_CV", "max_hour_diff", "grid_size", "date",
                     "satellite_source", "satellite_sensor", "satellite_atm_corr", "satellite_algorithm",
                     "where_to_save_MU_results"))
  
  # args$region_name <- ifelse(args$region_name == 'Manche orientale - Mer du Nord', 'Baie de Seine',
  #                     ifelse(args$region_name == 'Baie de Seine', 'Baie de Seine',
  #                     ifelse(args$region_name == 'Manche occidentale', 'Bretagne Sud',
  #                     ifelse(args$region_name == 'Bretagne Sud', 'Bretagne Sud',
  #                     ifelse(args$region_name == 'Pays de la Loire - Pertuis', 'Golfe de Gascogne',
  #                     ifelse(args$region_name == 'Sud Golfe de Gascogne', 'Golfe de Gascogne',
  #                     ifelse(args$region_name == 'Golfe du Lion', 'Golfe du Lion',
  #                     ifelse(args$region_name == 'Mer ligurienne - Corse', 'Golfe du Lion',
  #                     "Not found"))))))))
  
  unique(args$insitu_variable) %>% l_ply(function(insitu_var) {
    
    plot_title = paste("In situ", insitu_var, "Vs.", "Satellite",  
                       paste(args$satellite_algorithm, args$satellite_source, args$satellite_sensor, args$satellite_atm_corr, sep = " / "))
    plot_subtitle = paste("Grid size = ", args$grid_size, " x ", args$grid_size," / Variation Coefficient ≤ ", args$max_CV, "% / n ≥ ", args$min_n, ' / time difference ≤ ', args$max_hour_diff, "h", sep = "")
    
    satellite_CV = 100 * (args$satellite_sd / args$satellite_median)
    Year = lubridate::year(args$date)
    Month = lubridate::month(args$date)
    
    mask_CV <- case_when(args$max_CV %>% is.finite() ~ satellite_CV <= args$max_CV, .default = TRUE)
    mask_n <- case_when(args$min_n %>% is.finite() ~ args$satellite_n >= args$min_n, .default = TRUE)
    mask_hour_diff <- case_when(args$max_hour_diff %>% is.finite() ~ 
                                  abs( as.numeric(hms(args$satellite_times) - hms(args$insitu_time), unit = "hours") ) <= args$max_hour_diff,
                                .default = TRUE)
    mask_positive_values <- (args$satellite_median > 0) & (args$insitu_value > 0)
    mask_insitu_parameter <- insitu_variable == insitu_var
    
    final_mask = mask_CV & mask_n & mask_hour_diff & mask_positive_values & mask_insitu_parameter

    if ( final_mask %>% is.null()) { return() }
    final_mask_index <- which(final_mask)
    
    regions_to_process <- c('Global', args$region_name[final_mask_index]) %>% unique()
    
    statistics_values <- regions_to_process %>% llply(function(region_name) {
      
      if (region_name == 'Global') {
        mask <- final_mask
      } else {
        mask <- final_mask & (args$region_name == region_name)
      }
      
      compute_stats(sat_values = args$satellite_median[which(mask)], insitu_values = args$insitu_value[which(mask)])  
      
    }, .inform = TRUE) %>% 
      setNames(regions_to_process)
    
    Figures <- make_the_figures(args$insitu_value[final_mask_index], 
                                args$satellite_median[final_mask_index], 
                                Year[final_mask_index], Month[final_mask_index],
                                args$region_name[final_mask_index],
                                args$insitu_Data_source[final_mask_index],
                                plot_title, plot_subtitle, 
                                args$satellite_algorithm, 
                                statistics_values[['Global']])
    
    # if (variable %>% str_detect("SST")) { 
    save_plot_as_png(Figures$scatterplot_with_side_histograms, 
                     name = paste("Insitu", insitu_var, "Vs", "Satellite", args$satellite_algorithm, args$satellite_source, args$satellite_sensor, args$satellite_atm_corr, sep = "_"), 
                     path = file.path(args$where_to_save_MU_results, "SCATTERPLOTS", "Per_sat_product"), 
                     width = 22, height = 16)
    
    save_file_as_csv(statistics_values %>% convert_list_to_df() %>% dplyr::select_at(vars(-ends_with("_line"), -contains("_line_"))) %>% 
                       mutate("insitu_var" = insitu_var, 'satellite_var' = args$satellite_algorithm), 
                     file.path(args$where_to_save_MU_results, "STATISTICS", "Per_sat_product",
                               paste("Insitu", insitu_var, "Vs", "Satellite", 
                                     args$satellite_algorithm, args$satellite_source, 
                                     args$satellite_sensor, args$satellite_atm_corr, sep = "_")))
    
  }, .inform = TRUE)
  
}




# satellite_median = c(4,0.549999987706542,3,2,0.5999999865889549,1.1099999751895666,1,0.9699999783188104,3,5)
# satellite_n = c(0, 1, 0, 0, 1, 1, 0, 1, 0, 0)
# satellite_sd = c(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
# satellite_times = c(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
# insitu_variable = c('SPM', 'SPM', 'SPM', 'SPM', 'SPM', 'SPM', 'TURB', 'TURB', 'TURB', 'TURB')
# insitu_value = c(0.034, 0.331, 0.787, 0.148, 0.297, 0.217, 0.342, 0.342, 0.182, 0.354)
# insitu_time = c('11:15:00','10:10:00','11:40:00','14:00:00','10:30:00','11:05:00','11:35:00','10:45:00','10:10:00','11:20:00')
# insitu_Data_source = c('REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY','REPHY')
# min_n = 1
# max_CV = NaN
# max_hour_diff = NaN
# grid_size = 1
# date = c('2018-01-02','2018-01-02','2018-01-08','2018-01-16','2018-01-16','2018-01-15','2018-01-15','2018-01-29','2018-01-29','2018-01-29')
# satellite_source = 'SEXTANT'
# satellite_sensor = 'modis'
# satellite_atm_corr = 'Standard'
# satellite_algorithm = 'suspended_matters'
# site_name = c('Marseillan (a)','SÃ¨te mer','Parc Leucate 2','Marseillan (a)','SÃ¨te mer','Barcares','Parc Leucate 2','Barcares','Marseillan (a)','Parc Leucate 2')
# region_name = c('Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion',
#                 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion', 'Golfe du Lion')
# zones = c('GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY')
where_to_save_MU_results <- '/home/terrats/Desktop/RIOMAR/TEST/MATCH_UP_DATA/GULF_OF_LION_&_BAY_OF_SEINE_&_BAY_OF_BISCAY_&_SOUTHERN_BRITTANY'

Summarize_statistics_in_a_table <- function(where_to_save_MU_results) {
  
  ### Redo the match-up at first ! Because the zones are not good. 
  stat_files <- file.path(where_to_save_MU_results, "STATISTICS", "Per_sat_product") %>% 
                  list.files(pattern = "*.csv", full.names = TRUE) %>% 
                  ldply(read.csv) %>% 
                  dplyr::select(Region, Bias = bias_Pahlevan, Error = error_Pahlevan, n, insitu_var) %>% 
                  mutate_at(vars(-Region, -insitu_var), ~ round(., 1)) %>% 
                  mutate(insitu_var = paste("When compared with in situ", ifelse(insitu_var == 'TURB', 'Turbidity', 'SPM')))

  desired_colnames <- names(stat_files) %>% str_replace_all("Bias", "Bias (%)") %>% str_replace_all("Error", "Error (%)") %>% str_remove_all("insitu_var")
  names(desired_colnames) <- names(stat_files)
  
  the_Table <- stat_files %>% 
    rowwise %>% mutate(Region = paste0("**", Region, "**")) %>% 
    gt(rowname_col = 'Region', groupname_col = 'insitu_var', process_md = TRUE) %>% 
    cols_label(.list = desired_colnames) %>% 
    tab_spanner(label = md('**Metrics**'), columns = c("Error", 'Bias', "n")) %>% 
    tab_header(title = 'Performances of satellite SPM', subtitle = 'Compared with SPM and Turbidity in situ measurements.') %>% 
    sub_missing(missing_text = "-") %>% 
    
    tab_options(data_row.padding = px(2),
                summary_row.padding = px(3), # A bit more padding for summaries
                row_group.padding = px(4)    # And even more for our groups
    ) %>% 
    opt_stylize(style = 6, color = 'gray') %>%
    tab_style(
      style = cell_text(align = "center"),
      locations = cells_column_labels()
    )
  
  the_Table %>% gtsave(file.path(where_to_save_MU_results, "STATISTICS", "Per_sat_product", "Table.html"))
  
}



#### Utils ####
convert_nan <- function(...) {
  list(...) %>% llply(function(x) ifelse(x == "nan", NaN, x))
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

ggplot_theme <-   function() {
  theme(
        text = element_text(size=35, colour = "black"), #25
        plot.title = element_text(hjust = 0.5, size = 55),
        plot.subtitle = element_text(hjust = 0.5, size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x=element_text(angle=0),
        axis.ticks.length=unit(.25, "cm"))}

compute_stats <- function(sat_values, insitu_values, axis_limits = c(NA, NA)) {

  if (length(sat_values) < 3) {return(list())}
  
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


unit_color_and_axis_limits_for_the_scatterplot <- function(variable) {
  
  if ( grepl( 'CHL', variable) ) {
    unit <- expression(mg~m^-3)
    unit_text = "mg m-3"
    color = "green3"
    axis_limits <- c(10^-2, 10^2)
  } else if ( grepl( 'SPM|SPIM|suspended_matters', variable) ) {
    unit <- expression(g~m^-3)
    unit_text = "g m-3"
    axis_limits <- c(10^-2, 10^3)
    if (variable %in% c("SPM-G", "SPIM", 'SPM', 'suspended_matters')) {
      color = "red3"
    } else if (variable == "SPM-R") {
      color = "brown"
    }
  } else if ( grepl( 'SST', variable) ) {
    unit <- expression('"°C"')
    unit_text = "°C"
    color = "orange"
    axis_limits <- c(3, 30)
  } else {
      stop(paste( "The function unit_color_and_axis_limits_for_the_scatterplot is not suited for", variable) )
  }
  
  to_return <- list("unit" = unit, 
                    "unit_text" = unit_text, 
                    "color" = color, 
                    "axis_limits" = axis_limits)
  
  return(to_return)
  
}


colors_of_stations <- function() {
  
  # color_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
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
  
  # color_values = c('Manche orientale - Mer du Nord' = mako(n = 1,begin = 0.8,end = 0.8), # Manche
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
  
  color_values = c('BAY OF SEINE' = mako(n = 1,begin = 0.8,end = 0.8),
                   'SOUTHERN BRITTANY' = viridis(n = 1,begin = 0.9,end = 0.9),
                   'GULF OF BISCAY' = plasma(n = 1,begin = 0.05,end = 0.05),
                   'GULF OF LION' = rocket(n = 1,begin = 0.80,end = 0.80))
  
  return(color_values)

}

save_file_as_csv <- function(data, file_name) {
  path <- file_name %>% str_split(pattern = "/") %>% .[[1]] %>% head(-1) %>% paste(collapse = "/")
  if (dir.exists(path) == FALSE) {dir.create(path = path, recursive = TRUE)}
  if (grepl(".csv", file_name) == FALSE) {file_name <- paste(file_name, ".csv", sep = "")}
  data.table::fwrite(data %>% as.data.frame(), file = file_name, row.names = TRUE, col.names = TRUE)
}

convert_list_to_df <- function(lst) {
  # Convert list of lists into a dataframe
  
  lst <- lst %>% purrr::discard(~ length(.) == 0)
  
  df <- do.call(rbind, lapply(lst, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  
  # Add the names of the outer list as the first column
  df <- cbind(Region = names(lst), df)
  
  return(df)
}


make_the_figures <- function(insitu_value, satellite_median, Year, Month,
                             site_name, insitu_Data_source, plot_title, plot_subtitle, satellite_algorithm, statistics_values) {
  
  unit_color_and_axis_limits <- unit_color_and_axis_limits_for_the_scatterplot(variable = satellite_algorithm)
  
  identity_line <- data.frame(x = unit_color_and_axis_limits$axis_limits, 
                              y = unit_color_and_axis_limits$axis_limits)
  
  MU_data = data.frame(insitu_value, satellite_median, site_name, Year, Month, insitu_Data_source) %>% 
              mutate_at(vars(site_name, insitu_Data_source), . %>% as.factor())
  
  color_values <- colors_of_stations()
  
  if (str_detect(satellite_algorithm, 'SST')) {
    Error_value <- statistics_values$error_Pahlevan_linear
    Bias_value <- statistics_values$bias_Pahlevan_linear
    Slope_value <- statistics_values$slope
  } else {
    Error_value <- statistics_values$error_Pahlevan
    Bias_value <- statistics_values$bias_Pahlevan
    Slope_value <- statistics_values$slope_log
  }

  scatterplot <- ggplot(data = MU_data, aes(x = insitu_value, y = satellite_median)) + 
    
    geom_point(aes(color = site_name, shape = insitu_Data_source), size = 6, show.legend = TRUE) + 
    geom_line(data = identity_line, aes(x = x, y = y, linetype = "Identity line")) +
    
    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('In~situ~measurements~(', unit_color_and_axis_limits$unit, ')'))) +
    
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Satellite~estimates~(', unit_color_and_axis_limits$unit, ')'))) + 
    
    coord_cartesian(xlim = unit_color_and_axis_limits$axis_limits, 
                    ylim = unit_color_and_axis_limits$axis_limits) +
    
    annotate(geom = 'text', x = unit_color_and_axis_limits$axis_limits[1], y = unit_color_and_axis_limits$axis_limits[2], 
             hjust = 0, vjust = 1, color = "black", size = 12.5,
             label = paste('Error = ', round( ifelse(Error_value %>% is.numeric(),Error_value, NA), 1), "%\n",
                           'Bias = ', round(ifelse(Bias_value %>% is.numeric(),Bias_value, NA), 1), " %\n",
                           # 'r² (linearity) = ', round(statistics_values$r2_log, 2),"\n",
                           # 'Slope = ', round(ifelse(Slope_value %>% is.numeric(),Slope_value, NA), 2),"\n",
                           'n = ', nrow(MU_data), sep = "")) + 
    
    labs(title = plot_title, subtitle = plot_subtitle, shape = "") +
    
    scale_color_manual(name = "", values = color_values, drop = FALSE) +
    
    scale_linetype_manual(values = c("Identity line" = "dashed",
                                     "Linear regression" = "solid"), name = "") +
    
    guides(color=guide_legend(ncol=1, override.aes = list(size = 10),
                              order = 1),
           shape=guide_legend(ncol=1, override.aes = list(size = 10),
                              order = 3),
           linetype = guide_legend(override.aes = list(color = c("black"),
                                                       shape = c(NA),
                                                       linetype = c("dashed")),
                                   ncol = 2), order = 2) +
    
    ggplot_theme() + 
    
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.position = c(.8,.2),
          legend.text = element_text(size = 30),
          legend.margin=margin(c(-30,5,5,5)),
          plot.subtitle = element_text(hjust = 0.5, color = "black", face = "bold.italic"),
          plot.title = element_text(color = unit_color_and_axis_limits$color, face = "bold", size = 35))
  
  if (satellite_algorithm %>% str_detect("SST")) { 
    scatterplot <- scatterplot +  
      # geom_line(data = statistics_values$lm_line, aes(x = x, y = y_pred, linetype = "Linear regression"),
      #           color = "black", linewidth = 1.5) + 
      scale_x_continuous(name = parse(text = paste0('In~situ~measurements~(', unit_color_and_axis_limits$unit, ')'))) +
      scale_y_continuous(name = parse(text = paste0('Satellite~estimates~(', unit_color_and_axis_limits$unit, ')')))
  } else {
    scatterplot <- scatterplot + 
      # geom_line(data = statistics_values$lm_line_log, aes(x = x, y = y_pred, linetype = "Linear regression"),
      #           color = "black", linewidth = 1.5) +
      annotation_logticks()
  }
  
  scatterplot_with_side_hist <- ggMarginal(scatterplot, type = "histogram", fill = "white")
  
  bar_plot_freq_per_year <- MU_data %>% 
    dplyr::count(Year) %>% 
    ggplot(aes(x = Year)) + 
    geom_col(aes(y = n), fill = "white", color = unit_color_and_axis_limits$color, linewidth = 1.5) + 
    scale_x_continuous(breaks = 1998:2030, name = "", labels = function(x) substring(x, 3, 4)) +
    scale_y_continuous(expand = c(0,0), name = "n per year") +
    ggplot_theme() 
  
  bar_plot_freq_per_month <- MU_data %>% 
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
