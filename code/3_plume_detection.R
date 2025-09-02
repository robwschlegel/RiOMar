list_of_packages <- c("tidyverse", "plyr", "ggpubr", "viridis")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, require, character.only = TRUE)

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

save_file_as_csv <- function(data, file_name) {
  path <- file_name %>% str_split(pattern = "/") %>% .[[1]] %>% head(-1) %>% paste(collapse = "/")
  if (dir.exists(path) == FALSE) {dir.create(path = path, recursive = TRUE)}
  if (grepl(".csv", file_name) == FALSE) {file_name <- paste(file_name, ".csv", sep = "")}
  data.table::fwrite(data %>% as.data.frame(), file = file_name, row.names = TRUE, col.names = TRUE)
}

plot_function <- function(data, y_variable) {
  
  data %>% 
    ggplot(aes(x = date, y = .data[[y_variable]], group = Satellite_sensor, color = Satellite_sensor)) +
    geom_line(linewidth = 1.5, na.rm = TRUE) + geom_point(size = 2, na.rm = TRUE) + 
    labs(y = case_when(grepl("^SPM_threshold", y_variable) ~ y_variable %>% str_replace_all("_", " "),
                       y_variable == "area_of_the_plume_mask_in_km2" ~ "Plume area (km²)")) + 
    ggplot_theme() 
  
}

make_the_plot_from_the_df <- function(data) {
  
  SPM_thresholds_columns <- data %>% names() %>% grep(pattern = "^SPM_threshold", value = TRUE)
  
  time_series_plot <- ggarrange(plotlist = c("area_of_the_plume_mask_in_km2", SPM_thresholds_columns) %>% llply(plot_function, data = data),
                                ncol = 1, align = 'v', common.legend = TRUE)
  
  unique(data$path_to_file) %>% 
    l_ply( save_plot_as_png, name = paste('Time_series_of', tolower(data$Var10[1]), 'plume_area_and_SPM_threshold', sep = "_"), 
                                          width = 28, height = 16, plot = time_series_plot )
  
}

work_dir = '/home/terrats/Desktop/RIOMAR/'
Zone = 'BAY_OF_SEINE'
Data_source = 'SEXTANT'
Satellite_sensor = c('merged', 'modis')
atmospheric_correction = 'Standard'
Time_resolution = 'DAILY'
Years = 2018:2021

plot_time_series_of_plume_area_and_thresholds <- function(work_dir, Zone, Data_source, Satellite_sensor, atmospheric_correction, Time_resolution, Years) {
  
  cases_to_process <- expand.grid(work_dir, 'RESULTS', Zone, Data_source, "Satellite_sensor", atmospheric_correction, 
                                  Years %>% as.character(), 'PLUME_DETECTION', 'SPM-G', Time_resolution)
  
  ts_data <- cases_to_process %>% 
        adply(1, function(row) {

          file_path <- row %>% as_vector() %>% paste(collapse = "/")
          
          files_to_load <- Satellite_sensor %>% 
                            laply(function(x) file_path %>% str_replace('Satellite_sensor', x)) %>% 
                            list.files( pattern = ".csv", full.names = TRUE, recursive = TRUE ) 
          
          if (files_to_load %>% is_empty()) {return()}
          
          data <- files_to_load %>%  
            ldply( function(file_to_load) {
              read.csv(file_to_load, stringsAsFactors = FALSE, colClasses = c("date" = "character")) %>% 
                mutate(date = date %>% as.Date(format = "%Y-%m-%d"),
                       Satellite_sensor = file_to_load %>% 
                                            str_extract(pattern = paste(row$Var4, '/[A-Za-z]*/', sep = '')) %>% 
                                            str_remove_all(paste('[/', row$Var4 ,']', sep = '|')),
                       path_to_file = file_to_load %>% str_remove('/[A-Za-z]*.csv'))  
              
            }  ) %>% bind_cols(row)
          
          make_the_plot_from_the_df(data)       
          
          return(data)
            
        })
  
  global_cases_to_process <- expand.grid(work_dir, 'RESULTS', Zone, Data_source, "*", atmospheric_correction, 
                                         "*", 'PLUME_DETECTION', 'SPM-G', Time_resolution)
  
  global_cases_to_process %>% 
    
    a_ply(1, function(row) {
    
      file_path <- row %>% as_vector() %>% paste(collapse = "/")
      
      ts_data_of_the_case <- ts_data %>% 
                              filter(Var3 == row$Var3, Var4 == row$Var4, Var6 == row$Var6, Var10 == row$Var10) %>% 
                              mutate(path_to_file = file_path %>% str_remove_all('[*][A-Za-z0-9/_-]*'))
      
      save_file_as_csv(ts_data_of_the_case, 
                       file.path(ts_data_of_the_case$path_to_file[1], paste('Time_series_of', 
                                                                         tolower(ts_data_of_the_case$Var10[1]), 
                                                                         'plume_area_and_SPM_threshold', sep = "_")))

      make_the_plot_from_the_df(data = ts_data_of_the_case)                   
    
  })
  
}
