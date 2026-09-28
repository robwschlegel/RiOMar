# Plotting ----------------------------------------------------------------

# Create pretty plot labels from zone values
make_pretty_title <- function(df){
  df <- df |> 
    mutate(plot_title = case_when(zone == "BAY_OF_SEINE" ~ "Bay of Seine",
                                  zone == "SOUTHERN_BRITTANY" ~ "Southern Brittany",
                                  zone == "BAY_OF_BISCAY" ~ "Gironde shelf",
                                  zone == "GULF_OF_LION" ~ "Rhône shelf"), .after = "zone") |>
    mutate(plot_title = factor(plot_title, levels = c("Bay of Seine", "Southern Brittany",
                                                      "Gironde shelf", "Rhône shelf")))
  return(df)
}

# Scale one value to another for tidier double-y-axis plots
sec_axis_adjustement_factors <- function(var_to_scale, var_ref){
  
  index_to_keep <- which(is.finite(var_ref))
  var_ref <- var_ref[index_to_keep]
  
  index_to_keep <- which(is.finite(var_to_scale))
  var_to_scale <- var_to_scale[index_to_keep]
  
  max_var_to_scale <- max(var_to_scale, na.rm = T) 
  min_var_to_scale <- min(var_to_scale, na.rm = T) 
  max_var_ref <- max(var_ref, na.rm = T) 
  min_var_ref <- min(var_ref, na.rm = T) 
  
  diff_to_scale <- max_var_to_scale - min_var_to_scale
  diff_to_scale <- ifelse(diff_to_scale == 0, 1 , diff_to_scale)
  diff_ref <- max_var_ref - min_var_ref
  diff <- diff_ref / diff_to_scale
  
  adjust <- (max_var_ref - max_var_to_scale*diff) 
  
  return(data.frame(diff = diff, adjust = adjust, operation = "scaled var = (var_to_scale * diff) + adjust",
                    trans_axis_operation = "var_to_scale = {scaled_var - adjust} / diff)"))
}

# Consistent theme for project
ggplot_theme <-   function(){
  theme(text = element_text(size = 35, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 55),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(angle = 0),
        axis.ticks.length = unit(.25, "cm"))
}

# Convenience wrapper for saving png files
save_plot_as_png <- function(plot, name = c(), width = 14, height = 8.27, path, res = 150){
  
  graphics.off()
  
  if (name |> length() == 1) {
    if (dir.exists(path) == FALSE) {dir.create(path, recursive = TRUE)}
    path <- file.path(path, paste(name, ".png", sep = ""))
  } else {
    path <- paste(path, ".png", sep = "")
  }
  
  if (grepl(pattern = ".png.png", path)) {path <- path |> gsub(pattern = ".png.png", replacement = ".png", x = .)}
  
  png(path, width = width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
  
}

# Comparison plots
comparison_plot <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2){
  
  if(grepl("seas", var_1)){
    df_sub <- df[,c("plot_title", "month", var_1, var_2)]
    colnames(df_sub) <- c("plot_title", "month", "var_1", "var_2")
  } else {
    df_sub <- df[,c("plot_title", "date", var_1, var_2)]
    colnames(df_sub) <- c("plot_title", "date", "var_1", "var_2")
  }
  
  # Plot base
  if(grepl("seas", var_1)){
    
    # Scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_sub$var_2, df_sub$var_1)
    df_scaling <- summarise(df_sub, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_sub, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    
    # Get range for ribbon plot
    df_scale_sub <- df_scale |> 
      summarise(var_1_min = min(var_1, na.rm = TRUE),
                var_1_mean = mean(var_1, na.rm = TRUE),
                var_1_max = max(var_1, na.rm = TRUE),
                var_2_min = min(var_2_scaled, na.rm = TRUE),
                var_2_mean = mean(var_2_scaled, na.rm = TRUE),
                var_2_max = max(var_2_scaled, na.rm = TRUE), .by = c("plot_title", "month")) |> 
      mutate(month_int = as.integer(month))
    
    # Plot them
    pl_base <- ggplot(data = df_scale_sub, aes(x = month_int)) + 
      # Var 1
      geom_ribbon(aes(ymin = var_1_min, ymax = var_1_max), fill = colour_1, alpha = 0.3) +
      geom_path(aes(y = var_1_mean), color = colour_1, linewidth = 2) +
      # Var 2
      geom_ribbon(aes(ymin = var_2_min, ymax = var_2_max), fill = colour_2, alpha = 0.3) +
      geom_path(aes(y = var_2_mean), color = colour_2, linewidth = 2) +
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      scale_x_continuous(expand = c(0, 0), breaks = 1:12, labels = month.abb)
  } else {
    
    # Perform rolling mean
    df_roll_mean <- df_sub |>
      mutate(date = date - lubridate::days(lubridate::wday(date)-1)) |>
      # mutate(date = round_date(date, unit = "months")) |>
      filter(date >= min(df$date)) |>
      group_by(plot_title, date) |>
      summarise(var_1 = mean(var_1, na.rm = TRUE),
                var_2 = mean(var_2, na.rm = TRUE), .groups = "keep") |>
      group_by(plot_title) |>
      mutate(var_1 = roll_mean(var_1, n = 48, fill = NA, align = "center"),
             var_2 = roll_mean(var_2, n = 48, fill = NA, align = "center")) |>
      ungroup()
    
    # Then get the scaling factor
    scaling_factor <- sec_axis_adjustement_factors(df_roll_mean$var_2, df_roll_mean$var_1)
    df_scaling <- summarise(df_roll_mean, sec_axis_adjustement_factors(var_2, var_1), .by = plot_title)
    df_scale <- left_join(df_roll_mean, df_scaling, by = "plot_title") |>
      mutate(var_2_scaled = var_2 * diff + adjust, .after = "var_2")
    unique_years <- df_scale$date |> year() |> unique()
    
    pl_base <- ggplot(data = df_scale) +
      # Var 1 data
      geom_point(aes(x = date, y = var_1), color = colour_1) +
      geom_path(aes(x = date, y = var_1), color = colour_1) +
      # geom_smooth(aes(x = date, y = var_1), method = "lm", se = FALSE, color = colour_1) +
      # Var 2 data
      geom_point(aes(x = date, y = var_2_scaled), color = colour_2) +
      geom_path(aes(x = date, y = var_2_scaled), color = colour_2) +
      # geom_smooth(aes(x = date, y = var_2_scaled), method = "lm", se = FALSE, color = colour_2) +
      # Facet
      facet_wrap(~plot_title, ncol = 1, scales = "free_y") +
      # X-axis labels
      scale_x_date(name = "", expand = c(0, 0),
                   breaks = paste(unique_years, "01-01", sep = "-") %>% as.Date(),
                   labels = unique_years %>% str_extract_all('[0-9][0-9]$') %>% unlist()) 
  }
  
  # Finish up the comparison plot
  pl_comp <- pl_base +
    # Y-axis labels
    scale_y_continuous(name = label_1,
                       sec.axis = sec_axis(transform = ~ {. - scaling_factor$adjust} / scaling_factor$diff, 
                                           name = label_2)) +
    labs( x = NULL) +
    # Extra bits
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y.left = element_text(color = colour_1),
          axis.ticks.y.left = element_line(color = colour_1),
          axis.line.y.left = element_line(color = colour_1),
          axis.title.y.left = element_text(color = colour_1, margin = unit(c(0, 7.5, 0, 0), "mm")),
          axis.text.y.right = element_text(color = colour_2),
          axis.ticks.y.right = element_line(color = colour_2),
          axis.line.y.right = element_line(color = colour_2),
          axis.title.y.right = element_text(color = colour_2, margin = unit(c(0, 0, 0, 7.5), "mm")),
          panel.border = element_rect(linetype = "solid", fill = NA))
  return(pl_comp)
}

# Convenience wrapper to run and save comparison plots
# NB: this is hard coded to work with four plots
comparison_plot_save <- function(df, var_1, var_2, colour_1, colour_2, label_1, label_2, file_stub){
  comp_list <- plyr::dlply(df, c("zone"), comparison_plot, var_1 = var_1, var_2 = var_2, 
                           colour_1 = colour_1, colour_2 = colour_2, label_1 = label_1, label_2 = label_2)
  comp_fig <- comp_list[[2]] + comp_list[[4]] + comp_list[[1]] + comp_list[[3]] + patchwork::plot_layout(ncol = 1, axes = "collect")
  ggsave(filename = paste0("figures/",file_stub,".png"), plot = comp_fig, width = 24, height = 24, dpi = 300)
}

# It does what it says on the tin
var_labels <- function(var_name){
  
  if(grepl("CHL|RRS", var_name)){
    var_colour <- "green4"
    if(grepl("CHL", var_name)){
      axis_limits <- c(10^-2, 10^2)
      unit <- expression(mg~m^-3)
    } else {
      unit <- expression(sr-1)
      axis_limits <- c(10^-5, 10^-1)
    }
  } else if(grepl("TEMP|SST", var_name)) {
    unit <- expression('"°C"')
    axis_limits <- c(3, 30)
    var_colour <- "orange4"
    # NB: Must be run after SST because the T for turbidity is an issue
  } else if(grepl("SPM|TUR|T|POC|CDOM", var_name)){
    axis_limits <- c(10^-2, 10^3)
    var_colour<- "brown4"
    if(grepl("TUR|T", var_name)){
      unit <- expression(NTU)
    } else if(var_name == "CDOM"){
      unit <- expression(m-1)
    } else if(var_name == "POC"){
      unit <- expression(mg~m^-3)
    } else {
      unit <- expression(g~m^-3)
    }
  } else {
    stop(paste("Could not find ", var_name))
  }
  
  # List and exit
  to_return <- list("unit" = unit,
                    "var_colour" = var_colour, 
                    "axis_limits" = axis_limits)
  return(to_return)
}

# Convenience colour wrapper
colours_of_stations <- function(){
  
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
                    'Gironde shelf' = viridis::plasma(n = 1, begin = 0.05,end = 0.05),
                    'Rhône shelf' = viridis::rocket(n = 1, begin = 0.80,end = 0.80))
  
  return(colour_values)
}

# Plot linear trends and stats for matched data
# var_combi = zone_all_in_situ_monthly$variable_combi[1]; df = zone_all_in_situ_monthly; df_stats = zone_all_monthly_lm
validation_lm_plots <- function(var_combi, sat_name, median_base, df, df_stats){
  
  # Filter and pivot data
  df_var_sub <- df |> 
    filter(variable_combi == var_combi) |> 
    pivot_longer(cols = value_in_situ:value_satellite) |> 
    mutate(name = gsub("value|_", "", name),
           name = gsub("insitu", "in situ", name)) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Gironde shelf", "Rhône shelf")))
  
  # Get y-axis label and units
  if(df_var_sub$variable[1] == "TEMP"){
    y_lab <- "Temperature (°C)"
    y_lim <- c(4, 32)
    var_sat <- "SST"
  } else if(df_var_sub$variable[1] == "CHLA"){
    y_lab <- "Chlorophyll-a (mg m⁻³)"
    y_lim <- c(0, 20)
    var_sat <- "CHL"
  } else if(df_var_sub$variable[1] == "TUR"){
    y_lab <- "Turbidity (NTU)"
    y_lim <- c(0, 20)
  } else if(df_var_sub$variable[1] == "SPM"){
    y_lab <- "SPM (g m⁻³)"
    y_lim <- c(0, 20)
  } else {
    # Not worried about other variables for the moment
    y_lab <- NA
    y_lim <- c(0, 20)
  }
  
  # Filter stats labels
  df_stats_var_sub <- df_stats |> 
    filter(variable_combi == var_combi) |> 
    mutate(y = max(y_lim)-2) |> 
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Gironde shelf", "Rhône shelf")))
  
  # Create title and subtitle
  var_name_sat <- df_var_sub$variable_sat[1]
  var_name_is <- df_var_sub$variable[1]
  cor_name <- df_var_sub$correction[1]
  proc_name <- df_var_sub$processing[1]
  grid_name <- df_var_sub$grid_size[1]
  plot_title = paste(sat_name, var_name_sat, "vs.", "in situ", var_name_is)
  plot_subtitle = paste("Correction:", cor_name, "| Processing:", proc_name, "| Grid size:", grid_name)
  
  # Create plot
  plot_zone_var_TS <- ggplot(data = df_var_sub, aes(x = date, y = value)) +
    geom_point(aes(y = value, colour = name, shape = source), alpha = 0.4) +
    geom_line(aes(y = value, colour = name, linetype = source), alpha = 0.4) +
    geom_smooth(aes(colour = name, linetype = source), method = "lm", se = FALSE) +
    # In situ values
    geom_label(data = filter(df_stats_var_sub, source == "REPHY"), 
               aes(x = as.Date("2005-01-01"), y = y, label = paste0(source," : ", round(slope_is, 2)," / year"), vjust = 1.0)) +
    geom_label(data = filter(df_stats_var_sub, source == "SOMLIT"), 
               aes(x = as.Date("2005-01-01"), y = y, label = paste0(source," : ", round(slope_is, 2)," / year"), vjust = -0.5)) +
    # Satellite values
    geom_label(data = filter(df_stats_var_sub, source == "REPHY"), colour = "red",
               aes(x = as.Date("2020-01-01"), y = y, label = paste0(source," : ", round(slope_sat, 2)," / year"), vjust = 1.0)) +
    geom_label(data = filter(df_stats_var_sub, source == "SOMLIT"),  colour = "red",
               aes(x = as.Date("2020-01-01"), y = y, label = paste0(source," : ", round(slope_sat, 2)," / year"), vjust = -0.5)) +
    facet_wrap(~zone_pretty,nrow = 2) +
    scale_colour_manual(values = c("black", "red")) +
    # scale_y_continuous(limits = c(min(df_var_sub$value, na.rm = TRUE), max(df_var_sub$value, na.rm = TRUE)*1.1)) +
    coord_cartesian(ylim = y_lim) +
    labs(title = plot_title, subtitle = plot_subtitle,
         y = y_lab, x = NULL, colour = "Type", shape = "Source", linetype = "Source") +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 25, face = "bold"), 
          strip.text = element_text(size = 20),
          legend.title = element_text(size = 23),
          legend.text = element_text(size = 20),
          axis.title = element_text(size = 23),
          axis.text = element_text(size = 20),
          legend.position = "bottom")
  # plot_zone_var_TS
  ggsave(paste0("figures/validation/ts/ts_",sat_name,"_",var_combi,"_",median_base,".png"), plot_zone_var_TS, width = 14, height = 8)
  return()
}

# The figure code wrapper
validation_plots <- function(var_combi, sat_name, median_base, match_up_df, match_up_stats){

  # Subset datasets for chosen variable
  match_up_df_var <- match_up_df |> 
    filter(variable_combi == var_combi) |>
    mutate(zone_pretty = factor(zone,
                                levels = c("BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"),
                                labels = c("Bay of Seine", "S. Brittany", "Gironde shelf", "Rhône shelf")))
  match_up_stats_var <- match_up_stats |> 
    filter(variable_combi == var_combi,
           zone == "GLOBAL",
           source == "ALL",
           site == "ALL",
           season == "ALL")
  
  # Get the variable names etc. for plotting
  var_name_is <- match_up_df_var$variable[1]
  var_name_sat <- match_up_df_var$variable_sat[1]
  cor_name <- match_up_df_var$correction[1]
  proc_name <- match_up_df_var$processing[1]
  grid_name <- match_up_df_var$grid_size[1]
  
  # Get axis labels
  plot_meta_is <- var_labels(var_name_is)
  plot_meta_sat <- var_labels(var_name_sat)

  # Get 1:1 line limits
  identity_line <- data.frame(x = plot_meta_is$axis_limits, 
                              y = plot_meta_sat$axis_limits)

  # Get colour values
  colour_values <- colours_of_stations()
  
  # Get stats to plot
  if(var_name_sat %in% c("SST", "SST-NIGHT")){
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope
    # R2_value <- match_up_stats_var$r2 # Not currently calculated
  } else {
    Error_value <- match_up_stats_var$Error
    Bias_value <- match_up_stats_var$Bias
    Slope_value <- match_up_stats_var$Slope_log
    # R2_value <- match_up_stats_var$r2_log # Not currently calculated
  }
  
  # Create title and subtitle
  plot_title = paste(sat_name, var_name_sat, "vs.", "in situ", var_name_is)
  plot_subtitle = paste("Correction:", cor_name, "| Processing:", proc_name, "| Grid size:", grid_name)
  
  # Create the plot
  scatterplot <- ggplot(data = match_up_df_var, aes(x = value_in_situ, y = value_satellite)) + 
    
    geom_point(aes(colour = zone_pretty, shape = source), size = 6, show.legend = TRUE) + 
    # TODO: Add 2:1 and 1:2 dashed lines
    geom_line(data = identity_line, aes(x = x, y = y), linetype = "dashed", show.legend = FALSE) +
    
    scale_x_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('In~situ~measurements~(', plot_meta_is$unit, ')'))) +
    
    scale_y_continuous(
      trans = "log10",
      labels = trans_format("log10", math_format(10^.x)),
      name = parse(text = paste0('Satellite~estimates~(', plot_meta_sat$unit, ')'))) + 
    
    coord_equal(xlim = plot_meta_is$axis_limits, 
                ylim = plot_meta_sat$axis_limits) +
    
    annotate(geom = 'text', x = plot_meta_is$axis_limits[1], y = plot_meta_sat$axis_limits[2], 
             hjust = 0, vjust = 1, color = "black", size = 12.5,
             label = paste('Error = ', round(ifelse(Error_value |> is.numeric(), Error_value, NA), 1), "%\n",
                           'Bias = ', round(ifelse(Bias_value |> is.numeric(), Bias_value, NA), 1), "%\n",
                           # TODO: Change this to show Slope_log or Slope depending on the variable tested (e.g. SST or not)
                           # 'R²_log = ', round(R2_value, 2),"\n", # Not currently calculated
                           'Slope = ', round(ifelse(Slope_value |> is.numeric(), Slope_value, NA), 2),"\n",
                           'n = ', nrow(match_up_df_var), sep = "")) +
            #  label = paste("Slope = ", round(ifelse(Slope_value |> is.numeric(), Slope_value, NA), 2),"\n",
            #                # 'R² = ', round(statistics_values$r2_log, 2),"\n",
            #                "n = ", nrow(match_up_df_var), sep = "")) +
    
    labs(title = plot_title, subtitle = plot_subtitle) +
    
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
          plot.subtitle = element_text(size = 20, hjust = 0.5, color = "black", face = "bold.italic"),
          plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # scatterplot
  
  if(var_name_sat == "SST"){ 
    scatterplot <- scatterplot +  
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      scale_x_continuous(name = parse(text = paste0('In~situ~measurements~(', plot_meta_is$unit, ')'))) +
      scale_y_continuous(name = parse(text = paste0('Satellite~estimates~(', plot_meta_sat$unit, ')')))
  } else {
    scatterplot <- scatterplot + 
      geom_smooth(method = "lm", colour = "black", se = FALSE) +
      annotation_logticks()
  }
  # scatterplot
  
  # Add histograms to x and y axes
  scatterplot_with_side_hist <- ggMarginal(scatterplot, type = "histogram", groupFill = TRUE, alpha = 1)
  # scatterplot_with_side_hist
  ggsave(paste0("figures/validation/scatterplot/",sat_name,"_",var_combi,".png"), 
         plot = scatterplot_with_side_hist, height = 16, width = 16, bg = "white")
  
  # Barplot of frequency per year
  bar_plot_freq_per_year <- match_up_df_var |> 
    mutate(Year = year(date)) |> 
    dplyr::count(Year) |> 
    ggplot(aes(x = Year)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_is$var_colour, linewidth = 1.5) + 
    scale_x_continuous(breaks = seq(1998, 2025, by = 3), name = "", labels = function(x) substring(x, 3, 4)) +
    scale_y_continuous(expand = c(0,0), name = "n per year") +
    labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_year
  ggsave(paste0("figures/validation/barplot/annual_",sat_name,"_",var_combi,".png"), 
         plot = bar_plot_freq_per_year, height = 10, width = 14, bg = "white")
  
  # Barplot of monthly counts
  bar_plot_freq_per_month <- match_up_df_var |> 
    mutate(Month = month(date)) |> 
    dplyr::count(Month) |> 
    ggplot(aes(x = Month)) + 
    geom_col(aes(y = n), fill = "white", colour = plot_meta_is$var_colour, linewidth = 2) + 
    scale_x_continuous(breaks = 1:12, name = "", labels = function(x) month.abb[x]) +
    scale_y_continuous(expand = c(0,0), name = "n per month") +
    coord_cartesian(xlim = c(1,12)) +
    labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot_theme() +
    theme(plot.title = element_text(color = plot_meta_is$var_colour, face = "bold", size = 35))
  # bar_plot_freq_per_month
  ggsave(paste0("figures/validation/barplot/monthly_",sat_name,"_",var_combi,".png"), 
         plot = bar_plot_freq_per_month, height = 10, width = 14, bg = "white")
  return()
}

# Run all validation stats and produce the plots
# sat_name = "SEXTANT"; median_base = "all"
# sat_name = "MODIS"; median_base = "all"
# sat_name = "OLCI-B"; median_base = "all"
validate_sensor <- function(sat_name, median_base){
  
  # Get the pixel cutoff based on product
  if(median_base == "small"){
    if(sat_name == "SEXTANT"){
      pixel_n_cut <- 1
    } else {
      pixel_n_cut <- 3
    }
  } else {
    if(sat_name == "SEXTANT"){
      pixel_n_cut <- 3
    } else {
      pixel_n_cut <- 16
    }
  }
  
  # Load prepped data
  sat_files <- dir("output/MATCH_UP_DATA/FRANCE", full.names = TRUE,
                   pattern = paste0("zone_median_",sat_name))
  sat_files <- sat_files[grepl(paste0("_",median_base), sat_files)]

  # No zone_median files means process_pixels() was never run for this sensor
  # (currently true for MODIS/MERIS/OLCI-A/OLCI-B, whose ODATIS-MR source data
  # lives only on the old machine's local mount) -- skip rather than crash.
  if(length(sat_files) == 0){
    message("Skipping validate_sensor(\"",sat_name,"\", \"",median_base,"\"): no zone_median files found.")
    return(invisible(NULL))
  }

  zone_median <- map_dfr(sat_files, data.table::fread) |> mutate(date = as.Date(date)) |>
    # Remove all rows that are below the pixel cutoff and CV cutoff of 20%
    mutate(sd = case_when(n <= 2 ~ 0, TRUE ~ sd), # Necessary for CV for 1 pixel count matchups for SEXTANT 'small'
           cv = sd / median) |> 
    filter(n >= pixel_n_cut, cv <= 0.20) #|> 
    # Complete all dates
    # complete(date = seq(min(date), max(date), by = "day"), fill = list(median = NA), 
    #          nesting(zone, source, site, variable))    
  
  # Make variable name conversions as necessary
  if(sat_name == "SEXTANT"){
    # Create a TUR data.frame of SPM data and re-add to the main dataset
    zone_median_tur <- zone_median |> 
      filter(variable == "SPM") |> 
      mutate(variable = "TUR")
    zone_median <- bind_rows(zone_median, zone_median_tur) |> 
      distinct() |>
      mutate(variable_is = variable) |> 
      mutate(variable = case_when(variable == "CHLA" ~ "CHL",
                                  variable == "TUR" ~ "SPM",
                                  TRUE ~ variable)) |> 
      mutate(correction = "Standard",
             processing = "OC5",
             grid_size = case_when(median_base == "small" ~ "1x1",
                                   median_base == "all" ~ "3x3"))
    rm(zone_median_tur); gc()
  } else {
    zone_median <- zone_median |> 
      mutate(variable_is = case_when(grepl("CHL", variable) ~ "CHLA",
                                     grepl("SST", variable) ~ "TEMP",
                                     grepl("SPM", variable) ~ "SPM",
                                     grepl("T-FNU", variable) ~ "TUR",
                                     TRUE ~ variable)) |> 
      mutate(correction = case_when(grepl("-AC", variable) ~ "AC",
                                     grepl("-PO", variable) ~ "PO",
                                     grepl("-NS", variable) ~ "NS")) |> 
      mutate(processing = case_when(grepl("-GONS-", variable) ~ "GONS",
                                     grepl("-OC5-", variable) ~ "OC5",
                                     grepl("-G-", variable) ~ "G",
                                     grepl("-R-", variable) ~ "R",
                                     grepl("-FNU-", variable) ~ "FNU")) |> 
      mutate(grid_size = case_when(median_base == "small" ~ "3x3",
                                  median_base == "all" ~ "7x7")) |> 
      # Clean up variable names for further use
      mutate(variable = gsub("-AC|-PO|-NS", "", variable)) |> 
      mutate(variable = gsub("-GONS|-OC5|-G|-R|-FNU", "", variable))
  }
  
  # Load in situ data and complete the date column
  zone_data_in_situ <- read_csv("data/INSITU_data/zone_data_in_situ.csv", show_col_types = FALSE) |> 
    dplyr::select(-lon, -lat) #|> 
    # complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA), 
    #          nesting(zone, zone_pretty, source, site, time_class, variable))
  
  # Combine extracted sat data with in situ
  zone_all_in_situ_base <- zone_data_in_situ |> 
    left_join(zone_median, by = c("zone", "source", "site", "date", "variable" = "variable_is"),
              relationship = "many-to-many") |> # For day and night time classes for SST
    dplyr::rename(value_in_situ = value, value_satellite = median, variable_sat = variable.y) |> 
    filter(!is.na(variable_sat)) |> 
    filter(!(variable_sat == "SST" & time_class == "night")) |>
    filter(!(variable_sat == "SST-NIGHT" & time_class == "day")) |>
    mutate(season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter", month(date) %in% 3:5  ~ "Spring",
      month(date) %in% 6:8  ~ "Summer", month(date) %in% 9:11 ~ "Autumn"), .after = "date") |> 
    dplyr::select(zone, dplyr::everything(), -time_class) |> 
    mutate(variable_combi = paste(variable, variable_sat, correction, processing, grid_size, sep = "_"))

  # Create monthly average TS for lm analysis
  zone_all_in_situ_monthly <- zone_all_in_situ_base |> 
    mutate(date = floor_date(date, unit = "month")) |> 
    summarise(value_in_situ = mean(value_in_situ, na.rm = TRUE),
              value_satellite = mean(value_satellite, na.rm = TRUE),
              .by = c("variable_combi", "correction", "processing", "grid_size", "zone", "zone_pretty", "source", "season", "date", "variable", "variable_sat"))
  
  # Calculate linear model stats to look at change over time
  zone_in_situ_monthly_lm <- zone_all_in_situ_monthly |> 
    filter(value_in_situ > 0) |> 
    group_by(variable_combi, correction, processing, grid_size, zone, zone_pretty, source, variable, variable_sat) |> 
    do(broom::tidy(lm(value_in_situ ~ date, data = .))) |> 
    filter(term == "date") |> 
    dplyr::rename(slope_is = estimate, p_is = p.value) |> 
    dplyr::select(variable_combi:variable_sat, slope_is, p_is)
  zone_sat_monthly_lm <- zone_all_in_situ_monthly |> 
    filter(value_satellite > 0) |> 
    group_by(variable_combi, correction, processing, grid_size, zone, zone_pretty, source, variable, variable_sat) |> 
    do(broom::tidy(lm(value_satellite ~ date, data = .))) |> 
    filter(term == "date") |> 
    dplyr::rename(slope_sat = estimate, p_sat = p.value) |> 
    dplyr::select(variable_combi:variable_sat, slope_sat, p_sat)
  
  # Combine results
  zone_all_monthly_lm <- left_join(zone_in_situ_monthly_lm, zone_sat_monthly_lm,
                                   by = join_by(variable_combi, correction, processing, grid_size, zone, zone_pretty, 
                                                source, variable, variable_sat)) |> 
    # Convert to values / year
    mutate(slope_is = slope_is*365.25, slope_sat = slope_sat*365.25)
  
  # Save statistics
  write_csv(zone_all_monthly_lm, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/",sat_name,"_lm_stats_",median_base,".csv"))

  # Plot the linear TS matchups
  plan(multisession, workers = 10)
  future_walk(unique(zone_all_in_situ_monthly$variable_combi), validation_lm_plots, 
              sat_name = sat_name, median_base = median_base,
              df = zone_all_in_situ_monthly, df_stats = zone_all_monthly_lm)
  plan(sequential)
  
  # Remove any missing values for stats matchups
  zone_all_in_situ <- zone_all_in_situ_base |> 
    filter(value_in_situ > 0, value_satellite > 0)
  
  # Create big grid of sites to look for obvious outliers
  # ggplot(data = zone_in_situ, aes(x = value_in_situ, y = value_satellite)) +
  #   geom_point(aes(colour = source, shape = variable)) +
  #   facet_wrap(~site, scales = "free")
  
  # Stats for all groups together by variable
  zone_in_situ_stats_01 <- zone_all_in_situ |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all groups together by variable and season
  zone_in_situ_stats_02 <- zone_all_in_situ |> mutate(zone = "GLOBAL", source = "ALL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones by variable
  zone_in_situ_stats_03 <- zone_all_in_situ |> mutate(source = "ALL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones by variable and season
  zone_in_situ_stats_04 <- zone_all_in_situ |> mutate(source = "ALL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sources by variable
  zone_in_situ_stats_05 <- zone_all_in_situ |> mutate(zone = "GLOBAL", site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sources by variable and season
  zone_in_situ_stats_06 <- zone_all_in_situ |> mutate(zone = "GLOBAL", site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones and sources by variable
  zone_in_situ_stats_07 <- zone_all_in_situ |> mutate(site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all zones and sources by variable and season
  zone_in_situ_stats_08 <- zone_all_in_situ |> mutate(site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sites by variable
  zone_in_situ_stats_09 <- zone_all_in_situ |> mutate(site = "ALL", season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for all sites by variable and season
  zone_in_situ_stats_10 <- zone_all_in_situ |> mutate(site = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for each site by variable
  zone_in_situ_stats_11 <- zone_all_in_situ |> mutate(season = "ALL") |> 
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Stats for each site by variable and season
  zone_in_situ_stats_12 <- zone_all_in_situ |>
    summarise(compute_stats(value_in_situ, value_satellite), 
      .by = c("correction", "processing", "grid_size", "zone", "source", "site", "season", "variable", "variable_sat"))
  
  # Bind all together
  zone_all_in_situ_stats <- bind_rows(zone_in_situ_stats_01, zone_in_situ_stats_02, zone_in_situ_stats_03,
                                      zone_in_situ_stats_04, zone_in_situ_stats_05, zone_in_situ_stats_06,
                                      zone_in_situ_stats_07, zone_in_situ_stats_08, zone_in_situ_stats_09,
                                      zone_in_situ_stats_10, zone_in_situ_stats_11, zone_in_situ_stats_12) |> 
    mutate(sensor = sat_name, .before = "correction")
  
  # Save results
  write_csv(zone_all_in_situ_stats, paste0("output/MATCH_UP_DATA/FRANCE/STATISTICS/",sat_name,"_stats_",median_base,".csv"))
  rm(zone_in_situ_stats_01, zone_in_situ_stats_02, zone_in_situ_stats_03,
     zone_in_situ_stats_04, zone_in_situ_stats_05, zone_in_situ_stats_06,
     zone_in_situ_stats_07, zone_in_situ_stats_08, zone_in_situ_stats_09,
     zone_in_situ_stats_10, zone_in_situ_stats_11, zone_in_situ_stats_12); gc()
  
  # Run all of the plots per variable pairing
  zone_all_in_situ_stats$variable_combi <- paste0(zone_all_in_situ_stats$variable,"_",
                                                 zone_all_in_situ_stats$variable_sat, "_",
                                                 zone_all_in_situ_stats$correction, "_",
                                                 zone_all_in_situ_stats$processing, "_",
                                                 zone_all_in_situ_stats$grid_size)
  plan(multisession, workers = 10)
  future_walk(unique(zone_all_in_situ_stats$variable_combi), validation_plots,
                     sat_name = sat_name, median_base = median_base,
                     match_up_df = zone_all_in_situ, match_up_stats = zone_all_in_situ_stats)
  plan(sequential)
}


