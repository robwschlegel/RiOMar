# func/river_flow_prep.R

# This script combines and cleans output files from HydroPortail


# Setup -------------------------------------------------------------------

# Needed libraries
library(tidyverse)

# The zones
zones_list <- c("GULF_OF_LION", "BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY")

# Load HydroPortail meta-data
sites_HP <- read_csv("metadata/HydroPortail_station_list.csv")

# Log-log linear model prediction function
# pred_x <- rhone_g_4$grand; pred_y <- rhone_g_4$all
# pred_x <- lay_5$debit.y; pred_y <- lay_5$debit.x
# rm(pred_x, pred_y, df, df_res, lm_fit)
predict_loglog <- function(pred_x, pred_y){
  
  # Prep dataframe of matching values
  df <- data.frame(pred_x = pred_x, pred_y = pred_y) |> 
    filter(pred_x > 1, pred_y > 1) |> # Log transform fall over under 1
    mutate(log_x = log(pred_x), log_y = log(pred_y),
           time_step = 1:n())
  
  # Test plots
  ggplot(df, aes(x = time_step, y = pred_x)) + geom_line() + geom_smooth(method = "lm", colour = "black") +  
    geom_line(aes(y = pred_y), colour = "red") + geom_smooth(method = "lm", aes(y = pred_y), colour = "red")
  ggplot(df, aes(x = pred_x, y = pred_y)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  ggplot(df, aes(x = log_x, y = log_y)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  
  # Get slope and intercept
  lm_fit <- lm(log_y ~ log_x, data = df)
  # summary(lm_fit)
  
  # Test out the matchup
  # df$pred_pred <- exp((log(df$pred_y)/coef(lm_fit)[2])-coef(lm_fit)[1])
  # ggplot(df, aes(x = time_step, y = pred_x)) + geom_line() + geom_smooth(method = "lm", colour = "black") +
  #   geom_line(aes(y = pred_y), colour = "blue") + geom_smooth(method = "lm", aes(y = pred_y), colour = "blue") +
  #   geom_line(aes(y = pred_pred), colour = "red") + geom_smooth(method = "lm", aes(y = pred_pred), colour = "red")
  # ggplot(df, aes(x = pred_x, y = pred_pred)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  
  # Correct the longest series based on the ideal series
  df_res <- data.frame(pred_long = exp((log(pred_y)/coef(lm_fit)[2])-coef(lm_fit)[1]))
  
  # Exit
  return(df_res$pred_long)
}


# HydroPortail loading function
load_HP <- function(file_name){
  df <- read_csv(file_name) |> 
      mutate(date = as.Date(`Date (TU)`)) |> 
      filter(`Valeur (en m³/s)` >= 0) |> 
      summarise(debit = mean(`Valeur (en m³/s)`, na.rm = TRUE), .by = date)
  return(df)
}

# HydroPortail preparation function
# zone <- zones_list[3]
prep_HP <- function(zone){
  
  # Automagic directory names
  dir_main <- file.path("data/RIVER_FLOW", zone)
  dir_HP <- file.path("data/RIVER_FLOW", zone, "HydroPortail")
  
  if(zone == "GULF_OF_LION"){
    
    # Grand Rhone
    rhone_g_1 <- map_dfr(dir(dir_HP, pattern = "V730000302", full.names = TRUE), load_HP) |> mutate(site = "grand")
    rhone_g_2 <- map_dfr(dir(dir_HP, pattern = "V720001002", full.names = TRUE), load_HP) |>  mutate(site = "all")
      # filter(date >= "2023-01-01") |> mutate(debit = debit*0.9)
    
    # Compare the two 
    rhone_g_3 <- bind_rows(rhone_g_1, rhone_g_2) |>
      filter(date >= "2018-01-01", date <= "2022-12-31")
    rhone_g_4 <- bind_rows(rhone_g_1, rhone_g_2) |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_dif = grand/all)
    # mean(rhone_g_4$prop_dif, na.rm = TRUE)
    # NB: The average proportion is 0.98, not 0.90 as indicated in the literature
    # The proportion of 0.90 is accurate for periods of crue
    # But the Arles station is the same or greater than Tarascon during daily conditions
    # This appears to require that a more sophisticated correction is used
    # ggplot(rhone_g_4, aes(x = all, y = grand)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # ggplot(rhone_g_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(rhone_g_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(rhone_g_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Log-log transform data and run linear model to find a correction between gauges
    rhone_g_4 <- rhone_g_4 |> 
      mutate(grand_pred = predict_loglog(pred_x = grand, pred_y = all),
             grand_90 = all*0.9) |> 
      mutate(prop_dif_pred = grand_pred/grand,
             prop_dif_90 = grand_90/grand)
    # mean(rhone_g_4$prop_dif_pred, na.rm = TRUE)
    # mean(rhone_g_4$prop_dif_90, na.rm = TRUE)
    # ggplot(rhone_g_4, aes(x = grand, y = grand_pred)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # ggplot(rhone_g_4, aes(x = grand, y = grand_90)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # rhone_g_4 |> pivot_longer(cols = c(all, grand, grand_pred, grand_90), names_to = "site", values_to = "debit") |>
    #   ggplot(aes(x = date, y = debit)) + geom_line(aes(colour = site), alpha = 0.6) + geom_smooth(method = "lm", aes(colour = site))
    # cor(rhone_g_4$grand, rhone_g_4$grand_pred, method = "pearson", use = "complete.obs")
    # cor(rhone_g_4$grand, rhone_g_4$grand_90, method = "pearson", use = "complete.obs")
    # The linear model method seems to work better than a flat conversion rate, especially during periods of crue, but the average is too high...
    
    # Correct current data based on the linear model match
    # rhone_g_5 <- left_join(rhone_g_2, rhone_g_1, by = "date") |> mutate(debit = predict_loglog(debit.y, debit.x)) |> 
    #   filter(date >= "2023-01-01") |> dplyr::select(date, debit, site.x) |> rename(site = site.x) |> mutate(site = "grand_pred")
    rhone_g_5 <- rhone_g_2 |> mutate(debit = debit*0.9) |> 
      filter(date >= "2023-01-01") |> dplyr::select(date, debit, site) |> mutate(site = "grand_90")
    
    # Combine and save
    grand_rhone <- rbind(rhone_g_1, rhone_g_5) |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date) |> dplyr::select(date, debit)
    # ggplot(grand_rhone, aes(x = date, y = debit)) +
    # geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
      # geom_line() + geom_smooth(method = "lm")
    write_csv(grand_rhone, file.path(dir_main, "grand_rhone.csv"))
    
    # Petit Rhone
    rhone_p_1 <- map_dfr(dir(dir_HP, pattern = "V730000202", full.names = TRUE), load_HP) |> mutate(site = "petit")
    rhone_p_2 <- map_dfr(dir(dir_HP, pattern = "V720001002", full.names = TRUE), load_HP) |> mutate(site = "all")
      # filter(date >= "2023-01-01") |> 
      # mutate(debit = debit*0.1)
    
    # Compare the two 
    rhone_p_3 <- bind_rows(rhone_p_1, rhone_p_2) |>
      filter(date >= "2018-01-01", date <= "2022-12-31")
    rhone_p_4 <- rhone_p_3 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_dif = petit/all)
    # mean(rhone_p_4$prop_dif, na.rm = TRUE)
    # NB: The average proportion is 0.11, which is close to the literature of 0.10
    # ggplot(rhone_p_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(rhone_p_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(rhone_p_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm") + geom_smooth(method = "lm")
    
    # Correct current data based on the linear model match
    # rhone_p_5 <- left_join(rhone_p_2, rhone_p_1, by = "date") |> mutate(debit = predict_loglog(debit.y, debit.x)) |> 
    #   filter(date >= "2023-01-01") |> dplyr::select(date, debit, site.x) |> rename(site = site.x)
    rhone_p_5 <- rhone_p_2 |> mutate(debit = debit*0.1) |> 
      filter(date >= "2023-01-01") |> dplyr::select(date, debit, site) |> mutate(site = "grand_90")
    # ggplot(rbind(rhone_p_1, rhone_p_5), aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    
    # Combine and save
    petit_rhone <- rbind(rhone_p_1, rhone_p_5) |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date) |> dplyr::select(date, debit)
    write_csv(petit_rhone, file.path(dir_main, "petit_rhone.csv"))
    
  } else if(zone == "BAY_OF_SEINE"){
    
    # Seine
    seine_1 <- rbind(map_dfr(dir(dir_HP, pattern = "H320000101", full.names = TRUE), load_HP),
                     map_dfr(dir(dir_HP, pattern = "H320000104", full.names = TRUE), load_HP)) |> 
      # There is one day of overlap at 2006-08-16
      summarise(debit = mean(debit, na.rm = TRUE), .by = "date") |> mutate(site = "complete")# |> 
      # mutate(debit = debit*1.1) # Correction based on closest gauge to river mouth
    
    # Compare
    seine_2 <- map_dfr(dir(dir_HP, pattern = "H322011003", full.names = TRUE), load_HP) |> mutate(site = "closest")
    seine_3 <- rbind(seine_1, seine_2) |>
      filter(date >= "1998-01-01", date <= "2005-12-31")
    seine_4 <- seine_3 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_dif = closest/complete)
    # mean(seine_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 1.10 times the complete gauge.
    # The peaks are noticeably higher in the closest gauge
    # ggplot(seine_3, aes(x = date, y = debit)) + geom_line(aes(colour = site))
    # ggplot(seine_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(seine_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Log-log transform data and run linear model to find a correction between gauges
    seine_5 <- left_join(seine_1, seine_2, by = "date") |> mutate(debit = predict_loglog(debit.y, debit.x)) |> 
      dplyr::select(date, debit, site.x) |> rename(site = site.x) |> mutate(site = "closest_pred")
    # ggplot(filter(rbind(seine_2, seine_5), date <= "2005-12-31"),
    #        aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    
    # Combine and save
    seine <- seine_5 |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date) |> dplyr::select(date, debit)
    write_csv(seine, file.path(dir_main, "seine.csv"))
    
    # Orne
    orne_1 <- map_dfr(dir(dir_HP, pattern = "I362101001", full.names = TRUE), load_HP)
    
    # Plot
    # ggplot(orne_1, aes(x = date, y = debit)) + geom_line() + geom_smooth(method = "lm")
    
    # Combine and save
    orne <- orne_1 |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(orne, file.path(dir_main, "orne.csv"))
    
  } else if(zone == "SOUTHERN_BRITTANY"){
    
    # Loire
    loire_1 <- map_dfr(dir(dir_HP, pattern = "L800001020", full.names = TRUE), load_HP)
    
    # Plot
    # ggplot(loire_1, aes(x = date, y = debit)) + geom_line() + geom_smooth(method = "lm")
    
    # Combine and save
    loire <- loire_1 |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(loire, file.path(dir_main, "loire.csv"))
    
    ### Lay
    lay_1 <- map_dfr(dir(dir_HP, pattern = "N3511610", full.names = TRUE), load_HP) |> mutate(site = "closest") 
    lay_2 <- map_dfr(dir(dir_HP, pattern = "N330161010", full.names = TRUE), load_HP) |> mutate(site = "complete") #|> 
      # filter(date >= "2023-01-01") |> 
      # mutate(debit = debit*2.33)
    
    # Compare
    lay_3 <- rbind(lay_1, lay_2) |>
      filter(date >= "2004-01-01", date <= "2025-12-31")
    lay_4 <- lay_3 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_dif = closest/complete)
    # mean(lay_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 2.33 times the complete gauge.
    # This may warrant a more sophisticated correction than a flat conversion rate
    # For example, more accurately track the proportion ration based on how large the debit values are
    # This should be doable with a linear model and applied the correction based on the predicted proportion ratio for each day
    # ggplot(lay_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(lay_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(lay_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Log-log transform data and run linear model to find a correction between gauges
    lay_5 <- left_join(lay_2, lay_1, by = "date") |> mutate(debit = predict_loglog(debit.y, debit.x)) |> 
      dplyr::select(date, debit, site.x) |> rename(site = site.x) |> mutate(site = "complete_pred")
    # ggplot(filter(rbind(lay_1, lay_2, lay_5), date >= "2003-01-01", date <= "2005-12-31"), 
    # ggplot(rbind(lay_1, lay_2, lay_5),
    #        aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    
    # Combine and save
    lay <- left_join(lay_5, lay_1, by = "date") |> mutate(debit = ifelse(is.na(debit.y), debit.x, debit.y)) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date) |> dplyr::select(date, debit)
    write_csv(lay, file.path(dir_main, "lay.csv"))
    
    ### Vilaine
    vilaine_1 <- map_dfr(dir(dir_HP, pattern = "J930061101", full.names = TRUE), load_HP) |> mutate(site = "closest")
    vilaine_2 <- map_dfr(dir(dir_HP, pattern = "J770061002", full.names = TRUE), load_HP) |> mutate(site = "complete") #|> 
      # filter(date <= "2002-08-08") |> mutate(debit = debit*2.34)
    
    # Compare
    vilaine_3 <- rbind(vilaine_1, vilaine_2)
    vilaine_4 <- vilaine_3 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_dif = closest/complete)
    # mean(vilaine_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 2.34 times the complete gauge.
    # ggplot(vilaine_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(vilaine_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(vilaine_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Log-log transform data and run linear model to find a correction between gauges
    # vilaine_5 <- left_join(vilaine_2, vilaine_1, by = "date") |> mutate(debit = predict_loglog(debit.y, debit.x)) |> 
    #   dplyr::select(date, debit, site.x) |> rename(site = site.x) |> mutate(site = "complete_pred")
    vilaine_5 <- vilaine_2 |> mutate(debit = debit* 2.34, site = "complete_pred") # NB: Values are too small for loglog transform
    # ggplot(filter(rbind(lay_1, lay_2, lay_5), date >= "2003-01-01", date <= "2005-12-31"),
    # ggplot(rbind(vilaine_1, vilaine_2, vilaine_5),
    #        aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    
    # Combine and save
    vilaine <- left_join(vilaine_5, vilaine_1, by = "date") |> mutate(debit = ifelse(is.na(debit.y), debit.x, debit.y)) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date) |> dplyr::select(date, debit)
    write_csv(vilaine, file.path(dir_main, "vilaine.csv"))
    
  } else if(zone == "BAY_OF_BISCAY"){
    
    # Garonne
    garonne_1 <- map_dfr(dir(dir_HP, pattern = "O909001001", full.names = TRUE), load_HP) |> mutate(site = "closest")
    garonne_2 <- map_dfr(dir(dir_HP, pattern = "O900001002", full.names = TRUE), load_HP) |> mutate(site = "backup")
      # mutate(debit = debit*1.04) # Correction based on closest gauge to river mouth
    
    # Compare
    garonne_3 <- rbind(garonne_1, garonne_2)
    garonne_4 <- garonne_3 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_dif = closest/backup)
    mean(garonne_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 1.04 times the backup gauge, which is nearby
    # This seems fine, no need for any corrections, but can use them to average out any funny records
    # ggplot(garonne_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(garonne_4, aes(x = closest, y = backup)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # ggplot(garonne_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(garonne_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Combine and save
    garonne <- rbind(garonne_1, garonne_2) |> summarise(debit = mean(debit, na.rm = TRUE), .by = date) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(garonne, file.path(dir_main, "garonne.csv"))
    
    # Dordogne
    dordogne_1 <- map_dfr(dir(dir_HP, pattern = "P555001001", full.names = TRUE), load_HP)
    
    # Plot
    # ggplot(dordogne_1, aes(x = date, y = debit)) + geom_line() + geom_smooth(method = "lm")
    
    # Combine and save
    dordogne <- dordogne_1 |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(dordogne, file.path(dir_main, "dordogne.csv"))
    
    # Charente
    charente_1 <- map_dfr(dir(dir_HP, pattern = "R523001001", full.names = TRUE), load_HP) |> mutate(site = "dist_1")
    charente_2 <- map_dfr(dir(dir_HP, pattern = "R520001001", full.names = TRUE), load_HP) |> mutate(site = "dist_2")
    charente_3 <- map_dfr(dir(dir_HP, pattern = "R423001001", full.names = TRUE), load_HP) |> mutate(site = "dist_3")
    charente_4 <- map_dfr(dir(dir_HP, pattern = "R314001001", full.names = TRUE), load_HP) |> mutate(site = "dist_4")
    charente_5 <- map_dfr(dir(dir_HP, pattern = "R222001001", full.names = TRUE), load_HP) |> mutate(site = "dist_5") #|> 
      # mutate(debit = debit*2.80) # Correction based on closest gauge to river mouth, which is 0.31 of the furthest gauge
    
    # Compare
    charente_6 <- rbind(charente_1, charente_2, charente_3, charente_4, charente_5)
    charente_7 <- charente_6 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_1 = dist_1/dist_3, prop_2 = dist_2/dist_3, prop_3 = dist_3/dist_4, prop_4 = dist_4/dist_5, prop_5 = dist_2/dist_5)
    # mean(charente_7$prop_1, na.rm = TRUE) # dist_1 is 1.19 greater than dist_3 (no dist_2 overlap)
    # mean(charente_7$prop_2, na.rm = TRUE) # dist_2 is 1.54 greater than dist_3
    # mean(charente_7$prop_3, na.rm = TRUE) # dist_3 is 1.36 of dist_4
    # mean(charente_7$prop_4, na.rm = TRUE) # dist_4 is 1.53 of dist_5
    # mean(charente_7$prop_5, na.rm = TRUE) # dist_2 is 4.74 of dist_5
    # NB: dist_1 is too short to provide reliable comparisons, taking dist_2 as the benchmark
    # ggplot(charente_6, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(charente_7, aes(x = dist_2, y = dist_5)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # ggplot(charente_7, aes(x = date, y = prop_5)) + geom_col() + geom_smooth(method = "lm")
    # ggplot(charente_6, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    
    # Log-log transform data and run linear model to find a correction between gauges
    # charente_8 <- left_join(charente_5, charente_2, by = "date") |> mutate(debit = predict_loglog(debit.y, debit.x)) |>
    #   dplyr::select(date, debit, site.x) |> rename(site = site.x) |> mutate(site = "dist_2_pred")
    charente_8 <- charente_5 |> mutate(debit = debit*2.2+15, site = "dist_2_pred")
    # ggplot(rbind(charente_2, charente_5, charente_8),
    #        aes(x = date, y = debit)) + geom_line(aes(colour = site), alpha = 0.3) + geom_smooth(method = "lm", aes(colour = site))
    # rbind(charente_2, charente_5, charente_8) |> 
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   ggplot(aes(x = dist_2, y = dist_2_pred)) + geom_point() + 
    #   geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    
    # Combine and save
    charente <- left_join(charente_8, charente_2, by = "date") |> mutate(debit = ifelse(is.na(debit.y), debit.x, debit.y)) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date) |> dplyr::select(date, debit)
    # ggplot(charente, aes(x = date, y = debit)) + geom_line() + geom_smooth(method = "lm")
    write_csv(charente, file.path(dir_main, "charente.csv"))
    
    # Sevre_niortaise
    sevre_1 <- map_dfr(dir(dir_HP, pattern = "N611061001", full.names = TRUE), load_HP) |> mutate(site = "dist_1") |> 
      filter(debit <= 250) # Some bad outliers
    sevre_2 <- map_dfr(dir(dir_HP, pattern = "N430062201", full.names = TRUE), load_HP) |> mutate(site = "dist_2") |> 
      mutate(debit = debit*4.6) # Correction based on closest gauge to river mouth
    sevre_3 <- map_dfr(dir(dir_HP, pattern = "N401061001", full.names = TRUE), load_HP) |> mutate(site = "dist_3") |> 
      mutate(debit = (debit*5.3)) # Correction based on closest gauge to river mouth
    
    # Compare
    sevre_4 <- rbind(sevre_1, sevre_2, sevre_3)
    sevre_5 <- sevre_4 |>
      pivot_wider(values_from = debit, names_from = site) |>
      mutate(prop_1 = dist_1/dist_2, prop_2 = dist_1/dist_3, prop_3 = dist_2/dist_3)
    mean(sevre_5$prop_1, na.rm = TRUE) # dist_1 is 9.44 of dist_2
    mean(sevre_5$prop_2, na.rm = TRUE) # dist_1 is 4.93 of dist_3
    mean(sevre_5$prop_3, na.rm = TRUE) # dist_2 is 0.84 of dist_3
    # NB: dist_2 appears numerically smaller than dist_3, but it captures the peaks better
    ggplot(sevre_4, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    ggplot(sevre_5, aes(x = dist_1, y = dist_2)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    ggplot(sevre_5, aes(x = dist_1, y = dist_3)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # ggplot(sevre_5, aes(x = date, y = prop_2)) + geom_col() + geom_smooth(method = "lm")
    # ggplot(sevre_4, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    
    # Combine and save
    sevre <- rbind(sevre_1, sevre_2, sevre_3) |> summarise(debit = mean(debit, na.rm = TRUE), .by = date) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(sevre, file.path(dir_main, "sevre_niortaise.csv"))
    
  }
}

# Run it
walk(zones_list, prep_HP)
