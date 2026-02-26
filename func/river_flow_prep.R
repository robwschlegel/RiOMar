# func/river_flow_prep.R

# This script combines and cleans output files from HydroPortail


# Setup -------------------------------------------------------------------

# Needed libraries
library(tidyverse)
library(mgcv) # For GAM model fits between gauge stations

# The zones
zones_list <- c("GULF_OF_LION", "BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY")

# Load HydroPortail meta-data
sites_HP <- read_csv("metadata/HydroPortail_station_list.csv")

# HydroPortail loading function
load_HP <- function(file_name){
  
  # NB: This can be avoided by never useing the batch download option in HydroPortail
  # Detect file format
  # df_line_1 <- readLines(file_name, n = 1)
  
  # if(grepl("SiteHydro", df_line_1)){
  #   df <- read_delim(file_name, skip = 1, delim = ";")
  # } else {
    df <- read_csv(file_name) |> 
      mutate(date = as.Date(`Date (TU)`)) |> 
      filter(`Valeur (en m³/s)` >= 0) |> 
      summarise(debit = mean(`Valeur (en m³/s)`, na.rm = TRUE), .by = date)
  # }
  return(df)
}
prep_HP <- function(zone){
  
  # Automagic directory names
  dir_main <- file.path("data/RIVER_FLOW", zone)
  dir_HP <- file.path("data/RIVER_FLOW", zone, "HydroPortail")
  
  if(zone == "GULF_OF_LION"){
    
    # Grand Rhone
    rhone_g_1 <- map_dfr(dir(dir_HP, pattern = "V730000302", full.names = TRUE), load_HP)# |> mutate(site = "grand")
    rhone_g_2 <- map_dfr(dir(dir_HP, pattern = "V720001002", full.names = TRUE), load_HP) |>#  mutate(site = "all")
      filter(date >= "2023-01-01") |> mutate(debit = debit*0.9)
    
    # Compare the two 
    # rhone_g_3 <- bind_rows(rhone_g_1, rhone_g_2) |>
    #   filter(date >= "2018-01-01", date <= "2022-12-31")
    # rhone_g_4 <- rhone_g_3 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_dif = grand/all)
    # mean(rhone_g_4$prop_dif, na.rm = TRUE)
    # NB: The average proportion is 0.98, not 0.90 as indicated in the literature
    # The proportion of 0.90 is accurate for periods of crue
    # But the Arles station is the same or greater than Tarascon during daily conditions
    # This may require that a more sophisticated correction is used
    # ggplot(rhone_g_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(rhone_g_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(rhone_g_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Combine and save
    grand_rhone <- rbind(rhone_g_1, rhone_g_2) |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(grand_rhone, file.path(dir_main, "grand_rhone.csv"))
    
    # Petit Rhone
    rhone_p_1 <- map_dfr(dir(dir_HP, pattern = "V730000202", full.names = TRUE), load_HP)# |> mutate(site = "petit")
    rhone_p_2 <- map_dfr(dir(dir_HP, pattern = "V720001002", full.names = TRUE), load_HP) |># mutate(site = "all")
      # filter(date >= "2023-01-01") |> 
      mutate(debit = debit*0.1)
    
    # Compare the two 
    # rhone_p_3 <- bind_rows(rhone_p_1, rhone_p_2) |>
    #   filter(date >= "2018-01-01", date <= "2022-12-31")
    # rhone_p_4 <- rhone_p_3 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_dif = petit/all)
    # mean(rhone_p_4$prop_dif, na.rm = TRUE)
    # NB: The average proportion is 0.11, which is close to the literature of 0.10
    # ggplot(rhone_p_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(rhone_p_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(rhone_p_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm") + geom_smooth(method = "lm")
    
    # Combine and save
    petit_rhone <- rbind(rhone_p_1, rhone_p_2) |>  summarise(debit = mean(debit, na.rm = TRUE), .by = date) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(petit_rhone, file.path(dir_main, "petit_rhone.csv"))
  } else if(zone == "BAY_OF_SEINE"){
    
    # Seine
    seine_1 <- rbind(map_dfr(dir(dir_HP, pattern = "H320000101", full.names = TRUE), load_HP),
                     map_dfr(dir(dir_HP, pattern = "H320000104", full.names = TRUE), load_HP)) |> 
      # There is one day of overlap at 2006-08-16
      summarise(debit = mean(debit, na.rm = TRUE), .by = "date") |> #mutate(site = "complete") |> 
      mutate(debit = debit*1.1) # Correction based on closest gauge to river mouth
    
    # Compare
    # seine_2 <- map_dfr(dir(dir_HP, pattern = "H322011003", full.names = TRUE), load_HP) |>  mutate(site = "closest")
    # seine_3 <- rbind(seine_1, seine_2) |> 
    #   filter(date >= "1998-01-01", date <= "2005-12-31")
    # seine_4 <- seine_3 |>
    #   pivot_wider(values_from = debit, names_from = site) |> 
    #   mutate(prop_dif = closest/complete)
    # mean(seine_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 1.10 times the complete gauge.
    # The peaks are noticeably higher in the closest gauge
    # ggplot(seine_3, aes(x = date, y = debit)) + geom_line(aes(colour = site))
    # ggplot(seine_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(seine_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Combine and save
    seine <- seine_1 |> complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
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
    
    # Lay
    lay_1 <- map_dfr(dir(dir_HP, pattern = "N3511610", full.names = TRUE), load_HP)# |> mutate(site = "closest") 
    lay_2 <- map_dfr(dir(dir_HP, pattern = "N330161010", full.names = TRUE), load_HP) |> #mutate(site = "complete") |> 
      # filter(date >= "2023-01-01") |> 
      mutate(debit = debit*2.33)
    
    # Compare
    # lay_3 <- rbind(lay_1, lay_2) |> 
    #   filter(date >= "2004-01-01", date <= "2025-12-31")
    # lay_4 <- lay_3 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_dif = closest/complete)
    # mean(lay_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 2.33 times the complete gauge.
    # This may warrant a more sophisticated correction than a flat conversion rate
    # For example, more accurately track the proportion ration based on how large the debit values are
    # This should be doable with a linear model and applied the correction based on the predicted proportion ratio for each day
    # ggplot(lay_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(lay_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(lay_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Combine and save
    lay <- rbind(lay_1, lay_2) |> summarise(debit = mean(debit, na.rm = TRUE), .by = date) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(lay, file.path(dir_main, "lay.csv"))
    
    # Vilaine
    vilaine_1 <- map_dfr(dir(dir_HP, pattern = "J930061101", full.names = TRUE), load_HP)# |> mutate(site = "closest")
    vilaine_2 <- map_dfr(dir(dir_HP, pattern = "J770061002", full.names = TRUE), load_HP) |># mutate(site = "complete") |> 
      filter(date <= "2002-08-08") |> mutate(debit = debit*2.34)
    
    # Compare
    # vilaine_3 <- rbind(vilaine_1, vilaine_2)
    # vilaine_4 <- vilaine_3 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_dif = closest/complete)
    # mean(vilaine_4$prop_dif, na.rm = TRUE)
    # NB: The closest gauge is 2.34 times the complete gauge.
    # ggplot(vilaine_3, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(vilaine_3, aes(x = site, y = debit)) + geom_boxplot(aes(fill = site))
    # ggplot(vilaine_4, aes(x = date, y = prop_dif)) + geom_col() + geom_smooth(method = "lm")
    
    # Combine and save
    vilaine <- rbind(vilaine_1, vilaine_2) |> summarise(debit = mean(debit, na.rm = TRUE), .by = date) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(vilaine, file.path(dir_main, "vilaine.csv"))
    
  } else if(zone == "BAY_OF_BISCAY"){
    
    # Garonne
    garonne_1 <- map_dfr(dir(dir_HP, pattern = "O909001001", full.names = TRUE), load_HP) #|> mutate(site = "closest")
    garonne_2 <- map_dfr(dir(dir_HP, pattern = "O900001002", full.names = TRUE), load_HP) |> #mutate(site = "backup")
      mutate(debit = debit*1.04) # Correction based on closest gauge to river mouth
    
    # Compare
    # garonne_3 <- rbind(garonne_1, garonne_2)
    # garonne_4 <- garonne_3 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_dif = closest/backup)
    # mean(garonne_4$prop_dif, na.rm = TRUE)
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
    charente_1 <- map_dfr(dir(dir_HP, pattern = "R523001001", full.names = TRUE), load_HP)# |> mutate(site = "dist_1")
    charente_2 <- map_dfr(dir(dir_HP, pattern = "R520001001", full.names = TRUE), load_HP)# |> mutate(site = "dist_2")
    charente_3 <- map_dfr(dir(dir_HP, pattern = "R423001001", full.names = TRUE), load_HP)# |> mutate(site = "dist_3")
    charente_4 <- map_dfr(dir(dir_HP, pattern = "R314001001", full.names = TRUE), load_HP)# |> mutate(site = "dist_4")
    charente_5 <- map_dfr(dir(dir_HP, pattern = "R222001001", full.names = TRUE), load_HP) |># mutate(site = "dist_5") |> 
      mutate(debit = debit*2.80) # Correction based on closest gauge to river mouth, which is 0.31 of the furthest gauge
    
    # Compare
    # charente_6 <- rbind(charente_1, charente_2, charente_3, charente_4, charente_5)
    # charente_7 <- charente_6 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_1 = dist_1/dist_3, prop_2 = dist_2/dist_3, prop_3 = dist_3/dist_4, prop_4 = dist_4/dist_5, prop_5 = dist_2/dist_5)
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
    
    # Combine and save
    charente <- rbind(charente_2, charente_5) |> summarise(debit = mean(debit, na.rm = TRUE), .by = date) |>
      complete(date = seq.Date(min(date), max(date), by = "day")) |> arrange(date)
    write_csv(charente, file.path(dir_main, "charente.csv"))
    
    # Sevre_niortaise
    sevre_1 <- map_dfr(dir(dir_HP, pattern = "N611061001", full.names = TRUE), load_HP) |># mutate(site = "dist_1") |> 
      filter(debit <= 250) # Some bad outliers
    sevre_2 <- map_dfr(dir(dir_HP, pattern = "N430062201", full.names = TRUE), load_HP) |># mutate(site = "dist_2") |> 
      mutate(debit = debit*3.4) # Correction based on closest gauge to river mouth, which is 0.84 of the furthest gauge
    sevre_3 <- map_dfr(dir(dir_HP, pattern = "N401061001", full.names = TRUE), load_HP) |># mutate(site = "dist_3") |> 
      mutate(debit = debit*4.93) # Correction based on closest gauge to river mouth, which is 0.20 of the furthest gauge
    
    # GAM models
    # This is an interesting idea, but the lack of consistent overlap is problematic
    # sevre_1_2 <- inner_join(sevre_1, sevre_2, by = "date")
    # gam_2 <- gam(debit.y ~ s(debit.x), data = sevre_1_2, method = "REML")
    # summary(gam_2)
    # sevre_2_pred <- predict(gam_2, newdata = data.frame(debit.x = sevre_1_2$debit.x))
    # sevre_2$debit <- sevre_2$debit - (sevre_2$debit - sevre_2_pred)
    
    # Compare
    # sevre_4 <- rbind(sevre_1, sevre_2, sevre_3)
    # sevre_5 <- sevre_4 |>
    #   pivot_wider(values_from = debit, names_from = site) |>
    #   mutate(prop_1 = dist_1/dist_2, prop_2 = dist_1/dist_3, prop_3 = dist_2/dist_3)
    # mean(sevre_5$prop_1, na.rm = TRUE) # dist_1 is 9.44 of dist_2
    # mean(sevre_5$prop_2, na.rm = TRUE) # dist_1 is 4.93 of dist_3
    # mean(sevre_5$prop_3, na.rm = TRUE) # dist_2 is 0.84 of dist_3
    # NB: dist_2 appears numerically smaller than dist_3, but it captures the peaks better
    # ggplot(sevre_4, aes(x = date, y = debit)) + geom_line(aes(colour = site)) + geom_smooth(method = "lm", aes(colour = site))
    # ggplot(sevre_5, aes(x = dist_1, y = dist_2)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    # ggplot(sevre_5, aes(x = dist_1, y = dist_3)) + geom_point() + geom_smooth(method = "lm") + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
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
