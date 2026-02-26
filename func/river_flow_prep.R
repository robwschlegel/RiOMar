# func/river_flow_prep.R

# This script combines and cleans output files from HydroPortail


# Setup -------------------------------------------------------------------

# Needed libraries
library(tidyverse)

# Load HydroPortail meta-data
sites_HP <- read_csv("metadata/HydroPortail_station_list.csv")

# HydroPortail loading function
load_HP <- function(file_name){
  df <- read_csv(file_name) |> 
    mutate(date = as.Date(`Date (TU)`)) |> 
    filter(`Valeur (en m³/s)` >= 0) |> 
    summarise(debit = mean(`Valeur (en m³/s)`, na.rm = TRUE), .by = date)
}

# Gulf of Lion

## Grand Rhone
rhone_1 <- map_dfr(dir("data/RIVER_FLOW/GULF_OF_LION/HydroPortail", pattern = "V730000302", full.names = TRUE), load_HP) |> mutate(site = "grand")
rhone_2 <- map_dfr(dir("data/RIVER_FLOW/GULF_OF_LION/HydroPortail", pattern = "V720001002", full.names = TRUE), load_HP) |> mutate(site = "all")
rhone_3 <- bind_rows(rhone_1, rhone_2) |> 
  filter(date >= "2018-01-01", date <= "2022-12-31")
rhone_4 <- rhone_3 |> 
  pivot_wider(values_from = debit, names_from = site) |> 
  mutate(prop_dif = grand/all)
mean(rhone_4$prop_dif, na.rm = TRUE)

# Compare the two 
ggplot(rhone_3, aes(x = date, y = debit)) +
  geom_line(aes(colour = site))
ggplot(rhone_3, aes(x = site, y = debit)) +
  geom_boxplot(aes(fill = site))
ggplot(rhone_4, aes(x = date, y = prop_dif)) +
  geom_col()
