# One-off: compute satellite-vs-in-situ percent error and percent bias for
# manuscript Table 4, from the raw SOMLIT/REPHY match-up file produced by
# func/validate.py (output/MATCH_UP_DATA/FRANCE/summary.csv). Uses the 3x3
# pixel-window satellite mean (matching Figure 2's validation scatterplot
# convention).
#
# Covers both in-situ variables in that file: SPM (SOMLIT) and turbidity/TURB
# (REPHY) -- REPHY does not measure SPM at all, only turbidity, which is why
# a separate turbidity-vs-satellite-SPM comparison is needed to make any use
# of the REPHY match-ups (see Section 2.3.1's discussion of Doxaran et al.
# 2009, who calibrated in-situ turbidity sensors to give SPM concentration in
# the Gironde, the precedent for treating the two as comparable here).
#
# Two data-quality fixes are applied in-place below (also fixed at the
# source in func/validate.py's region_mapping/zone_mapping for future
# reruns, since Stage 1 -- code/1_validate.py -- currently fails downstream
# of where this file is written; see the pipeline map's "Still open" list):
# 'GULF OF BISCAY' -> 'BAY OF BISCAY' (typo, inconsistent with every other
# table/figure in the manuscript), and the SOMLIT station 'Lanveoc' (Rade de
# Brest), which validate.py's region_mapping didn't recognise and left as
# 'Region to fill', reassigned to Southern Brittany alongside the
# neighbouring Portzic station.
#
# Definitions (following standard ocean-colour validation practice, e.g.
# Pahlevan et al. 2021):
#   percent error = mean(|satellite - in situ| / in situ) x 100  (MAPE)
#   percent bias  = mean((satellite - in situ) / in situ) x 100
#
# Run from repo root: Rscript func/compute_validation_stats.R
library(dplyr)
library(readr)

raw <- read_csv("output/MATCH_UP_DATA/FRANCE/summary.csv", show_col_types = FALSE)

zone_display <- c(
  "BAY OF SEINE"      = "Bay of Seine",
  "SOUTHERN BRITTANY" = "Southern Brittany",
  "BAY OF BISCAY"     = "Bay of Biscay",
  "GULF OF LION"       = "Gulf of Lions"
)

validation_stats <- raw |>
  mutate(REGION = recode(REGION, "GULF OF BISCAY" = "BAY OF BISCAY"),
         REGION = if_else(SITE == "Lanvéoc", "SOUTHERN BRITTANY", REGION)) |>
  filter(Insitu_variable %in% c("SPM", "TURB"), REGION %in% names(zone_display),
         !is.na(`3x3_mean`), Insitu_value > 0) |>
  mutate(pct_error = abs(`3x3_mean` - Insitu_value) / Insitu_value * 100,
         pct_bias  = (`3x3_mean` - Insitu_value) / Insitu_value * 100) |>
  summarise(n = n(),
            error_pct = mean(pct_error, na.rm = TRUE),
            bias_pct  = mean(pct_bias, na.rm = TRUE),
            .by = c("REGION", "Insitu_variable")) |>
  mutate(zone = zone_display[REGION], .after = "REGION") |>
  arrange(Insitu_variable, zone)

write_csv(validation_stats, "output/STATS/validation_error_bias_summary.csv")
print(validation_stats, n = Inf)
