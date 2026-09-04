# One-off: compute satellite-vs-in-situ percent error and percent bias for
# the manuscript's "validation_summary_stats" table (see
# manuscript/figure_table_registry.csv for its current table number -- Table
# S1 as of the validation table/figure's move to the supplement), from the
# raw SOMLIT/REPHY match-up file produced by func/validate.R
#
# Covers both in-situ variables in that file: SPM (SOMLIT) and turbidity/TURB
# (REPHY) -- REPHY does not measure SPM at all, only turbidity, which is why
# a separate turbidity-vs-satellite-SPM comparison is needed to make any use
# of the REPHY match-ups (see Section 2.3.1's discussion of Doxaran et al.
# 2009, who calibrated in-situ turbidity sensors to give SPM concentration in
# the Gironde, the precedent for treating the two as comparable here).
#
# Definitions (following standard ocean-colour validation practice, e.g.
# Pahlevan et al. 2021):
# NB: These equations are not exactly correct.
#   percent error = median(|satellite - in situ| / in situ) x 100  (MAPE)
#   percent bias  = median((satellite - in situ) / in situ) x 100
#
# Run from repo root: Rscript func/compute_validation_stats.R
library(dplyr)
library(readr)

zone_display <- c(
  "BAY_OF_SEINE"      = "Bay of Seine",
  "SOUTHERN_BRITTANY" = "Southern Brittany",
  "BAY_OF_BISCAY"     = "Gironde shelf",
  "GULF_OF_LION"       = "Rhône shelf"
)

zone_table <- read_csv("output/MATCH_UP_DATA/FRANCE/STATISTICS/table_all.csv", show_col_types = FALSE) |> 
       filter(sensor == "SEXTANT", zone != "GLOBAL", season == "ALL", site == "ALL", source == "ALL") |> 
       rename(ZONE = zone) |> 
       mutate(zone = zone_display[ZONE], .after = "ZONE") |> 
       arrange(variable, zone)

write_csv(zone_table, "output/STATS/validation_error_bias_summary.csv")
print(zone_table, n = Inf)
