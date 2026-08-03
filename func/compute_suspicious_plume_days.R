# One-off: run the suspicious-plume-day flagging pipeline (func/surface.R)
# over every zone and both threshold runs, and write the flagged subset to
# output/STATS/suspicious_plume_days.csv for manual review (and, eventually,
# masking upstream of the plume-area trend/driver analyses in func/multi.R).
# Run from repo root: Rscript func/compute_suspicious_plume_days.R
source("func/surface.R")

all_days <- flag_suspicious_days_all()

suspicious <- all_days |>
  dplyr::filter(suspicious) |>
  dplyr::select(zone, threshold, date, n_pixel_in_the_plume_area,
                area_of_the_plume_mask_in_km2, mean_SPM_in_the_plume_area,
                confidence_index_in_perc, local_median, ratio) |>
  dplyr::arrange(zone, date)

readr::write_csv(suspicious, "output/STATS/suspicious_plume_days.csv")
message(nrow(suspicious), " suspicious plume-day(s) flagged; written to output/STATS/suspicious_plume_days.csv")
print(suspicious, n = Inf)

plot_suspicious_plume_days(all_days, "dynamic")
plot_suspicious_plume_days(all_days, "static")
