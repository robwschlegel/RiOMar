# func/generate_driver_x11_comparison_figures.R
#
# Rollback (2026-09-23) of the Aug-11 driver-vs-plume X11 migration's Python
# plotting -- replaces func/compute_driver_x11_figures.py (now a stub, kept
# for git history, matching the func/plume.py/func/validate.py convention)
# as the source of figures/driver_x11_comparison/Figure6_style_<driver>.png,
# one per driver (wind/tide/wave/current), each a 4-zone stacked dual-axis
# comparison of plume area vs. the driver's own X11 Interannual signal, with
# a Pearson-r annotation per panel -- same output path/filename convention
# as the script it replaces, so nothing else needs to change to keep finding
# these PNGs.
#
# Python stays limited to the X11 *calculation* (func/compute_x11_driver_signals.py,
# which persists each driver's weekly Seasonal_signal/Interannual_signal to
# output/panache/dynamic/<Zone>/X11_ANALYSIS/<driver>/<driver>_WEEKLY.csv --
# run that first, or here automatically if missing); this script and
# func/X11.R::plot_driver_x11_trend_comparison() do all the plotting, in R,
# matching func/X11.R::make_the_plot()'s house style for river flow's own
# Fig 6/7 rather than matplotlib's.
#
# Run from repo root: Rscript func/generate_driver_x11_comparison_figures.R

source("func/X11.R")  # sources multi.R (and, transitively, util.R) itself

OTHER_DRIVERS <- c("wind", "tide", "wave", "current")
X11_DIR <- "output/panache/dynamic"
OUT_DIR <- "figures/driver_x11_comparison"

driver_csv_missing <- function(driver_name){
  !all(file.exists(file.path(X11_DIR, zones, "X11_ANALYSIS", driver_name, paste0(driver_name, "_WEEKLY.csv"))))
}

if (any(vapply(OTHER_DRIVERS, driver_csv_missing, logical(1)))) {
  system("python func/compute_x11_driver_signals.py")
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (driver_name in OTHER_DRIVERS) {
  message("Building Figure-6-style stacked comparison for ", driver_name, " (Census X11, weekly)...")

  zone_plots <- purrr::map(zones, plot_driver_x11_trend_comparison,
                           where_are_saved_X11_results = X11_DIR, driver_name = driver_name)
  composite <- ggpubr::ggarrange(plotlist = zone_plots, ncol = 1, nrow = length(zones), align = "v")

  out_file <- file.path(OUT_DIR, paste0("Figure6_style_", driver_name, ".png"))
  ggplot2::ggsave(out_file, composite, width = 10, height = 12, dpi = 200)
  message("Wrote ", out_file)
}
