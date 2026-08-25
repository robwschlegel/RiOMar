# func/run_figure_3_zone_maps.R
# Standalone Rscript entry point for plot_methodology_zone_maps_panel() (func/figure.R),
# run as a separate process by func/figure.py::plot_methodology_zone_maps_panel() instead
# of via in-process rpy2 -- same sf/GEOS conflict and same workaround as
# run_figure_1.R (see that file's header for the full explanation).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript func/run_figure_3_zone_maps.R <where_to_save_the_figure_3>")
}
where_to_save_the_figure <- args[[1]]

source("func/figure.R")
plot_methodology_zone_maps_panel(where_to_save_the_figure = where_to_save_the_figure)
