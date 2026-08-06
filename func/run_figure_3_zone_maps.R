# func/run_figure_3_zone_maps.R
# Standalone Rscript entry point for Figure_3_zone_maps() (func/figure.R),
# run as a separate process by func/figure.py::Figure_3_zone_maps() instead
# of via in-process rpy2 -- same sf/GEOS conflict and same workaround as
# run_figure_1.R (see that file's header for the full explanation).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript func/run_figure_3_zone_maps.R <where_to_save_the_figure_3>")
}
where_to_save_the_figure <- args[[1]]

source("func/figure.R")
Figure_3_zone_maps(where_to_save_the_figure = where_to_save_the_figure)
