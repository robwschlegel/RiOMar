# func/run_figure_3_panels.R
# Standalone Rscript entry point for plot_methodology_worked_example_panel() (func/figure.R), run as
# a separate process by func/figure.py::Figure_3_panels() instead of via
# in-process rpy2 -- same sf/GEOS conflict and same workaround as
# run_figure_1.R (see that file's header for the full explanation).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript func/run_figure_3_panels.R <where_to_save_the_figure_3>")
}
where_to_save_the_figure <- args[[1]]

source("func/figure.R")
for (letter in c("A", "B", "C", "D", "E")) {
  plot_methodology_worked_example_panel(where_to_save_the_figure = where_to_save_the_figure, name_of_the_plot = letter)
}

plot_methodology_transect_panel(where_to_save_the_figure = where_to_save_the_figure)
