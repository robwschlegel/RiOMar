# func/run_figure_1.R
# Standalone Rscript entry point for Figure_1() (func/figure.R), run as a
# separate process by func/figure.py::Figure_1() instead of via in-process
# rpy2. Figure_1() is the only figure that needs the 'sf' package (for the
# high-resolution France coastline used in its zoomed insets), and 'sf'
# links against the system's GDAL/GEOS/PROJ, which conflict with the
# GDAL/GEOS/PROJ conda's Python geospatial stack (shapely/geopandas/pyproj,
# loaded in-process via panache) already has loaded by the time Figure_1()
# runs -- two incompatible copies of GEOS/PROJ/HDF5 in one process crashes.
# Running R as its own process sidesteps this: separate process, separate
# library address space, no conflict possible.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript func/run_figure_1.R <where_to_save_the_figure>")
}
where_to_save_the_figure <- args[[1]]

source("func/figure.R")
Figure_1(where_to_save_the_figure = where_to_save_the_figure)
