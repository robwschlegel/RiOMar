#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to create the figures used in the publication of this workflow.


# =============================================================================
#### Modules
# =============================================================================

import os, sys, subprocess
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, figure
from figure import (Figure_1_mean_spm_map, Figure_2_methodology_panels, Figure_2_methodology_zone_maps,
                    Figure_2_methodology, Figure_3_S2_timeseries, Figure_4_monthly_median_heatmap,
                    Figure_X11_weekly_results, Figure_8_driver_rose,
                    Figure_S1_validation, Figure_S8_gam_partial)

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')

# Function names carry the current manuscript figure number; slot names in
# the comments are the keys in metadata/figure_table_registry.csv, which
# remains the source of truth for numbers and output folders.
# Figure 5 (monthly_trend_pct_heatmap) is not built here -- it is written by
# func/analysis/generate_monthly_trend_pct_heatmap.R, run from code/4_time_series.py.


# =============================================================================
# ### Main-text figures
# =============================================================================

# study_zone_map (Fig. 1)
Figure_1_mean_spm_map(where_to_save_the_figure = "figures")

# plume_methodology_panel (Fig. 2): panels A-E, then the zone-maps panel,
# then the composite -- the composite must run after both
Figure_2_methodology_panels(where_are_saved_panache_outputs = "output",
                            where_to_save_the_figure = "figures")
Figure_2_methodology_zone_maps(where_are_saved_panache_outputs = "output",
                               where_to_save_the_figure = "figures")
Figure_2_methodology(where_to_save_the_figure = "figures")

# plume_area_timeseries (Fig. 3) + thresholds_comparison (Fig. S2)
Figure_3_S2_timeseries(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                       where_are_saved_plume_results_with_fixed_threshold = "output/panache/static",
                       where_to_save_the_figure = "figures")

# seasonal_boxplot_heatmap (Fig. 4): heatmap of each month's median plume
# property/driver value as a ratio to the zone's 1998-2025 median. Per-month
# trends are computed separately by func/analysis/compute_seasonal_trend.R.
Figure_4_monthly_median_heatmap(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
                                where_are_saved_plume_results_with_static_threshold = "output/panache/static",
                                where_to_save_the_figure = "figures")

# x11_seasonal_river_flow (Fig. 6), x11_interannual_river_flow (Fig. 7),
# x11_residual_river_flow (Fig. S6), and the dynamic-vs-static comparisons
# x11_seasonal (Fig. S3) / x11_interannual (Fig. S4) / x11_residual (Fig. S5)
Figure_X11_weekly_results(where_are_saved_X11_results_dynamic = "output/panache/dynamic",
                          where_are_saved_X11_results_static = "output/panache/static",
                          where_to_save_the_figure = "figures")

# driver_rose_diagram (Fig. 8): wind/wave/current direction roses
Figure_8_driver_rose(where_to_save_the_figure = "figures")


# =============================================================================
# ### Supplementary figures
# =============================================================================

# validation_scatterplot_panel (Fig. S1)
Figure_S1_validation(where_to_save_the_figure = "figures")

# x11_driver_correlation_heatmap (Fig. S7): Pearson r between plume area's
# X11 seasonal/interannual component and each driver's own component, per
# zone. Self-contained -- reads its own registry row for the output path,
# and regenerates output/panache/dynamic/<Zone>/X11_ANALYSIS/<driver>/*_WEEKLY.csv
# (func/analysis/compute_x11_driver_signals.py) itself if missing. A plain top-level
# script (unlike the other R figure code here, which defines functions
# sourced then called), so run via Rscript subprocess exactly as its own
# header documents -- was previously only ever run by hand, not wired into
# any code/ stage (see metadata/figure_table_registry.csv).
subprocess.run(
    ['Rscript', 'func/analysis/generate_x11_driver_correlation_heatmap.R'],
    cwd=proj_dir, check=True,
)

# gam_partial_effects (Fig. S8): GAM partial-dependence curves
Figure_S8_gam_partial(where_to_save_the_figure = "figures")

