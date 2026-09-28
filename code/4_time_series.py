#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to run the time series analyses.


# =============================================================================
#### Modules
# =============================================================================

import os
import sys
import subprocess
import concurrent.futures
import matplotlib as mpl
import rpy2.robjects as robjects

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import config

from X11 import Apply_X11_method_on_time_series, Apply_X11_method_on_time_series_per_river

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')

# The zones for mapping
zones_list = config.zones()

# Basic arguments to be used throughout the script
sextant_spm_all = config.satellite_dict('SPM')


# =============================================================================
# ### X11 analyses — dynamic threshold (main results)
# =============================================================================

# NB: X11 can only be used on weekly or monthly data, not daily
Apply_X11_method_on_time_series(sextant_spm_all,
                                Zones = zones_list,
                                plume_time_step = "WEEKLY",
                                plume_dir_in = "output/panache/dynamic",
                                X11_dir_out = "output/panache/dynamic",
                                include_river_flow = True)

# Per-river variant (metadata/river_discharge_mapping.csv): now that panache
# v5.0.0+ tracks individual river-mouth plume time series, also run the same
# X11 plume-vs-flow comparison once per river instead of only per zone.
Apply_X11_method_on_time_series_per_river(sextant_spm_all,
                                          plume_time_step = "WEEKLY",
                                          plume_dir_in = "output/panache/dynamic",
                                          X11_dir_out = "output/panache/dynamic")


# =============================================================================
# ### X11 analyses — static threshold (supplementary)
# =============================================================================

Apply_X11_method_on_time_series(sextant_spm_all,
                                Zones = zones_list,
                                plume_time_step = "WEEKLY",
                                plume_dir_in = "output/panache/static",
                                X11_dir_out = "output/panache/static",
                                include_river_flow = True)

Apply_X11_method_on_time_series_per_river(sextant_spm_all,
                                          plume_time_step = "WEEKLY",
                                          plume_dir_in = "output/panache/static",
                                          X11_dir_out = "output/panache/static")


# =============================================================================
# ### Multi-driver interaction analysis (GLM / GAM / RF, both thresholds)
# =============================================================================

# Run via a Rscript subprocess rather than rpy2's embedded R: ranger (loaded
# by driver_interactions.R) initialises its own OpenMP runtime, which
# collides with the one numpy/scipy already loaded into this Python process
# (macOS-only "OMP: Error #15: libomp.dylib already initialized" abort) --
# a separate process keeps the two OpenMP runtimes apart.
driver_interactions_R_path = os.path.join(func_dir, 'driver_interactions.R')

# Runs the dynamic-threshold main analysis; see func/driver_interactions.R
# subprocess.run(
#     ['Rscript', '-e', f"source('{driver_interactions_R_path}'); run_driver_interactions_analysis()"],
#     cwd=proj_dir, check=True
# )


# =============================================================================
# ### Monthly multi-driver interaction analysis (sec:seasonal_methods)
# =============================================================================

# Re-runs the same six-step GLM/GAM/RF sequence above independently within
# each calendar month's data subset, dynamic threshold only. Feeds the
# Supplementary monthly driver-dominance table (manuscript.tex).
#
# The 12 months are independent, so they're dispatched as separate Rscript
# subprocesses running concurrently (each via run_monthly_driver_interactions_analysis_for_month(),
# with ranger pinned to 1 thread there -- see func/driver_interactions.R for
# why) rather than looping over them sequentially in one process. A
# ThreadPoolExecutor is enough here, not multiprocess: each worker just
# launches and waits on its own Rscript subprocess, no Python-level
# computation happens in the workers themselves.
nb_of_cores_to_use = max(1, os.cpu_count() - 2)

def _run_month_driver_interactions(m):
    subprocess.run(
        ['Rscript', '-e',
         f"source('{driver_interactions_R_path}'); run_monthly_driver_interactions_analysis_for_month({m})"],
        cwd=proj_dir, check=True
    )

with concurrent.futures.ThreadPoolExecutor(max_workers=min(12, nb_of_cores_to_use)) as executor:
    list(executor.map(_run_month_driver_interactions, range(1, 13)))

# Aggregates the 12 months' output into the compact Supplementary table,
# now that every month above has finished.
subprocess.run(
    ['Rscript', '-e', f"source('{driver_interactions_R_path}'); write_monthly_driver_dominance_summary()"],
    cwd=proj_dir, check=True
)


# =============================================================================
# ### Plume shape (compactness), both thresholds
# =============================================================================

# Derives PlumeShape.csv per zone/threshold from panache's PlumeMasks.nc (see
# func/compute_plume_shape.py); read by func/figure.R's compactness panels
# and func/compute_shape_alongcoast_trend.R.
import compute_plume_shape  # noqa: F401


# =============================================================================
# ### Manuscript stats scripts
# =============================================================================

# These are one-off scripts that generate stats used in the manuscript; 
# wiring them in here keeps their outputs (output/STATS/*.csv) in sync 
# with the panache/X11/ driver-interactions data above. 
# Order matters: generate_monthly_trend_pct_heatmap.R and 
# generate_table_s_monthly_trends.R read compute_seasonal_trend.R's output; 
# generate_table_s_octant_trends.R reads compute_direction_octant_trend.R's output.
# compute_driver_correlation_matrices.R sources driver_interactions.R (for
# its zone/driver helpers), which loads ranger -- run it via Rscript
# subprocess rather than rpy2's embedded R for the same OpenMP-collision
# reason as driver_interactions.R itself, above. The other 9 don't load
# ranger, so they stay on rpy2.
stats_scripts = [
    'compute_area_trend.R',
    'compute_mass_spm_trend.R',
    'compute_shape_alongcoast_trend.R',
    'compute_driver_correlation_trend.R',
    'compute_seasonal_trend.R',
    'generate_monthly_trend_pct_heatmap.R',
    'generate_table_s_monthly_trends.R',
    'compute_direction_octant_trend.R',
    'generate_table_s_octant_trends.R',
]
for script in stats_scripts:
    robjects.r['source'](os.path.join(func_dir, script))

subprocess.run(
    ['Rscript', os.path.join(func_dir, 'compute_driver_correlation_matrices.R')],
    cwd=proj_dir, check=True
)

