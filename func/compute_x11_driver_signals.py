#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# func/compute_x11_driver_signals.py
#
# Persists the weekly Census-X11 Seasonal and Interannual components for
# wind/tide/wave/current to disk, closing a gap left by the 2026-08-11
# migration off monthly R/X-13ARIMA-SEATS (which rejects weekly-frequency
# input) onto Python's weekly Census-Pezzulli decomposition: only the
# correlation-table and plotting scripts were ported forward, not a
# persisted-CSV step, so unlike river flow (func/X11.py's
# apply_X11_method_and_save_results(), called with variable_name="river_flow"
# from Apply_X11_method_on_time_series()) these four drivers' decomposed
# series have only ever existed in-memory.
#
# Deliberately does NOT edit func/X11.py (frozen, per CLAUDE.md) -- calls
# temporal_decomp_V2_7_x11() directly, the same core function every existing
# X11.py wrapper (decompose_driver_series(), apply_X11_method_and_save_results())
# already calls, rather than routing through decompose_driver_series()
# (which only returns Interannual_signal) or apply_X11_method_and_save_results()
# (which expects a plume-processing "info" row from
# get_all_cases_to_process_for_regional_maps_or_plumes_or_X11() -- machinery
# this standalone script has no reason to depend on).
#
# Reads the daily plume+driver CSVs func/export_driver_x11_inputs.R writes
# (run automatically here if missing).
# Writes one CSV per zone x driver:
#   output/panache/dynamic/<Zone>/X11_ANALYSIS/<driver>/<driver>_WEEKLY.csv
#   (columns: date, Raw_signal, Seasonal_signal, Interannual_signal, Residual_signal)
# -- same directory convention and column names river flow's own X11 output
# already uses, so func/generate_x11_driver_correlation_heatmap.R can read
# all 5 drivers (flow + these 4) the same way.
#
# A driver/zone combination that fails X11's internal data-quality cutoff
# (flag == -100) is skipped (nothing written for it, a warning is printed)
# rather than writing an all-NaN file -- the R side checks file.exists()
# and renders a missing combination as a grey "failed" cell rather than
# erroring.
#
# Run from repo root: python func/compute_x11_driver_signals.py

import os
import subprocess
import sys

import pandas as pd

proj_dir = os.path.dirname(os.path.abspath('__file__'))
func_dir = os.path.join(proj_dir, 'func')
sys.path.append(func_dir)

import X11 as X11mod  # noqa: E402

ZONES = ["BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"]
OTHER_DRIVERS = ["wind", "tide", "wave", "current"]
INPUT_DIR = os.path.join(proj_dir, "output", "STATS", "driver_x11_inputs")
X11_DIR_OUT = os.path.join(proj_dir, "output", "panache", "dynamic")


def ensure_daily_inputs_exist():
    if os.path.isdir(INPUT_DIR) and len(os.listdir(INPUT_DIR)) > 0:
        return
    subprocess.run(["Rscript", os.path.join(func_dir, "export_driver_x11_inputs.R")],
                   cwd=proj_dir, check=True)


def compute_and_write(zone, driver_name):
    df = pd.read_csv(os.path.join(INPUT_DIR, f"{zone}_{driver_name}.csv"))
    binned = X11mod.bin_to_pseudo_weekly(df["date"], df["value"])

    results = X11mod.temporal_decomp_V2_7_x11(
        values=binned["value"].tolist(), dates=pd.to_datetime(binned["date"]),
        time_frequency="WEEKLY", filter_outlier=False, overall_cutoff=50,
        out_limit=3, perc_month_limit=50, var_stationary=False,
        lin_interpol=False, cutoff_fill=30, season_test=True)

    if results['24_flag'] == -100:
        print(f"[SKIP] {zone} / {driver_name}: X11 decomposition failed quality cutoff")
        return

    out = pd.DataFrame({'date': results['7_dates'],
                        'Raw_signal': results['8_values_ini'],
                        'Seasonal_signal': results['9_Seasonal_signal'],
                        'Interannual_signal': results['10_Interannual_signal'],
                        'Residual_signal': results['11_Residual_signal']})

    out_dir = os.path.join(X11_DIR_OUT, zone, 'X11_ANALYSIS', driver_name)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{driver_name}_WEEKLY.csv")
    out.to_csv(out_file, index=False)
    print(f"Wrote {out_file}")


if __name__ == "__main__":
    ensure_daily_inputs_exist()
    for zone in ZONES:
        for driver_name in OTHER_DRIVERS:
            compute_and_write(zone, driver_name)
