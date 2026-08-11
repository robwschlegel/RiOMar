#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# func/compute_driver_x11_correlation_table.py
#
# Python replacement for the X11-decomposition rows in
# func/compute_driver_x11_correlation_table.R -- computes, per zone x driver
# (wind/tide/wave/current), the correlation between plume area's and the
# driver's own weekly Census-X11 Interannual (trend-cycle) signal, using
# func/X11.py's existing temporal_decomp_V2_7_x11 (via the new
# bin_to_pseudo_weekly()/decompose_driver_series() wrappers added there),
# rather than R's seasonal::seas()/X-13ARIMA-SEATS. Same family of method
# (Cleveland 1976 Census-X11), but not numerically identical to X-13 -- see
# CLAUDE.md and manuscript/TODO.md for why this replacement was made
# (X-13 rejects weekly-frequency input outright).
#
# The wind/wave category rows and the river-flow reference column are NOT
# X11 decomposition and are left computed by R (func/multi.R's existing,
# unmodified driver-loading machinery, plus the category-correlation helpers
# still in compute_driver_x11_correlation_table.R) -- this script merges
# those in via an Rscript subprocess rather than re-implementing them.
#
# Reads the daily plume+driver CSVs func/export_driver_x11_inputs.R writes
# (Rscript func/export_driver_x11_inputs.R must be run first, or is run
# here automatically if the input CSVs are missing).
#
# Run from repo root: python func/compute_driver_x11_correlation_table.py

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
CATEGORY_ROWS_CSV = os.path.join(proj_dir, "output", "STATS", "driver_x11_correlation_summary_category_rows.csv")


def ensure_daily_inputs_exist():
    if os.path.isdir(INPUT_DIR) and len(os.listdir(INPUT_DIR)) > 0:
        return
    subprocess.run(["Rscript", os.path.join(func_dir, "export_driver_x11_inputs.R")],
                   cwd=proj_dir, check=True)


def x11_trend_r(zone, driver_name):
    """
    Weekly Census-X11 Interannual-signal correlation for one zone x driver,
    the Python/X11.py equivalent of compute_driver_x11_correlation_table.R's
    x11_trend_r() (monthly, seasonal::seas()-based).
    """
    df = pd.read_csv(os.path.join(INPUT_DIR, f"{zone}_{driver_name}.csv"))

    plume_binned = X11mod.bin_to_pseudo_weekly(df["date"], df["plume_area"])
    driver_binned = X11mod.bin_to_pseudo_weekly(df["date"], df["value"])

    plume_res = X11mod.decompose_driver_series(plume_binned["value"].tolist(), plume_binned["date"])
    driver_res = X11mod.decompose_driver_series(driver_binned["value"].tolist(), driver_binned["date"])
    if plume_res is None or driver_res is None:
        return None

    merged = plume_res.rename(columns={"Interannual_signal": "plume_trend"}).merge(
        driver_res.rename(columns={"Interannual_signal": "driver_trend"}), on="date")
    return merged["plume_trend"].corr(merged["driver_trend"])


def compute_x11_rows():
    rows = []
    for zone in ZONES:
        for driver_name in OTHER_DRIVERS:
            r = x11_trend_r(zone, driver_name)
            rows.append({"zone": zone, "driver": driver_name, "subset": "all (X11 weekly trend)",
                        "r": round(r, 2) if r is not None else None, "lag_days": pd.NA})
    return pd.DataFrame(rows)


def compute_category_rows():
    """
    Runs the trimmed compute_driver_x11_correlation_table.R (category
    correlation rows + river-flow reference column only, no X11
    decomposition) as a subprocess, writing CATEGORY_ROWS_CSV, then reads
    it back in.
    """
    subprocess.run(["Rscript", os.path.join(func_dir, "compute_driver_x11_correlation_table.R")],
                   cwd=proj_dir, check=True)
    return pd.read_csv(CATEGORY_ROWS_CSV)


if __name__ == "__main__":
    ensure_daily_inputs_exist()

    print("Computing X11 weekly trend rows (Census-Pezzulli, Python)...")
    x11_rows = compute_x11_rows()

    print("Computing category + river-flow-reference rows (R, unchanged logic)...")
    category_rows = compute_category_rows()

    # river_flow_daily_r/river_flow_lag_days are the same reference column
    # for every row of a given zone -- carry it from the category rows
    # (which already computed it) onto the X11 rows too.
    flow_ref = category_rows[["zone", "river_flow_daily_r", "river_flow_lag_days"]].drop_duplicates("zone")
    x11_rows = x11_rows.drop(columns=["lag_days"]).merge(flow_ref, on="zone").assign(lag_days=pd.NA)
    x11_rows = x11_rows[["zone", "driver", "subset", "r", "lag_days", "river_flow_daily_r", "river_flow_lag_days"]]

    results = pd.concat([x11_rows, category_rows], ignore_index=True)
    results["zone"] = pd.Categorical(results["zone"], categories=ZONES, ordered=True)
    results = results.sort_values(["zone", "driver", "subset"]).reset_index(drop=True)

    out_file = os.path.join(proj_dir, "output", "STATS", "driver_x11_correlation_summary.csv")
    results.to_csv(out_file, index=False)
    print(f"Wrote {out_file}")
    print(results.to_string())
