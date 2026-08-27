#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# func/compute_driver_x11_figures.py
#
# Python replacement for the X11 dual-axis comparison figures in Section 1
# of func/compute_driver_x11_figures.R (x11_dual_axis_df/plot_x11_dual_axis/
# the driver loop) -- one Figure-6-style stacked (4 zones) comparison PNG
# per driver (wind/tide/wave/current), using func/X11.py's weekly
# Census-Pezzulli decomposition instead of R's monthly-only
# seasonal::seas()/X-13ARIMA-SEATS. See compute_driver_x11_correlation_table.py
# for the matching correlation-table replacement.
#
# Section 2 of compute_driver_x11_figures.R (plot_by_category, raw daily
# wind/wave category-split scatter plots) is NOT X11 decomposition and is
# untouched -- run separately via Rscript func/compute_driver_x11_figures.R
# as before.
#
# Reads the daily plume+driver CSVs func/export_driver_x11_inputs.R writes
# (run automatically here if missing).
#
# Run from repo root: python func/compute_driver_x11_figures.py

import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

proj_dir = os.path.dirname(os.path.abspath('__file__'))
func_dir = os.path.join(proj_dir, 'func')
sys.path.append(func_dir)

import X11 as X11mod  # noqa: E402
from util import zone_title  # noqa: E402

ZONES = ["BAY_OF_SEINE", "SOUTHERN_BRITTANY", "BAY_OF_BISCAY", "GULF_OF_LION"]
OTHER_DRIVERS = ["wind", "tide", "wave", "current"]
INPUT_DIR = os.path.join(proj_dir, "output", "STATS", "driver_x11_inputs")
OUT_DIR = os.path.join(proj_dir, "figures", "driver_x11_comparison")

# Same driver_name -> (label, colour) mapping as func/multi.R's driver_display
# tribble, restricted to the 4 drivers used here.
DRIVER_DISPLAY = {
    "wind":    ("Wind speed (m s⁻¹)",     "purple"),
    "tide":    ("Tidal range (m)",        "darkgreen"),
    "wave":    ("Wave height (m)",        "steelblue"),
    "current": ("Current speed (m s⁻¹)",  "orchid"),
}


def ensure_daily_inputs_exist():
    if os.path.isdir(INPUT_DIR) and len(os.listdir(INPUT_DIR)) > 0:
        return
    subprocess.run(["Rscript", os.path.join(func_dir, "export_driver_x11_inputs.R")],
                   cwd=proj_dir, check=True)


def build_trend_df(zone, driver_name):
    df = pd.read_csv(os.path.join(INPUT_DIR, f"{zone}_{driver_name}.csv"))

    plume_binned = X11mod.bin_to_pseudo_weekly(df["date"], df["plume_area"])
    driver_binned = X11mod.bin_to_pseudo_weekly(df["date"], df["value"])

    plume_res = X11mod.decompose_driver_series(plume_binned["value"].tolist(), plume_binned["date"])
    driver_res = X11mod.decompose_driver_series(driver_binned["value"].tolist(), driver_binned["date"])
    if plume_res is None or driver_res is None:
        return None

    return plume_res.rename(columns={"Interannual_signal": "plume_trend"}).merge(
        driver_res.rename(columns={"Interannual_signal": "driver_trend"}), on="date")


if __name__ == "__main__":
    ensure_daily_inputs_exist()
    os.makedirs(OUT_DIR, exist_ok=True)

    for driver_name in OTHER_DRIVERS:
        driver_label, driver_colour = DRIVER_DISPLAY[driver_name]
        print(f"Building Figure-6-style stacked comparison for {driver_name} (Census X11, weekly)...")

        fig, axs = plt.subplots(nrows=len(ZONES), ncols=1, figsize=(10, 12))
        for ax_row, zone in zip(axs, ZONES):
            df = build_trend_df(zone, driver_name)
            if df is None:
                ax_row.set_title(f"{zone_title(zone)} (insufficient data)")
                continue

            r_value = df["plume_trend"].corr(df["driver_trend"])
            ax_row.plot(df["date"], df["plume_trend"], color="brown", marker="o", markersize=2)
            ax_row.set_ylabel("Plume area (km²)", color="brown")
            ax_row.tick_params(axis="y", labelcolor="brown")

            ax2 = ax_row.twinx()
            ax2.plot(df["date"], df["driver_trend"], color=driver_colour, marker="o", markersize=2)
            ax2.set_ylabel(driver_label, color=driver_colour)
            ax2.tick_params(axis="y", labelcolor=driver_colour)

            ax_row.set_title(zone_title(zone))
            ax_row.text(0.02, 0.9, f"r = {r_value:.2f}", transform=ax_row.transAxes,
                       ha="left", va="top", fontsize=11, color="black")

        fig.tight_layout()
        out_file = os.path.join(OUT_DIR, f"Figure6_style_{driver_name}.png")
        fig.savefig(out_file, dpi=200)
        plt.close(fig)
        print(f"Wrote {out_file}")
