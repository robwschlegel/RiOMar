#!/usr/bin/env python3
"""
Diagnostic plots for panache Results.csv outputs.

Produces, per zone:
  1. figures/qc/panache_timeseries_{zone}.png  — 4-panel stacked time series
     (plume area, mean SPM, lon centroid, lat centroid) for dynamic and static
     thresholds overlaid.
  2. figures/qc/panache_missing_{zone}.png     — monthly bar chart of missing days
     (days absent from Results.csv) for dynamic and static.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

zones = ["BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION"]
zone_labels = {
    "BAY_OF_SEINE":       "Bay of Seine (Seine)",
    "BAY_OF_BISCAY":      "Bay of Biscay (Gironde)",
    "SOUTHERN_BRITTANY":  "Southern Brittany (Loire)",
    "GULF_OF_LION":       "Gulf of Lion (Rhône)",
}

COLOURS = {"dynamic": "#1f77b4", "static": "#d62728"}
ALPHA   = {"dynamic": 0.8,       "static": 0.6}


def load_results(zone, threshold):
    path = os.path.join(proj_dir, "output", "panache", threshold, zone, "Results.csv")
    df = pd.read_csv(path, parse_dates=["date"])
    df["date"] = pd.to_datetime(df["date"]).dt.normalize()
    df = df.sort_values("date").reset_index(drop=True)
    # Treat area == 0 as no plume: set area and centroid columns to NaN so they
    # don't plot as spurious zeros.
    no_plume = df["area_of_the_plume_mask_in_km2"] == 0
    for col in ["area_of_the_plume_mask_in_km2", "mean_SPM_in_the_plume_area",
                "lon_weighted_centroid_of_the_plume_area",
                "lat_weighted_centroid_of_the_plume_area"]:
        if col in df.columns:
            df.loc[no_plume, col] = np.nan
    return df


def missing_by_month(df):
    """Return a Series indexed by Period('M') with count of missing days."""
    if df.empty:
        return pd.Series(dtype=int)
    full_range = pd.date_range(df["date"].min(), df["date"].max(), freq="D")
    present = set(df["date"])
    missing_dates = [d for d in full_range if d not in present]
    if not missing_dates:
        return pd.Series(dtype=int)
    return (pd.Series(1, index=pd.DatetimeIndex(missing_dates))
              .resample("MS").sum()
              .astype(int))


# ── 1. Time series plots ─────────────────────────────────────────────────────

metrics = [
    ("area_of_the_plume_mask_in_km2",           "Plume area (km²)"),
    ("mean_SPM_in_the_plume_area",               "Mean SPM (g m⁻³)"),
    ("lon_weighted_centroid_of_the_plume_area",  "Centroid longitude (°E)"),
    ("lat_weighted_centroid_of_the_plume_area",  "Centroid latitude (°N)"),
]

for zone in zones:
    fig, axes = plt.subplots(4, 1, figsize=(16, 12), sharex=True)
    fig.suptitle(zone_labels[zone], fontsize=13, fontweight="bold", y=0.995)

    data = {}
    for thresh in ["dynamic", "static"]:
        data[thresh] = load_results(zone, thresh)

    for ax, (col, ylabel) in zip(axes, metrics):
        for thresh in ["dynamic", "static"]:
            df = data[thresh]
            if col not in df.columns or df[col].isna().all():
                continue
            ax.plot(df["date"], df[col],
                    color=COLOURS[thresh], alpha=ALPHA[thresh],
                    linewidth=0.5, label=thresh.capitalize())

        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    # Legend on first panel only
    axes[0].legend(fontsize=8, frameon=False, loc="upper right")

    # x-axis: year ticks
    axes[-1].xaxis.set_major_locator(mdates.YearLocator(5))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes[-1].xaxis.set_minor_locator(mdates.YearLocator(1))
    plt.setp(axes[-1].xaxis.get_majorticklabels(), rotation=0, ha="center")

    fig.tight_layout(rect=[0, 0, 1, 0.995])
    out = os.path.join(proj_dir, "figures", "qc", f"panache_timeseries_{zone}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ── 2. Monthly missing-data bar plots ────────────────────────────────────────

for zone in zones:
    fig, axes = plt.subplots(2, 1, figsize=(18, 7), sharex=False)
    fig.suptitle(f"{zone_labels[zone]} — missing days per month", fontsize=12,
                 fontweight="bold")

    for ax, thresh in zip(axes, ["dynamic", "static"]):
        df = load_results(zone, thresh)
        miss = missing_by_month(df)

        if miss.empty:
            ax.text(0.5, 0.5, "No missing days", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
        else:
            ax.bar(miss.index, miss.values,
                   width=25, color=COLOURS[thresh], alpha=0.8)
            ax.set_ylabel("Missing days", fontsize=9)

        ax.set_title(thresh.capitalize(), fontsize=10, loc="left")
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

        # Annotate total
        total = int(miss.sum()) if not miss.empty else 0
        ax.text(0.99, 0.97, f"Total missing: {total:,} days",
                ha="right", va="top", transform=ax.transAxes,
                fontsize=8, color="grey")

    for ax in axes:
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_minor_locator(mdates.YearLocator(1))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=0, ha="center")

    fig.tight_layout()
    out = os.path.join(proj_dir, "figures", "qc", f"panache_missing_{zone}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
