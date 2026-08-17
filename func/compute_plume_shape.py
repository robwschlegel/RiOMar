"""
One-off: derive a daily plume-shape (compactness) time series per zone from
the panache PlumeMasks.nc daily masks, since no shape metric currently exists
anywhere in the pipeline (see manuscript/manuscript.tex Table 5, "Shape: TBD
metric"). Compactness = 4*pi*area / perimeter^2 (isoperimetric ratio; 1 for a
circle, lower for elongated/irregular shapes). Perimeter is estimated with
skimage's Crofton estimator (more accurate than naive edge-counting for
digitised shapes) and converted to km using the grid's near-isotropic pixel
size (dx/dy differ by <10% at these latitudes -- checked per zone before
writing this).

Run from repo root with the RiOMar conda env:
  /home/calanus/anaconda3/envs/RiOMar/bin/python3 func/compute_plume_shape.py

Covers both the dynamic and static thresholds (added 2026-08-07, per Robert,
so the seasonal-analysis Figure S3 dynamic-vs-static comparison can include
compactness -- previously only the dynamic threshold had been run).
"""
import numpy as np
import pandas as pd
import xarray as xr
from skimage.measure import perimeter_crofton

ZONES = ["BAY_OF_SEINE", "BAY_OF_BISCAY", "SOUTHERN_BRITTANY", "GULF_OF_LION"]
PLUME_DIRS = ["output/panache/dynamic", "output/panache/static"]

for PLUME_DIR in PLUME_DIRS:
    for zone in ZONES:
        ds = xr.open_dataset(f"{PLUME_DIR}/{zone}/PlumeMasks.nc")

        dlon = float(np.diff(ds.lon.values).mean())
        dlat = float(np.diff(ds.lat.values).mean())
        mean_lat = float(ds.lat.values.mean())
        dx_km = dlon * 111.32 * np.cos(np.radians(mean_lat))
        dy_km = dlat * 111.32
        pixel_area_km2 = dx_km * dy_km
        pixel_size_km = float(np.mean([dx_km, dy_km]))  # near-isotropic (checked: ratio 0.97-1.10 across zones)

        n_time = ds.sizes["time"]
        dates = pd.to_datetime(ds.time.values)
        n_pixel = np.full(n_time, np.nan)
        area_km2 = np.full(n_time, np.nan)
        perimeter_km = np.full(n_time, np.nan)
        compactness = np.full(n_time, np.nan)

        # panache v5.0.0+ adds a 'river' dim (one mask per river mouth plus
        # a combined 'ALL' union layer) -- keep only the zone-level union.
        masks = ds.plume_mask.sel(river="ALL").values.astype(bool)
        for i in range(n_time):
            mask = masks[i]
            n = int(mask.sum())
            n_pixel[i] = n
            if n < 4:
                continue
            a_km2 = n * pixel_area_km2
            p_px = perimeter_crofton(mask, directions=4)
            p_km = p_px * pixel_size_km
            area_km2[i] = a_km2
            perimeter_km[i] = p_km
            if p_km > 0:
                compactness[i] = 4 * np.pi * a_km2 / (p_km ** 2)
            if (i + 1) % 2000 == 0:
                print(f"  {PLUME_DIR}/{zone}: {i + 1}/{n_time}", flush=True)

        out = pd.DataFrame({
            "date": dates.date,
            "n_pixel": n_pixel,
            "area_km2": area_km2,
            "perimeter_km": perimeter_km,
            "compactness": compactness,
        })
        out_path = f"{PLUME_DIR}/{zone}/PlumeShape.csv"
        out.to_csv(out_path, index=False)
        print(f"{PLUME_DIR}/{zone}: wrote {out_path} "
              f"(n_valid_days={out['compactness'].notna().sum()}, "
              f"mean_compactness={out['compactness'].mean():.4f})", flush=True)
        ds.close()
