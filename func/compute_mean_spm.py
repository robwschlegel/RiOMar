# func/compute_mean_spm.py
#
# One-off (but reusable/re-runnable): computes the pixel-wise temporal mean
# SPM field across the entire SEXTANT daily archive (1998-2025), for use as
# manuscript Figure 1's national main-panel background layer (replacing the
# single hardcoded snapshot date, 2011/02/02, that func/figure.py::Figure_1()
# used previously). Streams one daily file at a time rather than loading the
# ~280 GB archive into memory at once, accumulating a running sum and a
# valid-pixel count per grid cell (not every pixel has data on every day --
# cloud, land, or otherwise masked), then divides at the end.
#
# Benchmarked at ~0.56 s/file (20-file sample) -> ~95 minutes for the full
# ~10,100-file archive.
#
# Run from repo root: python func/compute_mean_spm.py

import os
import sys
import glob
import time

import numpy as np
import xarray as xr

proj_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

DATA_DIR = os.path.expanduser("~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY")
OUTPUT_PATH = os.path.join(proj_dir, "output", "STATS", "mean_spm_national.nc")


def national_mean_spm(data_dir=DATA_DIR, output_path=OUTPUT_PATH):
    files = sorted(glob.glob(os.path.join(data_dir, "*", "*", "*", "*.nc")))
    print(f"Found {len(files)} daily SPM files under {data_dir}")
    if len(files) == 0:
        raise FileNotFoundError(f"No .nc files found under {data_dir}")

    sum_spm = None
    count_valid = None
    lat = lon = None
    n_failed = 0

    t0 = time.time()
    for i, f in enumerate(files):
        try:
            with xr.open_dataset(f) as ds:
                spm = ds['analysed_spim'].isel(time=0).values
                if sum_spm is None:
                    sum_spm = np.zeros_like(spm, dtype=np.float64)
                    count_valid = np.zeros_like(spm, dtype=np.int64)
                    lat = ds['lat'].values.copy()
                    lon = ds['lon'].values.copy()
                valid = ~np.isnan(spm)
                sum_spm[valid] += spm[valid]
                count_valid[valid] += 1
        except Exception as e:
            n_failed += 1
            print(f"  [skip] {f}: {e}")

        if (i + 1) % 500 == 0 or (i + 1) == len(files):
            elapsed = time.time() - t0
            print(f"  processed {i + 1}/{len(files)} ({elapsed/60:.1f} min elapsed)")

    mean_spm = np.divide(sum_spm, count_valid, out=np.full_like(sum_spm, np.nan), where=count_valid > 0)

    out = xr.Dataset(
        {"mean_analysed_spim": (("lat", "lon"), mean_spm),
         "n_valid_days": (("lat", "lon"), count_valid)},
        coords={"lat": lat, "lon": lon},
        attrs={"description": "Pixel-wise temporal mean of SEXTANT daily analysed_spim, "
                              f"{len(files)} candidate days ({n_failed} unreadable, skipped), "
                              "for manuscript Figure 1's national SPM background layer.",
              "source_dir": data_dir}
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    out.to_netcdf(output_path)
    print(f"Wrote {output_path} ({n_failed} files skipped)")
    return out


if __name__ == "__main__":
    national_mean_spm()
