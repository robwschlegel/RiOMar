# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

RiOMar is a scientific data-processing pipeline for analysing river plumes along the French coast. It downloads, validates, and analyses satellite-derived chlorophyll a (Chl a) and suspended particulate matter (SPM) data alongside river discharge, wind, tide, and oceanographic driver data. The primary output is publication-quality figures and processed datasets. This repository is a Python-first rework of an earlier R codebase by Louis Terrats (the original module is at [myRIOMAR_dev](https://github.com/louis-terrats/myRIOMAR_dev)).

## Running the pipeline

Each numbered script in [code/](code/) corresponds to a pipeline stage. Run them sequentially from the repo root with the project Python environment active:

```bash
python code/0_download_data.py   # Download satellite + driver data (~hours, ~280 GB)
python code/1_validate.py        # Satellite vs in situ match-up
python code/2_regional_maps.py   # Create & QC regional maps
# Step 3: plume detection via the external panache module (must be called from terminal directly).
# Two config variants per zone -- dynamic (main results) and static (supplementary) -- 8 calls total:
panache metadata/zone_config_dynamic_GULF_OF_LION.json
panache metadata/zone_config_dynamic_BAY_OF_BISCAY.json
panache metadata/zone_config_dynamic_SOUTHERN_BRITTANY.json
panache metadata/zone_config_dynamic_BAY_OF_SEINE.json
panache metadata/zone_config_static_GULF_OF_LION.json
panache metadata/zone_config_static_BAY_OF_BISCAY.json
panache metadata/zone_config_static_SOUTHERN_BRITTANY.json
panache metadata/zone_config_static_BAY_OF_SEINE.json
python code/4_time_series.py     # X11 decomposition + driver comparisons
python code/5_figures.py         # Publication figures
```

All scripts prepend `func/` to `sys.path` by setting `proj_dir` from `os.path.abspath('__file__')` — they must be executed from the repo root, not from inside `code/`.

## Pipeline map (living document)

A browsable map of which script produces which figure/table, where each output lands, and what's currently a known gap is maintained at https://claude.ai/code/artifact/fd8a00ad-f109-48cd-94c1-2cca329f657f. It is a living document, not a one-time snapshot — update it in place (same URL) whenever the pipeline's wiring changes (a figure function moves, an output path changes, a `\figplaceholder` gets repointed, a gap gets fixed or a new one is found).

## Data storage

Large datasets are stored **outside** this repo under `~/pCloudDrive/data/` and are never committed. The `.gitignore` also excludes most of `output/` and `data/SEXTANT`, `data/INSITU_data`, etc. Only shapefiles, metadata CSVs, and zone config JSONs are tracked.

## Architecture

### Study zones
Four coastal zones used throughout: `GULF_OF_LION`, `BAY_OF_SEINE`, `BAY_OF_BISCAY`, `SOUTHERN_BRITTANY`.

### func/ modules (Python)
- [func/dl.py](func/dl.py) — FTP/CMEMS downloads via `copernicusmarine`; `Download_satellite_data`, `download_cmems_subset`, `daily_integral`
- [func/util.py](func/util.py) — shared helpers: file discovery, parameter parsing, zone coordinates, path templating
- [func/validate.py](func/validate.py) — satellite vs in situ match-up (`Match_up_with_insitu_measurements`)
- [func/regmap.py](func/regmap.py) — regional map creation and QC (`create_regional_maps`, `QC_of_regional_maps`)
- [func/plume.py](func/plume.py) — plume geometry and pixel extraction (supports the panache workflow)
- [func/X11.py](func/X11.py) — X11 seasonal decomposition, calls R via `rpy2`; `Apply_X11_method_on_time_series`
- [func/figure.py](func/figure.py) — all publication figures (`Figure_1` through `Figure_8_9_10`)

### func/ modules (R)
Parallel R implementations exist for most modules (`util.R`, `validate.R`, `regmap.R`, `X11.R`, etc.). These are used for analyses that rely on R packages (e.g. X11 seasonal decomposition) and are called from Python via `rpy2`.

### metadata/
Zone configuration JSONs consumed directly by `panache` and zone-pixel CSVs (one per sensor × variable × atmospheric correction combination) used for plume pixel extraction.

### Satellite data dict convention
A Python dict like:
```python
{'Data_sources': ['SEXTANT'], 'Sensor_names': ['merged'], 'Satellite_variables': ['SPM'],
 'Atmospheric_corrections': ['Standard'], 'Temporal_resolution': ['DAILY'],
 'start_day': '1998/01/01', 'end_day': '2025/12/31'}
```
is the standard argument passed to every major pipeline function. `util.define_parameters` converts it into a named-tuple `info` object.

### Multiprocessing
`dl.py`, `regmap.py`, `validate.py`, and `plume.py` all use `multiprocess` (not the stdlib `multiprocessing`). The start method is forced to `'spawn'` for macOS compatibility — do not change this.

### Rhône apportionment
The Rhône splits near Arles into Grand Rhône and Petit Rhône, each with its own gauge (Arles, Fourques) through 2022. From 2023 onward, both branches are extended from the upstream Tarascon gauge via a MOVE.1 log-log regression fit separately per branch on the full 1998–2022 overlap (`func/river_flow_prep.R`, `river_config`) — not a flat fraction of the combined discharge. (An earlier fixed-fraction rule — 11% to Petit Rhône before 2012-05-28, 10% thereafter — was superseded by this; it was never actually implemented in code and has been removed from the manuscript.)