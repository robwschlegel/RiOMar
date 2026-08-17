# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

RiOMar is a scientific data-processing pipeline for analysing river plumes along the French coast. It downloads, validates, and analyses satellite-derived chlorophyll a (Chl a) and suspended particulate matter (SPM) data alongside river discharge, wind, tide, wave, and ocean current data. The primary output is publication-quality figures and processed datasets. This repository is a rework of an earlier codebase by Louis Terrats (the original module is at [myRIOMAR_dev](https://github.com/louis-terrats/myRIOMAR_dev)).

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
python code/4_time_series.py     # X11 decomposition + driver comparisons + monthly seasonal driver analysis (~2h, see below)
python code/5_figures.py         # Publication figures
```

`code/4_time_series.py`'s multi-driver stage now includes `func/driver_interactions.R::run_monthly_driver_interactions_analysis()`, which reruns the full six-step GLM/GAM/regime/RF sequence independently for each of the 12 calendar months (dynamic threshold only) — this is the slow part of the stage (~90 min on a 16-core machine), not the annual dynamic/static passes that used to be all this script did.

All scripts prepend `func/` to `sys.path` by setting `proj_dir` from `os.path.abspath('__file__')` — they must be executed from the repo root, not from inside `code/`.

## Pipeline map (living document)

A browsable map of which script produces which figure/table, where each output lands, and what's currently a known gap is maintained at https://claude.ai/code/artifact/fd8a00ad-f109-48cd-94c1-2cca329f657f. It is a living document, not a one-time snapshot — update it in place (same URL) whenever the pipeline's wiring changes (a figure function moves, an output path changes, a `\figplaceholder` gets repointed, a gap gets fixed or a new one is found).

## Roadmap to co-author review (living document)

The tracker of remaining manuscript work before it goes to co-authors — open decisions, a "what's left" Gantt, task detail, and a completed-work log — is published at https://claude.ai/code/artifact/04edb9a0-c5b2-48e3-a446-84f5d5e49b45. When asked to update "the roadmap," update this artifact in place (same URL). Resolve an open decision by moving it out of the decisions section into the completed-work log with today's date; add newly-discovered work the same way other entries are dated and grouped.

## French river plume literature review (gaps document)

When considering weaknesses or gaps in the current RiOMar plume-analysis methodology (e.g. for Discussion/Conclusion writing, or before proposing a methodological change), consult [manuscript/french_plume_literature_review.md](manuscript/french_plume_literature_review.md). It reviews every French river plume study (Seine, Loire, Gironde, Rhône, plus the Adour as an out-of-sample comparator) in `manuscript/references.bib`, grouped by zone, and synthesises seven recurring gap themes (no dynamical/process model, exclusion of the near-mouth/turbidity-maximum zone, no sub-daily/tidal-phase resolution, surface-only detection, plume detachment invisible to a threshold-and-flood-fill detector, ROFI used only as a static check, no compositional/biogeochemical SPM breakdown).

A live, formatted version of the same content is published at https://claude.ai/code/artifact/a282bc0a-b486-4e88-897c-76d7a94db34b. Like the pipeline map above, this is a living document: if RiOMar's methodology changes in a way that closes, changes, or adds to one of the gaps it identifies, update both the `.md` file and that artifact URL together, in place.

This document is a strong source of material for the manuscript's Discussion and Conclusion sections (which gaps are worth naming as limitations, which are natural future work).

## Target journal

The manuscript (`manuscript/manuscript.tex`) targets *Remote Sensing of Environment* (Elsevier). Guide for authors: https://www.sciencedirect.com/journal/remote-sensing-of-environment/publish/guide-for-authors. Consult it (or fetch it fresh, since Elsevier updates these pages) before drafting end-matter/boilerplate sections (CRediT, Declaration of competing interest, Data availability, Declaration of Generative AI use, Funding, Acknowledgements) or checking structural requirements (word limits, reference style, appendix numbering). Note: Research Articles are capped at 15,000 words including references and figure captions (Review Articles get 20,000) — the 20,000-word figure in `manuscript.tex`'s header comment is the wrong article type's limit and should be corrected.

## Google Doc sync (co-author review copy)

`manuscript.tex` is mirrored into a Google Doc for co-author comments/editing, in place (same URL, same open comment threads survive each push). To push the current `.tex` to the Doc:

```bash
manuscript/google_doc_sync/sync.sh
```

~8s to build the docx (figures downscaled to JPEG, citations resolved via `references.bib`, section/figure/table numbering reconstructed since pandoc doesn't do this on its own — see script comments) + ~30-40s for Drive's server-side conversion. Full mechanism, setup, and known limitations (payload-size ceiling on the relay, the one cosmetic equation-rendering gap) are in [manuscript/google_doc_sync/README.md](manuscript/google_doc_sync/README.md). Credentials live in `manuscript/google_doc_sync/.env` (gitignored, not in this repo's history — read it directly rather than asking the user to re-paste the token).

## Data storage

Large datasets are stored **outside** this repo under `~/pCloudDrive/data/` and are never committed. The `.gitignore` also excludes most of `output/` and `data/SEXTANT`, `data/INSITU_data`, etc. Only shapefiles, metadata CSVs, and zone config JSONs are tracked.

## Architecture

### Study zones
Four coastal zones used throughout: `GULF_OF_LION`, `BAY_OF_SEINE`, `BAY_OF_BISCAY`, `SOUTHERN_BRITTANY`.

### func/ modules (Python)
- [func/dl.py](func/dl.py) — FTP/CMEMS downloads via `copernicusmarine`; `Download_satellite_data`, `download_cmems_subset`, `daily_integral`
- [func/util.py](func/util.py) — shared helpers: file discovery, parameter parsing, zone coordinates, path templating
- [func/validate.py](func/validate.py) — placeholder, kept for git history; satellite vs in situ match-up now lives entirely in `func/validate.R` (see below)
- [func/regmap.py](func/regmap.py) — regional map creation and QC (`create_regional_maps`, `QC_of_regional_maps`)
- [func/plume.py](func/plume.py) — placeholder, kept for git history; plume detection is now handled entirely by the external `panache` package (`panache.plume_algorithm`, `panache.utils`); the one RiOMar-specific figure-prep helper this file used to hold (`preprocess_annual_dataset_and_compute_land_mask`) was removed from `func/figure.py` in favour of `panache.plume_algorithm.derive_masks_from_bathymetry`
- [func/X11.py](func/X11.py) — X11 seasonal decomposition, calls R via `rpy2`; `Apply_X11_method_on_time_series`. **Frozen**: this file is intentionally excluded from renaming/refactor passes — do not edit its existing functions, even for naming consistency. New wrapper functions may be added alongside them. `func/X11.R` is not frozen.
- [func/figure.py](func/figure.py) — all publication figures (`Figure_1`, `Figure_2`, `Figure_3`/`Figure_3_panels`/`Figure_3_zone_maps`, `Figure_4_S1_timeseries`, `Figure_5_seasonal_analysis`, `Figure_S_daily_flow`, `Figure_X11_weekly_results`, `Figure_7_driver_rose`, `Figure_8_gam_partial`, `Figure_S3_seasonal_boxplots`). `Figure_8_driver_category` is deprecated (uncalled, kept under a "Deprecated functions" heading in both `figure.py` and `figure.R` for git history)

### func/ modules (R)
Parallel R implementations exist for most modules (`util.R`, `validate.R`, `X11.R`, etc. — there is no `regmap.R`, regional-map creation is Python-only). These are used for analyses that rely on R packages (e.g. base stats and all plotting) and are called from Python via `rpy2`. `func/validate.R` is the authoritative satellite-vs-in-situ match-up pipeline (writes both `output/MATCH_UP_DATA/FRANCE/summary.csv`, feeding manuscript Table 4, and the SEXTANT/ODATIS-MR `STATISTICS/*.csv` tables feeding Figure 2).

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
`dl.py` and `regmap.py` use `multiprocess` (not the stdlib `multiprocessing`). The start method is forced to `'spawn'` for macOS compatibility — do not change this.

## Future work ideas

- An SST (sea surface temperature) analysis of river plumes — how the plume footprint appears in a high-resolution SST product — is still an idea worth pursuing, not yet scoped or started. A stub subsection for this in `manuscript.tex` (§ Results, "Sea surface temperature") was removed 2026-08-12 since it held no content; revisit as a possible future-work mention or a follow-up study, not a manuscript gap to fill now.

## Bug history

Bug-fix history (root cause, symptom, before/after) is intentionally not kept here — it bloats a file loaded every session regardless of relevance. Each fix is commented in place at its own function; check `git log`/`git blame` on the relevant file for the full story.

