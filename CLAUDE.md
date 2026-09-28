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
python code/3_plumes.py           # Plume detection via the external panache CLI (~60 min per zone x mode)
# Step 3 first writes one panache JSON per zone x threshold mode (dynamic = main results,
# static = supplementary) to output/panache/configs/ via tools/write_panache_configs.py,
# then calls `panache <json>` 8 times. To run a single one by hand:
#   python tools/write_panache_configs.py
#   panache output/panache/configs/zone_config_dynamic_GULF_OF_LION.json
python code/4_time_series.py     # X11 decomposition + driver comparisons + monthly seasonal driver analysis (~2h, see below)
python code/5_figures.py         # Publication figures
```

`code/4_time_series.py`'s multi-driver stage runs `func/driver_interactions.R::run_driver_interactions_analysis()` once (dynamic threshold): GLM comparison, per-metric GLM/GAM models and random-forest importance, with per-cell checkpoints under `output/STATS/.checkpoints/` (`overwrite = FALSE` resumes an interrupted run). The monthly rerun, zone-level GAM figures, regime GLMs and RF H-statistic were removed 2026-09-28 (dead ends with no manuscript consumer; their last CSVs stay on disk as frozen results the text cites, and the code is in git history at commit 533e42c). GLM/GAM/RF outputs are excluded from `tools/snapshot_outputs.py`'s golden-hash check (exploratory, and RF drifts by design).

All scripts prepend `func/` to `sys.path` by setting `proj_dir` from `os.path.abspath('__file__')` — they must be executed from the repo root, not from inside `code/`.

### Snakemake (in progress)
`Snakefile` at the repo root is being built up step by step to describe the same pipeline as rules, with declared inputs and outputs. It runs nothing unless asked. `snakemake -n` dry-runs it and `snakemake -n --summary` shows what's out of date. The numbered `code/*.py` scripts remain the way to run the pipeline until the Snakefile covers it. Its settings come from `metadata/riomar_config.yml`. Logs go to the gitignored `logs/`, and Snakemake's bookkeeping to the gitignored `.snakemake/`.

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

Large datasets are stored **outside** this repo under the pCloud data folder (`~/pCloud Drive/data/` on macOS, `~/pCloudDrive/data/` on Linux) and are never committed. Code never hardcodes this path: Python uses `func/config.py` (`config.data_root()`, `config.data_path('WIND', zone)`) and R uses `func/config.R` (`riomar_data_root()`, `riomar_data_path(...)`, sourced by `util.R`). They pick the folder by OS (falling back to the other spelling), unless `data_root` in `metadata/riomar_config.yml` or the `RIOMAR_DATA_ROOT` env var overrides it. The same YAML holds the zone list and the default satellite dict. Panache settings live in the same YAML's `panache` section (including `input_path`, the SEXTANT SPM folder panache reads — currently the external `/Volumes/Toshi` drive); `tools/write_panache_configs.py` renders them into machine-specific JSONs under the gitignored `output/panache/configs/`, and `func/figure.py` reads them via `config.panache_zone_config(zone, mode)`. Never edit the generated JSONs. The `.gitignore` also excludes most of `output/` and `data/SEXTANT`, `data/INSITU_data`, etc. Only shapefiles, metadata CSVs, and `metadata/riomar_config.yml` are tracked.

Two distinct kinds of "outside the repo" apply here, and only one of them is backed up automatically:

- **Pipeline-native pCloud data** — `SEXTANT`, `WIND`, `WAVE`, `GLORYS`, `DOWNLOAD_REPORTS` under `~/pCloudDrive/data/`. `code/0_download_data.py`, `code/2_regional_maps.py`, `code/5_figures.py`, and the R side (`func/util.R`, `func/multi.R`, `func/tools/VOG.R`, `func/analysis/compute_driver_spatial_variance.R`) read and write these paths directly — the data never lives inside the repo working tree, so it's inherently durable across machine migrations.
- **Repo-local gitignored content** — everything else the `.gitignore` excludes (`manuscript/`, `data/EUROPE_shapefile`, `data/HydroRIVERS_v10_eu_shp`, `data/INSITU_data`, `data/RIVER_FLOW`, `data/TIDES`, `output/MATCH_UP_DATA`, `output/STATS`, `output/panache`, `figures/ARTICLE/*/DATA`, `figures/ARTICLE/gam_monthly_breakdown`, and the diagnostic figure folders `figures/{driver_comparison,driver_octant_trend,driver_spatial_variance,driver_x11_comparison,qc,rhone_side_analyses,validation}` + `animations/`, untracked 2026-09-28) lives only in the repo working tree on whichever machine produced it. Nothing copies this automatically — a lost or wiped machine loses it (this is what happened to `data/ROFI` in the September 2026 migration; still unresolved, see memory). Back it up by hand into the pre-existing mirror at `~/pCloud Drive/Documents/OMTAB/RiOMar/`, which mirrors the repo's structure path-for-path:
  ```bash
  rsync -au <repo_relative_path>/ "~/pCloud Drive/Documents/OMTAB/RiOMar/<repo_relative_path>/"
  ```
  `-au` (`--archive --update`) only overwrites a pCloud file when the local copy is newer, so it's safe to re-run anytime — e.g. before a machine migration, or after a pipeline re-run that regenerates `output/`. `data/ROFI` currently has no local or pCloud copy — not yet restorable, but no longer needed: the ROFI analysis (`func/ROFI.R`) was removed from the pipeline 2026-09-28. `RIOMAR_old/` inside that same pCloud folder is an unrelated pre-reorg archive, not this mirror.

## Architecture

### Study zones
Four coastal zones used throughout: `GULF_OF_LION`, `BAY_OF_SEINE`, `BAY_OF_BISCAY`, `SOUTHERN_BRITTANY`.

### func/ modules (Python)
- [func/dl.py](func/dl.py) — FTP/CMEMS downloads via `copernicusmarine`; `Download_satellite_data`, `download_cmems_subset`, `daily_integral`
- [func/util.py](func/util.py) — shared helpers: file discovery, parameter parsing, zone coordinates, path templating
- [func/regmap.py](func/regmap.py) — regional map creation and QC (`create_regional_maps`, `QC_of_regional_maps`)
- [func/X11.py](func/X11.py) — X11 seasonal decomposition, calls R via `rpy2`; `Apply_X11_method_on_time_series`. **Frozen**: this file is intentionally excluded from renaming/refactor passes — do not edit its existing functions, even for naming consistency. New wrapper functions may be added alongside them. `func/X11.R` is not frozen.
- [func/figure.py](func/figure.py) — all publication figures, entry points named after their current manuscript number (renamed 2026-09-28; `metadata/figure_table_registry.csv` stays the source of truth if they drift again): `Figure_1_mean_spm_map`, `Figure_2_methodology_panels`/`Figure_2_methodology_zone_maps`/`Figure_2_methodology`, `Figure_3_S2_timeseries` (Fig. 3 + Fig. S2), `Figure_4_monthly_median_heatmap`, `Figure_X11_weekly_results` (Figs. 6, 7, S3–S6), `Figure_8_driver_rose`, `Figure_S1_validation`, `Figure_S8_gam_partial`. Fig. 5 is written by `func/analysis/generate_monthly_trend_pct_heatmap.R` from `code/4_time_series.py`. Output folders `figures/ARTICLE/FIGURE_<n>/` also match the manuscript number. Figures whose manuscript slot is marked `(removed)` in `metadata/figure_table_registry.csv` have no code left in the repo — each such row's notes name the commit to recover it from

### func/ layout
Shared libraries (`util`, `multi`, `figure`, `config`, `dl`, `regmap`, `driver_interactions`, `tide`, `surface`, `validate`, `river_flow_prep`, `X11`) sit directly in `func/`. Stand-alone scripts live in `func/analysis/` (manuscript statistics, several called from `code/4_time_series.py` as `analysis/<script>`) and `func/tools/` (one-off runners and diagnostics, e.g. `run_figure_1.R`, called by `figure.py`). All are still run from the repo root. After moving files, `python tools/check_paths.py` checks that every `.R`/`.py` reference still resolves, with no data needed. `util.R`, `multi.R` and `figure.R` are thin loaders (split 2026-09-28): each `source()`s its topic files from `func/sections/` (e.g. `util_4_loading.R`, `multi_5_rhone.R`, `figure_3_timeseries.R`) in a fixed order with `local = TRUE`, so keep sourcing the loader, and edit functions in the section files. Order matters: later sections use earlier ones, and some sections run code at load time (e.g. `multi_4_surface_missing.R` writes the missing-data CSVs). `python tools/check_split.py` checks that the sections exist and that no function is defined twice. Older comments citing `util.R:<line>`/`multi.R:<line>` refer to the pre-split line numbers.

### func/ modules (R)
Parallel R implementations exist for most modules (`util.R`, `validate.R`, `X11.R`, etc. — there is no `regmap.R`, regional-map creation is Python-only). These are used for analyses that rely on R packages (e.g. base stats and all plotting) and are called from Python via `rpy2`. `func/validate.R` is the authoritative satellite-vs-in-situ match-up pipeline (writes both `output/MATCH_UP_DATA/FRANCE/summary.csv`, feeding the manuscript's `validation_summary_stats` table, and the SEXTANT/ODATIS-MR `STATISTICS/*.csv` tables feeding the `validation_scatterplot_panel` figure — see `metadata/figure_table_registry.csv` for their current numbers).

### Per-river panache output
`panache`'s `Results.csv`/`PlumeMasks.nc` carry one row/mask layer per individual river mouth within a zone, plus an `'ALL'` union-mask row/layer (the zone total). RiOMar's loaders (`util.R::load_plume_ts()`, `figure.py::_load_results()`, `compute_plume_shape.py`, etc.) default to `river == "ALL"` for all zone-level stats. A real per-river analysis layer also exists (`func/analysis/compute_river_plume_correlation.R`, `X11.py::Apply_X11_method_on_time_series_per_river()`, `metadata/river_discharge_mapping.csv`) but is intentionally not surfaced in the manuscript — all published tables/figures stay at zone level, per project convention.

### metadata/
`riomar_config.yml` (project settings, including panache's — see Data storage) and zone-pixel CSVs (one per sensor × variable × atmospheric correction combination) used for plume pixel extraction.

Also tracked here (moved out of the gitignored `manuscript/` on 2026-09-28 so a fresh clone can run the pipeline): `figure_table_registry.csv` (slot → current figure/table number, output folder, rendering function; read by `util.py::get_registry_row()` and `util.R`), `paragraph_source_registry.csv` (read by `util.R`), `TODO.md` (the manuscript/pipeline to-do list), and `make_figures_tables.R` (the figure/table/paragraph-source checklist; run `Rscript metadata/make_figures_tables.R` from the repo root — it still reads `manuscript/manuscript.tex` and `references.bib`). There are no copies left in `manuscript/` (`google_doc_sync/` never read them).

### Satellite data dict convention
A Python dict like:
```python
{'Data_sources': ['SEXTANT'], 'Sensor_names': ['merged'], 'Satellite_variables': ['SPM'],
 'Atmospheric_corrections': ['Standard'], 'Temporal_resolution': ['DAILY'],
 'start_day': '1998/01/01', 'end_day': '2025/12/31'}
```
is the standard argument passed to every major pipeline function; build it with `config.satellite_dict('SPM')` / `('CHLA')` rather than retyping it. `util.define_parameters` converts it into a named-tuple `info` object.

### Multiprocessing
`dl.py` and `regmap.py` use `multiprocess` (not the stdlib `multiprocessing`). The start method is forced to `'spawn'` for macOS compatibility — do not change this.

## Future work ideas

- An SST (sea surface temperature) analysis of river plumes — how the plume footprint appears in a high-resolution SST product — is no longer a possible RiOMar analysis (2026-09-23: Robert no longer has access to the server hosting the hi-res SST data). A stub subsection for this in `manuscript.tex` (§ Results, "Sea surface temperature") was removed 2026-08-12 since it held no content; as of 2026-09-23 it's kept alive only as a one-sentence future-work mention at the end of the Conclusion (§`sec:conclusion`, third paragraph) — not a manuscript gap, and not scoped as a follow-up study for this project.

## Bug history

Bug-fix history (root cause, symptom, before/after) is intentionally not kept here — it bloats a file loaded every session regardless of relevance. Each fix is commented in place at its own function; check `git log`/`git blame` on the relevant file for the full story.


## Rules for automated changes

Added 2026-09-28 after an agent deleted code and outputs the manuscript still depended on (commits 321fc22/7cedb76, rolled back in 4017af8). These apply to every Claude session, local or cloud:

- **Never delete without per-item approval.** Don't delete or untrack a file, or remove a function, unless Robert has approved that exact item in writing. Present the list first (path or function name, what uses it, why it looks dead) and wait. A general instruction like "clean up dead code", or a registry row marked `(removed)`/`(orphaned)`, is not approval. Record each approved item in `metadata/approved_deletions.txt`, with the date and where the approval was given.
- **`manuscript/` is invisible to cloud sessions,** so "I found no reference" is not evidence that nothing uses a file. Treat anything named in `metadata/paragraph_source_registry.csv` or `metadata/figure_table_registry.csv` as in use, whatever its status column says.
- **Enforced by `tools/check_protected.py`.** The `.claude/settings.json` PreToolUse hook runs it before every `git commit`/`git push` issued through Bash, and it blocks deleted files and removed functions. Moves and renames pass: same file content at a new path, or same argument list under a new function name. Run it by hand with `python tools/check_protected.py` (uncommitted changes) or `--range A..B`. Don't bypass it (for example with `--no-verify` or by editing the hook) without Robert's say-so.
- **Only deliver work as pull requests.** Put changes on a new branch and open a PR. Never push to `main` or to a branch Robert is running the pipeline from; he merges when no run is in progress.
- **Prove refactors are result-neutral.** A refactor must leave `python tools/snapshot_outputs.py --check` unchanged on Robert's machine. Keep structural PRs to moves/renames, and never mix them with deletions.
