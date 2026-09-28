# Snakefile -- RiOMar pipeline as a Snakemake workflow (being built step by step).
#
# The numbered code/*.py scripts are still the normal way to run the pipeline;
# this file describes the same steps as rules and doesn't replace them yet.
# Nothing here runs unless you ask it to:
#
#   snakemake -n                  # dry run: list what WOULD run, and why
#   snakemake --dag | dot -Tpdf > dag.pdf   # draw the rules as a graph (needs graphviz)
#   snakemake -n --summary        # every output, its date, and whether it is out of date
#   snakemake -c1 <output file>   # actually (re)build one output, using 1 core
#
# The first dry run may also say "missing provenance/metadata": that only means
# the existing outputs were made outside Snakemake (by code/*.py). It goes away
# once Snakemake has built a file itself.
#
# Core idea: each rule declares the files it reads (input) and writes (output),
# plus the command that turns one into the other. To build a file, Snakemake
# finds the rule that outputs it, and runs it only if the output is missing
# or older than any of its inputs -- code files included, so editing a
# script also marks its outputs out of date.
# (Once Snakemake has built a file itself, it also stores the inputs'
# checksums: an input that is newer but has identical content, e.g. after
# `touch` or a git checkout, does not trigger a rerun. Real edits do.)


# Settings ------------------------------------------------------------------

import csv
import json
import glob  # a Snakefile is Python plus rule blocks, so ordinary Python works here

# The same YAML the Python/R code reads (func/config.py, func/config.R):
# it becomes the `config` dict here.
configfile: "metadata/riomar_config.yml"

ZONES = config["zones"]
PANACHE_MODES = ["dynamic", "static"]

# Wildcards ({zone}, {mode}) may only take these values, so a path like
# output/panache/dynamic/typo/Results.csv can't match a rule by accident.
wildcard_constraints:
    zone = "|".join(ZONES),
    mode = "|".join(PANACHE_MODES),

# Shared R code that most analysis scripts source(). util.R/multi.R are loaders
# for their func/sections/ parts, so the parts are listed too -- editing any
# of them marks every rule that uses SHARED_R as out of date.
SHARED_R = ["func/config.R", "func/util.R", "func/multi.R"] + sorted(
    glob.glob("func/sections/util_*.R") + glob.glob("func/sections/multi_*.R"))


def registry_figure(slot_key):
    """A figure's output path, from metadata/figure_table_registry.csv --
    the same lookup the R scripts do (get_registry_row()/registry_filename()),
    so renumbering a figure there moves it here too."""
    with open("metadata/figure_table_registry.csv", newline="") as f:
        row = next(r for r in csv.DictReader(f) if r["slot_key"] == slot_key)
    subdir = row["output_subdir"]
    return f"figures/ARTICLE/{subdir}/{subdir.replace('FIGURE_', 'Figure_', 1)}.png"


# Default target ---------------------------------------------------------------
# `snakemake` with no target builds the FIRST rule in the file, so by convention
# that is `all`: a rule with no command whose inputs are the final outputs you
# want. Snakemake then works backwards to every rule needed to make them.

rule all:
    input:
        "output/STATS/area_trend_summary.csv",
        # cited directly in the manuscript text, so final outputs too (a file
        # that is only an intermediate isn't rebuilt when missing if everything
        # downstream of it is up to date)
        "output/STATS/monthly_trend_summary.csv",
        "output/STATS/monthly_trend_compact_summary.csv",
        registry_figure("monthly_trend_pct_heatmap"),
        "output/STATS/driver_rf_importance.csv",
        "output/STATS/driver_correlation_matrix_seasonal.png",


# Step 1: one rule ------------------------------------------------------------
# Long-term linear trend in daily plume area, per zone (feeds the
# panache_stats_table "Surface area" row and the Fig. 3 trend lines).
# Runs the script exactly as its header documents (Rscript, from the repo root);
# code/4_time_series.py runs the same script via rpy2.

rule area_trend:
    input:
        script = "func/analysis/compute_area_trend.R",
        code = SHARED_R,
        zone_metadata = "metadata/panache_zone_metadata.csv",
        # expand() fills {zone} in for every zone: 4 Results.csv files
        plumes = expand("output/panache/dynamic/{zone}/Results.csv", zone=ZONES),
        # the river-flow folders (a folder's date changes when files are added to it)
        flow = expand("data/RIVER_FLOW/{zone}", zone=ZONES),
    output:
        "output/STATS/area_trend_summary.csv",
    # The script's printed output goes here instead of the terminal
    # (matters once several rules run in parallel).
    log:
        "logs/area_trend.log",
    shell:
        "Rscript {input.script} > {log} 2>&1"


# Step 2: chaining rules -----------------------------------------------------
# generate_monthly_trend_pct_heatmap.R (Fig. 5) reads the CSV that
# compute_seasonal_trend.R writes. Nothing here says "run seasonal_trend
# first": Snakemake sees that the heatmap rule's input
# output/STATS/monthly_trend_summary.csv is the seasonal_trend rule's output,
# and orders them itself -- the "order matters" comment in
# code/4_time_series.py becomes unnecessary.
#
# Choosing inputs is a judgement call: the in-repo data is declared, but the
# pCloud driver data (WIND/WAVE/GLORYS, outside the repo) is not -- it is
# treated as fixed, so re-downloading it would not by itself trigger a rerun.

rule seasonal_trend:
    input:
        script = "func/analysis/compute_seasonal_trend.R",
        code = SHARED_R + ["func/tide.R"],
        zone_metadata = "metadata/panache_zone_metadata.csv",
        plumes = expand("output/panache/{mode}/{zone}/Results.csv",
                        mode=["dynamic", "static"], zone=ZONES),
        shapes = expand("output/panache/{mode}/{zone}/PlumeShape.csv",
                        mode=["dynamic", "static"], zone=ZONES),
        flow = expand("data/RIVER_FLOW/{zone}", zone=ZONES),
        tides = "data/TIDES",
    output:
        detail = "output/STATS/monthly_trend_summary.csv",
        compact = "output/STATS/monthly_trend_compact_summary.csv",
    log:
        "logs/seasonal_trend.log",
    shell:
        "Rscript {input.script} > {log} 2>&1"


rule monthly_trend_heatmap:
    input:
        script = "func/analysis/generate_monthly_trend_pct_heatmap.R",
        code = SHARED_R + ["func/tide.R"],
        registry = "metadata/figure_table_registry.csv",
        trends = rules.seasonal_trend.output.detail,   # <- the link between the two rules
        plumes = expand("output/panache/dynamic/{zone}/Results.csv", zone=ZONES),
        shapes = expand("output/panache/dynamic/{zone}/PlumeShape.csv", zone=ZONES),
        flow = expand("data/RIVER_FLOW/{zone}", zone=ZONES),
        tides = "data/TIDES",
    output:
        registry_figure("monthly_trend_pct_heatmap"),
    log:
        "logs/monthly_trend_heatmap.log",
    shell:
        "Rscript {input.script} > {log} 2>&1"


# Step 3: wildcards -----------------------------------------------------------
# One rule stands for all 8 panache runs (4 zones x 2 threshold modes).
# {zone} and {mode} are wildcards: when another rule needs, say,
# output/panache/static/BAY_OF_SEINE/Results.csv, Snakemake matches it against
# this rule's output pattern, fills in zone=BAY_OF_SEINE and mode=static, and
# uses those values in the input path and command. This replaces the loop in
# code/3_plumes.py.
#
# SAFETY (panache takes ~1 h per zone x mode): a panache job reruns only
# when its own Results.csv/PlumeMasks.nc are missing, or -- once Snakemake
# has built them itself -- when that zone x mode's panache settings in
# metadata/riomar_config.yml change (the `params` below; Snakemake records
# them). Each job writes its own JSON config as part of its command instead
# of depending on a shared config file, so editing unrelated YAML settings,
# or one missing output, never cascades into all 8 runs. The SEXTANT input on
# the external drive is not declared. Always check that `snakemake -n` lists
# no panache jobs before running anything, unless you mean to rerun detection.

rule panache:
    input:
        writer = ancient("tools/write_panache_configs.py"),
    params:
        # this zone x mode's settings; a change here marks the run out of date
        settings = lambda wc: json.dumps(
            {**config["panache"][wc.mode], **config["panache"]["common"],
             "input_path": config["panache"]["input_path"]}, sort_keys=True),
    output:
        results = "output/panache/{mode}/{zone}/Results.csv",
        masks = "output/panache/{mode}/{zone}/PlumeMasks.nc",
    # panache itself uses nb_cores cores per run; declaring them here stops
    # `snakemake -c8` from starting 8 runs at once (8 x 8 cores)
    threads: config["panache"]["common"]["nb_cores"]
    log:
        "logs/panache_{mode}_{zone}.log",
    shell:
        "python {input.writer} --zone {wildcards.zone} --mode {wildcards.mode} > {log} 2>&1"
        " && panache output/panache/configs/zone_config_{wildcards.mode}_{wildcards.zone}.json >> {log} 2>&1"


# Step 4: one expensive rule, many outputs ----------------------------------
# The multi-driver GLM/GAM/random-forest analysis is the slowest analysis step
# (tens of minutes). As its own rule, it reruns only when one of its inputs
# really changes -- tweaking a figure or another stats script never triggers
# it. A rule can declare several outputs; Snakemake treats them as one unit,
# so asking for any one of them (rule all asks for the RF CSV) brings in all.
#
# Its per-cell checkpoints (output/STATS/.checkpoints/) are deliberately not
# declared: they are the script's own resume mechanism, not results.
# pCloud driver data (wind/wave/current) is not declared (as in step 2).

# Must match metric_responses in func/driver_interactions.R
METRIC_RESPONSES = ["plume_area", "mean_SPM_in_the_plume_area", "mass_SPM_in_the_plume_area_in_t",
                    "lon_weighted_centroid_of_the_plume_area", "lat_weighted_centroid_of_the_plume_area"]

rule driver_interactions:
    input:
        script = "func/driver_interactions.R",
        code = SHARED_R + ["func/tide.R"],
        zone_metadata = "metadata/panache_zone_metadata.csv",
        plumes = expand("output/panache/dynamic/{zone}/Results.csv", zone=ZONES),
        flow = expand("data/RIVER_FLOW/{zone}", zone=ZONES),
        tides = "data/TIDES",
    output:
        matrices = expand("output/STATS/daily_driver_matrix_{zone}.csv", zone=ZONES),
        glm = "output/STATS/driver_glm_comparison.csv",
        metric_models = expand("output/STATS/driver_metric_models_{resp}.csv", resp=METRIC_RESPONSES),
        rf = "output/STATS/driver_rf_importance.csv",
    log:
        "logs/driver_interactions.log",
    shell:
        # the same call code/4_time_series.py makes (a separate Rscript process,
        # which avoids the ranger/numpy OpenMP clash noted there)
        """Rscript -e "source('{input.script}'); run_driver_interactions_analysis()" > {log} 2>&1"""


rule driver_correlation_matrices:
    input:
        script = "func/analysis/compute_driver_correlation_matrices.R",
        code = SHARED_R + ["func/tide.R", "func/driver_interactions.R"],
        matrices = rules.driver_interactions.output.matrices,   # chained, as in step 2
    output:
        "output/STATS/driver_correlation_matrix_overall.png",
        "output/STATS/driver_correlation_matrix_seasonal.png",
    log:
        "logs/driver_correlation_matrices.log",
    shell:
        "Rscript {input.script} > {log} 2>&1"
