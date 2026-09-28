# Snakefile -- RiOMar pipeline as a Snakemake workflow (being built step by step).
#
# The numbered code/*.py scripts are still the normal way to run the pipeline;
# this file describes the same steps as rules and doesn't replace them yet.
# Nothing here runs unless you ask it to:
#
#   snakemake -n                  # dry run: list what WOULD run, and why
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


# Settings ------------------------------------------------------------------

import glob  # a Snakefile is Python plus rule blocks, so ordinary Python works here

# The same YAML the Python/R code reads (func/config.py, func/config.R):
# it becomes the `config` dict here.
configfile: "metadata/riomar_config.yml"

ZONES = config["zones"]

# Shared R code that most analysis scripts source(). util.R/multi.R are loaders
# for their func/sections/ parts, so the parts are listed too -- editing any
# of them marks every rule that uses SHARED_R as out of date.
SHARED_R = ["func/config.R", "func/util.R", "func/multi.R"] + sorted(
    glob.glob("func/sections/util_*.R") + glob.glob("func/sections/multi_*.R"))


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
