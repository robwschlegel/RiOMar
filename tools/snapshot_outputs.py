#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fingerprint the pipeline's numerical outputs, to prove a refactor changed no
results.

    python tools/snapshot_outputs.py            # write tests/golden_hashes.json
    python tools/snapshot_outputs.py --check    # compare against it

Each file matched by OUTPUT_GLOBS is hashed (SHA-256) by content. --check
lists every file whose hash changed, that disappeared, or that is new, and
exits non-zero if anything differs. Run from the repo root on a machine that
has run the pipeline (output/ is gitignored, so there is nothing to hash in a
fresh clone). PNGs are not hashed: plotting libraries can embed timestamps or
render slightly differently between versions, so the CSVs behind the figures
are the reliable signal.

Files matching IGNORED_GLOBS are skipped entirely: the multi-driver GLM/GAM
outputs are exploratory dead ends not used in the main text, and random
forest importances drift between runs by design, so neither is a useful
refactor signal.
"""

import argparse
import fnmatch
import glob
import hashlib
import json
import os
import sys

OUTPUT_GLOBS = [
    'output/STATS/**/*.csv',
    'figures/ARTICLE/**/DATA/*.csv',
]
# fnmatch patterns ('*' also matches '/'), applied to repo-relative paths
IGNORED_GLOBS = [
    'output/STATS/driver_glm_comparison.csv',
    'output/STATS/driver_gam_summary.csv',
    'output/STATS/driver_regime_glm.csv',
    'output/STATS/driver_metric_models_*.csv',
    'output/STATS/driver_rf_*.csv',
    'output/STATS/monthly/*',
]
GOLDEN_PATH = os.path.join('tests', 'golden_hashes.json')


def hash_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def is_ignored(path):
    return any(fnmatch.fnmatch(path, pattern) for pattern in IGNORED_GLOBS)


def current_hashes():
    paths = sorted({p for pattern in OUTPUT_GLOBS for p in glob.glob(pattern, recursive=True)})
    return {p: hash_file(p) for p in paths if not is_ignored(p)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check', action='store_true', help=f'compare against {GOLDEN_PATH} instead of writing it')
    args = parser.parse_args()

    hashes = current_hashes()
    if not hashes:
        sys.exit('No output files found -- run from the repo root on a machine that has run the pipeline.')

    if not args.check:
        os.makedirs(os.path.dirname(GOLDEN_PATH), exist_ok=True)
        with open(GOLDEN_PATH, 'w') as f:
            json.dump(hashes, f, indent=1, sort_keys=True)
            f.write('\n')
        print(f'Wrote {len(hashes)} hashes to {GOLDEN_PATH}')
        return

    with open(GOLDEN_PATH) as f:
        golden = {p: h for p, h in json.load(f).items() if not is_ignored(p)}
    changed = sorted(p for p in golden.keys() & hashes.keys() if golden[p] != hashes[p])
    missing = sorted(golden.keys() - hashes.keys())
    new = sorted(hashes.keys() - golden.keys())
    for label, paths in (('CHANGED', changed), ('MISSING', missing), ('NEW', new)):
        for p in paths:
            print(f'{label:8} {p}')
    if changed or missing or new:
        sys.exit(f'{len(changed)} changed, {len(missing)} missing, {len(new)} new (of {len(golden)} golden files)')
    print(f'All {len(golden)} files match {GOLDEN_PATH}')


if __name__ == '__main__':
    main()
