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
"""

import argparse
import glob
import hashlib
import json
import os
import sys

OUTPUT_GLOBS = [
    'output/STATS/**/*.csv',
    'figures/ARTICLE/**/DATA/*.csv',
]
GOLDEN_PATH = os.path.join('tests', 'golden_hashes.json')


def hash_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def current_hashes():
    paths = sorted({p for pattern in OUTPUT_GLOBS for p in glob.glob(pattern, recursive=True)})
    return {p: hash_file(p) for p in paths}


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
        golden = json.load(f)
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
