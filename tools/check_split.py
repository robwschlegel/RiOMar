#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the loader files that source their parts from func/sections/
(func/util.R, func/multi.R, func/figure.R).

    python tools/check_split.py              # parts exist; no function defined twice
    python tools/check_split.py --ref REF    # also: parts reassemble to REF's single file

The split (2026-09-28) cut each former single file into contiguous slices,
sourced in their original order with local = TRUE, so behaviour is unchanged.
--ref proves that for a given pre-split commit: concatenating the parts in
loader order must reproduce that commit's file byte for byte (only meaningful
until someone edits a part).
"""

import argparse
import collections
import os
import re
import subprocess
import sys

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOADERS = ['func/util.R', 'func/multi.R', 'func/figure.R']
SOURCE_LINE = re.compile(r'^source\("(func/sections/[^"]+)", local = TRUE\)')
R_DEF = re.compile(r'^([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*\(', re.M)


def loader_parts(loader):
    with open(os.path.join(proj_dir, loader)) as f:
        return [m.group(1) for m in map(SOURCE_LINE.match, f) if m]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ref', help='pre-split commit to compare the reassembled files against')
    args = parser.parse_args()

    problems = []
    definitions = collections.defaultdict(list)
    for loader in LOADERS:
        parts = loader_parts(loader)
        if not parts:
            problems.append(f'{loader}: no parts sourced')
        texts = []
        for part in parts:
            path = os.path.join(proj_dir, part)
            if not os.path.exists(path):
                problems.append(f'{loader}: missing part {part}')
                continue
            texts.append(open(path).read())
            for name in R_DEF.findall(texts[-1]):
                definitions[name].append(part)
        if args.ref:
            original = subprocess.run(['git', 'show', f'{args.ref}:{loader}'], cwd=proj_dir,
                                      capture_output=True, text=True, check=True).stdout
            if ''.join(texts) != original:
                problems.append(f'{loader}: parts do not reassemble to {args.ref}:{loader}')
            else:
                print(f'{loader}: {len(parts)} parts reassemble exactly to {args.ref}')
    for name, where in sorted(definitions.items()):
        if len(where) > 1:
            problems.append(f'{name}() defined in more than one part: {", ".join(where)}')
    for p in problems:
        print(f'PROBLEM  {p}')
    if problems:
        sys.exit(1)
    print(f'OK: {sum(len(loader_parts(l)) for l in LOADERS)} parts, {len(definitions)} functions, none duplicated')


if __name__ == '__main__':
    main()
