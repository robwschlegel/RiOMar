#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write one panache zone_config JSON per zone x threshold mode, from the
`panache` section of metadata/riomar_config.yml.

    python tools/write_panache_configs.py
    panache output/panache/configs/zone_config_dynamic_GULF_OF_LION.json

The JSONs go to output/panache/configs/ (gitignored): their paths are
absolute for the machine that wrote them, so they are regenerated rather
than committed. Edit settings in riomar_config.yml, never in these files.
"""

import argparse
import json
import os
import sys

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(proj_dir, 'func'))

import config  # noqa: E402

CONFIG_DIR = os.path.join(proj_dir, 'output', 'panache', 'configs')


def config_path(zone, mode):
    return os.path.join(CONFIG_DIR, f'zone_config_{mode}_{zone}.json')


def write_all(proj_dir_for_paths=proj_dir):
    os.makedirs(CONFIG_DIR, exist_ok=True)
    paths = []
    for mode in config.PANACHE_MODES:
        for zone in config.zones():
            path = config_path(zone, mode)
            with open(path, 'w') as f:
                json.dump(config.panache_zone_config(zone, mode, proj_dir_for_paths), f, indent=2)
                f.write('\n')
            paths.append(path)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--proj-dir', default=proj_dir,
                        help='repo root to write into the JSON paths (default: this checkout)')
    args = parser.parse_args()
    for path in write_all(args.proj_dir):
        print(os.path.relpath(path, proj_dir))


if __name__ == '__main__':
    main()
