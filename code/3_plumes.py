#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to analyse the river plumes used in the RiOMar project.

import os
import sys
import subprocess

proj_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(os.path.join(proj_dir, 'tools'))

from write_panache_configs import write_all, config_path

# NB: The plume detection is done via the panache module.
# NB: Takes roughly 60 minutes per zone.

# Write this machine's panache JSONs (output/panache/configs/) from the
# `panache` section of metadata/riomar_config.yml -- edit settings there.
write_all()

# Same order as the original hand-written calls
zones_in_run_order = ['GULF_OF_LION', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY', 'BAY_OF_SEINE']


# =============================================================================
#### Dynamic thresholds, then static thresholds
# =============================================================================

for mode in ['dynamic', 'static']:
    for zone in zones_in_run_order:
        subprocess.run(['panache', config_path(zone, mode)], cwd=proj_dir, check=True)
