#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to analyse the river plumes used in the RiOMar project.

import os
import subprocess

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# NB: The plume detection is done via the panache module.
# NB: Takes roughly 60 minutes per zone.


# =============================================================================
#### Dynamic thresholds
# =============================================================================

subprocess.run("panache metadata/zone_config_dynamic_GULF_OF_LION.json", shell=True, cwd=proj_dir, check=True)
subprocess.run("panache metadata/zone_config_dynamic_BAY_OF_BISCAY.json", shell=True, cwd=proj_dir, check=True)
subprocess.run("panache metadata/zone_config_dynamic_SOUTHERN_BRITTANY.json", shell=True, cwd=proj_dir, check=True)
subprocess.run("panache metadata/zone_config_dynamic_BAY_OF_SEINE.json", shell=True, cwd=proj_dir, check=True)


# =============================================================================
#### Static thresholds
# =============================================================================

subprocess.run("panache metadata/zone_config_static_GULF_OF_LION.json", shell=True, cwd=proj_dir, check=True)
subprocess.run("panache metadata/zone_config_static_BAY_OF_BISCAY.json", shell=True, cwd=proj_dir, check=True)
subprocess.run("panache metadata/zone_config_static_SOUTHERN_BRITTANY.json", shell=True, cwd=proj_dir, check=True)
subprocess.run("panache metadata/zone_config_static_BAY_OF_SEINE.json", shell=True, cwd=proj_dir, check=True)

