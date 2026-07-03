#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to analyse the river plumes used in the RiOMar project. 
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

# import os
# import sys
import subprocess

# NB: The plume detection is done via the panache module, which is called directly via terminal.


# =============================================================================
# ### Detect plumes
# =============================================================================

# Run panache module directly via terminal
# NB: Takes roughly 20 minutes per zone
# TODO: At the moment this does not run via Python. Needs to be called directly within a terminal.
subprocess.run("panache metadata/zone_config_GULF_OF_LION.json", shell=True)
subprocess.run("panache metadata/zone_config_BAY_OF_BISCAY.json", shell=True)
subprocess.run("panache metadata/zone_config_SOUTHERN_BRITTANY.json", shell=True)
subprocess.run("panache metadata/zone_config_BAY_OF_SEINE.json", shell=True)

    