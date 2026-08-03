#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to validate the full dataset used in the RiOMar project.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import rpy2.robjects as robjects

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )


# =============================================================================
# ### Satellite - in situ Match up for Chl a and SPM data
# =============================================================================

# The full validation pipeline (pixel extraction, match-up stats, scatterplots,
# and tables) is implemented as a sequential R script. Run it directly.
validate_R_path = os.path.join(func_dir, 'validate.R')
robjects.r['source'](validate_R_path)

