#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to create the figures used in the publication of this workflow.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, figure
from figure import Figure_1, Figure_4, Figure_5, Figure_6_7, Figure_8_9_10

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')


# =============================================================================
# ### Main figures (dynamic threshold)
# =============================================================================

Figure_1(where_are_saved_satellite_data = "../pCloudDrive/data",
         where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_4(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_5(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

# Figure 6-7: compare dynamic vs static threshold plume detection
Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold = "output/panache/dynamic",
           where_are_saved_plume_results_with_fixed_threshold = "output/panache/static",
           where_to_save_the_figure = "figures")

# Figures 8-10: X11 decomposition results — dynamic threshold (main)
Figure_8_9_10(where_are_saved_X11_results = "output/panache/dynamic",
              where_to_save_the_figure = "figures")


# =============================================================================
# ### Supplementary figures (static threshold)
# =============================================================================

# Figures 8-10 equivalent using static threshold results
Figure_8_9_10(where_are_saved_X11_results = "output/panache/static",
              where_to_save_the_figure = "figures/supplementary")
