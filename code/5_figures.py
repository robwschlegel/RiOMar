#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to create the figures used in the publication of this workflow. 
# It can be designed to be called by the Makefile.


# =============================================================================
#### Modules
# =============================================================================

import os, sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

import util, figure
from figure import Figure_1, Figure_2, Figure_4, Figure_5, Figure_6_7, Figure_8_9_10

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')


# =============================================================================
# ### Create figures
# =============================================================================

Figure_1(where_are_saved_satellite_data = "data",
         where_to_save_the_figure = "figures")

Figure_2(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_4(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

Figure_5(where_are_saved_regional_maps = "output",
         where_to_save_the_figure = "figures")

# TODO: Still need to run dynamic threshold plumes to get the correct figure
Figure_6_7(where_are_saved_plume_results_with_dynamic_threshold = "output/DYNAMIC_THRESHOLD",
           where_are_saved_plume_results_with_fixed_threshold = "output/FIXED_THRESHOLD",
           where_to_save_the_figure = "figures")

Figure_8_9_10(where_are_saved_X11_results = "output/FIXED_THRESHOLD",
              where_to_save_the_figure = "figures")

print("All figures have been created and saved in the 'figures' folder.")

