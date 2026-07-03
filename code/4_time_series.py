#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The code needed to run the time series analyses.


# =============================================================================
#### Modules
# =============================================================================

import os
import sys
import matplotlib as mpl

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

# import util
# import X11
from X11 import Apply_X11_method_on_time_series

# Set matplotlib backend to prevent plots from displaying
mpl.use('agg')

# The zones for mapping
zones_list = ['GULF_OF_LION', 'BAY_OF_SEINE', 'BAY_OF_BISCAY', 'SOUTHERN_BRITTANY']

# Basic arguments to be used throughout the script
sextant_spm_all = {'Data_sources':['SEXTANT'],
                   'Sensor_names':["merged"],
                   'Satellite_variables':['SPM'],
                   'Atmospheric_corrections':['Standard'],
                   'Temporal_resolution':['DAILY'],
                   'start_day':'1998/01/01',
                   'end_day':'2025/12/31'}


# =============================================================================
# ### X11 analyses
# =============================================================================

# NB: X11 can only be used on weekly or monthly data, not daily

# Full analysis in one go
Apply_X11_method_on_time_series(sextant_spm_all,
                                # Zones = ["SOUTHERN_BRITTANY"],
                                Zones = zones_list,
                                plume_time_step = "WEEKLY",
                                plume_dir_in = "output/panache",
                                X11_dir_out = "output/panache",
                                include_river_flow = True)


# =============================================================================
# ### Wind comparison
# =============================================================================


# =============================================================================
# ### Tide comparisons
# =============================================================================


# =============================================================================
# ### Other comparison
# =============================================================================


# =============================================================================
# ### Multivariate analysis
# =============================================================================

# The idea here is to use a GAM to understand what factors are driving plume surface
# variations. This will include river flow, wind speed and direction, precipitation,
# and perhaps others.

