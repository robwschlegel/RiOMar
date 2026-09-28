#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project-wide settings, read from metadata/riomar_config.yml. The R twin is
func/config.R -- keep the data-root resolution logic identical in both.
"""

import copy
import functools
import os
import platform

import yaml

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_PATH = os.path.join(proj_dir, 'metadata', 'riomar_config.yml')

# pCloud's mount folder is named differently per OS (macOS spelling first)
PCLOUD_FOLDERS = ['pCloud Drive', 'pCloudDrive']


@functools.lru_cache(maxsize=None)
def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def data_root():
    """
    Absolute path of the large-dataset root (SEXTANT, WIND, WAVE, GLORYS...).

    Resolution order: RIOMAR_DATA_ROOT env var > `data_root` in the YAML >
    the OS's pCloud folder (~/pCloud Drive/data on macOS, ~/pCloudDrive/data
    elsewhere), falling back to the other spelling if the OS default is
    missing. Raises FileNotFoundError naming every path tried.
    """
    explicit = os.environ.get('RIOMAR_DATA_ROOT') or load_config().get('data_root')
    if explicit:
        candidates = [explicit]
    else:
        folders = PCLOUD_FOLDERS if platform.system() == 'Darwin' else PCLOUD_FOLDERS[::-1]
        candidates = [os.path.join('~', folder, 'data') for folder in folders]
    for candidate in candidates:
        path = os.path.abspath(os.path.expanduser(candidate))
        if os.path.isdir(path):
            return path
    raise FileNotFoundError(
        'RiOMar data root not found; tried: ' + ', '.join(candidates) +
        '. Set RIOMAR_DATA_ROOT or data_root in metadata/riomar_config.yml.')


def data_path(*parts):
    """data_path('WIND', zone) -> <data_root>/WIND/<zone>"""
    return os.path.join(data_root(), *parts)


def zones():
    return list(load_config()['zones'])


def satellite_dict(variable):
    """The standard satellite-data dict (see CLAUDE.md) for 'SPM' or 'CHLA'."""
    sat = copy.deepcopy(load_config()['satellite'])
    return {'Data_sources': sat['Data_sources'],
            'Sensor_names': sat['Sensor_names'],
            'Satellite_variables': [variable],
            'Atmospheric_corrections': sat['Atmospheric_corrections'],
            'Temporal_resolution': sat['Temporal_resolution'],
            'start_day': sat['start_day'],
            'end_day': sat['end_day']}


PANACHE_MODES = ('dynamic', 'static')


def panache_zone_config(zone, mode, proj_dir=proj_dir):
    """
    The panache settings for one zone and threshold mode ('dynamic' or
    'static'), as the dict panache expects in its zone_config JSON (same key
    order as the hand-written JSONs this replaced). Paths are absolute for
    this machine; proj_dir can be overridden to render them for another.
    """
    if mode not in PANACHE_MODES:
        raise ValueError(f"mode must be one of {PANACHE_MODES}, got {mode!r}")
    panache = load_config()['panache']
    input_path = os.environ.get('RIOMAR_SEXTANT_SPM_PATH') or panache.get('input_path')
    if not input_path:
        input_path = data_path('SEXTANT', 'SPM') + os.sep
    output_dir = os.path.join(proj_dir, 'output', 'panache', mode, zone)
    return {'zone': zone,
            'input_path': input_path,
            'output_dir': output_dir,
            'bathymetry_path': os.path.join(output_dir, 'Bathy_data.pkl'),
            'coast_shapefile': os.path.join(proj_dir, panache['coast_shapefile']),
            **panache[mode],
            **panache['common']}
