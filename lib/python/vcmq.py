"""Module for quickly loading VACUMM essentials"""

# System
import os, sys, shutil
import re
import glob
import gc
import warnings
from collections import OrderedDict
import argparse
from traceback import format_exc

# Numeric

import numpy as N
npy = N
np = N


# CDAT

import MV2, cdms2,  cdtime, cdutil
from cdutil import averager
from genutil import minmax

# Matplotlib

import matplotlib.pyplot as plt
import pylab as P
mp = P
from matplotlib import use, rc, rcdefaults, rcParams
from matplotlib.dates import DayLocator, HourLocator, WeekdayLocator, MinuteLocator,  \
    YearLocator, MonthLocator, DateFormatter, SecondLocator
from matplotlib.ticker import AutoLocator, MaxNLocator, FixedLocator, IndexLocator, \
    LinearLocator, LogLocator, MultipleLocator, AutoMinorLocator

# VACUMM

from vacumm import VACUMMError, VACUMMWarning, help as vchelp


# - config

from vacumm.config import \
    edit_user_conf_file, print_config, get_config_value, data_sample

# - misc 

from vacumm.misc.misc import \
    lonlab, latlab, deplab, kwfilter, broadcast, ArgList, closeto, dict_merge, \
    dict_check_defaults, cp_atts, dict_filter,  cp_props, numod, set_atts, get_atts, \
    set_lang, geodir, auto_scale
    
from vacumm.misc.plot import \
    map2, section2, hov2, curve2, bar2, plot2d, stick2, make_movie, \
    savefigs, add_grid, xhide, yhide, xrotate, yrotate, add_key, taylor, dtaylor, dtarget, \
    add_shadow, add_glow, add_map_lines, add_map_line, add_map_point, minimap
    
from vacumm.misc.color import plot_cmap, show_cmap, get_cmap, simple_colors, \
    cmap_rs, cmap_srs, cmap_custom, Scalar2RGB, darken, whiten

from vacumm.misc.axes import \
    create_lon, create_lat, create_time, create_depth, \
    create_dep, islon, islat, isdep, isaxis, istime
    
from vacumm.misc.atime import comptime, strftime, strptime, Intervals, IterDates, now, \
    time_split, time_split_nmax, to_utc, tz_to_tz, utc_to_paris, paris_to_utc,  \
    lindates, ch_units, add_margin, round_date, midnight_date, midnight_interval, mpl


from vacumm.misc.bases import psinfo, code_base_name
    
from vacumm.misc.io import list_forecast_files, netcdf3, netcdf4, ncread_files, ncread_var, \
    ncread_axis, ncget_var, ncget_axis, ncfind_var, ncfind_axis, NcIterBestEstimate,  \
    ncget_fgrid
    
from vacumm.misc.filters import generic1d, generic2d, shapiro1d, shapiro2d, blackman1d, \
    bartlett1d, kaiser1d, hanning1d, blackman1d, gaussian1d, gaussian2d

from vacumm.misc.remote import InputWorkFiles, OutputWorkFile

from vacumm.misc.log import Logger

from vacumm.misc.config import ConfigManager, cfgargparse

from vacumm.misc.stats import StatAccum

from vacumm.misc.phys.units import convert_units
   
# - grid

from vacumm.misc.grid.misc import \
    get_grid, set_grid, create_grid, create_grid2d, resol, monotonic, isdepthup, depth2dz, \
    meshcells, meshbounds, meshgrid, bounds1d, bounds2d, coord2slice, isregular, \
    isrect, curv2rect, isgrid, get_xy, create_axes2d, gridsel, varsel, xshift, rotate_grid, \
    get_closest, makedepthup, get_axis, transect_specs
    
from vacumm.misc.grid.regridding import \
    regrid1d, regrid2d, interp1d, interp2d, cellave1d, cellave2d, \
    cargen, xy2xy, shift1d, shift2d, extend1d, extend2d, regrid_method, \
    transect, griddata, CDATRegridder, grid2xy, \
    nearest1d, shiftgrid, extendgrid, fill1d, fill2d

from vacumm.misc.grid.basemap import \
    merc, create_map, get_proj, reset_cache
    
from vacumm.misc.grid.masking import \
    erode_coast, polygon_mask, GetLakes, get_coast, polygons, polygon_select, zcompress, \
    envelop, get_coastal_indices, grid_envelop
    
# - data

from vacumm.data import setup_dataset, DS

from vacumm.data.cf import format_var, format_axis, format_grid, match_var, change_loc, \
    set_loc, get_loc

from vacumm.data.misc.sigma import NcSigma, sigma2depths, SigmaGeneralized, SigmaStandard

from vacumm.data.misc.arakawa import ArakawaGrid, locations as arakawa_locations, \
    CGrid, AGrid, ArakawaGridTransfer

# - diag

from vacumm.diag.thermdyn import density, mixed_layer_depth

from vacumm.diag.dynamics import \
    barotropic_geostrophic_velocity, coriolis_parameter, kinetic_energy

