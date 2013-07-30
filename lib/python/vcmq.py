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


# CDAT

import MV2, cdms2,  cdtime, cdutil
from cdutil import averager
from genutil import minmax

# Matplotlib

import matplotlib.pyplot as plt
import pylab as P
from matplotlib import use, rc, rcdefaults, rcParams
from matplotlib.dates import DayLocator, HourLocator, WeekdayLocator, MinuteLocator,  \
    YearLocator, MonthLocator, DateFormatter, SecondLocator
from matplotlib.ticker import AutoLocator, MaxNLocator, FixedLocator, IndexLocator, \
    LinearLocator, LogLocator, MultipleLocator, AutoMinorLocator

# VACUMM

from vacumm import VACUMMError, VACUMMWarning


# - config

from vacumm.config import \
    edit_user_conf_file, print_config, get_config_value, data_sample

# - misc 

from vacumm.misc.misc import \
    lonlab, latlab, deplab, kwfilter, broadcast, ArgList, \
    dict_check_defaults, cp_atts, dict_filter,  cp_props, numod
    
from vacumm.misc.plot import \
    map2, section2, hov2, curve2, bar2, plot2d, \
    savefigs, add_grid, xhide, yhide, xrotate, yrotate, add_key, \
    add_shadow, add_glow
    
from vacumm.misc.color import plot_cmap, show_cmap, get_cmap, simple_colors

from vacumm.misc.axes import \
    create_lon, create_lat, create_time, create_depth, \
    create_dep, islon, islat, isdep, isaxis, istime
    
from vacumm.misc.atime import comptime, strftime, strptime, Intervals, IterDates, now, \
    time_split, time_split_nmax, to_utc, tz_to_tz, utc_to_paris, paris_to_utc,  \
    lindates, ch_units, add_margin, round_date, midnight_date, midnight_interval


from vacumm.misc.bases import psinfo, code_base_name
    
from vacumm.misc.io import list_forecast_files, netcdf3, netcdf4, ncread_files, ncread_var, \
    ncread_axis, ncget_var, ncget_axis, ncfind_var, ncfind_axis, NcIterBestEstimate,  \
    ncget_fgrid
    
from vacumm.misc.filters import generic1d, generic2d, shapiro1d, shapiro2d
   
# - grid

from vacumm.misc.grid import \
    get_grid, set_grid, create_grid, create_grid2d, \
    meshcells, meshbounds, meshgrid, bounds1d, bounds2d, \
    isrect, curv2rect, isgrid, get_xy, create_axes2d
    
from vacumm.misc.grid.regridding import \
    regrid1d, regrid2d, transect, griddata, CDATRegridder, grid2xy, interp1d, interp1d, \
    nearest1d

from vacumm.misc.grid.basemap import \
    merc, create_map, get_proj, reset_cache
    
# - data

from vacumm.data import setup_dataset, DS

from vacumm.data.cf import format_var, format_axis

# - diag

from vacumm.diag.thermdyn import density, mixed_layer_depth
