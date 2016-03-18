"""Module for quickly loading VACUMM essentials"""

# System
import os
import sys
import shutil
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

import MV2
import cdms2
import cdtime
import cdutil
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

from vacumm.config import (
    edit_user_conf_file, print_config, get_config_value, data_sample
    )

# - misc

from vacumm.misc.misc import (
    lonlab, latlab, deplab, kwfilter, broadcast, ArgList, closeto, dict_merge,
    dict_check_defaults, cp_atts, dict_filter,  cp_props, numod, set_atts, get_atts,
    set_lang, geodir, auto_scale, MV2_concatenate, MV2_axisConcatenate,
    zoombox, squeeze_variable, set_lang_fr, scalebox, checkdir,
    )

from vacumm.misc.axes import (
    create_lon, create_lat, create_time, create_depth,
    create_dep, islon, islat, isdep, isaxis, istime, create_axis,
    )

from vacumm.misc.atime import (
    comptime, strftime, strptime, Intervals, IterDates, now,
    time_split, time_split_nmax, to_utc, tz_to_tz, utc_to_paris, paris_to_utc,
    lindates, ch_units, add_margin, round_date, midnight_date, midnight_interval,
    mpl, tic, toc, pat2freq, pat2glob, has_time_pattern,
    is_interval, is_time, is_axistime, is_cdtime, is_datetime, is_in_range,
    is_numtime, is_strtime, itv_intersect, itv_union, julday,
    are_valid_units,
    )


from vacumm.misc.bases import psinfo, code_base_name, code_file_name, code_dir_name

from vacumm.misc.io import (
    list_forecast_files, netcdf3, netcdf4, ncread_files, ncread_var,
    ncread_axis, ncget_var, ncget_axis, ncfind_var, ncfind_axis, NcIterBestEstimate,
    ncget_fgrid
    )

from vacumm.misc.filters import (
    generic1d, generic2d, shapiro1d, shapiro2d, blackman1d,
    bartlett1d, kaiser1d, hanning1d, blackman1d, gaussian1d, gaussian2d
    )

from vacumm.misc.remote import InputWorkFiles, OutputWorkFile

from vacumm.misc.log import Logger

from vacumm.misc.config import ConfigManager, cfgargparse

from vacumm.misc.stats import StatAccum, qtmax, qtmin, qtminmax

from vacumm.misc.phys.units import (
    convert_units, kt2ms, ms2kt, deg2m, m2deg, ms2bf, dms2deg, deg2dms, mph2ms,
    ms2mph, tometric, kel2degc, degc2kel, basic_proj,
    )


# - plot

from vacumm.misc.plot import (
    map2, section2, hov2, curve2, bar2, plot2d, stick2, make_movie,
    savefigs, add_grid, xhide, yhide, xrotate, yrotate, add_key, taylor, dtaylor, dtarget,
    add_shadow, add_glow, add_map_lines, add_map_line, add_map_point, minimap,
    add_logo, add_left_label, add_right_label, add_top_label, add_bottom_label
    )

from vacumm.misc.color import (
    plot_cmap, show_cmap, get_cmap, simple_colors,
    cmap_rs, cmap_srs, cmap_custom, Scalar2RGB, darken, whiten, cmaps_mpl,
    cmap_gmt, cmaps_registered, cmaps_vacumm, print_cmaps_gmt, plot_cmaps,
    anamorph_cmap, discretize_cmap, StepsNorm
    )

# - grid

from vacumm.misc.grid.misc import (
    get_grid, set_grid, create_grid, create_grid2d, resol, monotonic, isdepthup, depth2dz,
    meshcells, meshbounds, meshgrid, bounds1d, bounds2d, coord2slice, isregular,
    isrect, curv2rect, isgrid, get_xy, create_axes2d, gridsel, varsel, xshift, rotate_grid,
    get_closest, get_closest_depth, makedepthup, get_axis, transect_specs,
    get_axis_slices, get_distances, dz2depth, meshweights,
    )

from vacumm.misc.grid.regridding import (
    regrid1d, regrid2d, interp1d, interp2d, cellave1d, cellave2d, cellerr1d,
    cargen, xy2xy, shift1d, shift2d, extend1d, extend2d, regrid_method,
    transect, griddata, CDATRegridder, grid2xy, cubic1d,
    nearest1d, shiftgrid, extendgrid, fill1d, fill2d
    )

from vacumm.misc.grid.basemap import (
    merc, create_map, get_proj, reset_cache
    )

from vacumm.misc.grid.masking import (
    erode_coast, polygon_mask, GetLakes, get_coast, polygons, polygon_select, zcompress,
    envelop, get_coastal_indices, grid_envelop, create_polygon, clip_shape, proj_shape,
    plot_polygon, merge_masks
    )

# - data

from vacumm.data import setup_dataset, DS

from vacumm.data.cf import (
    format_var, format_axis, format_grid, match_var, change_loc, match_obj,
    set_loc, get_loc
    )

from vacumm.data.misc.sigma import NcSigma, sigma2depths, SigmaGeneralized, SigmaStandard

from vacumm.data.misc.arakawa import (
    ArakawaGrid, ARAKAWA_LOCATIONS as arakawa_locations,
    CGrid, AGrid, ArakawaGridTransfer
    )

# - diag

from vacumm.diag.thermdyn import density, mixed_layer_depth

from vacumm.diag.dynamics import (
    barotropic_geostrophic_velocity, coriolis_parameter, kinetic_energy
    )

