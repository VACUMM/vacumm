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

import numpy
npy = N = np = numpy
testing = N.testing
assert_allclose = testing.assert_allclose
assert_raises = testing.assert_raises

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

from vacumm import VACUMMError, VACUMMWarning, help as vchelp, vacumm_warning, vcwarn


# - config

from vacumm.config import (
    edit_user_conf_file, print_config, get_config_value, data_sample
    )

# - misc

from vacumm.misc.misc import (
    lonlab, latlab, deplab, kwfilter, broadcast, ArgList, closeto, dict_merge,
    dict_check_defaults, cp_atts, dict_filter,  cp_props, numod, set_atts, get_atts,
    set_lang, geodir, auto_scale, MV2_concatenate, MV2_axisConcatenate, phaselab,
    zoombox, squeeze_variable, set_lang_fr, scalebox, checkdir, create_selector,
    numod, deg2str, deg_from_dec, dict_filter_out, dict_copy_items, geo_scale,
    grow_depth, grow_lat, grow_variables, history, intersect, is_iterable,
    isempty, isnumber, lunique, main_geodir, selector2str, rm_html_tags,
    tunique, split_selector, xls_style, Att, splitidx, CaseChecker, check_case,
    indices2slices, filter_level_selector, filter_selector, match_atts,
    match_string)

from vacumm.misc.axes import (
    create_lon, create_lat, create_time, create_depth,
    create_dep, islon, islat, isdep, isaxis, istime, create_axis,
    check_axes, check_id, check_order, get_axis_type, get_checker,
    guess_timeid, is_geo_axis, is_geo_axis_type, merge_orders,
    order_match, set_order,
    )

from vacumm.misc.atime import (
    comptime, strftime, strptime, Intervals, IterDates, now, Gaps, DateSorter,
    time_split, time_split_nmax, to_utc, tz_to_tz, utc_to_paris, paris_to_utc,
    lindates, ch_units, add_margin, round_date, midnight_date, midnight_interval,
    mpl, tic, toc, pat2freq, pat2glob, has_time_pattern,
    is_interval, is_time, is_axistime, is_cdtime, is_datetime, is_in_time_interval,
    is_numtime, is_strtime, itv_intersect, itv_union, julday, reltime,
    are_valid_units, SpecialDateFormatter, check_range,
    daily, daily_bounds, day_of_the_year, unit_type, toc,
    hourly, hourly_bounds, hourly_exact, monthly, tsel2slice,
    reduce, yearly, are_same_units, ascii_to_num, detrend, from_utc,
    time_selector, selector, filter_time_selector,
    plot_dt, strtime, interp_clim,
    datetime as adatetime,
    add,
    )
datetime_ = adatetime
add_time = add
is_time_interval = is_interval

from vacumm.misc.bases import (
    psinfo, code_base_name, code_file_name, code_dir_name, describe, Class, Object,
    add_logging_proxies, func_name, get_noconflict_metaclass, get_prog_dir,
    remove_redundant, skip_redundant, stack_trace, classmaker,
    classinstancemethod, isfreezed, prog,
    )

from vacumm.misc.io import (
    list_forecast_files, netcdf3, netcdf4, ncread_files, ncread_var,
    ncread_axis, ncget_var, ncget_axis, ncfind_var, ncfind_axis, NcIterBestEstimate,
    ncget_fgrid, grib2nc, grib_get_names, ncfind_obj, ncget_grid,
    ncget_lon, ncget_lat, ncget_time, ncget_var, ncmatch_obj, ncread_best_estimate,
    XYZ, XYZMerger, Shapes, NcIterBestEstimate, NcIterBestEstimateError,
    NcFileObj, ColoredFormatter, TermColors,
    Logger as SimpleLogger,
    )

from vacumm.misc.filters import (
    generic1d, generic2d, shapiro1d, shapiro2d, blackman1d,
    bartlett1d, kaiser1d, hanning1d, blackman1d, gaussian1d, gaussian2d,
    deriv, deriv2d,
    )

from vacumm.misc.remote import InputWorkFiles, OutputWorkFile

from vacumm.misc.log import (Logger, ColoredStreamHandler, LoggerAdapter,
    StreamLogWrapper, get_str_level, get_int_level)

from vacumm.misc.config import (ConfigManager, cfgargparse, ConfigException,
    ValidationWarning, cfg2rst, cfgoptparse, filter_section, get_secnames,
    getspec, list_options, opt2rst, option2rst, print_short_help)

from vacumm.misc.stats import (StatAccum, qtmax, qtmin, qtminmax, ensrank,
    corr_proba, StatAccumError)

from vacumm.misc.phys.units import (
    convert_units, kt2ms, ms2kt, deg2m, m2deg, ms2bf, dms2deg, deg2dms, mph2ms,
    ms2mph, tometric, kel2degc, degc2kel, basic_proj, rad2deg, uuconvert,
    vect2dir, vect2mod, vect2moddir, moddir2vectx, moddir2vecty, moddir2vectxy,
    strfsize, strpsize, basic_proj,
    )

from vacumm.misc.phys.constants import (
    EARTH_RADIUS, R, G, GRAVITATIONAL_CONSTANT, M, EARTH_MASS, g, GRAVITY,
    )

from vacumm.misc.file import (
    mkdirs, mkfdirs, rollover, strfsize, strpsize, walk, xefind, xfind, find, efind,
    )

from vacumm.misc.namelist import Namelist

from vacumm.misc.xml import (XmlConfig, XmlConfigDict, XmlConfigList)

# - plot

from vacumm.misc.plot import (
    map2, section2, hov2, curve2, bar2, plot2d, stick2, make_movie,
    savefigs, add_grid, xhide, yhide, xrotate, yrotate, add_key, taylor, dtaylor, dtarget,
    add_shadow, add_glow, add_map_lines, add_map_line, add_map_point, minimap,
    add_logo, add_left_label, add_right_label, add_top_label, add_bottom_label,
    ellipsis, get_cls, hldays, rotate_xlabels, rotate_ylabels,
    scale_xlim, scale_ylim, wedge, set_major_locator, set_minor_locator,
    xdate, ydate, get_quiverkey_value, add_lightshading,
    )
curve = curve2
map = map2
hov = hov2
stick = stick2
section = section2

from vacumm.misc.color import (
    plot_cmap, show_cmap, get_cmap, simple_colors, RGB, RGBA,
    cmap_srs, Scalar2RGB, darken, whiten, cmaps_mpl,
    cmap_gmt, cmaps_registered, cmaps_vacumm, print_cmaps_gmt, plot_cmaps,
    anamorph_cmap, discretise_cmap, StepsNorm,
    cmap_ajete, cmap_ajets, cmap_bathy, cmap_br, cmap_bwr, cmap_bwre, cmap_custom,
    cmap_chla, cmap_currents, cmap_dynamic_cmyk_hymex, cmap_eke, cmap_jet,
    cmap_jete, cmap_jets, cmap_land, cmap_linear, cmap_magic, cmap_mg,
    cmap_ncview_rainbow, cmap_nice_gfdl, cmap_pe, cmap_previmer, cmap_previmer2,
    cmap_rainbow, cmap_rainbow_sst_hymex, cmap_rb, cmap_red_tau_hymex,
    cmap_regular_steps, cmap_rnb2_hymex, cmap_rs, cmap_smoothed_regular_steps,
    cmap_smoothed_steps, cmap_ss, cmap_ssec, cmap_steps, cmap_topo,
    cmap_white_centered_hymex, cmap_wjet, cmap_wjets, cmap_wr, cmap_wre,
    to_shadow, to_grey, saturate, desaturate, change_saturation,
    change_luminosity, change_value, pastelise,
    discretise_cmap as discretize_cmap,
    pastelise as pastelize,

    )

# - grid

from vacumm.misc.grid.misc import (
    get_grid, set_grid, create_grid, create_grid2d, resol, monotonic, isdepthup, depth2dz,
    meshcells, meshbounds, meshgrid, bounds1d, bounds2d, coord2slice, isregular,
    isrect, curv2rect, isgrid, get_xy, create_axes2d, gridsel, varsel, xshift, rotate_grid,
    get_closest, get_closest_depth, makedepthup, get_axis, transect_specs,
    get_axis_slices, get_distances, dz2depth, meshweights, axis1d_from_bounds,
    bounds2mesh, cells2grid, check_xy_shape, curv_grid, create_var2d,
    get_geo_area, get_grid_axes, get_resolution, get_zdim, isoslice, mask2ind,
    merge_axis_slice, merge_axis_slices, num2axes2d, t2uvgrids, xextend,
    xshift,
    )

from vacumm.misc.grid.regridding import (
    regrid1d, regrid2d, interp1d, interp2d, cellave1d, cellave2d, cellerr1d,
    cargen, xy2xy, shift1d, shift2d, extend1d, extend2d, regrid_method,
    transect, griddata, CDATRegridder, grid2xy, cubic1d,
    nearest1d, shiftgrid, extendgrid, fill1d, fill2d
    )

from vacumm.misc.grid.basemap import (
    merc, create_map, get_proj, reset_cache, gshhs_autores, GSHHS_RESLIST,
    RSHPERE_WGS84
    )

from vacumm.misc.grid.masking import (
    erode_coast, polygon_mask, GetLakes, get_coast, polygons, polygon_select, zcompress,
    envelop, get_coastal_indices, grid_envelop, create_polygon, clip_shape, proj_shape,
    plot_polygon, merge_masks
    )

from vacumm.misc.grid.kriging import (
    OrdinaryCloudKriger, variogram, variogram_fit, variogram_model,
    variogram_model_type, variogram_multifit, cloud_split,
    KrigingError, VariogramModel, VariogramModelError, SimpleCloudKriger,
    )

# - data

from vacumm.data import setup_dataset, DS

from vacumm.data.cf import (
    format_var, format_axis, format_grid, match_var, change_loc, match_obj,
    set_loc, get_loc, GENERIC_NAMES, GENERIC_AXIS_NAMES, GENERIC_VAR_NAMES,
    AXIS_SPECS, VAR_SPECS, GRID_SPECS, ARAKAWA_SUFFIXES, HIDDEN_CF_ATTS,
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

