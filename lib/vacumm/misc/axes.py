# -*- coding: utf8 -*-
"""Generic tools dealing with information about longitude,
latitude, depth and time axes


.. seealso::

    Tutorials: :ref:`user.tut.misc.variables.axes`
"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2018)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model of IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
from __future__ import absolute_import
import re
from re import match
from numbers import Number
from collections import Sequence
from warnings import warn
import six
from six.moves import range
from six.moves import zip


import numpy as N
import cdms2
import MV2
from cdms2.axis import AbstractAxis, FileAxis
from cdms2.coord import TransientAxis2D

from vacumm import VACUMMError
from .misc import (check_def_atts, dict_check_defaults, match_atts,
                   set_atts, get_atts)
from .cf import CF_AXIS_SPECS, cf2atts, cf2search

__all__ = ['isaxis', 'islon', 'islat', 'islev', 'islevel',
           'islongitude', 'islatitude',
           'isdep', 'isdepth', 'istime',
           'isforecast', 'isalt',
           'check_axes', 'is_geo_axis', 'check_axis',
           'get_axis_type', 'check_id',
           'get_checker', 'is_geo_axis_type', 'axis_type',
           'create_time', 'create_lon', 'create_lat',
           'create_dep', 'create_depth', 'create_forecast',
           'create_altitude', 'create_alt',
           'guess_timeid', 'get_order', 'set_order',
           'order_match', 'merge_orders',
           'check_order',  'create_axis',
           'create_axes2d', 'axes2d', 'num2axes2d']


def isaxis(obj):
    if hasattr(obj, 'isAbstractCoordinate') and obj.isAbstractCoordinate():
        return True
    return isinstance(obj, (AbstractAxis, FileAxis, TransientAxis2D))


def islongitude(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
                **attchecks):
    """Check if an object is of longitude type"""
    specs = CF_AXIS_SPECS['lon']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    return is_geo_axis_type(obj, 'x', defaults=defaults, ro=ro,
                            checkatts=checkatts,
                            checkaxis=checkaxis, **attchecks)


islon = islongitude


def islatitude(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
               **attchecks):
    """Check if an object is of latitude type"""
    specs = CF_AXIS_SPECS['lat']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    return is_geo_axis_type(obj, 'y', defaults=defaults, ro=ro,
                            checkatts=checkatts,
                            checkaxis=checkaxis, **attchecks)


islat = islatitude


def islevel(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
            **attchecks):
    """Check if an object is of level type"""
    specs = CF_AXIS_SPECS['level']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    if not checkaxis:
        if 'units' in attchecks:
            del attchecks['units']
        if 'long_name' in attchecks:
            del attchecks['long_name']
    return is_geo_axis_type(obj, 'z', defaults=defaults, ro=ro,
                            checkatts=checkatts,
                            checkaxis=checkaxis, **attchecks)


islev = islevel


def isdepth(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
            **attchecks):
    """Check if an object is of depth type"""
    specs = CF_AXIS_SPECS['depth']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    if not checkaxis:
        if 'units' in attchecks:
            del attchecks['units']
        if 'long_name' in attchecks:
            del attchecks['long_name']
    return is_geo_axis_type(obj, 'z', defaults=defaults, ro=ro,
                            checkatts=checkatts,
                            checkaxis=checkaxis, **attchecks)


isdep = isdepth


def isaltitude(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
               **attchecks):
    """Check if an object is of altitude type"""
    specs = CF_AXIS_SPECS['altitude']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    if not checkaxis:
        if 'units' in attchecks:
            del attchecks['units']
        if 'long_name' in attchecks:
            del attchecks['long_name']
    return is_geo_axis_type(obj, 'z', defaults=defaults, ro=ro,
                            checkatts=checkatts,
                            checkaxis=checkaxis, **attchecks)


isalt = isaltitude


def istime(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
           **attchecks):
    """Check if an object is of time type"""
    specs = CF_AXIS_SPECS['time']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    units = attchecks.setdefault('units', [])
    if not isinstance(units, (list, tuple)):
        units = [units]
    from .atime import are_good_units
    units.append(are_good_units)
    attchecks['units'] = units
    myistime = is_geo_axis_type(obj, 't', defaults=defaults, ro=ro,
                                checkatts=checkatts,
                                checkaxis=checkaxis, **attchecks)
    if myistime and not ro:
        try:
            obj.calendar = 'gregorian'
        except Exception:
            pass
    return myistime


def isforecast(obj, defaults=None, ro=False, checkaxis=True, checkatts=True,
               **attchecks):
    """Check if an object is of level type"""
    specs = CF_AXIS_SPECS['forecast']
    if defaults is None:
        defaults = cf2atts(specs)
    dict_check_defaults(attchecks, **cf2search(specs))
    return is_geo_axis_type(obj, 'f', defaults=defaults, ro=ro,
                            checkatts=checkatts,
                            checkaxis=checkaxis, **attchecks)


def get_checker(name):
    """Get the function that checks if an axis of required type

    Parameters
    ----------
    name:
        Generic name of the axis.

    Returns
    -------
    :func:`islon`, :func:`islat`, :func:`islevel`, :func:`istime`
    or raises :exc:`TypeError`

    Example
    -------
    >>> get_checker('x')
    >>> get_checker('lon')(myaxis)
    """
    errmsg = 'Input must be a generic name axis, like "x" or "lon"'
    if not isinstance(name, six.string_types):
        raise TypeError(errmsg)
    name = name.lower()
    if name in ('x', 'lon', 'longitude'):
        return islongitude
    if name in ('y', 'lat', 'latitude'):
        return islatitude
    if name in ('dep', 'depth'):
        return isdepth
    if name in ('alt', 'altitude'):
        return isaltitude
    if name in ('lev', 'level'):
        return islevel
    if name in ('z', ):
        return lambda *args, **kwargs: (
                islevel(*args, **kwargs) or
                isdepth(*args, **kwargs) or
                isaltitude(*args, **kwargs))
    if name in ('t', 'time'):
        return istime
    if name in ('f', 'forecast'):
        return isforecast
    raise TypeError(errmsg)


def is_geo_axis_type(obj, atype, defaults=None, ro=False, checkaxis=True,
                     checkatts=True, **attchecks):
    """Check if an object is of a specific type

    Parameters
    ----------
    axis:
        CDAT 1D or 2D axis.
    atype:
        Axis type as one of 'x', 'y', 'z', 't' or 'f'.
    ids: optional
        List od ids to check.
    ro: optional
        Read-only mode?
    checkatts: optional
        If False, do not check units and long_name attributes.
    attchecks: optional
        Extra keywords are attributes name and checklist that
        will checks using :func:`~vacumm.misc.misc.match_atts`.

    Return
    ------
    bool

    """
    if defaults is None:
        defaults = {}

    # Pure axis checks
    if checkaxis:

        # Is it an axis?
        if not isaxis(obj):
            return False

        # Use for instance obj.isLongitude()
        if atype in ['t', 'time']:
            name = 'Time'
        elif atype in ['f', 'forecast']:
            name = 'Forecast'
        elif atype in ['z', 'lev', 'dep', 'level', 'depth']:
            name = 'Level'
        elif atype in ['y', 'lat']:
            name = 'Latitude'
        elif atype in ['x', 'lon']:
            name = 'Longitude'
            atype = 'x'
        else:
            False
        if getattr(obj, 'axis', '-').lower() not in ['-', atype]:
            return False
        isfunc = getattr(obj, 'is' + name)

        def valfunc():
            designate = getattr(obj, 'designate'+name, None)
            if designate:
                designate()
            check_id(obj)
            check_def_atts(obj, **defaults)
        if isfunc():
            if not ro:
                valfunc()
            return True

        # Check axis attribute
        if getattr(obj, 'axis', '') == atype.upper():
            if not ro:
                valfunc()
            return True

    # Check from attributes
    if not checkatts:
        return False
    # TODO: merge with ncmatch_obj
    valid = match_atts(obj, attchecks, ignorecase=True,
                       transform=lambda ss: (re.compile(ss, re.I).match
                             if isinstance(ss, six.string_types) else None))
    if not valid:
        return False
    if not ro:
        if not ro and checkaxis:
            valfunc()

    return True


_isgeoaxis_ = is_geo_axis_type  # Backward compat


def check_axes(var, **kw):
    """Check the format of all axes of a cdms variable"""
    if cdms2.isVariable(var):
        for axis in var.getAxisList():
            check_axis(axis, **kw)


def is_geo_axis(axis, **kw):
    """Return True if axis is time, level, lat or lon"""
    return (istime(axis, **kw) or islev(axis, **kw)
            or islat(axis, **kw) or islon(axis, **kw))


def check_axis(axis, **kw):
    """Check the format an axis"""
    is_geo_axis(axis, ro=False, **kw)


def get_axis_type(axis, genname=False, **kw):
    """Return the axis type as a signle letter (CDAT standards): -, t, z, y, or x

    Parameters
    ----------
    axis:
        CDAT axis.
    genname: optional
        Return a generic name or None.
    **kw
        Other keywords are passed to checking functions (:func:`islon`...).

    Example
    -------
    >>> get_axis_type(create_time((5,),'days since 2000'))
    't'
    >>> get_axis_type(axis, genname=True, ro=True, checkatts=False)
    'time'
    """
    if islon(axis, **kw):
        at = "x"
    elif islat(axis, **kw):
        at = "y"
    elif islev(axis, **kw):
        at = "z"
    elif istime(axis, **kw):
        at = "t"
    else:
        at = '-'
    if not genname:
        return at
    if at == '-':
        return
    return {'x': 'lon', 'y': 'lat', 'z': 'level', 't': 'time'}[at]


axis_type = get_axis_type  # Backward compat


def check_id(axis, **kwargs):
    """Verify that an axis has a suitable id (not like 'axis_3' but 'lon')"""
    aliases = dict(X='lon', Y='lat', Z='depth', T='time')
    aliases.update(kwargs)
    if (hasattr(axis, 'axis') and match('^axis_\d+$', axis.id) is not None
            or match('^variable_\d+$', axis.id)):
        axis.id = aliases[axis.axis]


def create_axis(values, atype='-', **atts):
    """Quickly create a :mod:`cdms2` axis

    Parameters
    ----------
    values:
        Numerical values.
    atype: optional
        axis type within 'x','y','z','t','-' [default: '-']
    **atts
        Other keywords are passed as attributes to the axis.

    Example
    -------
    >>> lon = create_axis(N.arange(-10., 0, 2), 'x')
    >>> lon = create_axis((-10., 0, 2), 't', id='temps',
    ...     units='seconds since 2000')
    """
    if N.isscalar(values):
        values = [values]
    if isinstance(values, tuple) and len(values) < 4:
        values = N.arange(*values, **{'dtype': 'd'})
    if cdms2.isVariable(values):
        for item in values.attributes.items():
            atts.setdefault(*item)
        values = values.asma()
    if not isaxis(values):
        axis = cdms2.createAxis(values)
    else:
        axis = values
    if atype in CF_AXIS_SPECS:
        defaults = cf2atts(CF_AXIS_SPECS[atype])
        dict_check_defaults(atts, **defaults)
    for att, val in atts.items():
        setattr(axis, att, val)
    if atype in 'xyztf-':
        axis.axis = atype.upper()
        check_axis(axis)
        if axis.axis == '-':
            del axis.axis
    return axis


create = create_axis


def create_time(values, units=None, **atts):
    """Create a time axis

    Parameters
    ----------
    values:
        Numeric values, or list of date objects
        (:class:`~datetime.datetime`, :func:`~cdtime.comptime`,
        :func:`~cdtime.reltime`).
    units: optional
        Time units like 'days since 2000-01-01'.
    **atts
        Other keywords are passed as attributes to the axis.

    Note
    ----
    Units must be provided explicitly if no date are passed.

    Example
    -------
    >>> from vacumm.misc.atime import create_time
    >>> from datetime import datetime
    >>> import cdtime
    >>> taxis = create_time([1,2],units='months since 2000',
                            long_name='My time axis')
    >>> taxis = create_time(taxis)
    >>> create_time([datetime(2000,1,1),'2000-2-1'],units='months since 2000')
    >>> create_time([cdtime.reltime(1,'months since 2000'),
                     cdtime.comptime(2000,1)])
    """
    if istime(values):
        return values
    from .atime import are_valid_units, comptime, strftime
    for var in values, units:
        if hasattr(var, 'units'):
            units = var.units
            break
    if units is not None and not are_valid_units(units):
        raise AttributeError('Bad time units: "%s"' % units)

    istuple = isinstance(values, tuple)
    if not istuple or (istuple and len(values) > 3):
        if isinstance(values, str) or not isinstance(values, Sequence):
            if hasattr(values, 'next') and hasattr(values, '__iter__'):
                values = [v for v in values]
            else:
                values = [values]
        newvalues = []
        for value in values:
            if isinstance(value, Number):
                newvalues.append(value)
            else:
                if units is None and hasattr(value, 'units'):
                    units = value.units
                value = comptime(value)
                if units is None:
                    units = strftime('hours since %Y-%m-%d %H:%M:%S', value)
                newvalues.append(value.torel(units).value)
    else:
        newvalues = values

    if units is None:
        raise ValueError('Unable to guess units. You must specify them.')

    return create_axis(newvalues, 'time', units=units, **atts)


def create_forecast(values, **atts):
    """Create a forecast axis

    Parameters
    ----------
    values:
        Numeric values
    **atts
        Keywords are passed as attributes to the axis.

    Example
    -------
    >>> create_forceast(numpy.arange(-18., -5.), units='hours')

    """
    return create_axis(values, 'forecast', **atts)


def create_lon(values, **atts):
    """Create a longitude axis

    Parameters
    ----------
    values:
        Numeric values
    **atts
        Keywords are passed as attributes to the axis.

    Example
    -------
    >>> create_lon(numpy.arange(-18., -5.))
    >>> create_lon(numpy.arange(-18., -5.), long_name='original_longitude')

    """
    if (isinstance(values, N.ndarray) and len(values.shape) == 2 and
            not isaxis(values)):
        atts.setdefault('long_name', 'Longitude')
        return create_axes2d(x=values, lonid=atts.pop('id', None), xatts=atts)
    return create_axis(values, 'lon', **atts)


create_longitude = create_lon


def create_lat(values, **atts):
    """Create a latitude axis

    Parameters
    ----------
    values:
        Numeric values
    **atts
        Keywords are passed as attributes to the axis.

    Example
    -------
    >>> create_lat(numpy.arange(40., 48., 1.5))
    >>> create_lat(numpy.arange(40., 48., 1.5), long_name='strange_latitude')

    """
    if (isinstance(values, N.ndarray) and len(values.shape) == 2
            and not isaxis(values)):
        atts.setdefault('long_name', 'Latitude')
        return create_axes2d(y=values, latid=atts.pop('id', None), yatts=atts)
    return create_axis(values, 'lat', **atts)


create_latitude = create_lat


def create_depth(values, **atts):
    """Create a depth axis

     Parameters
     ----------
     values:
        Numeric values
    **atts
        Keywords are passed as attributes to the axis.

    Example
    -------
    >>> create_dep(numpy.arange(-1000., -500., 10.))
    >>> create_dep(numpy.arange(-1000., -500., 10.), long_name='deep depths')

    """
    return create_axis(values, 'depth', **atts)


create_dep = create_depth


def create_altitude(values, **atts):
    """Create an altitude axis

     Parameters
     ----------
     values:
        Numeric values
    **atts
        Keywords are passed as attributes to the axis.

    Example
    -------
    >>> create_altitude(numpy.arange(0., 1000.))

    """
    return create_axis(values, 'altitude', **atts)


create_alt = create_altitude


def create_axes2d(x=None, y=None, bounds=False, numeric=False,
                  lonid=None, latid=None, iid='ni', jid='nj',
                  xatts=None, yatts=None,
                  xbounds2d=None, ybounds2d=None, nobounds=False):
    """Create 2D numerical of complete axes

    Example
    -------
    >>> lon2d, lat2d = create_axes2d(x2d, y2d)

    Parameters
    ----------
    xaxis: optional
        1D or 2D X axis
    xaxis: optional
        1D or 2D Y axis
    xatts: optional
        Attributes for output 2D X axis
    yatts: optional
        Attributes for output 2D Y axis
    xbounds2d: optional
        2D bounds of input xaxis
    ybounds2d: optional
        2D bounds of input yaxis
    nobounds: optional
        create (True) or not (False - default) bounds of axis
    numeric: optional
        Only return numerical values instead of complete axes
    bounds: optional
        Return extended axes positioned on bounds (useful for pcolor).
        Deprecated: use :func:`meshbounds` instead.
    lonid: optional
        Id of longitude axis [defaut='lon'].
    latid: optional
        Id of latitude axis [defaut='lat'].
    iid: optional
        Id of i axis [defaut='iid'].
    jid: optional
        Id of j axis [defaut='jid'].


    Return
    ------
    array(ny,nx), array(ny,nx)
        x and y 2D axes
    """
    if x is None:
        hasx = 0
    elif isinstance(x, TransientAxis2D):
        hasx = 2
    else:
        hasx = 1
    if y is None:
        hasy = 0
    elif isinstance(y, TransientAxis2D):
        hasy = 2
    else:
        hasy = 1

    # - ids
    xaxis1d = yaxis1d = None
    if hasx:
        lonid = getattr(x, 'id', lonid)
        if x[:].ndim == 2 and hasattr(x, 'getAxis'):
            xaxis1d = x.getAxis(-1)
            yaxis1d = x.getAxis(-2)
    if hasy:
        latid = getattr(y, 'id', latid)
        if y[:].ndim == 2 and hasattr(x, 'getAxis') and xaxis1d is None:
            xaxis1d = y.getAxis(-1)
            yaxis1d = y.getAxis(-2)

    # Numeric part
    if hasx:
        xn = N.asarray(x[:])
        if xn.ndim == 1 and y is None:
            raise VACUMMError("Can't create 2D from a single 1D X axis")
    if hasy:
        yn = N.asarray(y[:])
        if yn.ndim == 1 and x is None:
            raise VACUMMError("Can't create 2D from a single 1D Y axis")
    if hasx and hasy:
        from .grid import meshgrid
        xx, yy = meshgrid(xn[:], yn[:])
    else:
        xx = xn if hasx else None
        yy = yn if hasy else None
    if hasx:
        del xn
    if hasy:
        del yn
#    if hasx and hasy and xx.shape != yy.shape:
#        raise VACUMMError('Incomptible shape between 2D X and Y coordinates: %s != %s'%(xx.shape, yy.shape))
#    if bounds:
#        if hasx: xx = meshbounds(xx)
#        if hasy: yy = meshbounds(yy)
    if numeric:
        return xx, yy

    # 2D axes
    if hasx == 2:
        xaxis2d = x
    elif hasx:
        xaxis2d = TransientAxis2D(xx)
    if hasy == 2:
        yaxis2d = y
    else:
        yaxis2d = TransientAxis2D(yy)

    # 2D bounds
    if not nobounds:
        from .grid import bounds2d
        if hasx:
            if xbounds2d is not None:
                xaxis2d.setBounds(xbounds2d)
            elif xaxis2d.getBounds() is None:
                xaxis2d.setBounds(bounds2d(xx))
        if hasy:
            if ybounds2d is not None:
                yaxis2d.setBounds(ybounds2d)
            elif yaxis2d.getBounds() is None:
                yaxis2d.setBounds(bounds2d(yy))

    # 1D axes
    for axis2d in ([xaxis2d] if hasx else [])+([yaxis2d] if hasy else []):
        if xaxis1d is not None:
            axis2d.setAxis(-1, xaxis1d)
        elif iid is not None:
            axis2d.getAxis(-1).id = iid
        if yaxis1d is not None:
            axis2d.setAxis(-2, yaxis1d)
        elif jid is not None:
            axis2d.getAxis(-2).id = jid

    # Format
    if hasx:
        xaxis2d.id = lonid or 'lon'
#        xaxis2d = create_lon(xaxis2d, id=lonid or 'lon')
        xaxis2d.getAxis(1).designateLongitude()
        xaxis2d.getAxis(0).designateLatitude()
        if hasattr(xaxis2d, 'axis'):
            del xaxis2d.axis
        if xatts is not None:
            xa = get_atts(x)
            xa.update(xatts)
            set_atts(xaxis2d, xatts)
    if hasy:
        yaxis2d.id = latid or 'lat'
#        yaxis2d = create_lat(yaxis2d, id=latid or 'lat')
        yaxis2d.getAxis(1).designateLongitude()
        yaxis2d.getAxis(0).designateLatitude()
        if hasattr(yaxis2d, 'axis'):
            del yaxis2d.axis
        if yatts is not None:
            ya = get_atts(y)
            ya.update(yatts)
            set_atts(yaxis2d, yatts)

    # Output
    if not hasx and not hasy:
        return
    if not hasx:
        return yaxis2d
    if not hasy:
        return xaxis2d
    return xaxis2d, yaxis2d


def axes2d(*args, **kwargs):
    """Alias for :func:`create_axes2d`"""
    return create_axes2d(*args, **kwargs)


def num2axes2d(*args, **kwargs):
    """Alias for :func:`axes2d`"""
    return axes2d(*args, **kwargs)


def guess_timeid(ncfile, vnames=None):
    """Guess the id of the time axis in a netcdf file

    Parameters
    ----------
    ncfile:
        Netcdf file name or descriptor.
    vnames: optional
        List of variables to look for
        a time axis (defaults to all variables)

    Return
    ------
    The id as a string, or ``None`` if not found.

    Example
    -------
    >>> tid = guess_timeid('file.nc')
    >>> f = cdms2.open('file.nc')
    >>> tid = guess_timeid(f)
    >>> f.close()
    """
    from .io import NcFileObj
    nfo = NcFileObj(ncfile)
    if vnames is None:
        vnames = nfo.f.listvariables()
    for vv in vnames:
        time = nfo.f[vv].getTime()
        if time is not None:
            break
    nfo.close()
    return time.id if time is not None else None


def set_order(var, order, replace=False):
    """Restore axis order of cdms variable

    Parameters
    ----------
    var:
        A cdms array.
    order:
        A cdms order  string(like 'tx')
    replace:
        Erase existing axis types?

    Example
    -------
    >>> set_order(temp, 'tyx')
    """
    current_order = var.getOrder()
    assert len(current_order) == len(order), \
        'Specified order must have length %i' % len(current_order)
    for i, (co, o) in enumerate(zip(current_order, order.lower())):
        if co == o or (co != '-' and not replace):
            continue
        axis = var.getAxis(i)
        if o == 'x':
            axis.designateLongitude()
        elif o == 'y':
            axis.designateLatitude()
        elif o == 'z':
            axis.designateLevel()
        elif o == 't':
            axis.designateTime()
        elif o == '-' and hasattr(axis, 'axis'):
            del axis.axis

    return var


def get_order(var):
    """Enhanced version of getOrder() method that handles 2D axes

    Parameters
    ----------
    var:
        axis or cdms variable.

    Return
    ------
    string
        containing letters x, y, z, t or -

    Example
    -------
    >>> get_order(var)
    "-yx"
    """
    # Already an order
    if isinstance(var, six.string_types):
        return var.lower()

    # Axis
    if isaxis(var):
        order = get_axis_type(var)
        if len(var.shape) == 2 and order in 'xy':
            return 'yx'
        return order

    # Variable
    if not cdms2.isVariable(var):
        return '-'*len(var.shape)
    order = var.getOrder()
    if (getattr(var,  '_nogridorder', False) or
            '-' not in order[-2:]):
        return order
    if (var.getGrid() is not None and
            'z' not in order[-2:] and 't' not in order[-2:]):
        if order[-1] == '-' and 'x' not in order:
            #            lon = var.getLongitude()
            #            if len(lon.shape) == 2:
            order = order[:-1]+'x'
        if order[-2] == '-' and 'y' not in order:
            #            lat = var.getLatitude()
            #            if len(lat.shape) == 2:
            order = order[:-2]+'y'+order[-1]
    return order


def order_match(order1, order2, asscore=False, strict=False):
    """Check that to axis orders are compatible

    Parameters
    ----------
    order1/2:
        Order strings containing ``x``, ``y``, ``z``, ``t`` or ``-`` chars.
    asscore: optional
        Return the total level of matching, where, for one char:

        - 0: no matching,
        - 1: matching with ``-``,
        - 2: letters equal.

    strict: optional
        Be more strict.

            - ``False``: Not strict.
            - ``True`` or ``"both"``: Fail even with ``-``.
            - ``"left"`` or ``"right"``: Designate the reference order, where
              the other one is not allowed to be different,
              except when the former has a ``-``.

    Examples
    --------
    >>> order_match('y', 'x')
    False
    >>> order_match('x-', 'xy')
    True
    >>> order_match('x-', 'xy', strict="right")
    False

    """
    order1 = get_order(order1)
    order2 = get_order(order2)
    if len(order1) != len(order2):
        if asscore:
            return 0
        assert False, ('Both orders must have the same length '
                       '({}, {})'.format(order1, order2))
    score = 1
    if strict is True:
        strict = "both"
    for ic in range(len(order1)):
        o1 = order1[ic]
        o2 = order2[ic]
        if '-' in o1+o2:
            if o1 != o2 and (
                strict == "both" or
                (strict == 'left' and o1 != '-') or
                    (strict == 'right' and o2 != '-')):
                return 0 if asscore else False
        elif o1 == o2:
            score *= 2
        else:
            return 0 if asscore else False
    return score if asscore else True


def merge_orders(order1, order2, raiseerr=True):
    """Merge two axis orders

    When two orders doesn't have the same length,
    they are right adjusted.

    Examples
    --------
    >>> merge_orders('t-x', 'y-')
    'tyx', 'yx'
    >>> merge_orders('yx', 'tz')
    'yx', 'tz'
    >>> merge_orders(myvar, zaxis)
    """
    order1 = get_order(order1)
    order2 = get_order(order2)
    rev = slice(None, None, 1-2*int(len(order2) < len(order1)))
    order1, order2 = (order1, order2)[rev]

    # Inner loop
    ishift = 0
    n1 = len(order1)
    n2 = len(order2)
    n12 = n2-n1
    for i in range(n12+1):
        j = n12-i
        if order_match(order1, order2[j:j+n1]):
            i1 = 0
            i2 = j
            ll = n1
            break

    else:  # Outerloops

        for ishift in range(1, min(n1, n2)):
            ll = min(n1, n2)-ishift
            if order_match(order1[:ll], order2[ishift:ishift+ll]):
                i1 = 0
                i2 = ishift
                break
            if order_match(order2[:ll], order1[ishift:ishift+ll]):
                i1 = ishift
                i2 = 0
                break
        else:
            if raiseerr:
                raise VACUMMError('orders are incompatible and cannot be '
                                  'safely merged: {} {}'.format(
                                          *(order1, order2)[rev]))
            return (order1, order2)[rev]

    # Merge
    neworder1 = order1[:i1]
    neworder2 = order2[:i2]
    for i in range(ll):
        c1 = order1[i1+i]
        c2 = order2[i2+i]
        if c1 == c2 or c2 == '-':
            neworder1 += c1
            neworder2 += c1
        elif c1 == '-':
            neworder1 += c2
            neworder2 += c2
        else:
            if raiseerr:
                raise VACUMMError('orders are incompatible and cannot be '
                                  'safely merged: {} {}'.format(
                                          *(order1, order2)[rev]))
            return (order1, order2)[rev]
    neworder1 += order1[i1+ll:]
    neworder2 += order2[i2+ll:]

    # Check multiples
    for c in 'xyztdf':
        if neworder1.count(c) > 2 or neworder2.count(c) > 2:
            warn('Merging of orders ({} and {}) may have not '
                 'properly worked (multiple axes are '
                 'of the same type)'.format(order1, order2))

    return (neworder1, neworder2)[rev]


def check_order(var, allowed, vertical=None, copy=False, reorder=False,
                extended=None, getorder=False):
    """Check that the axis order of a variable is matches
    at least one the specifed valid orders

    Parameters
    ----------
    var:
        MV2 array.
    allowed:
        A single order string or a list. It should contain one or
        several of these letters:

        - ``x``: longitudes,
        - ``y``: latitudes,
        - ``z``: vertical levels,
        - ``t``: time axis,
        - ``d``: data values (ignored),
        - ``-``: any kind of axis.

    Return
    ------
    array
        ``var``, or ``var, order, reordered`` if **reorder** is True.
    """

    # Check allowed orders
    # - consistency
    if not isinstance(allowed, (list, tuple)):
        allowed = [allowed]
    else:
        allowed = list(allowed)
    withd = 'd' in allowed[0]

    def get_rank(o): return len(o.replace('d', ''))
    rank = get_rank(allowed[0])
    for order in allowed:
        try:
            cdms2.orderparse(order.lower().replace('d', ''))
        except Exception:
            raise VACUMMError("Wrong allowed order: "+order)
        if ('d' in order and not withd) or ('d' not in order and withd):
            raise VACUMMError(
                "'d' only partially present in allowed order: %s" % allowed)
        if get_rank(order) != rank:
            raise VACUMMError("Inconsistent ranks between allowed orders: %s" % [
                              get_rank(o) for o in allowed])
    # - check extended mode
    if extended is None:  # extended?
        re_upper = re.compile('[XYZT]').search
        for order in allowed:
            if re_upper(order) is not None:
                extended = True  # force extended mode
                break
        else:
            extended = False
    if extended is False:  # lower
        allowed = [order.lower() for order in allowed]
    else:  # add tolerance for lower case orders
        re_sub = re.compile('[xyzt]').sub
        allowed = allowed+[re_sub('-', order) for order in allowed]
    # - unique and lower case
    _, idx = N.unique(allowed, return_index=True)
    idx = N.sort(idx)
    allowed = N.array(allowed)[idx].tolist()
    # - restrict to vertical or horizontal (1D data only)
    if vertical is not None and len(allowed[0]) == 2:
        allowed = [oo for oo in allowed if oo[int(vertical)] == 'd']

    # Data order
    data_cdms_order = get_order(var)

    # Loop on allowed orders
    reordered = False
    for allowed_order in allowed:

        # Handle data case
        d = allowed_order.find('d')
        if d != -1:
            allowed_order = allowed_order.replace('d', '')  # pure axis

        # Check cdms order
        allowed_cdms_order = allowed_order.lower()  # lower case
        if order_match(data_cdms_order, allowed_cdms_order, strict='right'):
            break  # It is already good

        # Try to reorder
        if reorder:
            from .grid import get_axis, var2d
            try:
                reordered = cdms2.order2index(
                    var.getAxisList(), allowed_cdms_order)
                new_var = var.reorder(allowed_cdms_order)
                # 2D axes
                if (allowed_cdms_order[-1] == 'x' and
                        len(get_axis(new_var, -1).shape) == 2):
                    del var
                    var = var2d(new_var,
                                MV2.transpose(get_axis(new_var, 0)),
                                MV2.transpose(get_axis(new_var, -1)), copy=0)
                    set_order(new_var, allowed_cdms_order)
                else:
                    del var
                    var = new_var
                data_cdms_order = get_order(var)
                break  # No error so it worked and we leave

            except Exception:
                continue

    else:
        raise VACUMMError(
            'Wrong type of axes. Possible forms are: %s' % ', '.join(allowed))

    if not getorder:
        return var
    if d != -1:
        data_cdms_order = cdms2.orderparse(data_cdms_order)
        data_cdms_order.insert(d, 'd')
        data_cdms_order = ''.join(data_cdms_order)
    return var, data_cdms_order, reordered
