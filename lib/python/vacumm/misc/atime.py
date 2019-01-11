# -*- coding: utf8 -*-
"""
Time utilities
"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
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
import os
import re
import time, datetime as DT
from operator import isNumberType, gt, ge, lt, le
from re import split as resplit,match, compile as recompile
import math

import numpy as N, MV2, numpy.ma as MA, cdtime, cdms2 as cdms, cdutil, cdms2
from cdms2.axis import CdtimeTypes,AbstractAxis,ComptimeType,ReltimeType,FileAxis
from dateutil.rrule import rrule, MO, TU, WE, TH, FR, SA, SU, YEARLY, \
     MONTHLY, WEEKLY, DAILY, HOURLY, MINUTELY, SECONDLY
from matplotlib.dates import DateFormatter, num2date, date2num
#from pytz import timezone, all_timezones



# Attributes
STR_UNIT_TYPES = ['years','months','days','hours','minutes','seconds']
RE_SPLIT_DATE = recompile(r'[ Z:T\-]')
RE_MATCH_TIME_PATTERN = r'[0-2]?\d?\d?\d(-[01]?\d(-[0-3]?\d([ TZ][0-2]?\d(:[0-6]?\d(:[0-6]?\d(\.\d+)?)?)?)?)?)?'
RE_MATCH_TIME = recompile(RE_MATCH_TIME_PATTERN+'$', re.I).match
RE_MATCH_UNITS_PATTERN = r'(%s) since '%'|'.join(STR_UNIT_TYPES)+RE_MATCH_TIME_PATTERN+'[ \w]*'
RE_MATCH_UNITS = recompile(RE_MATCH_UNITS_PATTERN+'$', re.I).match

#: Time units for CNES julian days
JULIANDAY_TIME_UNITS_CNES = 'days since 1950-01-01'

#: Time units for NASA julian days
JULIANDAY_TIME_UNITS_NASA = 'days since 1958-01-01'


__all__ = ['STR_UNIT_TYPES','RE_SPLIT_DATE','now', 'add', 'axis_add', 'add_time',
    'mpl', 'are_same_units', 'are_valid_units', 'ch_units', 'comptime', 'reltime', 'datetime',
    'is_cdtime', 'is_reltime', 'is_comptime', 'is_datetime', 'check_range', 'is_in_time_interval',
    'num_to_ascii', 'Gaps', 'unit_type', 'get_dt', 'compress', 'plot_dt', 'reduce', 'yearly',
    'monthly', 'hourly', 'daily', 'hourly_exact', 'trend', 'detrend', 'strftime', 'strptime',
    'tz_to_tz', 'from_utc', 'to_utc', 'paris_to_utc', 'DateSorter', 'SpecialDateFormatter',
    'interp', 'is_time', 'round_date', 'Intervals', 'utc_to_paris', 'ascii_to_num',
    'lindates', 'itv_intersect', 'itv_union','day_of_the_year', 'pat2freq',  'strtime',
    'is_strtime', 'time_type', 'is_axistime', 'notz',  'IterDates', 'numtime',  'is_numtime',
    'pat2glob', 'midnight_date', 'midnight_interval','reduce_old', 'daily_bounds',
    'hourly_bounds', 'time_split', 'time_split_nmax', 'add_margin', 'fixcomptime',
    'is_interval', 'has_time_pattern', 'tsel2slice', 'tic', 'toc',
    'filter_time_selector','time_selector', 'selector',
    'julday', 'interp_clim', 'round_interval']

__all__.sort()

def pat2freq(pattern):
    """Get the maximal frequency ("days", "hours", etc) associated with a date pattern

    :Params:

        - **pattern**: A string containg date patterns (like '%d' or '%%d')

    :Return:

        One of ``None``, ``"seconds"``, ``"minutes"``, ``"hours"``,
        ``"days"``, ``"months"``, ``"yaers"``.

    :Example:

    >>> print pat2freq('%d/%m/%Y')
    days
    """

    # Patterns
    pats = (
        ("seconds", ['S', 'T', 'X']),
        ('minutes', ['M', 'R']),
        ('hours', ['H']),
        ('days', ['d', 'D', 'F', 'a', 'A', 'e', 'j', 'u', 'w']),
        ('months', ['b', 'B', 'h', 'm']),
        ('years', ['C', 'g', 'G', 'y', 'Y']),
    )

    # Search
    for freq, pp in pats:
        for p in pp:
            if '%'+p in pattern:
                return freq

def pat2glob(pattern, esc=False):
    """Convert a date pattern to a rough global matching pattern

    .. note::

        A global matching pattern is a UNIX file pattern.
        See :func:`glob.glob` for more information.

    .. warning::

        The output global pattern is NOT strict, and thus
        may embed more than requested by the date pattern.

    :Params:

        - **pattern**: A string containg date patterns (like '%d' or '%%d')
        - **esc**, optional: If ``True``, date patterns are escaped (like '%%d')
    """
    pats = {
        'a':    '???',
        'A':    '*',
        'b':    '*',
        'B':    '???',
        'c':    '*',
        'd':    '[0-3][0-9]',
        'H':    '[012][0-9]',
        'I':    '[01][0-9]',
        'j':    '[0-3][0-9][0-9]',
        'm':    '[0-1][0-9]',
        'M':    '[0-5][0-9]',
        'p':    '[AP]M',
        'S':    '[0-6][0-9]',
        'U':    '[0-5][0-9]',
        'w':    '[0-6]',
        'W':    '[0-5][0-9]',
        'x':    '*',
        'X':    '*',
        'y':    '[0-9][0-9]',
        'Y':    '[0-2][0-9][0-9][0-9]',
    }
    pc = '%%' if esc else '%'
    for d, g in pats.items():
        pattern = pattern.replace(pc+d, g)
    return pattern



def day_of_the_year(mytime):
    """Compute the day of the year of a specified date"""
    ctime = comptime(mytime)
    ctyear = cdtime.comptime(ctime.year)
    tunits = 'days since %s'%ctyear
    return ctime.torel(tunits).value+1

def lindates(first, last, incr, units=None):
    """Create a list of linearly incrementing dates

    :Params:

        - **first**: first date
        - **last**: first date
        - **incr**: increment step OR number of steps if units is None
        - **units**: units like "days" (see :func:`unit_type`) or None

    :Example:

        >>> dates = lindates('2000', '2002-05', 3, 'months')
        >>> dates = lindates('2000', '2002-05', 8)

    :Return:

        - list of :func:`cdtime.comptime` dates
    """
    ct0 = comptime(first)
    ct1 = comptime(last)
    tunits = 'hours since 2000'
    rt0 = ct0.torel(tunits).value
    rt1 = ct1.torel(tunits).value
    if rt1 < rt0:
        return []

    # Fixed number of steps
    if units is None:
        tt = N.linspace(rt0,  rt1, int(incr))
        return [cdtime.reltime(t, tunits).tocomp() for t in tt]

    # Fixed step
    units = unit_type(units)
    dates = [ct0]
    while dates[-1].torel(tunits).value < rt1:
        dates.append(add_time(dates[-1], incr, units))
    if dates[-1].torel(tunits).value > rt1:
        del dates[-1]
    return dates

def now(utc=False):
    """Return current time as :func:`cdtime.comptime`

    :Params:

        - **utc**, optional: Return UTC time [default: False]"""
    if utc:
        return comptime(DT.datetime.utcnow())
    return comptime(DT.datetime.now())



def add_time(mytime, amount, units=None, copy=True):
    """Add value to time

    :Params:

        - **mytime**: Target time, list of times or time axis.
        - **amount**: Value to add.
        - **units**, optional: Units, needed if ``mytime`` is not a time axis.
        - **copy**, optional: Copy list or time axis?

    Example:

        >>> add_time('2008-02', 3, 'days')
        '2008-2-4 0:0:0.0'
    """

    # Time axis
    if istime(mytime):
        if copy:
            mytime = mytime.clone()
        if units is None:
            mytime[:] = mytime[:] + amount
            return mytime
        assert hasattr(mytime,'units'), 'Your time axis must have an "units" attribute'
        mytime[:] = N.array([add_time(ct, amount, units).value for ct in mytime.asRelativeTime()])
        return mytime

    # Handle lists
    if not isinstance(mytime, list): copy = False
    LH = _LH_(mytime)
    intimes = LH.get()
    outtimes = list(intimes) if copy else intimes
    for i, mytime in enumerate(intimes):

        # Guess type
        typeback = time_type(mytime,  out='func')
        if typeback == reltime:
            myUnits = mytime.units
        mytime = comptime(mytime)

        # Units
        units = unit_type(units)
        assert units is not None, 'You must provide valid time units'

        # Add
        newtime = mytime.add(amount, units)

        # Correction
        if units != cdtime.Second:
            mytimes = [mytime.year,mytime.month,mytime.day,
                mytime.hour,mytime.minute,mytime.second]
            newtimes = [newtime.year,newtime.month,newtime.day,
                newtime.hour,newtime.minute,newtime.second]
            myunits = [cdtime.Year,cdtime.Month,cdtime.Day,
                cdtime.Hour,cdtime.Minute,cdtime.Second]
            iunits = myunits.index(units)
            if iunits < 2:
                newtimes[iunits+1:] = mytimes[iunits+1:]
            newtime = cdtime.comptime(*newtimes)

        # Back to correct type
        if typeback == reltime:
            mytime = typeback(newtime, myUnits)
        else:
            mytime = typeback(newtime)
        outtimes[i] = mytime

    return LH.put(outtimes)


# Obsolete :
add = add_time
axis_add = add_time



def mpl(mytimes, copy=True):
    """Convert to Matplotlib units.

    Time objects (or list of them) are converted to
    pure numerical values.

    :Params:

        - **mytimes**: Time object (or list), or time axis.
        - **copy**, optional: Copy time axis or not?

    .. note::

        This function does not use :func:`ch_units`.
        It uses Matplotlib internal function
        :func:`matplotlib.dates.date2num` instead.

    :Example:

        >>> taxis = create_time((0, 10.), 'hours since 2000-01-01')
        >>> maxis = mpl(taxis)
        >>> print maxis[0]
        730120.0
        >>> print maxis.units
        days since 0001
    """
    # Variable
    if cdms2.isVariable(mytimes):
        taxis = mytimes.getTime()
        if taxis is not None:
            taxis = mpl(taxis, copy=copy)
            mytimes.setAxis(mytimes.getOrder().index('t'), taxis)
        return mytimes

    # Time axis
    if istime(mytimes):
        if copy: mytimes = mytimes.clone()
        mytimes._data_ = mytimes._data_ .astype('d')
        if hasattr(mytimes, 'matplotlib_compatible') and mytimes.matplotlib_compatible:
            return mytimes
        # Values
        mytimes.assignValue(mytimes.getValue().astype('d'))
        datetimes = datetime(mytimes)
        mytimes[:] = [date2num(dt) for dt in datetimes]
        # Units
        ct0 = comptime(datetimes[0])
        del datetimes
        rt0 = ct0.torel('days since %s'%ct0)
        ctref = ct0.add(rt0.value-mytimes[0], cdtime.Days)
        mytimes.units = 'days since '+str(ctref)
        mytimes.matplotlib_compatible = True
        return mytimes

    # Time or list
    mytimes = datetime(mytimes)
    LH = _LH_(mytimes)
    intimes = LH.get()
    return LH.put([date2num(dt) for dt in datetime(intimes)])



# Obsolete :
axis_to_mpl = mpl
to_mpl = mpl

def are_same_units(units1, units2):
    """Compare time units

    >>> are_same_units('days since 1900-1', 'days since 1900-01-01 00:00:0.0')
    True

    """
    return cdtime.reltime(100,str(units1).lower())==cdtime.reltime(100, str(units2).lower())

def are_valid_units(units):
    """Check that units are well formatted

    :Example:

    >>> from vacumm.misc.atime import are_good_units
    >>> are_good_units('months since 2000')
    True
    >>> are_good_units('months since 2000-01-01 10h20')
    False
    """
    if not isinstance(units, basestring): return False
#    try:
#        cdtime.reltime(100,  units)
#        return True
#    except:
#        return False
#    for rem in ' UTC.', ' UTC':
#        units = units.replace(rem, '')
    return RE_MATCH_UNITS(units) is not None
are_good_units = are_valid_units
check_units = are_valid_units
def ch_units(mytimes, newunits, copy=True):
    """Change units of a CDAT time axis or a list (or single element) or cdtime times

    :Params:

        - **mytimes**: CDAT axis time object, a CDAT variable with a time axis,
          a valid time ot list of times.
        - **newunits**: New time units.
        - **copy**, optional: Create a new axis instead of updating
          the current one.

           .. note:: If ``mytimes is True``, then ``copy`` is set to ``False``.

    :Return:

        New time axis, or relatives times.

    :Example:

        >>> mytimes = ch_units(time_axis,relative=False,value=False,copy=True)

    .. warning::

        In the case of a single time object or a list of them, all object
        that is not :func:`cdtime.reltime` object will not be converted.

    """

    # If a variable with time axis
    if cdms.isVariable(mytimes):
        check_axes(mytimes)
        torder = mytimes.getOrder().find('t')
        if torder < 0:
            raise TypeError('Variable must have a valid time axis')
        return ch_units(mytimes.getTime(), newunits, copy=False)

    # Time axis
    if istime(mytimes):
        if copy: mytimes = mytimes.clone()
        mytimes._data_ = mytimes._data_ .astype('d')
        mytimes.toRelativeTime(newunits)
        return mytimes

    # Single or list of times
    if not isinstance(mytimes, list): copy = False
    LH = _LH_(mytimes)
    intimes = LH.get()
    outtimes = list(intimes) if copy else intimes
    for i, mytime in enumerate(intimes):
        if not is_reltime(mytime): continue
        outtimes[i] = mytime.tocomp().torel(newunits)
    return LH.put(outtimes)


def time_type(mytime, out='string', check=False):
    """Get the type of time

    :Params:

        - **mytime**: A time of one of the following types:
          :func:`cdtime.comptime`, :func:`cdtime.reltime`,
          :func:`datetime.datetime` or a date string.
        - **out**, optional: Output format, one of:

            - ``"string"``: Return one of ``"comptime"``, ``"reltime"``
              ``"datetime"``, ``"strtime"``, ``"numtime"``.
            - ``"type"``: Simply return result of ``type(mytime)``
            - ``"func"``: Return the function used to convert
              to this type : :func:`comptime`, :func:`reltime`,
              :func:`datetime`, :func:`numtime`.

    :Return: ``None`` if not a time, or see **out**.

    :Example:

            >>> time_type('2000-10')
            'str'
            >>> time_type(cdtime.comptime(2000,10), out='func')
            <function comptime at 0x31c4aa0>
    """
    for stype in 'strtime', 'comptime',  'reltime',  'datetime',  'numtime':
        if eval('is_'+stype)(mytime):
            if out=='string': return stype
            if out=='type': return type(mytime)
            return eval(stype)
    if check:
        raise TypeError, 'Time of wrong type: %s'%mytime


def _nummode_(nummode):
    if nummode is None:
        nummode = 'mpl'
    valid_nummodes = ['mpl', 'julday', 'cnes', 'nasa']
    nummode = str(nummode).lower()
    if nummode in valid_nummodes:
        if nummode=='mpl':
            return nummode
        if nummode=='julday':
            nummode = 'cnes'
        if nummode=='cnes':
            return JULIANDAY_TIME_UNITS_CNES
        else:
            return JULIANDAY_TIME_UNITS_NASA
    if are_valid_units(nummode):
        return nummode
    raise VACUMMError('nummode must be either a valid time units'
        ' string, or within: '+', '.join(valid_nummodes))

def comptime(mytime, nummode='mpl'):
    """Convert to :func:`cdtime.comptime` format

    :Params:

        - **mytime**: Time as string, :class:`~datetime.datetime`,
          :func:`cdtime.comptime`, :func:`cdtime.reltime`
          or a: mod:`cdms2` time axis.
        - **nummod**, optional: Numeric case mode.

            - ``"mpl"``: Converted using :func:`matplotlib.dates.num2date`
              and :class:`cdtime.comptime`.
            - ``"julday"``, ``"cnes"``, ``"nasa"``: Considered as juldian days.
            - Valid time units string: Converted using :class:`cdtime.reldate`.

    .. note::

        If not an :mod:`cdms2` time axis, first argument may be a list.
        In this case, a list is also returned.

    :Example:

    >>> from datetime import datetime ; import cdtime
    >>> from vacumm.misc.atime import comptime
    >>> comptime(datetime(2000,1,1))
    2000-1-1 0:0:0.0
    >>> comptime(cdtime.reltime(10,'days since 2008')).day
    10
    >>> comptime('1900-01-01').year
    1900

    :Sea also:

        :func:`reltime()` :func:`datetime()`
    """

    # Numeric mode
    tunits = _nummode_(nummode)

    # Time axis
    if is_axistime(mytime):
        mytime = mytime.asComponentTime()

    # Argument
    LH = _LH_(mytime)
    mytimes = LH.get()
    res = []

    for mytime in mytimes:

        # Component time
        if is_comptime(mytime):
            res.append(mytime)
            continue

        # Float or int
        if is_numtime(mytime):
            if tunits=='mpl':
                mytime = num2date(mytime)
            else:
                mytime = cdtime.reltime(mytime, tunits).tocomp()

        # Datetime
        if is_datetime(mytime):
            mytime = str(mytime)
        elif isinstance(mytime, time.struct_time):
            res.append(cdtime.comptime(*mytime[:6]))
            continue

        # String
        if isinstance(mytime, basestring):
            for c in 'TZ_':
                mytime = mytime.replace(c,' ')
            mytime = mytime.replace('/', '-').replace('h', ':')
            res.append(cdtime.s2c(mytime))
            continue

        # Relative time
        if is_reltime(mytime):
            res.append(mytime.tocomp())
            continue

        res.append(mytime)

    # Fix 60s
    res = fixcomptime(res)

    return LH.put(res)

def fixcomptime(mytime, decimals=3, copy=False):
    """Fix the 60s bug of :class:`cdtime.comptime` objects

    :Params:

        - **mytime**: Comptime or list of comptimes"""
    LH = _LH_(mytime)
    mytimes = LH.get()
    if copy and LH.listtype is list:
        mytimes = list(mytimes)
    for it, ct in enumerate(mytimes):
        if type(mytime) is not ComptimeType:
            continue
        if N.round(ct.second, decimals=decimals)==60:
            ct = cdtime.comptime(ct.year, ct.month, ct.day, ct.hour, ct.minute, 0)
            ct = ct.add(1, cdtime.Minute)
            mytimes[it] = ct
    return LH.put(mytimes)

def reltime(mytime, units):
    """Convert to func:`cdtime.reltime` format

    :Params:

        - **mytime**: Time as string, :class:`~datetime.datetime`,
          :func:`cdtime.comptime`, :func:`cdtime.reltime`, a number
          or a: mod:`cdms2` time axis.

    .. note::

        If not an :mod:`cdms2` time axis, first argument may be a list.
        In this case, a list is also returned.

    :Sea also:

        :func:`comptime` :func:`datetime` :func:`strtime`   :func:`numtime`
    """
    # Time axis
    if istime(mytime):
        if are_same_units(units, mytime.units):
            return mytime
        return mytime.toRelativeTime(mytime)

    # Other
    mytime = comptime(mytime)
    LH = _LH_(mytime)
    return LH.put([ct.torel(units) for ct in LH.get()])

def datetime(mytimes, nummode='mpl'):
    """Convert to :class:`datetime.datetime` format

    :Params:

        - **mytimes**: Time as string, :class:`~datetime.datetime`,
          :func:`cdtime.comptime`, :func:`cdtime.reltime`, a number,
          or a: mod:`cdms2` time axis.

    .. note::

        If not an :mod:`cdms2` time axis, first argument may be a list.
        In this case, a list is also returned.

    :Sea also:

        :func:`comptime` :func:`reltime`   :func:`strtime`    :func:`numtime`
    """
    # Numeric mode
    tunits = _nummode_(nummode)

    if istime(mytimes):
        mytimes = comptime(mytimes, nummode=tunits)
    LH = _LH_(mytimes)
    mytimes = LH.get()
    res = []
    for mytime in mytimes:

        # Numeric
        if is_numtime(mytime):
            if tunits=='mpl':
                res.append(num2date(mytime))
                continue
            else:
                mytime = cdtime.reltime(mytime, tunits).tocomp()

        # Tuple
        if isinstance(mytime, tuple):
            if not mytime: continue
            if len(mytime)==1: mytime += (1, )
            if len(mytime)==2: mytime += (1, )
            res.append(DT.datetime(*mytime))
            continue

        # Others
        ct = comptime(mytime)
        ct_seconds = int(math.floor(ct.second))
        ct_microseconds = int((ct.second-ct_seconds)*1e6)

        # Quick fix for stupid seconds at 60
        if ct_seconds == 60:
            ct_seconds = 0
            ct = ct.add(1,cdtime.Minute)
        res.append(DT.datetime(ct.year, ct.month, ct.day,
            ct.hour, ct.minute, ct_seconds, ct_microseconds))
    return LH.put(res)

def strtime(mytime):
    """Convert to valid string date

    :Params:

        - **mytime**: Time as string, :class:`~datetime.datetime`,
          :func:`cdtime.comptime`, :func:`cdtime.reltime`, a number,
          or a: mod:`cdms2` time axis.

    .. note::

        If not an :mod:`cdms2` time axis, first argument may be a list.
        In this case, a list is also returned.

    :Sea also:

        :func:`comptime` :func:`reltime`  :func:`datetime`  :func:`numtime`
    """
    ctimes = comptime(mytime)
    LH = _LH_(ctimes)
    ctimes = LH.get()
    return LH.put([str(ct) for ct in ctimes])

def numtime(mytime):
    """Convert to a numeric time using :func:`~matplotlib.dates.date2num`

    :Params:

        - **mytime**: Time as string, :class:`~datetime.datetime`,
          :func:`cdtime.comptime`, :func:`cdtime.reltime`, a number
          or a: mod:`cdms2` time axis.

    .. note::

        If not an :mod:`cdms2` time axis, first argument may be a list.
        In this case, a list is also returned.

    :Sea also:

        :func:`comptime` :func:`reltime`  :func:`datetime`  :func:`strtime`
    """
    dtimes = datetime(mytime)
    LH = _LH_(dtimes)
    dtimes = LH.get()
    return LH.put([date2num(dt) for dt in dtimes])

def julday(mytime, mode='cnes'):
    """Convert to a CNES julian days, i.e. days from 1950 (CNES) or 1958 (NASA)

    :Params:

        - **mytime**: Time as string, :class:`~datetime.datetime`,
          :func:`cdtime.comptime`, :func:`cdtime.reltime`, a number
          or a: mod:`cdms2` time axis.
        - **ref**, optional: computation mode

            - ``"cnes"``: days from 1958-01-01
            - ``"nasa"``: days from 1958-01-01

    .. note::

        If not an :mod:`cdms2` time axis, first argument may be a list.
        In this case, a list is also returned.

    :Sea also:

        :func:`comptime` :func:`reltime`  :func:`datetime`  :func:`strtime`
    """
    if mode is None:
        mode = 'cnes'
    else:
        mode = str(mode)
        assert mode in ['cnes', 'nasa']
    if mode=='cnes':
        tunits = JULIANDAY_TIME_UNITS_CNES
    else:
        tunits = JULIANDAY_TIME_UNITS_NASA
    rtimes = reltime(mytime, units=tunits)
    LH = _LH_(rtimes)
    rtimes = LH.get()
    return LH.put([rt.value for rt in rtimes])


def notz(mytime):
    """Suppres time zone

    :Params:

        - A :class:`datetime.datetime`

    :Return:

        - A :class:`datetime.datetime` instance with not TZ
    """
    return DT.datetime(*mytime.timetuple()[:6])

def is_cdtime(mytime):
    """Check if a mytime is a cdat time (from cdtime)

    Equivalent to::

        is_reltime(mytime) or is_comptime()

    :Params:

        - **mytime**: object to check

    :Example:

        >>> import cdtime
        >>> from datetime import datetime
        >>> from vacumm.misc.atime import is_cdtime
        >>> is_cdtime(cdtime.reltime(2,'days since 2000'))
        True
        >>> is_cdtime(cdtime.comptime(2000,2))
        True
        >>> is_cdtime(datetime(2000,2,1)
        False

    :See also:

        :func:`is_comptime()`:func:`is_reltime()`  :func:`is_time()`   :func:`is_datetime()`

    """

    return type(mytime) in CdtimeTypes

def is_reltime(mytime):
    """Check if a time is a cdat reltime (from :mod:`cdtime`)

    :Params:

        - **mytime**: object to check

    :Sea also:

        :func:`is_comptime()` :func:`is_reltime()` :func:`is_cdtime()`   :func:`is_time()`
    """

    return isinstance(mytime, ReltimeType)

def is_comptime(mytime):
    """Check if a time is a cdat comptime (from :mod:`cdtime`)

    :Params:

        - **mytime**: object to check

    :Sea also:

        :func:`is_datetime()` :func:`is_reltime()` :func:`is_cdtime()`  :func:`is_time()`
    """

    return isinstance(mytime, ComptimeType)

def is_datetime(mytime):
    """Check if a time is a :class:`datetime.datetime` time

    :Params:

        - **mytime**: object to check

    :Sea also:

        :func:`is_comptime()` :func:`is_reltime()` :func:`is_cdtime()`   :func:`is_time()`
    """

    return isinstance(mytime, DT.datetime)

def is_axistime(mytime):
    """Check if a time is a :mod:`cdms2` time axis

    Simple shortcut to :func:`~vacumm.misc.axes.istime`.

    :Params:

        - **mytime**: object to check

    :Sea also:

        :func:`~vacumm.misc.axes.istime` :func:`is_comptime()`
        :func:`is_reltime()` :func:`is_cdtime()`   :func:`is_time()`
    """

    return istime(mytime)

def is_strtime(mytime):
    """Check if a time is a valid string date

    :Params:

        - **mytime**: object to check

    :Sea also:

        :func:`is_datetime()` :func:`is_comptime()` :func:`is_reltime()`
        :func:`is_cdtime()`   :func:`is_time()`
    """
    if not isinstance(mytime, basestring): return False
    m = RE_MATCH_TIME(mytime)
    if m: return True
    #try:
        #cdtime.s2c(mytime)
        #return True
    #except:
        #return False


def is_numtime(mytime):
    """Simply check if a mytime is a number !

    :Params:

        - **mytime**: object to check

    :Sea also:

        :func:`is_datetime()` :func:`is_comptime()` :func:`is_reltime()`
        :func:`is_cdtime()`   :func:`is_time()`
    """
    return isinstance(mytime, (int, float)) or (N.isscalar(mytime) and N.isreal(mytime))



def check_range(this_time,time_range):
    """Check wether a time is before, within or after a range

    :Params:

    - **this_time**: time to check (string or cdat time)
    - **time_range**: 2(or 3)-element range (strings or cdat times) like ('1975',1980-10-01','co)

    :Example:

        >>> from vacumm.misc.atime import comptime, reltime, check_range
        >>> check_range('2000-12', ('2000-11', '2001'))
        0
        >>> check_range(comptime(2000), (comptime(2000), ('2000','20001','oc')))
    -1

    :Returns: -1 is before, 0 if within, 1 if after
    """

    # Range type
    time_range = list(time_range)
    if len(time_range) < 3:
        range_type = 'co'
    else:
        range_type = time_range[2].lower()

    # Convert to component time
    ctime = comptime(this_time)
    for i in 0,1:
        if type(time_range[i]) == type('s'):
            time_range[i] = comptime(time_range[i])

    # Before
    if ctime.cmp(time_range[0]) < int(range_type[0] is 'o'):
        return -1

    # After
    if ctime.cmp(time_range[1]) > -int(range_type[1] is 'o'):
        return 1

    # Within
    return 0




def is_in_time_interval(this_time,time_range):
    """Check if a time is in specified closed/open range

    :Params:

         - **this_time**: time to check (string or cdat time)
         - **time_range**: 2(or 3)-element range (strings or cdat times) like ('1975',1980-10-01','co)

    :Example:

        >>> is_in_range('2000-12', ('2000-11', '2001'))
        True

    :Sea also:

        :func:`check_range()`
    """

    return not check_range(this_time,time_range)

is_in_range = is_in_time_interval

def num_to_ascii(yyyy=1,mm=1,dd=1,hh=0,mn=0,ss=0):
    """Convert from [yyyy,mm,dd,hh,mn,ss] or component time  or relative time to 'yyyy-mm-dd hh:mn:ss'

    :Params:

        - *yyyy*: int year OR component OR relative time (cdms) [default: 1]
        - *mm*: month [default: 1]
        - *dd*: day [default: 1]
        - *hh*: hour [default: 0]
        - *mn*: minute [default: 0]
        - *ss*: second [default: 0]

    :Example:

        >>> num_to_ascii(month=2)
        0001-02-01 00:00:00
        >>> num_to_ascii(comptime(2000,10))
        2000-10-01 00:00:00
    """
    if is_cdtime(yyyy):
        yyyy = comptime(yyyy)
        ss = yyyy.second
        mn = yyyy.minute
        hh = yyyy.hour
        dd = yyyy.day
        mm = yyyy.month
        yyyy = yyyy.year
    return '%04i-%02i-%02i %02i:%02i:%02g' % \
           (yyyy,mm,dd,hh,mn,ss)



def ascii_to_num(ss):
    """ Convert from 'yyyy-mm-dd hh:mn:ss' to [yyyy,mm,dd,hh,mn,ss]

    :Example:

        >>> ascii_to_num('2000-01')
        2000, 1, 1, 0, 0, 0
    """
    sstr =  ss.split()

    yyyy = 1
    mm = 1
    dd = 1
    hh = 0
    mn = 0
    ss = 0
    ret = []

    date = ss.split()[0]
    for xx in sstr[0].split('-'):
        ret.append(int(xx))
    if len(ret) < 2:
        ret.append(mm)
    if len(ret) < 3:
        ret.append(dd)

    if len(sstr) > 1:
        for xx in sstr[1].split(':'):
            ret.append(int(xx))
    if len(ret) < 4: ret.append(hh)
    if len(ret) < 5: ret.append(mn)
    if len(ret) < 6: ret.append(ss)

    return ret



class Gaps(cdms.tvariable.TransientVariable):
    """Find time gaps in a variable

    A gap is defined as a missing time step or data.
    With this class, you can:

    - get the gaps as a variable (the object itself)
    - where 1 refers to the start of gap and -1 to the end
    - print them (:meth:`show()`)
    - plot them (:meth:`plot()`)
    - save them in a netcdf or an ascii file (:meth:`save()`)

    :Parameters:

        - **var**: A cdms variable with a cdms time axis
        - *dt*: Time step [default: minimal time step]
        - *tolerance*: There is a gap when dt varies  of more than tolerance*min(dt) [default: 0.1]
        - *keyparam*: verbose Verbose mode

    :Example:

        >>> import vacumm.misc as M
        >>> import MV2
        >>> time = M.axes.create_time([1,2,4,5],units='days since 2000')
        >>> var = MV2.arange

    """

    def __init__(self, var, tolerance=0.1, verbose=True, dt=None, **kwargs):

        # Check that we have a valid time axis
        if not cdms.isVariable(var):
            raise TypeError,'Your variable is not a cdms variable'
        mytime = var.getTime()
        if mytime is None:
            raise TypeError,'Your variable has no time axis'
        try:
            time_units =  mytime.units
        except:
            raise RuntimeError, 'Your time axis needs units'
        kwdt = kwfilter(kwargs, 'dt')

        if dt is None: dt = get_dt(mytime, **kwdt)

        # Time compression according to mask
        if len(var.shape) > 1:
            var = var(order='...t')
            mask = MV2.getmaskarray(var)
            while mask.ndim > 1:
                mask = N.logical_and.reduce(mask, axis=0)
        else:
            mask = MV2.getmaskarray(var)
        del var
        ctime = mytime.asComponentTime()
        self._ctbounds = [ctime[0], ctime[-1]]
        self._bounds = mpl(self._ctbounds)
        tt = mytime.getValue()[~mask].astype('d')

        # Derivatives
        dtf = N.diff(dt)
        dtr =  dtf/ dt # Relative time steps
        dtr = N.where((dtr % 1.)<tolerance,N.floor(dtr),dtr)
        dtr = N.where((dtr % 1.)>1.-tolerance,N.ceil(dtr),dtr)
        gg = MA.masked_less(dtr,1.+0.5*tolerance)-1.
        self.ngap = MA.count(gg)
        if not self.ngap:
            self.total_gap_len = 0.
            gaps = N.array([0,0])
            gapstime = tt[[0,-1]]
            self._gaps = 0.
        else:
            mask = MA.getmaskarray(gg)
            self._gaps = gg.compressed()
            self.total_gap_len = self._gaps.sum()
            ibefore = N.arange(len(tt)-1,dtype='l')[~mask]
            gaps = N.ones(2*self.ngap)
            gaps[1::2] = -1
            gapstime = N.zeros(2*self.ngap, dtype='d')
            ii = N.arange(self.ngap)*2
            gapstime[ii] = tt[ibefore] + 0.5*dt
            gapstime[ii+1] = tt[ibefore+1] - 0.5*dt
#            N.put(gapstime,ii,N.take(tt,ibefore)+0.5*dt)
#            N.put(gapstime,ii+1,N.take(tt,ibefore+1)-0.5*dt)


        # Time axis of gaps
        gapstime = cdms.createAxis(gapstime)
        gapstime.id = 'time'
        gapstime.long_name = 'Time'
        gapstime.units = mytime.units
        gapstime.designateTime(calendar=cdtime.DefaultCalendar)
        gapstime.dt = dt
        gapstime.tolerance = tolerance
        self._ctime = gapstime.asComponentTime()

        # Create cdms variable
        cdms.tvariable.TransientVariable.__init__(self, gaps, axes=(gapstime,), id='gaps')
        self.name = self.id
        self.long_name = 'Time gaps'
        self.setAxis(0,gapstime)

        if verbose:
            print 'Found %i gaps' % self.ngap
            self.show(headsep=None)

        del gaps,tt
        self._toplot = None


    def show(self,**kwargs):
        """Print out gaps

        Parameters are passed to col_printer()
        """
        if not self.ngap:
            print 'No gaps to print'
            return
        from .io import col_printer
        cp = col_printer([['LENGTH',7,'%i'],['START',22,'%s'],['END',22,'%s']],**kwargs)
        for igap in xrange(self.ngap):
            cp(self._gaps[igap],self._ctime[igap*2],self._ctime[igap*2+1])
        del cp

    def plot(self,show=True,savefig=None,figsize=(8,2.),title=None,color='red',
        subplots_adjust = dict(bottom=0.55,left=0.05,top=0.8),show_time_range=True,**kwargs):
        """Plot the gaps

        :Params:

            - *color*: Color of gaps [default: 'red']
            - *figure*: Show results on a figure if there are gaps [default: True]
            - *show*: Show the figure plotted if gaps are found [default: True]
            - *title*: Use this title for the figure
        """
        if not self.ngap:
            print 'No gaps to plot'
            return
        from plot import xdate,P
        if self._toplot is None:
            self._toplot = [0,]
            self._times = [self._bounds[0]]
            times = mpl(self.getTime())
            for igap in xrange(self.ngap):
                self._toplot.extend([0,1,1,0])
                self._times.extend([times[igap*2]]*2)
                self._times.extend([times[igap*2+1]]*2)
            self._toplot.append(0)
            self._times.append(self._bounds[1])
        if figsize is not None:
            P.figure(figsize=figsize)
        P.subplots_adjust(**subplots_adjust)
        P.fill(self._times,self._toplot,facecolor=color)
        P.gca().xaxis_date(None)
        P.gca().autoscale_view()
        P.legend(('Gaps',),numpoints=2,loc=0)
        xdate(**kwfilter(kwargs,'xdate'))
        P.xlim(min(self._times),max(self._times))
        P.ylim(0,1)
        P.yticks([])
        if title is None:
            title = 'Time gaps'
        P.title(title)
        if show_time_range:
            P.figtext(0,0,'Start: %s  /  End: %s'%tuple(self._ctbounds),color='#888888')
        P.grid(True)
        if savefig is not None:
            P.savefig(savefig)
        if show:
            P.show()
        else:
            P.close()

    def save(self,file):
        """Save gaps to a netcdf or an ascii file

        If the file name does not end with 'cdf' or 'nc',
        gaps are printed to an ascii file with the save format
        as displayed by (show())

        :Params:

            - **file**: Netcdf file name
        """
        if file.split('.')[-1] in ['nc','cdf']:
            f = cdms.open(file,'w')
            f.write(self)
            f.close()
        else:
            if not self.ngap:
                f = open(file,'w')
                f.write('NO GAPS\n')
                f.close()
                return
            self.show(file=file)


def unit_type(units, string_type=False, s=True, raiseerr=True):
    """Returns a type of units in a suitable form with checkings.

    :Params:

        - **units**: A valid string or cdtime type of units.
        - *string_type*: Returns a string type instead of a cdtime type [default: False]

    :Return: A valid string or cdtime type of units

    :Example:

        >>> unit_type('minutes')
        >>> unit_type(cdtime.Minutes,True)
    """
    if hasattr(units, 'units'):
        units = units.units.split()[0]
    if isinstance(units, basestring):
        units = units.lower()
        if not units.endswith('s'): units += 's'
    for tt in STR_UNIT_TYPES:
        this_cdt_type = eval('cdtime.'+tt.title())
        if units in [tt,this_cdt_type]:
            if string_type:
                return tt[:len(tt)-1+s]
            else:
                return this_cdt_type
    if raiseerr:  raise TypeError('Wrong type of units: "%s". Valid units are: %s'%(units, STR_UNIT_TYPES))




def get_dt(axis, units=None):
    """Returns the time steps of an axis.
    Value is computed according to the median time step,
    and returned in original or specified units.

    :Params:

        - **axis**: A time axis.
        - **units**, optional: Another valid unit type (cdtime or string)

    :Return: Time step

    :Example:

        >>> get_dt(var.getTime()
        >>> get_dt(var.getTime(), cdtime.Months)
        >>> get_dt(var.getTime(), 'months')

    :See also: :func:`unit_type()`
    """

    assert istime(axis), 'You must specify a valid time axis'
    assert len(axis) > 1, 'your time axis must have at least 2 elements'
    if units is not None:
        units = unit_type(units, True)

    # Same units
    time_units = axis.units.split()
    dt = N.median(N.diff(axis[:]))
    if units is None or time_units[0] == units:
        return dt

    # Other units
    time_units[0] = units
    time_units = ' '.join(time_units)
    t0 = cdtime.r2r(cdtime.reltime(axis[0], axis.units), time_units).value
    t1 = cdtime.r2r(cdtime.reltime(axis[0]+dt, axis.units), time_units).value
    return t1-t0


def compress(data):
    """Compress along time"""
    data = MV2.asarray(data)

    # Find time axis
    if data.getTime() is None:
        data.getAxis(0).designateTime(calendar=cdtime.DefaultCalendar)
    tt = data.getTime()
    nt = len(tt)

    # Time is first axis
    order = data.getOrder()
    data = data(order='t...')

    # Reference slab
    if data.ndim == 1:
        ref = data
    else:
        snap = data[0]
        ii = N.arange(len(snap.ravel()))
        mask = MA.getmaskarray(snap)
        ii0 = MA.masked_array(ii,mask=mask).compressed()[0]
        ref = MV2.reshape(data,(nt,N.size(snap)))[:,ii0]
    slab = N.arange(nt).compress(~MA.getmaskarray(ref)).tolist()

    # Select data
    res = MV2.take(data, slab)
    res.getAxis(0)[:] = tt[slab]
    cp_atts(data,res,id=True)
    cp_atts(tt,res.getAxis(0),id=True)
    if tt.getBounds() is not None:
        res.getAxis(0).setBounds(tt.getBounds()[slab])
    return res(order=order)


def plot_dt(file, time_axis=None, nice=False):
    """Plot axis time step of a netcdf file"""
    from pylab import plot_date,show,plot,ylim
    if time_axis is None:
        time_axis = 'time'
    f = cdms.open(file)
    tt = f.getAxis(time_axis).clone()
    tt[:] = tt[:].astype('d')
    f.close()
    dt = tt[1:]-tt[:-1]
    if nice:
        mm = axis_to_mpl(tt)
        plot_date(mm[:-1],dt)
    else:
        plot(dt)
    ylim(ymin=0.,ymax=(1.1*max(dt)))
    show()


def reduce_old(data, comp=True, fast=True):
    """Reduce a variable in time by performing time average on time step that have the same time or time bounds.

    - **data**: A cdms variable with first axis supposed to be time.
    - *comp*: Call to :meth:`compress()` before reducing [default: True]
    - *fast*: Convert to pure numpy before processing, then convert back to cdms variable [default: True]

    Return: The new variable on its new time axis.
    """

    from grid.misc import set_grid

    assert data.getTime() is not None, 'Your data must have a valid time axis'

    # Time compression
    if fast: comp = True
    if comp:
        data = compress(data)

    # Find interval for averages
    tt = data.getTime()
    nt = len(tt)
    tv = tt.getValue().astype('d')
    if tt.getBounds() is None:
        tt.setBounds(G.bounds1d(tt[:]))
    bb = tt.getBounds().astype('d')
    dt = tv[1:]-tv[:-1]
    dt = N.concatenate((dt,[1.]))
    dbb = bb[1:]-bb[:-1]
    dbb = N.concatenate((dbb,[[1.,1.]]))
    itv = MA.masked_where(MV2.equal(dt,0.),MA.arange(nt))
    if bb is not None:
        itv = MA.masked_where(N.logical_and(N.equal(dbb[:,0],0.),
          N.equal(dbb[:,1],0.)),itv)
    itv = itv.compressed()
    nt = len(itv)
    gg = data.getGrid()

    # Fast algo = use Numeric
    if fast:
        mydata = data.filled()
        missing_value = data.getMissing()
        AV = N
    else:
        mydata = data
        AV = MV
    sh = list(data.shape)
    sh[0] = nt
    adata = AV.zeros(sh,typecode=data.dtype.char,savespace=True)

    # Compute averages
    atime = N.zeros(nt).astype('d')
    if bb is not None:
        abounds = N.zeros((nt,2)).astype(bb.dtype.char)
    else:
        abounds = None
    it_first = 0
    for it,it_last in enumerate(itv):
        adata[it] = AV.average(mydata[it_first:it_last+1].astype(adata.dtype.char),0)
        if bb is None:
            atime[it] = tv[it_last]
        else:
            abounds[it] = bb[it_last]
            atime[it] = N.average(abounds[it],0)
        it_first = it_last+1

    # New time axis
    atime = cdms.createAxis(atime,bounds=abounds)
    cp_atts(tt,atime,id=True)
    atime.designateTime(calendar=cdtime.DefaultCalendar)

    # New variable
    if fast:
        adata = MV2.masked_object(adata,missing_value)
    adata.setAxis(0,atime)
    for i in xrange(1,data.rank()):
        adata.setAxis(i,data.getAxis(i))
    if data.getGrid() is not None:
        set_grid(adata,data.getGrid())
        #adata.setGrid(data.getGrid())
    cp_atts(data,adata,id=True)
    del data
    return adata

def yearly(data,**kwargs):
    """Convert to yearly means

    :Params:

        - **data**: A cdms variable.

    :Return: A cdms variable on a hourly time axis
    """
    assert data.getTime() is not None, 'Your data must have a valid time axis'

    hdata = MV2.array(data)
    cdutil.setTimeBoundsYearly(hdata)
    new_data = reduce(hdata,**kwargs)
    del hdata
    cdutil.setTimeBoundsYearly(new_data)
    return new_data

def monthly(data,**kwargs):
    """Convert to monthly means

    :Params:

        - **data**: A cdms variable.

    :Return: A cdms variable on a hourly time axis
    """
    assert data.getTime() is not None, 'Your data must have a valid time axis'

    hdata = MV2.array(data)
    cdutil.setTimeBoundsMonthly(hdata)
    new_data = reduce(hdata,**kwargs)
    del hdata
    cdutil.setTimeBoundsMonthly(new_data)
    return new_data

def hourly(data,frequency=24,**kwargs):
    """Convert to hourly means

    :Params:

        - **data**: A cdms variable.
        - *frequency* Used when different from hourly is requested [default: 24]

    :Return: A cdms variable on a hourly time axis
    """
    return daily(data,frequency=24,**kwargs)


def daily(data, hstart=0, **kwargs):
    """Convert to daily means

    :Params:

        - **data**: A cdms variable with a time axis.
        - *hstart*: First hour of daily intervals.

    :Example:

        >>> dsst = daily(sst, hstart=12)

    :See also: :func:`daily_bounds` :func:`reduce`

    :Return: A cdms variable on a daily time axis
    """
    taxis = data.getTime()
    assert taxis is not None, 'Your data must have a valid time axis'
    oldbounds = taxis.getBounds()
    taxis.setBounds(daily_bounds(taxis, hstart=hstart))
    varo = reduce(data)
    taxis.setBounds(oldbounds)
    return varo


def hourly_exact(data,time_units=None,maxgap=None, ctlims=None):
    """Linearly interpolate data at exact beginning of hours

    :Params:

        - **data**: Cdms variable.
        - *time_units*: Change time units.
        - *maxgap*: Maximal gap in hours to enable interpolation [default: None].
    """

    from grid.misc import set_grid

    # Check time
    taxis = data.getTime()
    assert taxis is not None, 'Your data must have a valid time axis'
    old_order = data.getOrder()
    data = data(order='t...')
    tbounds = G.bounds1d(taxis)
    if maxgap is not None:
        dt = get_dt(taxis,'hour')
        if dt > maxgap:
            warn('maxgap (%gh) is greater than your time step (%gh)'%(maxgap,dt))
            maxgap = None

    # Create the new hourly time axis
    if ctlims is None:
        ctlims = taxis.subAxis(0,len(taxis),len(taxis)-1).asComponentTime()
    else:
        ctlims = [comptime(ct) for ct in ctlims]
    if ctlims[0].minute == 0 and ctlims[0].second == 0:
        hctlim0 = ctlims[0]
    else:
        hctlim0 = cdtime.comptime(ctlims[0].year,ctlims[0].month,ctlims[0].day,ctlims[0].hour).add(1,cdtime.Hour)
    hunits = 'hours since '+str(hctlim0)
    nt = int(ctlims[1].torel(hunits).value-ctlims[0].torel(hunits).value)
    htaxis = cdms.createAxis(N.arange(nt,typecode='d'))
    htaxis.id = 'time'
    htaxis.long_name = 'Time'
    htaxis.units = hunits
    htaxis.designateTime(calendar=cdtime.DefaultCalendar)
    htvalue = htaxis.getValue()
    htaxis.setBounds(G.bounds1d(htvalue))
    cts = htaxis.asComponentTime()

    # Create the new variable
    sh = list(data.shape)
    sh[0] = nt
    axes = data.getAxisList()
    axes[0] = htaxis
    hdata = MV2.zeros(sh,typecode=data.dtype.char,savespace=1)
    hdata[:] *= MV2.masked
    cp_atts(data,hdata,id=True)
    hdata.setAxisList(axes)
    set_grid(hdata,data.getGrid())
    #hdata.setGrid(data.getGrid())

    # Convert to hunits
    def hvalue(value):
        return cdtime.r2r(cdtime.reltime(value,taxis.units),hunits).value

    # Interpolate to hours
    ith = 0
    t1 = t0 = hvalue(taxis[0])
    it0 = it1 = 0
    while ith < len(htaxis):
        th = htvalue[ith]
        # Find the right input interval
        while t1 < th:
            it1 += 1
            t1 = hvalue(taxis[it1])
            it0 = it1-1
            t0 = hvalue(taxis[max(it0,0)])
        if maxgap is not None and (t1-t0) > maxgap: # Too large interpolation
            dith = int(t1-t0)+1
##          hdata[ith:ith+dith]
            ith += dith # We skip these hours
        elif it1 == it0: # Strict equality
            hdata[ith] = data[it0]
            ith += 1
        else: # Linear interpolation
            v0 = data[it0]
            v1 = data[it1]
            hdata[ith] = v0 + (v1-v0) * (th-t0) / (t1-t0)
            ith += 1
        it0 = it1

    # Change time units
    if time_units is not None:
        htaxis = ch_units(htaxis,time_units)
        hdata.setAxis(0,htaxis)

    return hdata(order=old_order)



def trend(var):
    """Get linear trend"""

    from genutil.statistics import linearregression
    from grid.misc import set_grid

    # Expansion of first axis
    tt = MV2.array(var.getAxis(0).getValue(),typecode=var.dtype.char)
    if var.rank() > 1:
        sh = var.shape
        tt = MV2.reshape(MV2.repeat(tt,N.multiply.reduce(sh[1:])),sh)

    # Coeffs
    coefs = linearregression(var)

    # Trend
    sh = var.shape
    var_trend = MV2.masked_array(MV2.resize(coefs[0],sh) * tt + MV2.resize(coefs[1],sh))
    cp_atts(var,var_trend)
    var_trend.id = var.id+'_trend'
    var_trend.name = var_trend.id
    if var_trend.attributes.has_key('long_name'):
        var_trend.long_name = 'Linear trend of '+var_trend.long_name
    var_trend.setAxisList(var.getAxisList())
    set_grid(var_trend,var.getGrid())
    #var_trend.setGrid(var.getGrid())

    return var_trend

def detrend(var):
    """Linear detrend"""

    from grid.misc import set_grid

    var_detrend = var - trend(var)
    cp_atts(var,var_detrend)
    var_detrend.id = 'detrended_'+var.id
    var_detrend.name = var_detrend.id
    if var_detrend.attributes.has_key('long_name'):
        var_detrend.long_name = 'Detrended '+var_detrend.long_name
    var_detrend.setAxisList(var.getAxisList())
    set_grid(var_detrend,var.getGrid())
    #var_detrend.setGrid(var.getGrid())

    return var_detrend


def strftime(fmt,mytime=None):
    """Convert current time, datetime, cdtime or string time to strftime

    :Params:

        - **fmt**: Time format with date patterns, like ``"%Y-%m-%d"``.
        - **mytime**, optional: If None, takes current time using
         (:meth:`~datetime.datetime.now()`)

    :Examples:

        >>> print strftime('%Y-%m-%d')
        2014-02-25
        >>> ctime = strftime('%Hh%M', '2020')
        00h00

    :Sea also:

        :meth:`datetime.datetime.strftime` and
        `this link <http://docs.python.org/dev/library/datetime.html#strftime-strptime-behavior>`_.
    """

    if mytime is None:
        mytime = DT.datetime.now()
    else:
        mytime = datetime(mytime)


    return mytime.strftime(str(fmt))

def strptime(mytime,fmt):
    """Parse a string according to a format to retreive a component time

    :Params:

        - **fmt**: Time format with date patterns, like ``"%Y-%m-%d"``.
        - **mytime**: Date string.

    :Example:

        >>> print strptime('25 Jan 2000, '%d %b %Y').month
        1

    :Sea also:

        :meth:`datetime.datetime.strptime` and
        `this link <http://docs.python.org/dev/library/datetime.html#strftime-strptime-behavior>`_.
    """
    return comptime(time.strptime(mytime,fmt))


_re_has_time_pattern = recompile('%[aAbBcdfHIjmpSUwWxXyYzZ]').search
def has_time_pattern(ss):
    """Does ss string contains date pattern like %Y?"""
    if not isinstance(ss, basestring): return False
    if not _re_has_time_pattern(ss): return False
    try:
        DT.datetime.strftime(DT.datetime.now(), ss)
        return True
    except:
        return False


def tz_to_tz(mytime, old_tz, new_tz, copy=True):
    """Convert time from one time zone to another one"""
    try:
        from pytz import timezone
    except ImportError:
        raise VACUMMError('You must install pytz package to add time zone support'
            ' to vacumm')
    if isinstance(old_tz,basestring): old_tz = timezone(old_tz)
    if isinstance(new_tz,basestring): new_tz = timezone(new_tz)

    # Variable
    if cdms.isVariable(mytime):
        check_axes(mytime)
        axis = mytime.getTime()
        assert axis is not None, 'Your variable must have a valid time axis'
        return tz_to_tz(axis, old_tz, new_tz, copy=copy)

    # Time axis case
    if istime(mytime):
        if copy: mytime = mytime.clone()
        ctimes = mytime.asComponentTime()
        for it,ct in enumerate(ctimes):
            new_ct = comptime(old_tz.localize(datetime(ct)).astimezone(new_tz))
            mytime[it] = new_ct.torel(mytime.units).value
        return mytime

    # Loop
    if not isinstance(mytime, list): copy = False
    LH = _LH_(mytime)
    intimes = LH.get()
    outtimes = list(intimes) if copy else intimes
    for i, mytime in enumerate(intimes):

        # Guess type
        typeback = time_type(mytime,  out='func')
        if typeback is None:
            raise TypeError, 'Time of wrong type: %s'%mytime
        if typeback == reltime:
            myUnits = mytime.units
        dt = datetime(mytime)

        # Change zone
        new_dt = old_tz.localize(dt).astimezone(new_tz)

        # Back to correct type
        if typeback == reltime:
            mytime = typeback(new_dt, myUnits)
        else:
            mytime = typeback(new_dt)
        outtimes[i] = mytime

    return LH.put(outtimes)


def from_utc(mytime,new_tz):
    return tz_to_tz(mytime,'UTC',new_tz)
def to_utc(mytime,old_tz):
    return tz_to_tz(mytime,old_tz,'UTC')
def paris_to_utc(mytime):
    return tz_to_tz(mytime,'Europe/Paris','UTC')
def utc_to_paris(mytime):
    return tz_to_tz(mytime,'UTC','Europe/Paris')

def tzname_to_tz(mytzname):
    """Get first time zone found from time zone short name

    :Example:

        >>> tz = tzname_to_tz('CEST')
    """
    try:
        from pytz import all_timezones
    except ImportError:
        raise VACUMMError('You must install pytz package to add time zone support'
            ' to vacumm')
    for name in all_timezones:
        tz = timezone(name)
        if not hasattr(tz, '_tzinfos'): continue
        for (utcoffset, daylight, tzname), _ in tz._tzinfos.iteritems():
            if tzname == mytzname:
                return tz



class DateSorter(object):
    """Sort a list of date string, optionally using a date pattern

    Example:

        >>> from vacumm.misc.atime import DateSorter
        >>> dates = ['annee 2006', 'annee 2002']
        >>> ds = DateSorter('annee %Y')
        >>> dates.sort(ds)
        >>> print dates
        ['2002', '2006']

    """
    def __init__(self,pattern=None,basename=True):
        self.pattern = pattern
        self.basename = basename
    def __call__(self,arg1,arg2):
        if isinstance(self.pattern, basestring):
            if self.basename:
                arg1 = os.path.basename(arg1)
                arg2 = os.path.basename(arg2)
            date1 = strptime(arg1,self.pattern)
            date2 = strptime(arg2,self.pattern)
        else:
            date1 = comptime(arg1)
            date2 = comptime(arg2)
        tu = 'seconds since 2000'
        date1 = date1.torel(tu).value
        date2 = date2.torel(tu).value
        if date1 == date2:return 0
        if date1 < date2:return -1
        return 1

class SpecialDateFormatter(DateFormatter):
    """Special formatter for dates
    Example: ['00h 01/10/2000' '02h' .... '23h'  '00h 01/10/2000'] to mark days and keep hours
    Here the 'phase' is 0 (00h) and the level is 3 (='day')
    """

    def __init__(self,level,fmt=None,special_fmt=None,join=None,phase=None, **kwargs):


        # Which level?
        slevels = ['year', 'month', 'day', 'hour', 'minute']
        if not isNumberType(level):
            if level.lower() not in slevels:
                level = -1
            else:
                level = slevels.index(level.lower().replace('s', ''))
        self.level = level
        # Autoformat
        if level == 0: # Year
            if fmt is None: fmt = '%b'
            if special_fmt is None: special_fmt = '%Y'
        elif level == 1: # Month
            if fmt is None: fmt = '%e'
            if special_fmt is None: special_fmt = '%B'
        elif level == 2: # Day
            if fmt is None: fmt = '%Hh'
            if special_fmt is None: special_fmt = '%d/%m/%Y'
        elif level == 3: # Hour
            if fmt is None: fmt = "%M'"
            if special_fmt is None: special_fmt = '%Hh'
            if join is None:
                join = False
        elif level == 4: # Minute
            if fmt is None: fmt = "%S''"
            if special_fmt is None: special_fmt = "%M'"
            if join is None:
                join = False
        else:
            if level == -2:
                if fmt is None: fmt = '%Y'
            elif level == -3:
                if fmt is None: self.fmt = "%M'%Ss''"
            else:
                if fmt is None: self.fmt = '%Y-%d-%m %H:%M'
                self.level = level = -1
            self.fmt = fmt

        DateFormatter.__init__(self,fmt,**kwargs)
        if level<0: return
        self._phase = [1, 1, 0, 0, 0]
        if phase is not None: self._phase[level] = phase

        # Joined format
        if join in [True, None]:
            join = '%n'
        if join is not False:
            self.special_fmt = join.join([fmt,special_fmt])
        else:
            self.special_fmt = special_fmt

    def __call__(self,x, pos=0):
        dt = num2date(x, self.tz)
        # Main label
        if self.level>=0 and dt.timetuple()[self.level+1:6] == tuple(self._phase[self.level:]):#(0,)*(4-self.level):
            return self.strftime(dt, self.special_fmt)
        # Intermediate label
        return self.strftime(dt, self.fmt)



def interp(vari,outtimes,squeeze=1, **kwargs):
    """Linear interpolation in time

    :Params:

        - **vari****: A cdms variable with a time axis
        - **outtimes**: A time axis, or a list (or single element) of date strings, comptime, reltime, datetime times.
        - *squeeze*: Remove uneeded output dimensions [default: 1]
        - all other keywords are passed to :func:`~vacumm.misc.grid.regridding.interp1d`.
    """
    # Convert to time axis
    if not istime(outtimes):
        outtimes = create_time(outtimes)

    # Call to interp1d
    varo = G.regridding.interp1d(vari, outtimes, **kwargs)

    # Squeeze?
    if squeeze:
        varo = varo(squeeze=1)
    return varo



def is_time(mytime, alsonum=False):
    """Check if mytime is a CDAT time, a :class:`datetime.datetime` or date valid string.

    Equivalent to::

        timetype(mytime) is not None

    :Sea also:

        :func:`is_datetime()` :func:`is_reltime()` :func:`is_cdtime()`
        :func:`is_numtime()` :func:`is_time()`
    """
    ttype = time_type(mytime)
    return ttype and (alsonum or ttype!='numtime')

def is_interval(interval):
    """Check if interval is a valid time interval

    :Params:

        - **interval**: It should be in the following generic form
        ``(<time0>,<time1>[,<bounds>])``, where ``<time?>`` is a valid time
        (see :func:`time_type` and :func:`is_valid`) and ``<bounds>`` are
        interval bounds such as ``"co"``.
    """
    if not isinstance(interval, (list, tuple)) or len(interval)<2 or len(interval)>3:
        return False
    if not is_time(interval[0]) or not is_time(interval[1]):
        return False
    if len(interval)==3 and (not isinstance(interval[2], basestring)
        or not len(interval[2])>1):
        return False
    return True


def round_date(mydate, round_type, mode='round'):
    """Round a date to a step in year, month, day, hour, minute or second

    :Params:

        - **mydate**: A date compatible with :func:`comptime`.
        - **round_type**: A string like "year" or a tuple like (3, "year").
        - **mode**, optional: Rounding mode

            - ``"ceil"``: Choose the upper time.
            - ``"floor"``: Choose the lower time.
            - Else choose the nearest.

    """
    if isinstance(round_type, tuple):
        step, round_type = round_type
    else:
        step = 1
    step = max(int(step), 1)
    if not isNumberType(round_type):
        round_type = STR_UNIT_TYPES.index(unit_type(round_type, string_type=True))
    stype = STR_UNIT_TYPES[round_type]
    ct = comptime(mydate)
    sdate = str(ct)
    values = [int(float(ss)) for ss in RE_SPLIT_DATE.split(sdate) if ss != ''][:round_type+1]
    base = [0, 1, 1, 0, 0, 0][round_type] # [year, month, ...]
    rectif = (values[-1]-base) % step
    ct0 = add_time(cdtime.comptime(*values), -rectif, stype)
    ct1 = ct0 if ct==ct0 else add_time(ct0, step, stype)
    units = 'days since '+sdate
    t = ct.torel(units).value
    if mode not in ['floor', 'ceil']:
        mode = 'round'
    if mode=='round':
        mode = 'ceil' if (t-ct0.torel(units).value > ct1.torel(units).value-t) else 'floor'
    if mode=='ceil':
        return ct1
    return ct0

def round_interval(time, round_type, mode='inner'):
    """Round an time interval using :func:`round_date`

    :Params:

        - **time**: A time interval like ``(time0, time1, 'cce')``
        - **round_type**: A string like "year" or a tuple like (3, "year").
        - **mode**, optional: Rounding mode

            - ``"inner"``: Upper for the lower bound and lower for the upper bound.
            - ``"outer"``: Opposite to "inner"
            - tuple: first is applied to lower bound, second to upper bound.
            - Else, passed to :func:`round_date` and applied to both bounds.

    """
    # Mode
    if mode=='inner':
        modes = 'upper', 'lower'
    elif mode=='outer':
        modes = 'lower', 'upper'
    elif isinstance(mode, tuple):
        modes = mode
    else:
        mode = mode, mode

    return (round_date(time[0], round_type, modes[0]),
        round_date(time[1], round_type, modes[1])) + time[2:]

def midnight_date(mydate):
    """Round a date to the closest midnight

    :Example:

        >>> print midnight_date('2010-11-29 23:10')
        2010-11-30 0:0:0.0


    :Return: :class:`cdtime.comptime` date at midnight
    """
    ctime = comptime(mydate)
    hour = ctime.hour
    ctime = cdtime.comptime(ctime.year, ctime.month, ctime.day)
    if hour>=12:
        ctime = ctime.add(1, cdtime.Day)
    return ctime

def midnight_interval(date):
    """Round dates of a closed interval to the closest midnight

    :Example:

        >>> print midnight_interval('2010-11-29 23:10')
        (2010-11-30 0:0:0.0, 2010-11-30 0:0:0.0, 'ccb')
        >>> print midnight_interval(('2010-11-29 23:10','2010-11-30 04'))
        (2010-11-30 0:0:0.0, 2010-11-30 0:0:0.0, 'ccb')

    :Return: :class:`cdtime.comptime` interval
    """
    if not isinstance(date, tuple):
        return (midnight_date(date), midnight_date(date), 'ccb')
    mydate = [midnight_date(date[0]), midnight_date(date[1])]
    if len(date)>2:
        if date[2][0] == 'o':
            mydate[0] = mydate[0].add(1, cdtime.Day)
        if date[2][1] == 'o':
            mydate[1] = mydate[1].sub(1, cdtime.Day)
    mydate.append('ccb')
    return tuple(mydate)


def daily_bounds(taxis, hstart=0):
    """Create a daily time bounds array from a time axis

    :Params:

        - **taxis**: A time axis or a variable with a time axis.
        - **hstart**, optional: First hour of each daily bounds.
    """
    if cdms2.isVariable(taxis):
        taxis = taxis.getTime()
        if taxis is None:
            raise ValueError('Input variable has no valid time axis')
    elif not istime(taxis):
        raise ValueError('Input axis in not a valid time axis')
    bb = N.zeros((len(taxis),2))
    tu = taxis.units
    for it, ct in enumerate(taxis.asComponentTime()):
        bref = cdtime.comptime(ct.year, ct.month, ct.day, hstart)
        if ct.hour >= hstart:
            bb[it,0] = bref.torel(tu).value
            bb[it,1] = bref.add(1, cdtime.Day).torel(tu).value
        else:
            bb[it,0] = bref.sub(1, cdtime.Day).torel(tu).value
            bb[it,1] = bref.torel(tu).value
    return bb

def hourly_bounds(taxis, mstart=0):
    """Create a hourly time bounds array from a time axis

    :Params:

        - **taxis**: A time axis or a variable with a time axis.
        - **mstart**, optional: First minute of each daily bounds.
    """
    if cdms2.isVariable(taxis):
        taxis = taxis.getTime()
        if taxis is None:
            raise ValueError('Input variable has no valid time axis')
    elif not istime(taxis):
        raise ValueError('Input axis in not a valid time axis')
    bb = N.zeros((len(taxis),2))
    tu = taxis.units
    for it, ct in enumerate(taxis.asComponentTime()):
        bref = cdtime.comptime(ct.year, ct.month, ct.day, ct.hour, mstart)
        if ct.minute >= mstart:
            bb[it,0] = bref.torel(tu).value
            bb[it,1] = bref.add(1, cdtime.Hour).torel(tu).value
        else:
            bb[it,0] = bref.sub(1, cdtime.Hour).torel(tu).value
            bb[it,1] = bref.torel(tu).value
    return bb


def reduce(vari, geterr=False, **kwargs):
    """Average time steps that have the same bounds or time

    :Params:

        - **vari**: Aray with a valid time axis.

    :Example:

        >>> taxis = sst.getTime()
        >>> tbounds = daily_bounds(taxis, hstart=12)
        >>> taxis.setBounds(tbounds)
        >>> nightly_sst = reduce(sst)

    """

    from grid.misc import set_grid

    # Inits
    taxis = vari.getTime()
    if taxis is None:
        raise ValueError('Input variable has no valid time axis')
    order = vari.getOrder()
    vari = vari.reorder('t...')
    bb = taxis.getBounds()
    tt = taxis.getValue().astype('d')
    nti = len(tt)

    # Compression intervals
    ii = N.arange(nti)
    if bb is None:
        dt = N.diff(tt)
        it_lasts = ii[dt!=0.].tolist()
        tv.append(nti-1)
    else:
        ddt = N.diff(bb, axis=0)
        it_lasts = ii[(ddt[:,0]!=0.)&(ddt[:,1]!=0.)].tolist()
    it_lasts.append(nti-1)

    # Init output
    nto = len(it_lasts)
    varo = MV2.zeros((nto,)+vari.shape[1:], vari.dtype)
    cp_atts(vari, varo)
    if vari.getGrid() is not None: set_grid(varo,vari.getGrid())
    #varo.setGrid(vari.getGrid())
    for i,ax in enumerate(vari.getAxisList()[1:]):
        varo.setAxis(i+1, ax)
    varo[:] = MV2.masked
    dtime = N.arange(nto, dtype='d')
    if bb is not None:
        dbounds = N.zeros((nto, 2), dtype='d')
    if geterr:
        err = varo.clone()
        err.id += '_error'
        if hasattr(err, 'long_name'):
            err.long_name += ' average error'

    # Loop on intervals
    mvari = vari.asma()
    it_first = 0
    for it,it_last in enumerate(it_lasts):
        varo[it] = mvari[it_first:it_last+1].mean(axis=0)
        if geterr:
            err[it] = mvari[it_first:it_last+1].std(axis=0)
        if bb is None:
            dtime[it] = tt[it_last]
        else:
            dbounds[it] = bb[it_last]
            dtime[it] = 0.5*(bb[it_last,0]+bb[it_last,1])
        it_first = it_last+1

    # Finish
    dtime = create_time(dtime, taxis.units)
    if bb is not None:
        dtime.setBounds(dbounds)
    varo.setAxis(0, dtime)
    varo = varo.reorder(order)
    if geterr:
        err.setAxis(0, dtime).reorder(order)
        err = err.reorder(order)
        return varo, err
    return varo


def add_margin(interval, lmargin, rmargin=None):
    """Add a margin to an interval

    :Params:

        - **interval**: A time interval or selector such as ``('2000', '20001', 'co')``
          Specified times may be of any valid type.
        - **lmargin**: Left margin. Examples:

            - ``(3,'months')``: explicit.
            - ``'hours'``: Same as ````(1,'hours')``.
            - ``4``: Relative margin of size ``(interval[1]-interval[0)/4``.
            - ``False``: no margin.

          A negative margin decrease the size of the interval.

        - **rmargin**, optional: Right margin (defaults to the left).

    :Example:

        >>> print add_margin(('2000', '2001', 'co'), (5, 'days'))
        ('1999-12-27 0:0:0.0', '2001-1-6 0:0:0.0', 'co')
    """
    # Interval
    if not is_interval(interval):
        raise VACUMMError('Wrong time interval')
    bb = interval[2] if len(interval)>2 else None

    # Format margin
    if lmargin is None: margin = 0
    margins = [lmargin, rmargin]
    dt = None
    for i, margin in enumerate(margins):
        if margin is not None:
            if margin is False: margin = 0
            if isinstance(margin, basestring): # only units
                newmargin =  1, unit_type(margin, string_type=True)
            elif isinstance(margin, (int, float)): # fraction of interval
                if dt is None:
                    rtimes = reltime(interval[:2], 'seconds since 2000')
                    dt = rtimes[1].value-rtimes[0].value
                if margin==0:
                    newmargin = 0, 'seconds'
                else:
                    newmargin = (dt/margin, 'seconds')
            elif not isinstance(margin, (tuple, list)) and len(margin)<2 and \
                (not isinstance(margin[0], (int, float) or not isinstance(margin[1], basestring))):
                raise VACUMMError('Wrong time margin')
            else:
                newmargin = margin
        margins[i] = list(newmargin)
    margins[0][0] = -margins[0][0]

    # Apply
    return type(interval)([add_time(interval[0], *margins[0]),
        add_time(interval[1], *margins[1])]+[bb]*int(bb is not None))


class Intervals(object):
    """Iterator on intervals

    :Params:

        - **time_range**: Time range (optionally with bounds).
        - **dt**: Size of intervals.
        - **reverse**, optional: Reverse the iterator.
        - **roundto**, optional: Round times  to this time units (like 'hours', etc).
        - **bounds**, optional: Add bounds specs to the intervals (like "co").
        - **l/rmargin**, optional: Add a margin to the left and/or right of each interval.
          See :func:`add_margin`.
        - **innerbounds**, optional: Add bounds specs to inner intervals (like "cc").
        - **roundmode**, optional:
            - "all": Round all dates.
            - "inner": Round only inner dates, not first and last.

    :Example:

        >>> for itv in Intervals(('2000','2001','co'),(2,'month')): print itv
        >>> Intervals(('2000','2001','co'),12).tolist()
        >>> Intervals((cdtime.comptime(2000), '2001', 'month',
        ... reverse=True, lmargin=(3,'hours'))

    """
    def __init__(self, time_range, dt, reverse=False, roundto=None, bounds=True,
        lmargin=0, rmargin=None, innerbounds='co', roundmode='all'):

        # Global range
        if not is_interval(time_range):
            raise VACUMMError('Wrong time interval')
        start_date = comptime(time_range[0])
        end_date = comptime(time_range[1])

        # dt
        if isNumberType(dt):
            units = 'seconds since 2000'
            dt = (
                (end_date.torel(units).value-start_date.torel(units).value)/dt,
                'second')
        elif isinstance(dt, list):
            dt = tuple(dt)
        elif not isinstance(dt, tuple):
            dt = (1, dt)

        # Round this range
        if roundto is True:
            roundto = dt[1]
        assert roundmode in ['all', 'inner']
        self._roundmode = roundmode
        self._roundto = roundto
        if roundmode=='all':
            start_date = self.round(start_date)
            end_date = self.round(end_date)
        self.lmargin = lmargin
        self.rmargin = rmargin

        # Closed or open?
        self._lastbounds = None
        self._innerbounds = innerbounds
        if len(time_range) == 3 and bounds is not False:
            self._lastbounds = time_range[2]
        elif bounds in [True, None]:
            self._lastbounds = 'co'
        elif isinstance(bounds, basestring):
            self._lastbounds = bounds
        else:
            self._lastbounds = self._innerbounds = False

        # Inits the iterator
        self._current_date = [start_date, end_date][reverse]
        self._first_date, self._last_date = [start_date, end_date][::1-2*reverse]
        self._reverse = reverse

        # Interval as (value,units)
        self._dt = ([1, -1][reverse]*dt[0], dt[1])

    def round(self, mydate):
#        if self._roundto and (self._roundmode=='all' or
#                (mydate!=self._first_date and mydate!=self._last_date)):
        if self._roundto:
            return round_date(mydate, self._roundto)
        return mydate

    def __iter__(self):
        return self

    def next(self):
        # Iterator is consumed
        if self._current_date == self._last_date:
            raise StopIteration

        # Compute end of interval
        next_date = self.round(add_time(self._current_date, *self._dt))

        # We just passed the end
        if ((self._reverse and next_date.cmp(self._last_date) < 0) or
                (not self._reverse and next_date.cmp(self._last_date) > 0)):
            next_date = self._last_date

        # Save new state
        current_date = self._current_date
        self._current_date = next_date

        # Return interval
        if self._reverse:
            out = next_date, current_date
            bounds = (self._lastbounds
                      if current_date.cmp(self._first_date) == 0
                      else self._innerbounds)
        else:
            out = current_date, next_date
            bounds = (self._lastbounds
                      if next_date.cmp(self._last_date)
                      else self._innerbounds)
        if bounds is not False:
            out += (bounds, )
        if self.lmargin!=0 and self.rmargin is not None:
            out = add_margin(out, self.lmargin, self.rmargin)
        return out

    def tolist(self):
        return list(self)

class IterDates(object):
    """Iterator on dates

    Example:

        >>> from vacumm.misc.atime import IterDates
        >>> for date in IterDates(('2000','2001'),(1,'month')): print date
        >>> for date in IterDates(('2000','2001'),12,closed=False): print date

    """
    def __init__(self, time_range, dt, reverse=False, roundto=None, closed=True):
        # Global range
        start_date = comptime(time_range[0])
        end_date = comptime(time_range[1])

        # dt
        if isNumberType(dt):
            units = 'seconds since 2000'
            dt = (
                (end_date.torel(units).value-start_date.torel(units).value)*1./dt,
                'second')
        elif not isinstance(dt, tuple):
            dt = (1, dt)

        # Round this range
        if roundto is True:
            roundto = dt[1]
        self._roundto = roundto
        start_date = self.round(start_date)
        end_date = self.round(end_date)

        # Closed or open?
        self._closed = closed

        # Inits the iterator
        self._first_date = [start_date, end_date][reverse]
        self._last_date = [start_date, end_date][1-reverse]
        self._reverse = reverse
        self._current_date = None
        oper = 'l' if self._reverse else 'g'
        oper += 't' if self._closed else 'e'
        self._oper = eval(oper)

        # Interval as (value,units)
        self._dt = ([1, -1][reverse]*dt[0], dt[1])

    def round(self, mydate):
        if self._roundto:
            return round_date(mydate, self._roundto)
        return mydate

    def __iter__(self):
        return self

    def next(self):
        # Next date
        if self._current_date is None:
            next_date = self._first_date
        else:
            next_date = self.round(add_time(self._current_date, *self._dt))

        # Iterator is consumed
        tu = 'seconds since 2000'
        if self._oper(next_date.torel(tu).value, self._last_date.torel(tu).value):
            raise StopIteration

        # Save and return
        self._current_date = next_date
        return next_date

    def tolist(self):
        return [date for date in self]


def time_split(what, how, roundit=None, bb='co'):
    """Generic function to split an interval into subintervals

    :Params:

        - **what**: Time interval, dates or time axis
          that can be converted using :func:`comptime`.
        - **how**: Splitting specifications.

            #. A number: divide the interval in equal subintervals.
            #. A single or a list of dates: build subintervals
               with theses dates.
            #. A :class:`IterDates` instance: generate a list of
               dates and make as in 2.
            #. A :class:`Intervals` instance: directly generate
               a list of intervals.

        - **roundit**, optional: Round interval. Valid only if
          ``how`` is an time step specification such as
          ``(1,'year')`` or ``"year"`` (see :class:`Intervals`).
    """
    # Convert to comptime
    bbw = 'cc'
    if isinstance(what, tuple):
        if len(what)==3:
            bbw = what[2]
        what = list(what[:2])
    what = comptime(what)
    if not isinstance(what, list):
        raise TypeError, "Can't split a single date"
    tmin = min(what)
    tmax = max(what)
    tminmax = (tmin, tmax, bbw)
    if not how:
        return tminmax

    # Single date
    if is_time(how): how = [how]

    # dt
    if (isNumberType(how) or isinstance(how, basestring) or
            isinstance(how, tuple)):
        how = Intervals(tminmax[:2], how, roundto=roundit)

    # Interval interator
    if isinstance(how, Intervals):
        how = [itv for itv in how]

    # Date iterator
    elif isinstance(how, IterDates):
        how = [date for date in IterDates]
    tu = 'hours since 2000'
    if isinstance(how, list):
        if not len(how): return [tminmax]
        if not isinstance(how[0], tuple):
            how = comptime(how)
            if len(how) == 1:  # single date
                if (tmin.torel(tu).value <= how[0].torel(tu).value and
                        tmax.torel(tu).value >= how[0].torel(tu).value):
                    how = [tmin, how, tmax]
            how = [(how[i], how[i+1], 'co') for i in xrange(len(how)-1)]
    else:
        raise TypeError, "Can't recognize the split specification"

    # Checks
    itvs = []
    for itv in how:
        itc = itv_intersect(itv, tminmax)
        if itc is not None and itc[0]!=itc[1]:
            itvs.append(itc)
    itvs[-1] = itvs[-1][:2]+(bb,)
    return itvs


def time_split_nmax(what, nmax, roundit=True):
    """Smartly split an interval with a length into subintervals of max length"""
    if isinstance(what, tuple):
        raise TypeError, "Please provide a time axis, a list of dates, or IterDates or Intervals iterators."
    ctimes = comptime(what)
    if not isinstance(ctimes, list):
        raise TypeError, "Can't split a single date"
    if len(ctimes) < nmax:
        return ctimes[0], ctimes[-1]
    taxis = create_time(ctimes)
    specs = ['year', (3, 'month'), (1, 'month'), (10, 'day'), (1, 'day')]
    for i,how in enumerate(specs):
        itvs = time_split(what, how, roundit=roundit)
        nn = []
        for itv in itvs:
            ijk = taxis.mapIntervalExt(itv)
            if ijk is None:
                continue
            nn.append(ijk[2]*(ijk[1]-ijk[0]))
        n = max(nn)
        if n<nmax: break
    return itvs



def itv_intersect(itv1, itv2, bb=None, aslogical=False):
    """Return the intersection of 2 time intervals


    :Return: The interval or ``False`` if not intersection is found
    """
    # Check bounds
    b1 = 'cc' if len(itv1)==2 else itv1[2]
    b2 = 'cc' if len(itv2)==2 else itv2[2]

    # Comptime
    itv1 = [comptime(d) for d in itv1[:2]]
    itv2 = [comptime(d) for d in itv2[:2]]
    tu = 'hours since 2000'

    # Min
    bbi = b1[0]+b2[0]
    if itv1[0]==itv2[0]:
        itv = itv1[0],
        if not bb: bbo = 'c' if 'c' in bbi else 'o'
    else:
        imin = itv1[0].torel(tu).value < itv2[0].torel(tu).value
        itv = (itv1[0], itv2[0])[imin],
        if not bb: bbo = bbi[imin]

    # Max
    bbi = b1[1]+b2[1]
    if itv1[1]==itv2[1]:
        itv += itv1[1],
        if not bb: bbo = 'c' if 'c' in bbi else 'o'
    else:
        imax = itv1[1].torel(tu).value > itv2[1].torel(tu).value
        itv += (itv1[1], itv2[1])[imax],
        if not bb: bbo += bbi[imax]
    if bb is None:
        bb = bbo

    # Check
    if itv[1].torel(tu).value < itv[0].torel(tu).value:
        return False
    if (itv[0].torel(tu).value == itv[1].torel(tu).value and bb is not False
            and 'o' in bb):
        return False

    if aslogical: return True

    return (itv+(bb,)) if bb is not False else itv


def itv_union(itv1, itv2, bb=None):
    """Return the union of 2 time intervals"""

    # Check bounds
    b1 = 'cc' if len(itv1) == 2 else itv1[2]
    b2 = 'cc' if len(itv2) == 2 else itv2[2]
    tu = 'hours since 2000'

    # Comptime
    itv1 = [comptime(d) for d in itv1[:2]]
    itv2 = [comptime(d) for d in itv2[:2]]


    # Min
    bbi = b1[0]+b2[0]
    if itv1[0]==itv2[0]:
        itv = itv1[0],
        if not bb is None: bbo = 'c' if 'c' in bbi else 'o'
    else:
        imin = itv1[0].torel(tu).value > itv2[0].torel(tu).value
        itv = (itv1[0], itv2[0])[imin],
        if not bb: bbo = bbi[imin]

    # Max
    bbi = b1[1]+b2[1]
    if itv1[1].torel(tu).value == itv2[1].torel(tu).value:
        itv += itv1[1],
        if not bb: bbo = 'c' if 'c' in bbi else 'o'
    else:
        imax = itv1[1].torel(tu).value < itv2[1].torel(tu).value
        itv += (itv1[1], itv2[1])[imax],
        if not bb: bbo += bbi[imax]
    if bb is None:
        bb = bbo

    return (itv+(bb,)) if bb is not False else itv


class _LH_(object):
    """To handle argument and result possibily as a list or tuple"""
    def __init__(self, myarg):
        if isinstance(myarg, N.ndarray):
            myarg = myarg.tolist()
        if isinstance(myarg, list):
            self.listtype = list
        elif isinstance(myarg, tuple):
            self.listtype = tuple
        elif hasattr(myarg, 'next') and hasattr(myarg, '__iter__'):
            myarg = [t for t in myarg]
            self.listtype = list
        else:
            self.listtype = None
        self.arg = myarg
    def get(self):
        """Get input argument as a list"""
        if self.listtype is list: return self.arg
        if self.listtype is tuple: return list(self.arg)
        return [self.arg]
    def put(self, result):
        """Get output argument as original type"""
        if self.listtype is None:
            return result[0]
        for func in list, tuple:
            if self.listtype is func:
                if isinstance(result, func):
                    return result
                return func(result)

def _d2sel_(taxis, dsel):
    """Extract time selections from a dict"""
    return [sel for key,sel in dsel.items()
        if key in (['time', taxis.id]+cdms2.convention.time_aliases)
            and sel is not None]

def tsel2slice_old(taxis, *args, **kwargs):
    """Convert time selections on a time axis to a valid slice or None

    :Params:

        - **asind**, optional: Return indices instead of a slice.
        - **nonone**, optional: Return the full slice instead of ``None`` if everything is selected.
        - Positional argument can be coordinates intervals, slices,
          dictionaries or cdms2 selectors.
        - Optional arguments must have a key identified as time or be the axis id,
          and a value as coordinates or slice.

    :Return:

        - A :class:`slice` or ``(i,j,k)`` when possible.
        - ``None`` or the full slice if no slice needed (everything is selected).
        - ``False`` if no intersection is found.

    :Example:

        >>> myslice = tsel2slice(taxis, ('2000', '2002', 'co'), time=slice(4,6))
        >>> myslice = tsel2slice(taxis, cdms2.selectors.Selector(lon=(5,6), time=('2000','2002'))
        >>> myslice = tsel2slice(taxis, slice(10,12), dict(time=slice(4,5), time=('2000','2002'))
    """
    # Inits
    if not istime(taxis): raise VACUMMError('taxis must be a valid time axis')
    asind = kwargs.pop('asind', False)
    nonone = kwargs.pop('nonone', False)
    fullslice = slice(*slice(None).indices(len(taxis)))

    # Convert to list of tuples or slices
    selects = []
    ntup = 0
    for arg in args:
        if arg is None: continue
        if isinstance(arg, cdms2.selectors.Selector):
            psel, asel = split_selector(arg)
            selects.extend(_d2sel_(taxis, asel))
        elif isinstance(arg, dict):
            selects.extend(_d2sel_(taxis,arg))
        else:
            if arg==':':
                args = slice(None)
            elif isinstance(arg,tuple):
                ntup +=1
            elif not isinstance(arg, slice):
                arg = (arg,arg,'ccb')
            selects.append(arg)
    selects.extend(_d2sel_(taxis, kwargs))
    if len(selects)==0:
        if asind: fullslice = fullslice.start,fullslice.stop,fullslice.step
        return fullslice if nonone else None

    # Apply successive selections
    ii = N.arange(len(taxis))
    for isel, sel in enumerate(selects):

        # From coordinates to slice
        if isinstance(sel, tuple):
            ijk = taxis.mapIntervalExt(sel)
            if ijk is None: return False
            sel = slice(*ijk)

        # Work on indices
        ii = ii[sel]
        if ii.size==0: return False

        # Subaxis for future use
        if ntup>1:
            i,j,k = sel.indices(len(ii))
            taxis = taxis.subAxis(i,j,k)

    # Deduce final indices
    i = ii[0]
    j = ii[-1]
    j += N.sign(j-i) if i!=j else 1
    if j<0: j = None
    k = 1 if ii.size==1 else (ii[1]-ii[0])

    # Return
    if not nonone and slice(i,j,k)==fullslice: return
    if asind: return i,j,k
    return slice(i,j,k)

def time_selector(arg0, arg1=None, bounds=None, round=False, utc=True):
    """Time selector formatter that returns start date and end date as component times

    :Example:

        >>> selector('2006','2007') # between two dates
        >>> selector(comptime(1950)) # from a date to now
        >>> selector(1,'month','co') # from now into the past
    """

    # Interval
    if arg1 is None: # from a date to now
        selection = (comptime(arg0), now(utc))
    elif isNumberType(arg0): # from now into the past
        nn = now(utc)
        selection = (add_time(nn, -arg0, arg1), nn)
    else: # between two dates
        selection = (comptime(arg0), comptime(arg1))

    # Bounds
    if isinstance(bounds, basestring):
        selection += (bounds, )

    return selection



def filter_time_selector(*args, **kwargs):
    """Create a pure time selector from all arguments

    All components that are not recognized as a time selection are not kept.

    :Params:

        - **ids**, optional: Special keyword to specify allowed time ids in addition
          to generic ones defined by :attr:`cdms2.convention.time_aliases`.
        - **out**, optional: Inverse the process by removing all time selections
          (see :func:`~vacumm.misc.misc.filter_selector`)?
        - **keeppos**, optional: Remove positional components
          (see :func:`~vacumm.misc.misc.filter_selector`)?
        - **noslice**, optional: Remove slices
          (see :func:`~vacumm.misc.misc.filter_selector`)?
        - Positional argument can be coordinates intervals, slices,
          dictionaries or cdms2 selectors.
        - Optional arguments must have a key identified as time,
          and a value as coordinates or slice.
    """
    # Valid ids for filtering
    ids = kwargs.pop('ids', None)
    out = kwargs.pop('out', False)
    keeppos = kwargs.pop('keeppos', False)
    if ids is None: ids = []
    if isinstance(ids, basestring): ids = [ids]
    ids = list(ids)+['time']+cdms2.convention.time_aliases

    # Build the selector using refinements
    newargs = []
    selector = cdms2.selectors.Selector()
    for arg in args:
        if arg is None: continue
        if isinstance(arg, dict):
            selector.refine(*arg)
        else:
            selector.refine(arg)
    selector.refine(**kwargs)

    # Filter
    filter_selector(selector, ids=ids, copy=False, out=out, keeppos=keeppos)

    return selector

selector = filter_time_selector

def tsel2slice(taxis, *args, **kwargs):
    """Convert time selections on a time axis to a valid slice or None

    :Params:

        - **asind**, optional: Return indices instead of a slice.
        - **nonone**, optional: Return the full slice instead of ``None`` if everything is selected.
        - Positional argument can be coordinates intervals, slices,
          dictionaries or cdms2 selectors.
        - Optional arguments must have a key identified as time or be the axis id,
          and a value as coordinates or slice.

    :Return:

        - A :class:`slice` or ``(i,j,k)`` when possible.
        - ``None`` or the full slice if no slice needed (everything is selected).
        - ``False`` if no intersection is found.

    :Examples:

        >>> myslice = tsel2slice(taxis, ('2000', '2002', 'co'), time=slice(4,6))
        >>> myslice = tsel2slice(taxis, cdms2.selectors.Selector(lon=(5,6), time=('2000','2002'))
        >>> myslice = tsel2slice(taxis, slice(10,12), dict(time=slice(4,5), time=('2000','2002'))
    """
    # Inits
    if not istime(taxis): raise VACUMMError('taxis must be a valid time axis')
    asind = kwargs.pop('asind', False)
    nonone = kwargs.pop('nonone', False)
    fullslice = slice(*slice(None).indices(len(taxis)))

    # Convert to list valid time selector
    kwargs['ids'] = [taxis.id]
    selector = filter_time_selector(*args, **kwargs)

    # No selection
    if len(selector.components())==0:
        if asind: fullslice = fullslice.start,fullslice.stop,fullslice.step
        return fullslice if nonone else None

    # Select
    ii = MV2.arange(len(taxis))
    taxis = taxis.clone()
    taxis.designateTime()
    ii.setAxis(0, taxis)
    try:
        ii = ii(selector).filled()
    except:
        return False
    if ii.size==0:
        return False

    # Deduce final indices
    i = ii[0]
    j = ii[-1]
    j += N.sign(j-i) if i!=j else 1
    if j<0: j = None
    k = 1 if ii.size==1 else (ii[1]-ii[0])

    # Return
    if not nonone and slice(i,j,k)==fullslice: return
    if asind: return i,j,k
    return slice(i,j,k)

def tic():
    """Launch a time counter at the begining of your program.


    :Return:

        - A time to be used with the toc() function.

    :Examples:

        >>> stime = tic()
    """
    import time as tc
    stime = tc.clock()
    print tc.asctime()
    return stime

def toc(stime=0.):
    """Compute the cost of the computation and display in an adapted format.

    :Params:

        - **stime**, optional: The initial time given by the tic() function.

    :Return:

        - Display the time spent in the program.

    :Examples:

        >>> stime = tic()
        >>>
        >>> toc(stime=stime)
    """
    import time as tc
    # print tc.asctime()
    r = tc.clock()-stime
    if r > 60:
        if r > 3600:
            print (r/3600), "hours"
        else:
            print (r/60), " minutes"
    else:
        print tc.clock()-stime, " seconds"


def interp_clim(clim, times, method='linear', day=15):
    """Interpolate a climatology at specified dates"""
    from .grid.regridding import regrid1d, extend1d
    assert method in ('linear', 'cubic')
    taxis = create_time(times)
    assert not (N.diff(taxis[:])<0).any(), 'Output times must be monotonically increasing'
    assert clim.shape[0]==12
    ctimes = taxis.asComponentTime()
    months = N.array([ct.month for ct in ctimes])
    years = N.array([ct.year for ct in ctimes])
    atts = get_atts(clim)
    cmonths = range(1, 13)
    cyears = [0]*12

    left_extent = 0
    right_extent = 0
    if method=='linear':
        left_extent = int((months==1).any())
        right_extent = int((months==12).any())
    else:
        if (months==1).any():
            left_extent = 2
        elif (months==2).any():
            left_extent = 1
        if (months==12).any():
            right_extent = 2
        elif (months==11).any():
            right_extent = 1

    if left_extent or right_extent:
        clim = extend1d(clim, ext=(left_extent, right_extent), axis=0, mode='cylic')
        cmonths = cmonths[12-left_extent:] + cmonths + cmonths[:right_extent]
        cyears = [-1]*left_extent + cyears + [1]*right_extent
    else:
        old_clim_taxis = clim.getAxis(0)

    climo = MV2.zeros((len(taxis), )+clim.shape[1:]) + MV2.masked
    for year in N.unique(years):
        i, j, k = taxis.mapIntervalExt((cdtime.comptime(year, 1, 1),
            cdtime.comptime(year+1, 1, 1), 'co'))
        year_axis = create_time([cdtime.comptime(year+y, m, day)
            for y, m in zip(cyears, cmonths)])
        clim.setAxis(0, year_axis)
        climo[i:j] = regrid1d(clim, taxis.subaxis(i, j), method=method, axis=0)
    if not left_extent and not right_extent:
        clim.setAxis(0, old_clim_taxis)
    climo.setAxisList([taxis]+clim.getAxisList()[1:])
    grid = clim.getGrid()
    if grid is not None:
        climo.setGrid(grid)
    set_atts(climo, atts)
    return climo





#####################################################################
######################################################################

import grid as G
from .axes import istime,check_axes,check_axis,create_time, isaxis
from .misc import (cp_atts,isnumber,kwfilter,split_selector, filter_selector,
    get_atts, set_atts)
from vacumm import VACUMMError

