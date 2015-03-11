# -*- coding: utf8 -*-
"""
Filters for tide signals

Each filter is defined by a half-list (second half) of coefficients and
the time units between each coefficients.

.. seealso::

    :ref:`user.tut.tide.filters`
"""
# Copyright or © or Copr. Actimar/IFREMER (2011-2015)
#
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
import cdtime, numpy as N,MV2,cdms2,genutil.filters as F

from vacumm.misc.axes import check_axis, check_axes, create_time
from vacumm.misc.atime import unit_type,  compress,strftime
from vacumm.misc.filters import running_average
from vacumm.misc.grid import isregular
from vacumm.misc.grid.regridding import interp1d
from vacumm.misc import cp_atts
from scipy.ndimage.filters import convolve1d

import datetime
from operator import isNumberType


_filter_doc = """- *units*: Units of coefficient spacings (like 'hours' or cdms2.Hour) [default: 'hours']
    - *get_tide*: Also return the tide signal [default: False].
    - *only_tide*: Only return the tide signal [default: False].

    :Returns: cdms2 variables with proper time axis.

    - ``fvar``: Filtered version of var (tide signal removed)
    - OR ``fvar,tvar``: Same + tide signal (``get_tide=True``)
    - OR ``tvar``: Only tide signal (``only_tide=True``)
"""

def demerliac(var,**kwargs):
    """Apply a Demerliac filter to signal (cdms2 variable).  First axis is supposed to be time and regular.

    - **var**: Data to filter with first axis  supposed to be time.

    %s
    """

    coefs = N.array([768,766,762,752,738,726,704,678,
      658,624,586,558,512,465,435,392,351,325,
      288,253,231,200,171,153,128,105,91,72 ,
      55, 45, 32, 21, 15,  8,  3,  1],dtype='f')


    return generic(var,coefs,filter_name='Demerliac',**kwargs)
if demerliac.__doc__ is not None:
    demerliac.__doc__ = demerliac.__doc__ % _filter_doc


def godin(var,**kwargs):
    """Apply a Godin filter to signal (cdms2 variable). First axis is supposed to be time and regular.

    - **var**: Data to filter with first axis supposed to be time.

    %s
    """

    coefs = N.array([0.0308,0.0308,0.0306,0.0302,0.0297,0.0291,0.0283,
      0.0274,0.0264,0.0252,0.0239,0.0224,0.0208,0.0192,
      0.0176,0.0160,0.0146,0.0132,0.0119,0.0106,0.0094,
      0.0083,0.0073,0.0063,0.0054,0.0046,0.0038,0.0031,
      0.0025,0.0019,0.0015,0.0010,0.0007,0.0004,0.0002,
      0.0001],'f')


    return generic(var,coefs,filter_name='Godin',**kwargs)
if godin.__doc__ is not None:
    godin.__doc__ = godin.__doc__ % _filter_doc


def generic(var,coefs,units='hours',get_tide=False,only_tide=False,filter_name='Generic',
    check_regular=False,strict=True, interp_method='cubic', axis=None):
    """Apply a filter to signal (cdms2 variable) using symetric coefficients. First axis is supposed to be time and regular.

    - **var**: Data to filter with first axis supposed to be time.
    - **coefs**: Second half of coefficients.
    - *strict*: Mask values where at least one value is missing in the filtering interval.
    - *filter_name*: Name to give to the filter [default: 'Generic']
    %s
    """
    # Return?
    if only_tide:
        get_tide = True

    # Coefficients
    units = unit_type(units)
    sunits = unit_type(units, string_type=True)
    coefs = coefs.astype('f')
    if coefs[-1] != 0.:
        coefs = N.concatenate((coefs,[0.,]))
    ncoefs = (len(coefs)-1)*2+1

    # Input variable
    var = MV2.asarray(var)
    check_axes(var)

    # Input time axis
    if axis is None:
        try:
            axis = var.getOrder().index('t')
        except:
            axis = 0
    time = var.getAxis(axis)
    if check_regular:
        assert isregular(time), 'Your time axis seems not regular'
    dt = N.diff(time[:]).mean()
    if not time.attributes.has_key('units'):
           time.units = 'hours since 1970-01-01'
           print 'The time axis of your variable has no units. Assuming : '+time.units
    if not time.isTime():
        time.designateTime(calendar=cdtime.DefaultCalendar)

    # Check time range
    time_range = time[-1]-time[0]
    coefs_units = sunits+' '+' '.join(time.units.split()[1:])
    coefs_range0 = cdtime.r2r(cdtime.reltime(time[0],coefs_units),time.units)
    coefs_range1 = coefs_range0.add(ncoefs-1,units)
    coefs_range = coefs_range1.value - coefs_range0.value
    coefs_dt = coefs_range0.add(1,units).value - coefs_range0.value
    if time_range < coefs_range:
        raise Exception,'Your sample is shorter than the time range of your coefficients'
    nt = len(var)

    # Remap coefficients to data time axis
    tc = var.dtype.char
#   print 'dt', dt
#   print 'coefs_dt', coefs_dt
    if abs((coefs_dt-dt)/dt) > 0.01:
        oldx = cdms2.createAxis(N.arange(len(coefs))*coefs_dt)
        old_coefs = MV2.array(coefs, axes=[oldx])
        newx = cdms2.createAxis(N.arange(0.,max(oldx),dt))
        coefs = interp1d(old_coefs, newx, interp_method).filled()
    coefs = N.concatenate((coefs[:0:-1],coefs))
    coefs /= N.sum(coefs)
    ncoefs = len(coefs)

    # Filtering using convolution
    fvar = convolve1d(var.filled(0.), coefs, axis=axis, mode='constant')
    if var.mask is MV2.nomask:
        mask = N.zeros(nt)
    else:
        mask = var.mask.astype('f')
    fmask = convolve1d(mask, coefs, axis=axis, mode='constant', cval=1.)!=0.
    fvar = MV2.asarray(fvar)
    fvar = MV2.masked_where(fmask, fvar, copy=0)

    # Output variables
    if not only_tide:
        fvar.id = var.id+'_cotes'
        fvar.name = fvar.id
        if var.attributes.has_key('long_name'):
            fvar.long_name = 'Tide removed signal from '+var.long_name.lower()
        fvar.long_name_fr = 'Signal sans la maree'
        if var.attributes.has_key('units'):
            fvar.units = var.units
        fvar.setAxis(0,time)
        fvar.filter_coefficients = coefs
        fvar.filter_name = filter_name
        res = [fvar,]
    if get_tide:
        tvar = var-fvar
        tvar.id = 'tidal_'+var.id
        tvar.name = tvar.id
        if var.attributes.has_key('long_name'):
            tvar.long_name = 'Tidal signal from '+var.long_name.lower()
        if var.attributes.has_key('units'):
            tvar.units = var.units
        tvar.setAxis(0,time)
        tvar.filter_coefficients = coefs
        tvar.filter_name = filter_name
        if only_tide:
            return tvar
        res.append(tvar)

    if isinstance(res, list) and len(res) == 1:
        return res[0]
    else:
        return res
if generic.__doc__ is not None:
    generic.__doc__ = generic.__doc__ % _filter_doc

#########################################################################

def _get_anomaly_(var,ref='mean',mean=None):

    # Basic checks
    assert var.ndim==1, 'Input variable must be 1D'
    assert cdms2.isVariable(var) and var.getTime() is not None, 'Input variable must have a proper time axis'

    # Get reference
    if ref is None: ref = 'mean'
    nt = len(var)
    if isinstance(ref,int):
        var_ref = running_average(var,ref)
    elif ref in ['demerliac','godin']:
        var_ref = eval(ref)(var,get_tide=False)
    elif ref == 'mean':
        if mean is None:
            var_ref = float(var.mean())
        else:
            var_ref = mean
    else:
        var_ref = ref
    if not isinstance(var_ref, (float, int)) and nt != len(var_ref):
        ntref = len(var_ref)
        ct = var_ref.getAxis(0).subAxis(0,ntref,ntref-1).asComponentTime()
        var = var((ct[0],ct[-1],'cc'))
        nt = ntref
    if var.getMissing() is None: var.setMissing(1.e20)

    # Departure
    vara = var - var_ref

    # Deal with mask
    if vara.mask is not MV2.nomask:
        vara = compress(vara)
        if not  isinstance(var_ref, (float, int)):
            var_ref = compress(MV2.masked_where(vara.mask, var_ref, copy=0))

    return vara, var_ref


def zeros(var, ref='mean',mean=None, getref=True, **kwargs):
    """Get the zeros of a tidal signal

    :Returns: A :mod:`cdms2` variable of signs (-1,1) with a time axis

    :Usage:

    >>> tidal_zeros = zeros(sea_level,ref='demerliac')
    >>> print tidal_zeros[0:1]
    >>> print tidal_zeros[0:1].getTime().asComponentTime()
    """
    # Get anomaly
    ref = kwargs.pop('reference', ref)
    vara, varref = _get_anomaly_(var, ref=ref,mean=mean)
    taxis = vara.getTime()
    vara = vara.filled()
    longref = hasattr(varref, '__len__')

    # Find indices
    sign = N.sign(vara)
    izeros =  N.arange(len(vara)-1).compress(sign[:-1]!=sign[1:])

    # Interpolate
    units = taxis.units
    times = taxis.getValue()
    zeros = N.zeros((len(izeros), ))
    if getref:
        ret = MV2.zeros(len(zeros), id='zeros')
        if not longref:
            ret[:] = varref
    for i, i0 in enumerate(izeros):
        dv = vara[i0+1]-vara[i0]
        zeros[i] = times[i0]*vara[i0+1]/dv - times[i0+1]*vara[i0]/dv
        if getref and longref:
            dt = times[i0+1]-times[i0]
            ret[i] = var_ref[i0]*vara[i0+1]/dv - var_ref[i0+1]*vara[i0]/dv

    # Format
    if not getref:
        ret = MV2.array(sign[izeros], id='zeros')
        ret.units = '1 up and -1 down'
    else:
        cp_atts(var, ret)
    ret.long_name = 'Zeros'
    zeros = create_time(zeros, units)
    ret.setAxis(0, zeros)
    return ret


def extrema(var,ref='mean',mean=None,getmax=True,getmin=True,getsign=False,getidx=False,spline=True,**kwargs):
    """Find extrema of 1D array using a reference. This is suited for detecting low and high tides.

    - **var**: 1D :mod:`cdms2` array of sea level with a time axis.
    - *ref*: Zero reference to compute the sea level anomaly:

        - a filter name (``'demerliac'``, ``'godin'``) where the residual as used as the reference (useful only on long time series),
        - ``'mean'``: the reference is the averaged signal,
        - an integer, so that the reference is a running average using a window of length ``ref``,
        - a float that is used as the reference.

    - *mean*: If ``ref='mean'``, ``mean`` can be substituted to the real mean value.
    - *spline*: If True, extremes are computed using spline interpolation (precision=minute).
    - *getmin*: Return an array of low tides.
    - *getmax*: Return an array of tides.
    - *getsign*: Return signs of tide.
    - *getidx*: Return array of index of low and high tides

    :Returns: ``[lows,][highs,][signs]`` depending on options.
    """

    # Get anomaly
    ref = kwargs.pop('reference', ref)
    vara, varref = _get_anomaly_(var, ref=ref,mean=mean)
    ctimes = vara.getTime().asComponentTime()
    varn = (vara+varref).filled()
    vara = vara.filled()

    # Extremes between zeros
    maxima = []
    minima = []
    sign = N.sign(vara)
    nt = len(var)
    izeros = N.arange(nt-1).compress(sign[:-1]!=sign[1:])
#   izeros = N.compress((sign[:-1]!=sign[1:])&(sign[:-1]!=0)&(sign[1:]!=0),N.arange(nt-1))
    if izeros[-1] != (nt-1):
        izeros = N.concatenate((izeros,[nt-1]))
    i0 = 0
    for i1 in izeros:
##      if i1 == 0: continue
        subvar = varn[i0:i1+1]
        # Type of extremum
        if sign[i0] == -1:
            func = N.argmin
            extrem = minima
        else:
            func = N.argmax
            extrem = maxima
        # Index of extremum
        i = func(subvar)
        # No extremum at limits
        if i != 0 and i != (i1-i0) :
            # Estimation with spline or not
            val,this_ctime = _extremum_(func,ctimes,i0,i,subvar,spline)
            # Append info
            extrem.append((this_ctime,val,i0+i-1,))
        i0 = i1+1

    # Returned data
    if getidx:
        ctime,var,idxmin = zip(*minima)
        idxmin=N.array(idxmin)+1
        ctime,var,idxmax = zip(*maxima)
        idxmax=N.array(idxmax)+1

    res = []
    if getmin: res.append(_extrema_var_(minima,long_name='Minima',id='minima', **kwargs))
    if getmax: res.append(_extrema_var_(maxima,long_name='Maxima',id='maxima', **kwargs))
    if getsign: res.append(sign)

    if getidx:
        res.append(idxmin)
        res.append(idxmax)
    if isinstance(res,list) and len(res) == 1:
        res = res[0]
    return res

def _extremum_(func,ctime,i0,i,var,spline):
    """Extremum possibly using splines"""
    nt = len(var)
    if spline and nt >= 4: # and i != 0 and i != (nt-1)
        if i == 0:
            ifirst, ilast = 0, 4
        elif i == nt-1:
            ifirst, ilast = nt-4, nt
        else:
            icenter = i - int(var[i-1] > var[i+1])
            ifirst = max(icenter-1, 0)
            ilast = ifirst + 4
            if ilast > nt:
                ilast -= 1
                ifirst -= 1
        mn_units = 'minutes since %s'%ctime[i0+ifirst]
        old_rts = cdms2.createAxis(N.asarray([ct.torel(mn_units).value for ct in ctime[i0+ifirst:i0+ilast]],dtype='d'))
        old_var = MV2.array(var[ifirst:ilast], axes=[old_rts], copyaxes=0)
        mn_rts =  cdms2.createAxis(N.arange(int(old_rts[-1]+1),dtype='d'))
        mn_var = interp1d(old_var, mn_rts, method='cubic')
        del old_var, old_rts
#       mn_var = spline_interpolate(old_rts,var[i-1:i+2],mn_rts)
#       mn_var = splev(mn_rts, splrep(old_rts,var[ifirst:ilast]))
        mn_i = func(mn_var)
        val = mn_var[mn_i]
        del mn_var
        this_ctime = cdtime.reltime(mn_i,mn_units).tocomp()
    else:
        this_ctime = ctime[i0+i]
        val = var[i]
    return val,this_ctime

def _extrema_var_(extrem,units=None,indices=False,**kwargs):
    ctime,var,idx = zip(*extrem)
    if indices:
        mytime = cdms2.createAxis(idx)
        mytime.id = 'time_index'
        mytime.long_name = 'Time index'
    else:
        if units is None:
            units = 'minutes since %s'%strftime('%Y-%m-%d %H:%M:%S',ctime[0])
        mytime = create_time(list(ctime), units)
    var = cdms2.createVariable(var,copy=0)
    var.setMissing(1.e20)
    var.setAxis(0,mytime)
    for att,val in kwargs.items():
        setattr(var,att,val)
    return var



