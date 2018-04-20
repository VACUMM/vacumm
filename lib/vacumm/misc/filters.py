# -*- coding: utf8 -*-
"""Various 1D and 2D filters"""

# Copyright or © or Copr. Actimar/IFREMER (2010-2016)
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
from __future__ import absolute_import
import warnings
from six.moves import range

import numpy as N,MV2, cdms2
from genutil.filters import *
import scipy.signal
from scipy.signal import convolve2d
from pylab import meshgrid

from .units import deg2m
from .axes import islon, islat

__all__ = ['generic1d', 'shapiro1d', 'gaussian1d', 'hamming1d','generic2d',
    'shapiro2d', 'gaussian2d', 'deriv', 'deriv2d',
    'norm_atan','running_average',
    'bartlett1d', 'kaiser1d', 'hanning1d', 'blackman1d']



try:
    from numpy import ComplexWarning
    warnings.filterwarnings('ignore', 'Casting complex values', ComplexWarning)
except ImportError: pass

MA = N.ma
MV = MV2
cdms = cdms2

def generic1d(data, weights, axis=0, mask='same', copy=True, cyclic=False):
    """Generic 1D filter applied to :mod:`MV2` variables using convolution.

    Parameters
    ----------
    data:
        Atleast 1D :mod:`MV2` variable.
    weights:
        integer, 1D weights.
        They are expected to be symmetric and of odd sizes.
    axis: optional
        axis on which to operate
    mask: optional
        mode of masking.
        The mask is also filtered, and its value helps
        defining the new mask, depending on this parameter:

        - If a float is provided, data are masked if the mask value
          is greater than this parameter.
        - ``"minimal"``: Equivalent to ``1.``. Data is not masked
          if a valid value was used to compute data.
        - ``"maximal"``: Equivalent to ``0.``. Data is masked
          if a invalid value was used to compute data.
        - ``"same"``: Mask is the same as input mask.

    copy: optional
        Copy variable before filtering it.

    Return
    ------
    :class:`MV2.array`
        The filtered variable.

    Example
    -------
    >>> generic1d(data, 3.) # running mean using a 3-points block
    >>> generic1d(data, [1.,2.,1], axis=2) # shapiro filter on third axis

    See also
    --------
    :func:`scipy.signal.convolve2d`

    """

    # Setup
    # - data
    data = MV2.asarray(data)
    assert data.ndim, 'Input data array must be at least 1D'
    if axis<0: axis += data.ndim
    if axis!=data.ndim-1:
        init_order = data.getOrder(ids=1)
        data = data.reorder('...%i'%axis)
    datan = data.filled(0.).copy()
    nx = datan.shape[-1]
    datan.shape = -1, nx
    # - weights
    if isinstance(weights, int):
        weights = N.ones(weights)
    else:
        weights = N.asarray(weights)
    assert weights.ndim, 'Input weights array must be at least 1D'
    assert weights.shape[-1] % 2 == 1, 'Shape of weights must be of odd size'
    ww = (~N.ma.getmaskarray(data)).astype('i') #  = good 2D
    nw2 = weights.shape[-1] / 2
    weights.shape = 1, -1
    ww.shape = datan.shape
    if data.mask is not MV2.nomask:
        one2d = N.ones(datan.shape, 'i')

    # Cyclic case
    if cyclic:
        mode = 'valid'
        fdatan = N.concatenate((datan[:, -nw2:], datan, datan[:, :nw2]), axis=1)
        fww = N.concatenate((ww[:, -nw2:], ww, ww[:, :nw2]), axis=1)
        one1d = N.ones(datan.shape[-1]+nw2*2, 'i')
    else:
        mode = 'same'
        fdatan = datan
        fww = ww
        one1d = N.ones(datan.shape[-1], 'i')

    # Filter
    kwf = dict(mode=mode)
    datan = convolve2d(fdatan, weights, **kwf)
    one2d = N.resize(N.convolve(one1d, weights[0], **kwf), datan.shape)
    del one1d
    if data.mask is MV2.nomask:
        ww = one2d
    else:
        ww = convolve2d(fww, weights, **kwf)
    del fdatan, fww
    bad = ww==0
    ww[bad] = 1
    datan[:] = N.where(bad, datan, datan/ww.astype('d'))
    ww[bad] = 0

    # Set
    if copy:
        datao = data.clone()
    else:
        datao = data
    datan.shape = data.shape
    bad.shape = data.shape
    if data.mask is not MV2.nomask:
        ww.shape = one2d.shape = data.shape
        if mask is 'same':
            bad = data.mask
        else:
            if mask == 'minimal':
                mask = 1.
            elif mask == 'maximal':
                mask = 0.
            else:
                mask = float(mask)
            bad |= (ww/one2d)<(1-mask) # threshold
        datao[:] = N.ma.masked_where(bad, datan, copy=False)
        del ww, one2d, bad
    else:
        datao[:] = datan

    if axis!=data.ndim-1:
        init_order = cdms2.order2index(datao.getAxisList(), init_order)
        return datao.reorder(init_order)
    return datao

def shapiro1d(data, **kwargs):
    """Shapiro (121) 1D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.array([1.,2.,1.],data.dtype.char)
    return generic1d(data, weights, **kwargs)

def hamming1d(data, M, **kwargs):
    """Hamming 1D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    M: int
        Size of the Hamming window.
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.hamming(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def hanning1d(data, M, **kwargs):
    """Hanning 1D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    M: int
        Size of the Hanning window.
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.hanning(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def bartlett1d(data, M, **kwargs):
    """Bartlett 1D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    M: int
        Size of the Bartlett window.
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.bartlett(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def blackman1d(data, M, **kwargs):
    """Blackman 1D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    M: int
        Size of the Blackman window.
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.blackman(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def kaiser1d(data, M, beta, **kwargs):
    """Kaiser 1D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    M: int
        Size of the Kaiser window.
    beta:
        Shape of the window.
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.kaiser(M, beta).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def gaussian1d(data,nxw,**kwargs):
    """Gaussian 1D filter

    Parameters
    ----------
    data:
        Data array
    nxw:
        Size of gaussian weights array along X
    **kwargs
        Keywords are passed to :func:`generic1d`.

    Return
    ------
    :class:`MV2.array`
    """
    assert nxw % 2 == 1 , 'nxw must be an odd number'
    xx = N.arange(nxw)-nxw/2.
    return generic1d(data, N.exp(-xx**2/nxw**2), **kwargs)



def generic2d(data, weights, mask='same', copy=True):
    """Generic 2D filter applied to 2D (or more) :mod:`MV2` variables using convolution.

    Parameters
    ----------
    data: :class:`MV2.array`
    weights:
        integer, 2D weights.
        They are expected to be symmetric and of odd sizes.
        If an integer is provided, a ``(weights,weights)``
        array of ones is used.
    mask: optional
        mode of masking.
        The mask is also filtered, and its value helps
        defining the new mask, depending on this parameter:

        - If a float is provided, data are masked if the mask value
          is greater than this parameter.
        - ``"minimal"``: Equivalent to ``1.``. Data is not masked
          if a valid value was used to compute data.
        - ``"maximal"``: Equivalent to ``0.``. Data is masked
          if a invalid value was used to compute data.
        - ``"same"``: Mask is the same as input mask.

    copy: optional
        Copy variable before filtering it.

    Return
    ------
    :class:`MV2.array`

    Example
    -------
    >>> generic2d(data, 3.) # running mean using a 3x3 block
    >>> generic2d(data, N.ones(3,3)) # the same
    >>> generic2d(data, N.ones(3,3), weights, mode='minimal') # crop coasts

    See also
    --------
    :func:`scipy.signal.convolve2d`

    """

    # Setup
    data = MV2.asarray(data)
    assert data.ndim>1, 'Input data array must be at least 2D'
    if isinstance(weights, int):
        weights = N.ones((weights, weights))
    else:
        weights = N.asarray(weights)
    assert weights.ndim, 'Input weights array must be at least 2D'
    for i in -2,-1:
        assert weights.shape[i] % 2 == 1, \
            'Shape of weights must be of odd size in the two directions'
    ww = -N.ma.getmaskarray(data).astype('f')+1
    datan = data.filled(0.).copy()
    if datan.ndim > 2:
        ny, nx = datan.shape[-2:]
        datan.shape = datan.size/(nx*ny), ny, nx
    elif datan.ndim == 2:
        datan.shape = (1,) + datan.shape[-2:]
    ww.shape = datan.shape
    one2d = N.ones(datan.shape[-2:])

    # Filter
    kwf = dict(mode='same', boundary='fill', fillvalue=0.)
    for i in range(datan.shape[0]): # TODO: multiproc filter2d
        datan[i] = scipy.signal.convolve2d(datan[i], weights, **kwf)
        if data.mask is MV2.nomask:
            ww[i] = scipy.signal.convolve2d(one2d, weights, **kwf)
        else:
            ww[i] = scipy.signal.convolve2d(ww[i], weights, **kwf)
    if data.mask is not  MV2.nomask:
        one3d = scipy.signal.convolve2d(one2d, weights, **kwf)
        one3d = N.resize(one3d, datan.shape)
    bad = ww==0
    ww[bad]=1.
    datan[:] = N.where(bad, datan, datan/ww)
    ww[bad]=0
#    del bad

    # Set
    if copy:
        datao = data.clone()
    else:
        datao = data
    datan.shape = data.shape
    if data.mask is not MV2.nomask:
        ww.shape = bad.shape = one3d.shape = data.shape
#       print 'mask crit', ww/one3d
        if mask is 'same':
            bad = data.mask
        else:
            if mask == 'minimal':
                mask = 1.
            elif mask == 'maximal':
                mask = 0.
            else:
                mask = float(mask)
            bad |= (ww/one3d) < (1-mask)
        datao[:] = N.ma.masked_where(bad, datan, copy=False)
        del ww, one3d
    else:
        datao[:] = datan
        del one2d
    return datao


def shapiro2d(data, corners=.5, **kwargs):
    """Shapiro (121) 2D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    corner: optional, float
        Value in (4) corners.
    **kwargs
        Keywords are passed to :func:`generic2d`.

    Return
    ------
    :class:`MV2.array`
    """
    weights = N.empty((3,3),data.dtype.char)
    weights[:] = corners
    weights[1,:] = [1.,2.,1.]
    weights[:,1] = [1.,2.,1.]
    return generic2d(data, weights, **kwargs)



def gaussian2d(data, nxw, nyw=None, sxw=1/3., syw=1/3., rmax=3., **kwargs):
    """Gaussian 2D filter

    Parameters
    ----------
    data: :class:`MV2.array`
    nxw:
        Size of gaussian weights array along X (and Y if nyw not given)
    nyw:
        Size of gaussian weights array along Y [default: nxw]
    sxw:
        Standard deviation of the gaussian distribution along X.
        If <1, its size is relative to nxw. If > 1, it is directly expressed in grid
        steps.
    syw:
        Same as sxw for Y direction.
    rmax:
        Distance relative to sqrt(sxw**2+syw**2) after with weights are
        nullified.
    **kwargs
        Keywords are passed to :func:`generic2d`.
    """
    if nyw is None: nyw = nxw
    assert nxw % 2 == 1 and nyw % 2 == 1, 'nxw and nyw must be odd numbers'
    assert sxw > 0 and syw > 0,  'sxw and syw must be positive'

    xx,yy = meshgrid(N.arange(nxw)-nxw/2, N.arange(nyw)-nyw/2)

    if sxw < 1:
        sxw *= nxw/2
    if syw < 1:
        syw *= nyw/2

    weights = N.exp(-(xx**2/sxw**2 + yy**2/syw**2))

    weights[weights<N.exp(-rmax**2)] = 0.

    return generic2d(data, weights, **kwargs)

def deriv(data, axis=0, fast=True, fill_value=None, physical=True, lat=None):
    """Derivative along a given axis

    Parameters
    ----------
    data:
        Data array (converted to MV array if needed)

    axis:
        Axis on which the derivative is performed [default: 0]
    fast:
        Filled masked array before derivating, so use Numeric which is faster than MA or MV [*WARNING* default: True]
    physical:
        Try physical derivative, taking axis units into account [default: True]
    lat:
        Latitude for geographical deriviative to convert positions in degrees to meters

    """
##  print 'deriv2d'
    data = MV.masked_array(data)

    # Reordering if needed
    if axis:
        init_order = data.getOrder(ids=1)
        data = data.reorder(str(axis)+'...')
    data_deriv = data.clone()
    data_deriv.id += '_deriv%i'%axis

    # cdms or Numeric variable  to work on?
    if data.mask is MV.nomask:
        fast = True
    if fast:
        if fill_value is None:
            data_to_use = data.filled()
        else:
            data_to_use = data.filled(fill_value)
    else:
        data_to_use = data

    # Derivative
    data_deriv[1:-1] = data_to_use[2:]-data_to_use[:-2]
    data_deriv[0] = data_to_use[1]-data_to_use[0]
    data_deriv[-1] = data_to_use[-2]-data_to_use[-1]
#   if not fast:
#       for i in 0,-1: data_deriv[i] = MV.masked

    # Physical derivative
    if physical:
        pos = N.resize(data.getAxis(0)[:],data.shape[::-1]).transpose()
        if islon(data.getAxis(0)):
            if lat is None:
                for i,lataxis in enumerate(data.getAxisList()):
                    if islat(lataxis):
                        sh = list(data.shape)
                        if i != data.ndim-1:
                            olen = sh[-1]
                            sh[-1] = len(lataxis)
                            sh[i] = olen
                            lat = N.swapaxes(N.resize(lataxis[:],sh),i,-1)
                        else:
                            lat = N.resize(lataxis[:],sh)
                        break
        if islon(data.getAxis(0)) or islat(data.getAxis(0)):
            pos = deg2m(pos,lat)
            units = 'm-1'
        else:
            units = getattr(data.getAxis(0),'units',None)
        data_deriv[1:-1] /= (pos[2:]-pos[:-2])
        data_deriv[0] /= (pos[1]-pos[0])
        data_deriv[-1] /= (pos[-2]-pos[-1])
    else:
        data_deriv[1:-1] = data_deriv[1:-1] * .5
        units = None

    # Mask
    if fast:
        mask = MV.getmaskarray(data)
        dmask = mask.copy()
        dmask[1:-1] = N.logical_or(mask[2:],mask[:-2])
#       for i in 0,-1: dmask[i] = 1
        data_deriv[:] = MV.masked_array(data_deriv,mask=dmask)
    else:
        for i in 0,-1: data_deriv[i] = MV.masked

    # Reordering back
    if axis:
        init_order = cdms2.order2index(data_deriv.getAxisList(), init_order)
        data_deriv = data_deriv.reorder(init_order)

    # Units
    if units is not None:
        data_units = getattr(data_deriv,'units',None)
        if data_units is None:
            data_deriv.units = units
        else:
            data_deriv.units += ' '+units
            if data_deriv.units == 'm m-1':
                del data_deriv.units

    return data_deriv

def deriv2d(data, direction=None, **kwargs):
    """Derivative in a 2D space

    Parameters
    ----------
    data:
        2D variable

    direction:
        If not None, derivative is computed in this direction, else the module is returned [default: None]
    **kwargs
        Keywords are passed to :func:`deriv`.
    """
    data = MV.masked_array(data)
    data_deriv = data.clone()
    data_deriv.id += '_deriv'
    assert data.ndim == 2, 'You need a 2D data set'
    xdata_deriv = deriv(data,1,**kwargs)
    ydata_deriv = deriv(data,0,**kwargs)
##  print '2d der'
    if direction is None:
        data_deriv[:] = MV.sqrt(xdata_deriv**2+ydata_deriv**2)
    else:
        aa = direction*MV.pi/180.
        data_deriv[:] = xdata_deriv*MV.cos(aa) + ydata_deriv*MV.sin(aa)
    if hasattr(xdata_deriv,'units'):
        data_deriv.units = xdata_deriv.units
    elif hasattr(data_deriv,'units'):
        del data_deriv.units
    return data_deriv


def norm_atan(var,stretch=1.):
    """Normalize using arctan (arctan(strecth*var/std(var))

    Parameters
    ----------
    stretch:
        If stretch close to 1, saturates values [default: 1]

    Return
    ------
    float
        Value in [-1,1]
    """
    if cdms2.isVariable(var):
        var_norm = var.clone()
        var_norm.id  += '_norm'
        mm = MV
    elif MA.isMA(var):
        var_norm = var.copy()
        mm = MA
    else:
        var_norm = N.array(var)
        mm = N
    var_norm[:] = mm.arctan(stretch*var/var.std())*2./N.pi
    return var_norm



