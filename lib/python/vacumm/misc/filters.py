# -*- coding: utf8 -*-
"""Various 1d and 2D filters"""

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
__all__ = ['generic1d', 'shapiro1d', 'gaussian1d', 'hamming1d','generic2d', 'shapiro2d', 'gaussian2d', 'deriv', 'deriv2d',
    'norm_atan','running_average', 'bartlett1d', 'kaiser1d', 'hanning1d', 'blackman1d']
__all__.sort()

import numpy as N,MV2, cdms2
from genutil.filters import *
MA = N.ma
MV = MV2
cdms = cdms2
import scipy.signal
from scipy.signal import convolve2d

from misc import cp_atts
from phys.units import deg2m
from pylab import meshgrid
from axes import islon,islat
import warnings

try:
    from numpy import ComplexWarning
    warnings.filterwarnings('ignore', 'Casting complex values', ComplexWarning)
except ImportError: pass


def generic1d(data, weights, axis=0, mask='same', copy=True, cyclic=False):
    """Generic 1D filter applied to :mod:`MV2` variables using convolution.

    :Params:

        - **data**: Atleast 1D :mod:`MV2` variable.
        - **weights**: integer, 1D weights.
          They are expected to be symmetric and of odd sizes.
        - **axis**, optional: axis on which to operate
        - **mask**, optional: mode of masking.
          The mask is also filtered, and its value helps
          defining the new mask, depending on this parameter:

            - If a float is provided, data are masked if the mask value
              is greater than this parameter.
            - ``"minimal"``: Equivalent to ``1.``. Data is not masked
              if a valid value was used to compute data.
            - ``"maximal"``: Equivalent to ``0.``. Data is masked
              if a invalid value was used to compute data.
            - ``"same"``: Mask is the same as input mask.

        - **copy**, optional: Copy variable before filtering it.

    :Return:

        - The filtered :mod:`MV2` variable.

    :Example:

        >>> generic1d(data, 3.) # running mean using a 3-points block
        >>> generic1d(data, [1.,2.,1], axis=2) # shapiro filter on third axis

    :See also: :func:`scipy.signal.convolve2d`

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
    nw = weights.shape[-1]
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

    :Params:

        - **data**: A :mod:`MV2` variable.
        - Keywords are passed to :func:`generic1d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.array([1.,2.,1.],data.dtype.char)
    return generic1d(data, weights, **kwargs)

def hamming1d(data, M, **kwargs):
    """Hamming 1D filter

    :Params:

        - **data**: A :mod:`MV2` variable.
        - **M**: Size of the Hamming window.
        - Keywords are passed to :func:`generic1d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.hamming(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def hanning1d(data, M, **kwargs):
    """Hanning 1D filter

    :Params:

        - **data**: A :mod:`MV2` variable.
        - **M**: Size of the Hanning window.
        - Keywords are passed to :func:`generic1d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.hanning(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def bartlett1d(data, M, **kwargs):
    """Bartlett 1D filter

    :Params:

        - **data**: A :mod:`MV2` variable.
        - **M**: Size of the Bartlett window.
        - Keywords are passed to :func:`generic1d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.bartlett(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def blackman1d(data, M, **kwargs):
    """Blackman 1D filter

    :Params:

        - **data**: A :mod:`MV2` variable.
        - **M**: Size of the Blackman window.
        - Keywords are passed to :func:`generic1d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.blackman(M).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def kaiser1d(data, M, beta, **kwargs):
    """Kaiser 1D filter

    :Params:

        - **data**: A :mod:`MV2` variable.
        - **M**: Size of the Kaiser window.
        - **beta**: Shape of the window.
        - Keywords are passed to :func:`generic1d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.kaiser(M, beta).astype(data.dtype.char)
    return generic1d(data, weights, **kwargs)

def gaussian1d(data,nxw,**kwargs):
    """Gaussian 1D filter

    - **data**: Data array
    - **nxw**: Size of gaussian weights array along X
    - Other keywords are passed to :func:`generic1d`
    """
    assert nxw % 2 == 1 , 'nxw must be an odd number'
    tc = data.dtype.char
    xx = N.arange(nxw)-nxw/2.
    return generic1d(data, N.exp(-xx**2/nxw**2),**kwargs)



def generic2d(data, weights, mask='same', copy=True):
    """Generic 2D filter applied to 2D (or more) :mod:`MV2` variables using convolution.

    :Params:

        - **data**: Atleast 2D :mod:`MV2` variable.
        - **weights**: integer, 2D weights.
          They are expected to be symmetric and of odd sizes.
          If an integer is provided, a ``(weights,weights)``
          array of ones is used.
        - **mask**, optional: mode of masking.
          The mask is also filtered, and its value helps
          defining the new mask, depending on this parameter:

            - If a float is provided, data are masked if the mask value
              is greater than this parameter.
            - ``"minimal"``: Equivalent to ``1.``. Data is not masked
              if a valid value was used to compute data.
            - ``"maximal"``: Equivalent to ``0.``. Data is masked
              if a invalid value was used to compute data.
            - ``"same"``: Mask is the same as input mask.

        - **copy**, optional: Copy variable before filtering it.

    :Return:

        - The filtered :mod:`MV2` variable.

    :Example:

        >>> generic2d(data, 3.) # running mean using a 3x3 block
        >>> generic2d(data, N.ones(3,3)) # the same
        >>> generic2d(data, N.ones(3,3), weights, mode='minimal') # crop coasts

    :See also: :func:`scipy.signal.convolve2d`

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
    for i in xrange(datan.shape[0]): # TODO: multiproc filter2d
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

def generic2d_old(data,weights,fast=False,fill_value=None,min_valid=0):
    """Generic 2D filter

    - **data**: 2D variable.
    - **weights**: Weights of the filter as 2D array of odd sizes.
    """
    if not cdms.isVariable(data):
        data = MV.masked_array(data, data.mask)
    assert data.rank() == 2, 'You need a 2D data set'
    assert len(weights.shape) == 2, 'You need to apply 2D weights'
    for i in 0,1:
        assert weights.shape[i] % 2 == 1, \
            'Shape of weights must be of odd size in the two directions'
    weights = weights.astype(data.dtype.char)
    nyw,nxw = weights.shape
    dxw = nxw/2
    dyw = nyw/2
    ny,nx = data.shape

    data_filt = data.clone()

    if data.mask is MV.nomask:
        fast = True
    data_min = data.min()
    data_max = data.max()
    if fast:
        if fill_value is None:
            if data.getMissing() is None:
                data.setMissing(data.max()*1.e20)
                fill_value = data.getMissing()
            data_to_use = data.filled()
        else:
            data_to_use = data.filled(fill_value)
        mm = N
    else:
        data_to_use = data
        mm = MV

    for j in xrange(ny):
        yyrange = N.clip([j-dyw,j+dyw],0,ny-1)
        jminw = min(yyrange)-j+dyw+1
        jmaxw = max(yyrange)-j+dyw+1
        for i in xrange(nx):
            xxrange = N.clip([i-dxw,i+dxw],0,nx-1)
            this_data = data_to_use[yyrange[0]:yyrange[1],xxrange[0]:xxrange[1]]
            if min_valid and MV.count(this_data) < min_valid:
                if not fast:
                    data_filt[j,i] = MV.masked
                else:
                    data_filt[j,i] = fill_value
                continue
            iminw = min(xxrange)-i+dxw+1
            imaxw = max(xxrange)-i+dxw+1
            data_filt[j,i] = mm.average(this_data,
                weights = weights[jminw:jmaxw,iminw:imaxw])

    if fast and data.mask is not MV.masked:
        data_filt[:] = MV.masked_greater(data_filt,data_max)
        data_filt[:] = MV.masked_less(data_filt,data_min)
    if data.mask is not MV.nomask:
        data_filt[:] = MV.masked_where(data.mask, data_filt, copy=0)
    return data_filt

def shapiro2d(data, corners=.5, **kwargs):
    """Shapiro (121) 2D filter

    :Params:

        - **data**: A :mod:`MV2` variable.
        - **corner**, optional: Value in (4) corners.
        - Keywords are passed to :func:`generic2d`.

    :Return:

        - A :mod:`MV2` variable
    """
    weights = N.empty((3,3),data.dtype.char)
    weights[:] = corners
    weights[1,:] = [1.,2.,1.]
    weights[:,1] = [1.,2.,1.]
    return generic2d(data, weights, **kwargs)


def shapiro2d_old(data,**kwargs):
    """Shapiro (121) 2D filter

    - **data**: Data array
    - Keywords are passed to :func:`generic2d`
    """
##  print 'shap2d'
    weights = N.zeros((3,3),data.dtype.char)
    weights[1,:] = [1.,2.,1.]
    weights[:,1] = [1.,2.,1.]
    for j in 1,-1:
        for i in 1,-1:
            weights[j,i] = 0.5
    shap = generic2d(data,weights,**kwargs)
    for i in 0, -1:
        shap[i] = data[i]
        shap[:, i] = data[:, i]
    return shap

def gaussian2d(data, nxw, nyw=None, sxw=1/3., syw=1/3., rmax=3., **kwargs):
    """Gaussian 2D filter

    - **data**: Data array
    - **nxw**: Size of gaussian weights array along X (and Y if nyw not given)
    - *nyw*: Size of gaussian weights array along Y [default: nxw]
    - *sxw*: Standard deviation of the gaussian distribution along X.
      If <1, its size is relative to nxw. If > 1, it is directly expressed in grid
      steps.
    - *syw*: Same as sxw for Y direction.
    - *rmax*: Distance relative to sqrt(sxw**2+syw**2) after with weights are
      nullified.
    - Other keywords are passed to :func:`generic2d`
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

    - **data**: Data array (converted to MV array if needed)

    - *axis*: Axis on which the derivative is performed [default: 0]
    - *fast*: Filled masked array before derivating, so use Numeric which is faster than MA or MV [*WARNING* default: True]
    - *physical*: Try physical derivative, taking axis units into account [default: True]
    - *lat*: Latitude for geographical deriviative to convert positions in degrees to meters

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

def deriv2d(data,direction=None,**kwargs):
    """Derivative in a 2D space

    - **data**: 2D variable

    - *direction*: If not None, derivative is computed in this direction, else the module is returned [default: None]
    - Other keywords are passed to deriv()
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

    - *stretch*: If stretch close to 1, saturates values [default: 1]

    Return: Value in [-1,1]
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
    std = var.std()
    var_norm[:] = mm.arctan(stretch*var/var.std())*2./N.pi
    return var_norm



def running_average(x, l, d = 0, w = None, keep_mask = True):
    """Perform a running average on an masked array. Average is linearly reduced near bounds, so that the input and output have the same size.

    :Params:

        - **x**: Masked array
        - **l**: Window size

        - *d*: Dimension over which the average is performed (0)
        - *w*: Weights (1...)
        - *keep_mask*: Apply mask from x to output array

    :Example:

        >>> running_average(x, l, d = 0, w = None, keep_mask = 1)

        Returns Average array (same size as x)

    .. warning:: This function deprecated. Please use :func:`generic1d` instead.
    """

    s= x.shape
    if cdms.isVariable(x):
        AV = MV
        out = x.clone()
    else:
        AV = N
        out = N.zeros(s,x.dtype.char)
        keep_mask = False
    ns = len(s)

    if d > (ns-1) or d < 0:
        raise '[running_average] bad dimension', d

    n = s[d]
    #if w is None:
    #   w = N.ones(n,'f')

    ii = []
    io = []
    for id in range(ns):
        if id == d:
            ii.append('i')
            io.append('imin:imax')
        else:
            ii.append(':')
            io.append(':')

    for i in range(n):
        imin = max(0,  i-l/2)
        imax = min(n-1,i+l/2)+1
        exec 'out['+','.join(ii)+'] = AV.average(x['+','.join(io)+'], d)'

    if keep_mask:
        out = MV.masked_array(out,mask=MV.getmask(x))

    return out

#print 'importing filters 1'


#def smooth1d(var):
#       # Smooth
#   varan = vara.filled()
#   try:
#       len(smooth)
#   except:
#       smooth = [1., 2., 1.]
#   vara[:] = N.convolve(varan, smooth,  mode='same')
#   vara[:] /= N.convolve(N.ones(len(vara)), smooth, mode='same')
#   if vara.mask is not MV2.nomask:
#       mask = N.convolve(vara.mask.astype('f'), smooth, mode='same')
#       vara[:] = MV2.masked_where(mask!=0., vara, copy=0)
