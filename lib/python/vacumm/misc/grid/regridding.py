# -*- coding: utf8 -*-
"""Regridding utilities

.. seealso::

    Tutorials: :ref:`user.tut.misc.grid.regridding`
"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2018)
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
import gc, os, subprocess
import re
import warnings
from collections import OrderedDict
import tempfile, shutil
from copy import deepcopy

import numpy as N, cdms2,  MV2,  regrid2
from cdms2.axis import TransientAxis
import cdtime
import genutil
from _geoslib import Point, Polygon

from ...__init__ import VACUMMError, VACUMMWarning, vacumm_warn, vcwarn
from .kriging import krig as _krig_
from .misc import (axis1d_from_bounds, get_xy, isgrid, t2uvgrids, get_grid,
    set_grid, bounds1d, bounds2d, get_axis,
    meshgrid, create_grid, resol, meshcells, curv2rect, merge_axis_slice,
    get_axis_slices, get_axis, transect_specs, create_axes2d,
    get_distances,
    get_tri, get_grid_type)
from .. import axes as A
from ...misc.misc import (cp_atts, intersect, kwfilter, get_atts, set_atts, closeto,
    splitidx, MV2_concatenate)
from ...misc.atime import are_same_units, ch_units, time_split, lindates
from .basemap import get_proj

MV=MV2
cdms=cdms2

# Python functions
__all__ = ['fill1d', 'regular', 'regular_fill1d', 'cellave1d', 'spline_interp1d',
    'refine', 'GridData', 'griddata', 'cargen', 'fill2d', 'regrid2d',
    'regrid1d', 'interp1d', 'nearest1d', 'cubic1d', 'regrid2dold',
    'xy2grid', 'grid2xy', 'fill1d', 'GriddedMerger', 'regrid_method',
    'cellave2d', 'interp2d', 'xy2xy', 'shift1d', 'shift2d',
    'shiftgrid',  'transect', 'CDATRegridder', 'extend1d', 'extend2d',
    'extendgrid', 'regrid2d_method_name', 'fill1d2', 'krig', 'CurvedInterpolator',
    'regrid1dold', 'regrid2d_tool_name', 'regrid2dnew', 'regrid1d_method_name',
    'cellerr1d']
__all__.sort()

# Fortran functions
_interp_funcs = ['interp1d', 'interp1dx', 'interp1dxx',
    'remap1d', 'remap1dx', 'remap1dxx', 'nearest2d', 'bilin', 'dstwgt',
    'mbilin2d', 'mixt2dx', 'cargen', 'extrap1d', 'curv2rel',
    'nearest2dto1d', 'nearest2dto1dc', 'nearest2dto1dc_reduc',
    'bilin2dto1d', 'bilin2dto1dc', 'bilin2dto1dc_reduc',
    'dstwgt2dto1d', 'dstwgt2dto1dc', 'dstwgt2dto1dc_reduc',
    'cellerr1d', 'cellerr1dx', 'cellerr1dxx', 'linear4dto1dxx',
    'nearest4dto1dxx',
    ]

# Load fortran
_interp_funcs = ['%s as _%s_'%(ff, ff) for ff in _interp_funcs]
_interp_funcs = ', '.join(_interp_funcs)
import_interp = "from _interp_ import %s" % (_interp_funcs, )
try:
    exec  import_interp
except Exception, e:
    print e
    print 'Trying to build it...'
    import subprocess
    cmd = ["make"] # Compilation of all vacumm extensions from the root of sources
    out = subprocess.Popen(cmd, cwd=os.path.dirname(__file__), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE).communicate()
    if out[1]!='':
        raise ImportError("Can't build _interp_ for importation:\n%s"%('\n'.join(out)))
    exec import_interp


# Interpolation methods
_griddata_methods = ['nearest', 'linear', 'cubic', 'krig', 'carg']
_cellave_methods = ['conservative', 'remap', 'cellave', 'conserv']
_cdat_methods = _cellave_methods+['bilinear', 'patch']
_regrid2d_methods = ['nearest', 'mixt', 'interp', 'bining'] + _cdat_methods
_interp1d_methods = ['nearest', 'linear', 'cubic', 'hermit']
_regrid1d_methods = _interp1d_methods+_cellave_methods+['cellerr']

def regrid1d_method_name(method, raiseerr=True):
    """Check the method name and return its generic name"""
    if method is None or method.lower()=='auto': return 'auto'
    method = method.lower()
    if method.startswith('cons'): return 'conserv'
    if method.startswith('cella') or method.startswith('remap'): return 'cellave'
    if method.startswith('celle'): return 'cellerr'
    if 'lin' in method or method=='interp': return 'linear'
#    if method.startswith('bin'): return 'bining'
    if method.startswith('near'): return 'nearest'
    if method.startswith('cub'): return 'cubic'
    if method.startswith('her'): return 'hermit'
    if raiseerr:
        raise VACUMMError('Invalid regrid1d method. Please use for example one of these: '+
            ', '.join(_regrid1d_methods))
    else:
        return method


def regrid1dold(vari, axo, method='auto', axis=None, xmap=None, xmapper=None, mask_thres=.5,
    extrap=0):
    """Interpolation along one axis

    :Params:

        - **vari**: Input cdms array.
        - **axo**: Output cdms axis.
        - **method**:

            - ``"nearest"``: Nearest neighbor
            - ``"linear"``: Linear interpolation
            - ``"cubic"``: Cubic interpolation
            - ``"cellave"``: Cell averaging
            - ``"conserv"``: Conservative cel averaging (like ``cellave`` but with integral preserved)

        - **axis**, optional: Axis (int) on which to operate. If not specified, it is guessed from the
          input and output axis types, or set to ``0``.
        - **xmap**, optional: Integer or tuple that specify on which axes input axis is varying.
        - **xmapper**, optional: Array that specify values of input axis along axes specified by ``xmap``.
          It is an array of size ``(...,len(var.getAxis(xmap[-2])), len(var.getAxis(xmap[-1])), len(var.getAxis(axis))]``.
        - **mask_thres**, optional: Time steps when interpolated mask is greater than
          this value are masked.
        - **extrap**, optional: Extrapolate outside input grid when the "nearest" method
          is used:

            - ``0`` or ``False``: No extrapolation.
            - ``-1`` or ``"min"``, or ``"bottom"``, or ``"lower"``, or ``"first"``:
              Extrapolate toward first values of the axis.
            - ``1`` or ``"max"``, or ``"top"``, or ``"upper"``, or ``"last"``:
              Extrapolate toward last values of the axis.
            - ``2`` or ``"both"``: Extrapolate toward both first and last values.

    .. note::

        Cubic method, use "linear" interpolation when less than 4 valid points are available.
        Linear interpolation uses "nearest" interpolation when less than 2 points are
        available.
    """

    # Input specs
    assert cdms.isVariable(vari), 'Works only with cdms variables'
    A.check_axes(vari)
    order = vari.getOrder()
    axes = vari.getAxisList()
    grid = vari.getGrid()
    missing_value = vari.getMissing()
    if missing_value is None:
        vari.setMissing(1.e20)
    missing_value = vari.getMissing()

    # On which axis?
    if not A.isaxis(axo):
        axo = N.asarray(axo)
        if axo.ndim == 2:
            axo = cdms.createAxis(axo)
#   if A.isaxis(axo): A.check_axis(axo)
    if axis is None:# Guess it
        if A.isaxis(axo): # From axis type
            axis = order.find(A.axis_type(axo))
        else:             # From axis length
            try: axis = list(vari.shape).index(axo.shape[0])
            except ValueError: axis = -1
        # Not found
        if axis == -1: axis = 0
    elif axis == -1:      # Last axis
        axis = vari.ndim-1
    else:                  # Failed
        assert axis >= 0 and axis < vari.ndim, 'Wrong axis'
    if axis in [vari.ndim-2, vari.ndim-1]:
        grid = None

    # Method
    if method == 'auto' or method == None:
        method = regrid_method(vari.getAxis(axis), axo)
    method = regrid1d_method_name(method)
    if method in _cellave_methods:
        conserv = method=='conserv'
#       conserv = method == 'conservative'
        method  = -1
    elif isinstance(method, basestring):
        if method=='cellerr':
            method = 4
        elif method in _interp1d_methods:
            method = _interp1d_methods.index(method)
        else:
            method = -2 # Wrong method
    else:
        method = int(method)
    assert method >= -1 and method < 5, 'Wrong method'
    if method == 2: method = 3 # Hermit is always better
    if xmap is None:
        interp_func = _interp1d_
    else:
        interp_func = _interp1dx_
    if method == -1:
        if xmap is None:
            regrid_func = _remap1d_
        else:
            regrid_func = _remap1dx_
        args = [conserv]
    else:
        regrid_func = interp_func
        args = [method]

    # Extrapolation
    if isinstance(extrap, basestring):
        extrap = extrap.lower()
        if extrap in ['both', 'all']:
            extrap = 2
        elif extrap in ['bottom', 'min', 'lower', 'first']:
            extrap = -1
        elif extrap in ['top', 'max', 'upper', 'last']:
            extrap = 1
        else:
            extrap = 0
    elif not isinstance(extrap, int):
        extrap = int(bool(extrap))

    # Convert to numeric
    varin = vari.filled()

    # Translate axes
    if xmap is None: # Simple case
        if axis != vari.ndim-1:
            varis = varin.swapaxes(axis, -1)
            del varin
        else:
            varis = varin
    else: # Extended case

        # Check xmap form
        if isinstance(xmap, int):
            xmap = [xmap]
        else:
            xmap = list(xmap)
        for ix in xrange(len(xmap)):
            if xmap[ix] == -1:
                xmap[ix] = vari.ndim-1
            else:
                assert xmap[ix] >= 0 and xmap[ix] < vari.ndim, 'Wrong xmap axis'

        # Roll axes
        oldmap = range(varin.ndim)
        newmap = xmap+[axis]
        for iax in oldmap[::-1]:
            if iax not in newmap:
                newmap = [iax]+newmap
        if oldmap != newmap:
#           backmap = [oldmap.index(iax) for iax in newmap]
            backmap = [newmap.index(iax) for iax in oldmap]
            varis = varin.transpose(*newmap)
            del varin
        else:
            backmap = None
            varis = varin

        # Size of extension block
        nxb = N.multiply.reduce(varis.shape[:-1])
        xshape = (nxb, varis.shape[-1])
    nyi = varis.shape[-1]

    # Axes content
    if xmap is None: # Simple case

        # Input axis
        axi = axes[axis]

        # Special case for time axis: must have same units !
        if (A.axis_type(axo), order[axis]) == ('t', 't') and \
            hasattr(axi, 'units') and hasattr(axo, 'units') and \
            not are_same_units(axi.units, axo.units):
            axi = ch_units(axi, axo.units, copy=True)

    else: # Extended case

        axi = xmapper.reshape(xshape, order='F')


    # Reshape var to get a 2D array
    if varis.ndim != 2:
        vari2d = varis.reshape(varis.size/nyi, nyi, order='F')
    else:
        vari2d = varis

    # First guess
    varo2d = regrid_func(vari2d, axi[:], axo[:], missing_value, *args)

    # Mask
    maski2d = vari2d==missing_value
    if method and N.any(maski2d):

        # Float mask
        maski2df = maski2d.astype('f')
        masko2d = regrid_func(maski2df, axi[:], axo[:], 1.e20, *args)
        masko2d[masko2d==1.e20] = 1.

        # Masking
        if method == -1:

            # Cell case: threshold
            varo2d[:] = N.where(masko2d>mask_thres, missing_value, varo2d)

        else:

            # Lower method
            varolow = interp_func(vari2d, axi[:], axo[:], missing_value,
                min(abs(method), 2)-1, extrap=extrap)

            # Cubic case: lower order again
            if method >= 2:
                masko2dc = interp_func(maski2df, axi[:], axo[:], 1.e20, 1)
                masko2dc[masko2dc==1.e20] = 1.
                varolowc = interp_func(vari2d, axi[:], axo[:], missing_value, 0)
                varolow[:] = N.where(masko2dc!=0., varolowc, varolow)
                del masko2dc, varolowc

            # Select between nearest and linear, or linear and cubic
            varo2d[:] = N.where(masko2d!=0., varolow, varo2d)
            del varolow

        del maski2df, masko2d

    # Reshape back
    if varis.ndim != 2:
        varos = varo2d.reshape(varis.shape[:-1]+(len(axo), ), order='F')
        del varo2d, vari2d
    else:
        varos = varo2d
    del varis

    # Retranslate axes back
    if xmap is not None and backmap is not None:
        varon = varos.transpose(*backmap)
        del varos
    elif axis != vari.ndim-1:
        varon = varos.swapaxes(-1, axis)
        del varos
    else:
        varon = varos

    # Convert to cdms
    varo = MV.masked_values(varon, missing_value)
    cp_atts(vari, varo, id=True)
    varo.setMissing(missing_value)
    axes[axis] = axo
    varo.setAxisList(axes)
    varo.setGrid(grid)
    gc.collect()
    return varo

def _subshape_(bigshape, subshape, axis=None):
    """Get the index of a unique of occurence of subshape in bigshape or None

    :Params:

        - **bigshape**: Tuple or big array.
        - **subshape**: Tuple or smaller array.
        - **axis**: Axis index of bigshape that is allowed to differ in subshape.
    """
    if hasattr(bigshape, 'shape'): bigshape = bigshape.shape
    if hasattr(subshape, 'shape'): subshape = subshape.shape
    bigshape = list(bigshape)
    subshape = list(subshape)
    nb = len(bigshape)
    ns = len(subshape)
    assert nb>=ns # TODO: make it more generic
    istart = None
    for i in xrange(nb-ns+1):
        subbigshape = bigshape[i:i+ns]
        if axis is not None and axis>=i and axis<i+ns:
            i0 = axis-i
            subbigshape[i0] = subshape[i0]
        if subbigshape == subshape:
            if istart is not None: return # two solutions => no valid solution
            istart = i
    return istart

def _getiax_(vari, ax, axis):
    """Guess the index of the target axis of vari in the ax array"""
    # 1D case
    if ax.ndim==1: return 0

    # A unique sub-tuple
    iax = _subshape_(vari, ax, axis)
    if iax is not None:
        iax = axis - iax
        if iax>=0: return iax
        iax = None

    # A unique axis type
    if cdms2.isVariable(ax):
        l = A.axis_type(vari.getAxis(axis))
        if l!='-' and l in ax.getOrder():
            iax = ax.getOrder().index(l)

    return iax

def _syncshapes_(axi, iaxi, axo, iaxo):
    """Resize two arrays to have the same shapes except when the second one is 1D

    :Params:

        - **axi**: First numpy array.
        - **iaxi**: Index of pivot (target axis) in first array.
        - **axo**: Second numpy array.
        - **iaxo**: Index of pivot (target axis) in second array.
    """
    # Nothing to resize
    if axo.ndim==1:
        return axi, iaxi, axo, iaxo

    # Pivot as a reference
    nil = iaxi # left, first
    nir = axi.ndim-iaxi-1 # right, second
    nol = iaxo # left, first
    nor = axo.ndim-iaxo-1 # right, second

    # Right adjusment
    if nir>nor: # expand the second to the right
        for ir in xrange(nir-nor):
            axo = axo.reshape(axo.shape+(1,))
            axo = N.repeat(axo, axi.shape[iaxi+nor+ir+1], axis=-1)
    elif nor>nir: # expand the first to the right
        for ir in xrange(nor-nir):
            axi = axi.reshape(axi.shape+(1,))
            axi = N.repeat(axi, axo.shape[iaxo+nir+ir+1], axis=-1)

    # Left adjustment
    if nil>nol: # expand the second to the left
        axo = N.resize(axo, axi.shape[:nil-nol]+axo.shape)
        iaxo += nil-nol
    elif nol>nil:
        axi = N.resize(axi, axo.shape[:nol-nil]+axi.shape)
        iaxi += nol-nil

    return axi, iaxi, axo, iaxo

def _toright_(ar, iax):
    """Reform an array so the iax dim becomes the last

    :Params:

        - **ar**: A numpy array.
        - **iax**: Index of the dim.

    :Return:

        ``newar,bakmap`` such as::

            newar.shape[-1] == ar.shape[iax]
            ar == newarr.translate(*bakmap)
    """
    oldmap = range(ar.ndim)
    newmap = list(oldmap)
    newmap.remove(iax)
    newmap.append(iax)
    bakmap = [newmap.index(iax) for iax in oldmap]
    newar = ar.transpose(*newmap)
    return newar, bakmap


def regrid1d(vari, axo, method='auto', axis=None, axi=None, iaxo=None, iaxi=None,
        xmap=None, xmapper=None, mask_thres=.5, extrap=0,
        erri=None, errl=None, geterr=False):
    """Interpolation along one axis

    :Params:

        - **vari**: Input cdms array.
        - **axo**: Output cdms axis or array. It can be of any dimensions.
        - **method**:

            - ``"nearest"|0``: Nearest neighbor
            - ``"linear"|2``: Linear interpolation
            - ``"cubic"|2``: Cubic interpolation (not used, switched to ``3``)
            - ``"hermit"|3``: Cubic hermit interpolation
            - ``"cellerr"|4``: Cell averaging based on errors
            - ``"cellave"|-1``: Cell averaging
            - ``"conserv"|-1``: Conservative cel averaging
              (like ``cellave`` but with integral preserved)

        - **axis**, optional: Dimension (int) on which the interpolation is performed.
          If not specified, it is guessed from the
          input and output axis types, or set to ``0``.
        - **axi**, optional: Input axis. It defaults to the axis-th axis of ``vari``.
          Like ``axo``, it can be of any dimensions.
        - **iaxo**, optional: Dimension of ``axo`` on which the interpolation is performed
          when ``axo`` has more than one dimension.
        - **iaxi**, optional: Same as ``iaxo`` but for ``axi``.
        - **mask_thres**, optional: Time steps when interpolated mask is greater than
          this value are masked.
        - **extrap**, optional: Extrapolate outside input grid when the "nearest" method
          is used:

            - ``0`` or ``False``: No extrapolation.
            - ``-1`` or ``"min"``, or ``"bottom"``, or ``"lower"``, or ``"first"``:
              Extrapolate toward first values of the axis.
            - ``1`` or ``"max"``, or ``"top"``, or ``"upper"``, or ``"last"``:
              Extrapolate toward last values of the axis.
            - ``2`` or ``"both"``: Extrapolate toward both first and last values.

        - **erri**, optional: Input "measurement" errors with the same shape as input
          variable.
        - **errl**, optional: Derivative of lag error with respect to lag.
          Note that the lag must expressed in **days** for time axes.
          If positive, it based on quadratic errors, else on error itself.
          Estimate for instance it using the slope of a linear regression.
          It is usually varying in space and constant in time.
        - **geterr**, optional: When method is "cellerr", also return the error
          along with the variable.

    :Examples:

        >>> varo = regrid1d(vari, taxis, method='linear') # interpolation in time
        >>> varo = regrid1d(vari, zo, axis=1) # Z interpolation on second axis
        >>> varo = regrid1d(vari, zzo, iaxo=1, axi=zzi, iaxi=1) # sigma to sigma

    .. note::

        Cubic method, use "linear" interpolation when less than 4 valid points are available.
        Linear interpolation uses "nearest" interpolation when less than 2 points are
        available.
    """

    # Input specs
    assert cdms.isVariable(vari), 'Works only with cdms variables'
    if xmap is not None or xmapper is not None:
        warnings.warn('xmap and xmapper keywords of regrid1d are deprecated. '
            'Please use axi/axo and iaxi/iaxo keywords instead.', VACUMMWarning)
    A.check_axes(vari)
    order = vari.getOrder()
    axes = vari.getAxisList()
    grid = vari.getGrid()
    missing_value = vari.getMissing()
    if missing_value is None or N.isnan(missing_value):
        vari.setMissing(1.e20)
    missing_value = vari.getMissing()

    # Output axis
    axon = axo[:]
    if N.ma.isMA(axon): axon = axon.filled(missing_value)

    # Working data axis
    if vari.ndim==1: # 1D
            axis = 0

    elif axis is None:# Guess it

        if A.isaxis(axo): # From axis type
            axis = order.find(A.axis_type(axo))

        elif axon.ndim==1:  # From axis length
            axis = _subshape_(vari, axo)
            if axis is None:
                raise VACUMMError('Please, specify the "axis" parameter for interpolation')

        # Not found, so 0 without verification now (later)
        if axis == -1 or axis is None: axis = 0

    else:
        if axis<0:
            axis += vari.ndim
        if axis < 0 or axis >= vari.ndim:
            raise VACUMMError('Wrong "axis" parameter for interpolation')
    if axis in [vari.ndim-2, vari.ndim-1]:
        grid = None

    # Input axis
    # - get it
    if axi is None: axi = vari.getAxis(axis)
    # - special case for time axis: must have same units !
    dxi2o = 1.
    if ((A.axis_type(axo), order[axis]) == ('t', 't') and
            hasattr(axi, 'units') and hasattr(axo, 'units') and
            not are_same_units(axi.units, axo.units)):
        axi = ch_units(axi, axo.units, copy=True)
        if errl is not None:
            axou = axo.units  #.split(1)[0] + axi.units.split(1)[1]
            dxi2o = (cdtime.reltime(1, axou).value -
                     cdtime.reltime(0, axou).value)
            dxi2o /= (cdtime.reltime(1, axi.units).value -
                      cdtime.reltime(0, axi.units).value)
    # - numeric version
    axin = axi[:]
    if N.ma.isMA(axin):
        axin = axin.filled(missing_value)

    # Subaxes
    nxo = axon.ndim
    if nxo>1:
        if iaxo is None:
            iaxo = _getiax_(vari, axo, axis)
        if iaxo is None:
            raise VACUMMError("Please, specifiy the 'iaxo' parameter")
    else:
        iaxo = 0
    if iaxo<0:
        iaxo = nxo+iaxo
    nxi = axin.ndim
    if nxi>1:
        if iaxi is None:
            iaxi = _getiax_(vari, axi, axis)
        if iaxi is None:
            raise VACUMMError("Please, specifiy the 'iaxi' parameter")
    else:
        iaxi = 0
    if iaxi<0:
        iaxi = nxi+iaxi

    # Verifications
    if vari.shape[axis-iaxi:axis-iaxi+nxi]!=axi.shape:
        raise VACUMMError("Input axis has a invalid shape: %s (!=%s)"%(
            axi.shape, vari.shape[axis-iaxi:axis-iaxi+nxi]))
    vs = list(vari.shape[axis-iaxo:axis-iaxo+nxo])
    vs[iaxo] = axo.shape[iaxo]
    vs = tuple(vs)
    if vs!=axo.shape:
        raise VACUMMError("Output axis has a invalid shape: %s (!=%s)"%(axo.shape, vs))

    # Homogeneize axi and axo shapes (except for z dim)
    iaxo_bak = iaxo
    axin, iaxi, axon, iaxo = _syncshapes_(axin, iaxi, axon, iaxo)

    # Convert to numeric
    varin = vari.filled()

    # Push interpolation dimension to the right and convert to 1D or 2D arrays
    varis, bakmapv = _toright_(varin, axis)
    vari2d = varis.reshape(-1, varis.shape[-1])
    axind, bakmapi = _toright_(axin, iaxi)

    if axind.ndim>2:
        axind = axind.reshape((-1, axin.shape[iaxi]))
    axond, bakmapo = _toright_(axon, iaxo)
    if axond.ndim>2:
        axond = axond.reshape((-1, axon.shape[iaxo]))
    nxb = vari2d.size/max(axind.size, axond.size)
    nxi = axind.ndim
    nxo = axond.ndim


    # Method
    if method == 'auto' or method == None:
        method = regrid_method(axin, axon, iaxi=-1, iaxo=-1)
    method = regrid1d_method_name(method)
    if method in _cellave_methods:
        conserv = int(method.startswith('conserv'))
#       conserv = method == 'conservative'
        method  = -1
    elif isinstance(method, basestring):
        if method == 'interp': method = 'linear'
        if method=='cellerr':
            method = 4
        elif method in _interp1d_methods:
            method = _interp1d_methods.index(method)
        else:
            method = -2 # Wrong method
    else:
        method = int(method)
    assert method >= -1 and method < 5, 'Wrong method'
    if method == 2: method = 3 # Hermit is always better

    # Cellerr
    if method==4:

        # Measurement errors
        if erri is None:
            raise VACUMMError('You must provide also erri for the "cellerr" method')
        if vari.shape!=erri.shape:
            raise VACUMMError('vari and erro must have the same shape: %s != %s'%
                (vari.shape, erri.shape))
        erris, bakmapv = _toright_(N.ma.asarray(erri).filled(missing_value), axis)
        errm = erris.reshape(-1, erris.shape[-1])

        # Lag errors
        if errl is None:
            errl = 0.
        errl = N.atleast_1d(errl)
        if N.ma.isMA(errl):
            errl = errl.filled(missing_value)
        errl = errl.ravel() * dxi2o
        if errl.size%nxi:
            raise VACUMMError('The size of errl (%i) in not a multiple nxi (%i).'%
                (errl.size, nxi))


    # Routine and arguments
    if nxi==1 and nxo==1: # 1D->1D
        interp_func = _interp1d_
        remap_func = _remap1d_
        cellerr_func = _cellerr1d_
    elif nxo==1: # 1D->ND
        interp_func = _interp1dx_
        remap_func = _remap1dx_
        cellerr_func = _cellerr1dx_
    else: # ND->ND
        interp_func = _interp1dxx_
        remap_func = _remap1dxx_
        cellerr_func = _cellerr1dxx_
    if method == -1: # Cellave
        regrid_func = remap_func
        kwargs=dict(conserv=conserv,extrap=0)
    elif method==4:
        regrid_func = cellerr_func
        kwargs=dict(errm=errm, errl=errl)
    else: # Interp
        regrid_func = interp_func
        kwargs=dict(method=method,extrap=0)


    # Regrid
    varo2d = regrid_func(vari2d, axind, axond, missing_value, **kwargs)
    if method==4:
        varo2d, erro2d = varo2d


    # Extrapolation
    if isinstance(extrap, basestring):
        extrap = extrap.lower()
        if extrap in ['both', 'all']:
            extrap = 2
        elif extrap in ['bottom', 'min', 'lower', 'first']:
            extrap = -1
        elif extrap in ['top', 'max', 'upper', 'last']:
            extrap = 1
        else:
            extrap = 0
    elif not isinstance(extrap, int):
        extrap = int(bool(extrap))
    if extrap:
        varo2d = _extrap1d_(varo2d, missing_value, extrap)

    # Reshape back
    varos = varo2d.reshape(varis.shape[:-1]+axond.shape[-1:]) ; del varo2d
    varon = varos.transpose(*bakmapv) ; del varos
    if geterr and method==4:
        erros = erro2d.reshape(varis.shape[:-1]+axond.shape[-1:]) ; del erro2d
        erron = erros.transpose(*bakmapv) ; del erros

    # Convert to cdms
    varo = MV2.masked_values(varon, missing_value)
    cp_atts(vari, varo, id=True)
    varo.setMissing(missing_value)
    if geterr and method==4:
        erro = MV2.masked_values(erron, missing_value)
        cp_atts(vari, erro, id=True)
        erro.id = erro.id+'_error'
        if hasattr(erro, 'long_name'):
            erro.long_name = erro.long_name+' error'
        else:
            erro.long_name = 'Error'
        erro.setMissing(missing_value)
    if A.isaxis(axo):
        axes[axis] = axo
    elif cdms2.isVariable(axo):
        axes[axis] = axo.getAxis(iaxo_bak)
    else:
        axes[axis] = varo.getAxis(axis)
        if nxo==1:
            odlaxis = axes[axis]
            axes[axis][:] = axo[:]
            cp_atts(oldaxis, axes[axis])
    varo.setAxisList(axes)
    varo.setGrid(grid)
    if geterr and method==4:
        erro.setAxisList(axes)
        erro.setGrid(grid)
        return varo, erro
    return varo


def nearest1d(vari, axo, **kwargs):
    """Interpolation along an axes

    :Params:

        - **vari**: Input cdms array
        - **axo**: Output cdms axis
        - **axis**, optional: Axis on wich to operate
        - Other keywords are passed to :func:`regrid1d`

    .. note::

        This is an wrapper to :func:`regrid1d` using ``nearest`` as a default method.
        See its help for more information.
    """
    return regrid1d(vari, axo, 'nearest', **kwargs)

def interp1d(vari, axo, method='linear', **kwargs):
    """Linear or cubic interpolation along an axes


    :Params:

        - **vari**: Input cdms array
        - **axo**: Output cdms axis
        - **axis**, optional: Axis on wich to operate
        - Other keywords are passed to :func:`regrid1d`

    .. note::

        This is an wrapper to :func:`regrid1d` using ``linear`` as a default method.
        See its help for more information.
    """
    return regrid1d(vari, axo, method, **kwargs)

def cubic1d(vari, axo, **kwargs):
    """Cubic interpolation along an axes


    :Params:

        - **vari**: Input cdms array
        - **axo**: Output cdms axis
        - **axis**, optional: Axis on wich to operate
        - Other keywords are passed to :func:`regrid1d`

    .. note::

        This is an wrapper to :func:`regrid1d` using ``cubic`` as a default method.
        See its help for more information.
    """
    return regrid1d(vari, axo, 'cubic', **kwargs)

def cellave1d(vari, axo, conserv=False, **kwargs):
    """Cell averaging  or conservative regridding along an axis

    :Params:

        - **vari**: Input cdms array
        - **axo**: Output cdms axis
        - **axis**, optional: Axis on wich to operate
        - **conservative**, optional: If True, regridding is conservative
        - Other keywords are passed to :func:`regrid1d`

    .. note::

        This is an wrapper to :func:`regrid1d` using ``cellave`` or ``conservative`` as a default method.
        See its help for more information.
    """
    if kwargs.pop('conservative', conserv):
        method = 'conservative'
    else:
        method = 'cellave'
    return regrid1d(vari, axo, method, **kwargs)

remap1d = cellave1d

def cellerr1d(vari, axo, erri, errl=None, **kwargs):
    """Cell averaing with weights based on errors


    :Params:

        - **vari**: Input cdms array
        - **axo**: Output cdms axis
        - **erri**: Input measurement errors
        - **errl**, optional: Input lag error relative to lag
        - **axis**, optional: Axis on wich to operate
        - Other keywords are passed to :func:`regrid1d`

    .. note::

        This is an wrapper to :func:`regrid1d` using ``cellerr`` method.
        See its help for more information.
    """
    return regrid1d(vari, axo, method, **kwargs)

def fill1d2(vi,axis=0, k=1,padding=None,clip=False,min_padding=None, method='linear'):
    """Fill missing values of a 1D array using spline interpolation.

    - **vi**: Input cdms variable

    - *k*: Order of splines [default: 1 = linear]
    - *padding*: Padding around an gap defining where on which part of the sample we must fit splines [default: max([min_padding,len(gap)*5])]
    - *min_padding*: See padding [default: k]
    - *method*: See :func:`interp1d` [default: linear]

    :Return:
        Filled :mod:`cdms2` variable
    """

    assert vi.rank() == 1,'Input variable must be of rank 1'
    import scipy.signal as S


    if min_padding is None:
        min_padding = k
    tc = vi.dtype.char

    xi = vi.getAxis(0).getValue().astype('d')

    # Removes missing points
    mask = MV.getmaskarray(vi)
    cxi = N.ma.masked_array(xi,mask=mask).compressed()
    cvi = vi.compressed()
    cii = N.ma.masked_array(N.arange(len(xi)),mask=mask).compressed()
    nc = len(cii)
    cjj = N.arange(nc)

    # Identify gaps
    gaps = cii[1:]-cii[:-1]
    igaps = N.ma.masked_where(gaps == 1,cjj[:-1]).compressed()

    # Output variable to refill
    vo = vi.clone()

    # Loop on gaps
    for igap in igaps:
        gap = gaps[igap]
        if padding is None:
            pad = max([min_padding,gap*5])
        else:
            pad = max([min_padding,padding])
        istart = max([0,igap-pad+1])
        iend = min([nc,igap+1+pad])
        splines = S.interpolate.splrep(cxi[istart:iend],cvi[istart:iend],k=k)
        vo[cii[igap]+1:cii[igap+1]] = I.interpolate.splev(xi[cii[igap]+1:cii[igap+1]],splines)


    return vo


def fill1d(vari, axis=0, method='linear', maxgap=0):
    """Fill missing values of a 1D array using interpolation.

    :Params:

        - **vari**: Input :mod:`cdms2` variable
        - *axis*: Axis number on which filling is performed
        - *method*: Interpolation method (see :func:`interp1d`)
        - *maxgap*: Maximal size of filled gaps (in steps)

    :Example:

        >>> fill1d(vari, axis=2, method='cubic', maxgap=5)

    :Return:
        Filled :mod:`cdms2` variable similar to input one
    """

    if vari.mask is MV2.nomask: return vari.clone()
    iaxis = vari.getAxis(axis)
    nx = len(iaxis)
    ny = vari.size/nx

    if vari.getMissing() is None:
        vari.setMissing(1.e20)

    # We need a 2D variable
    if vari.ndim == 1:
        vari2d = MV2.reshape(vari, (1, nx))
        vari2d.setAxis(1, vari.getAxis(0))
    else:
        # Numeric
        varin = vari.filled()
        # Order
        if axis != vari.ndim-1:
            varis = N.rollaxis(varin, axis, vari.ndim)
        else:
            varis = varin
        del varin
        # 2D
        ndshape = varis.shape
        if vari.ndim == 2:
            vari2d = varis
        else:
            vari2d = varis.reshape((ny, nx))
        # MV2
        vari2d = MV2.masked_values(vari2d, vari.getMissing(), copy=0)
        vari2d.setAxis(-1, iaxis)
        vari2d.setMissing(vari.getMissing())
        del varis


    # Loop on extra dim
    vari2dn = vari2d.filled()
    varo2d = vari2d.clone()
    if maxgap >=1:
        ii = N.arange(nx)
        dd = N.ones(nx-1)
    keep = N.ones(nx, '?')
    vari2dm = N.ma.ones(nx)
    for iy in xrange(ny):

        # Mask
        # - base mask
        xmask = vari2d[iy].mask
        if xmask is MV2.nomask or xmask.all() or ~xmask.any(): continue

        # - check gaps
        keep[:] = ~xmask
        if maxgap >=1:
            dd[:] = N.diff(xmask.astype('i'))
            gapstart = ii[:-1][dd==1]+1
            gapend = ii[1:][dd==-1]
            ie0 = int(gapend[0]<gapstart[0])
            is1 = (gapstart[-1]>gapend[-1]) and -1 or nx
            for istart, iend in zip(gapstart[:is1], gapend[ie0:]):
                if (iend-istart+1) > maxgap:
                    keep[istart:iend+1] = True
            del gapstart, gapend

        # Compressed axis
        caxis = cdms2.createAxis(iaxis[keep])
        cp_atts(iaxis, caxis)

        # Compressed variable
        vari2dm[:] = vari2d[iy].asma()
        vari2dc = MV2.asarray(vari2dm[keep])
        vari2dc.setAxis(0, caxis)

        # Interpolate
        varo2d[iy] = regrid1d(vari2dc, iaxis, axis=0, method=method)
#       iaxisc = cdms2.createAxis(iaxis.getValue()[~keep])
#       if hasattr(iaxis, 'units'): iaxisc.units = iaxis.units
#       varo2d[iy][~keep] = regrid1d(vari2dc, iaxisc, axis=0, method=method)
        del caxis, vari2dc

    if maxgap >=1: del ii, dd
    del keep, vari2dm

    # Back to correct dims
    varo = vari.clone()
    if vari.ndim == 1:
        varo[:] = MV2.where(vari.mask, varo2d[0], varo)
        del varo2d
    else:
        varo2dn = varo2d.filled()
        del varo2d
        # 2D
        if vari.ndim == 2:
            varon = varo2dn
        else:
            varon = varo2dn.reshape(ndshape)
        del varo2dn
        # Order
        if axis != vari.ndim-1:
            varos = N.rollaxis(varon, vari.ndim-1, axis)
        else:
            varos = varon
        del varon
        # MV2
        varo[:] = MV2.where(vari.mask, MV2.masked_values(varos, vari.getMissing(), copy=0), varo)
        del varos
    return varo


def fill2d(var, xx=None, yy=None, mask=None, copy=True, **kwargs):
    """Fill missing value of 2D variable using inter/extrapolation

    - **var**: A cdms 2D variable.
    - *xx/yy*: Substitutes for axis coordinates [default: None]
    - Other keywords are passe to :func:`griddata`
    """
    # Checkings
    if not cdms.isVariable(var): var = cdms.createVariable(var)
    if copy: var = var.clone()
    assert var.ndim >= 2, 'Input var must be at least 2D'

    # Get Axes
    xo = get_axis(var, -1)
    yo = get_axis(var, -2)
    if xx is None: xx = xo.getValue()
    if yy is None: yy = yo.getValue()
    if xx.ndim != 2 or yy.ndim !=2:
        xx, yy = N.meshgrid(xx, yy)
    xx = xx.astype('d')
    yy = yy.astype('d')

    # Var
    assert xx.shape == var.shape[-2:], '2D axes and variable are not compatible in shape (%s against %s)' %(xx.shape, var.shape)
    nex = var.size/(var.shape[-1]*var.shape[-2])
    if var.ndim != 3:
        var3d = var.reshape(nex, var.shape[-2],  var.shape[-1])
    else:
        var3d = var

    # Mask
    if mask is None:
        mask = var3d.mask
    elif mask is not MV.nomask:
        assert mask.shape[-2:] == var3d.shape[-2:]
        mask = mask.copy()
        if mask.ndim!=2:
            mask.shape = var3d.shape

    # Loop on extra dimensions
    for iex in xrange(nex):

        var2d = var3d[iex].asma()
        if mask is not MV.nomask:
            if mask.shape==2:
                mask2d = mask
            else:
                mask2d = mask[iex]

        # Unmasking
        if mask is MV.nomask or mask2d.all() or ~mask2d.any(): continue
        good = ~mask2d
        xi = xx[good]
        yi = yy[good]
        zi = var2d[good]

        # Gridding and filling
        filled = griddata(xi, yi, zi, (xo, yo), compress=True, **kwargs)
        var3d[iex] = MV2.where(good, var2d, filled)
        del filled, xi, yi,  good

    # Shape back
    if var3d.ndim == 3:
        var[:] = var3d[:]
    else:
        var[:] = var3d.reshape(var.shape[:-2]+(var.shape[-2],  var.shape[-1]))
    return var


def regular(vi,dx=None,verbose=True,auto_bounds=False):
    """Fill a variable with missing values when step of first axis is increasing

    - **vi**: Input array on almost regular axis

    - *dx*: Force grid step to this. Else, auto evaluated.

    """
    from vacumm.misc import cp_atts

    # Guess info from input var
    axes = vi.getAxisList()
    xi = axes[0]
    if len(xi) < 2: return vi
    xid = xi.getValue().astype('d')
    xtc = xi.dtype.char
    missing_value = vi.getMissing()
    sh = vi.shape

    # Get dx and gaps
    ddx = (xid[1:]-xid[:-1])
    if dx is None:
        dx = N.minimum.reduce(ddx)
    gaps = ddx / dx - 1.
    gaps = N.where(N.less((gaps+1.) % 1.,0.1),N.floor(gaps),gaps)
    gaps = N.where(N.greater((gaps+1.) % 1.,0.9),N.ceil(gaps),gaps).astype('l')
    ngaps = N.sum(gaps)
    gaps = MV.masked_object(gaps,0)
    mask = MV.getmaskarray(gaps)
    gaps = gaps.compressed()
    igaps = MV.masked_array(N.arange(len(xi)-1,typecode='l'),mask=mask).compressed()
    ng = len(igaps)

    if not ngaps:
        return vi
    if verbose:
        print 'Filled %i gaps with missing values' % ngaps

    # Init output var
    nxi = len(xi)
    nxo = nxi+ngaps
    sho = list(sh) ;  sho[0] = nxo
    xo = N.zeros(nxo,typecode='d')
    vo = MV.resize(vi,sho)
    #vo.id = vi.id+'_regular' ; vo.name = vo.id

    # Loop on first axis
    j = igaps[0]+1
    xo[0:j] = xi[0:j]
    vi[0:j] = vi[0:j]
    for ig,i in enumerate(igaps):
        # Fill the gap
        gap = gaps[ig]
        xo[j:j+gap] = N.arange(1.,gap+1).astype(xtc)*dx + float(xi[i])
        vo[j:j+gap] = missing_value
        # Fill values after the gap
        i += 1
        j += gap
        if ig == ng-1:
            full = nxi-i
        else:
            full = igaps[ig+1]-i+1
        xo[j:j+full] = xi[i:i+full]
        vo[j:j+full] = vi[i:i+full]
        j += full

##  xo[-1] = xi[-1]
##  vo[-1] = vi[-1]

    # Axes
    vo = MV.masked_object(vo,missing_value)
    xo = cdms.createAxis(xo)
    cp_atts(xi,xo,id=True)
    axeso = list(axes)
    axeso[0] = xo
    vo.setAxisList(axeso)
    vo.setGrid(vi.getGrid())
    if auto_bounds:
        vo.getAxis(0).setBounds(bounds1d(vo.getAxis(0)))
    else:
        vo.getAxis(0).setBounds(None)

    return vo



def regular_fill1d(var,k=1,dx=None):
    """Combination: fill1d(regular_fill)) (with their parameter)"""
    return fill1d(regular(var,dx=dx),k=k)


def spline_interp1d(old_var,new_axis,check_missing=True,k=3,**kwargs):
    """Backward compatibility function

    See :func:`regrid1d`
    """
    return regrid1d(old_var,new_axis,method=['nearest', 'linear', 'cubic'][k], **kwargs)


class GridData(object):
    """2D interpolator from a randomly spaced sample on a regular grid

    Possible algorithms:

    - Natural Neighbor using nat.Natgrid:

        - http://www.ncarg.ucar.edu//ngmath/natgrid/nnhome.html
        - http://www.cisl.ucar.edu/zine/98/spring/text/3.natgrid.html
        - http://dilbert.engr.ucdavis.edu/~suku/nem/nem_intro/node3.html
        - http://www.ems-i.com/gmshelp/Interpolation/Interpolation_Schemes/Natural_Neighbor_Interpolation.htm

    - 2D splines using css.css.Cssgrid

    :Parameters:

        - **xi**: Input 1D X positions.
        - **yi**: Input 1D Y positions.
        - **ggo**: Output grid. Can either (xo,yo), a cdms grid or a cdms variable with a grid.
        - *method*: Interpolator type, either 'nat' (Natural Neighbors) or 'css' (='splines' using splines) [default: 'nat']
        - *nl*: Nonlinear interpolator (usually gives better results) [default: False]
        - *ext*: Extrapolate value outsite convex hull [**'nat' only**, default: False]
        - *mask*: Mask to apply to output data [default: None]
        - *compress*: If ``True``, interolate only unmasked data, and thus does not try guess the best mask (that's more efficient but very bad if data are masked!).
        - *sub*: Size of blocks for subblocking [**"nat" only**]
        - *margin*: Margin around ouput grid (or block) to select input data.
          Value is relative to X and Y extent.
        - Other keywords are set as attribute to the interpolator instance ; to get the list of parameters: ::

            >>> import nat ; nat.printParameterTable()
            >>> import css ; css.printParameterTable()

    :Example:

    >>> r = GridData(gridi, grido, method='nat', ext=False, margin=.7)
    >>> varo1 = r(vari1)
    >>> varo2 = r(vari2)
    """

    def __init__(self, xi, yi, ggo, nl=False, ext=False, geo=None, method='nat',
        sub=30, margin=.5, compress=False, **kwargs):

        vcdwarn('GridData will no longer be supported in future version')

        # Helper
        self._GDH = _GridDataHelper_(xi, yi, ggo, geo=geo, compress=compress)

        # Init interpolator
        if method in ['css', 'splines']:
            raise VACUMMError('css/spline method is no longer supported')
        elif method in ['nat', None, 'natgrid']:
            from nat import Natgrid
            self.r = Natgrid
            self.sub = sub
            self.method = 'nat'

        self.margin = margin

        # Attributes of interpolator
        kwargs['igr'] = 0
        kwargs['ext'] = ext
        kwargs.setdefault('hor', -1)
        kwargs.setdefault('ver', -1)
        self.ratts = kwargs
        if self.sub is not None:
            self._GDH.xo.shape = self._GDH.grid.shape
            self._GDH.yo.shape = self._GDH.grid.shape

    def __call__(self, zi, missing_value=None, **kwargs):
        """Interpolate zi on output grid

        - **zi**: At least a 1D array.
        """

        # Init
        zi2d, zo3d, mo3d, nex = self._GDH.init_data(zi, missing_value)

        # Prepare sub-blocks
        if not self.sub:
            jbs = xrange(1)
            ibs = xrange(1)
        else:
            jbs = xrange((self._GDH.ny-1)/self.sub+1)
            ibs = xrange((self._GDH.nx-1)/self.sub+1)

        # Loop on supplementary dims
        for iex in xrange(nex):

            # Get data and mask
            get = self._GDH.get(zi2d, iex)
            if get is None: continue
            xi, yi, zzi, mmi = get
            constant = N.allclose(zzi, zzi[0])
            if len(xi)<4: break
            # Interpolate by blocks
            for ib in ibs:
                for jb in jbs:

                    # Block
                    if self.sub is None:
                        xslice = slice(None)
                        yslice = slice(None)
                    else:
                        xslice = slice(ib*self.sub,(ib+1)*self.sub)
                        yslice = slice(jb*self.sub,(jb+1)*self.sub)
                    xo = self._GDH.xo[yslice, xslice].ravel().astype('d')
                    yo = self._GDH.yo[yslice, xslice].ravel().astype('d')

                    # Restriction of input data around output grid
                    if self.margin > 0:
                        margin = self.margin
                        ngoodi = 0
                        igoodi = 0
                        enlarge = .5
                        while ngoodi < 5: # Loop to get a sufficient number of pts
                            margin *= (1+enlarge*igoodi)
                            dxo = xo.ptp()*margin
                            dyo = yo.ptp()*margin
                            ximin = xo.min()-dxo
                            ximax = xo.max()+dxo
                            yimin = yo.min()-dyo
                            yimax = yo.max()+dyo
                            goodi = xi>ximin
                            goodi &= xi<ximax
                            goodi &= yi>yimin
                            goodi &= yi<yimax
                            ngoodi = goodi.sum()
                            igoodi += 1
                            if igoodi>100: break
                        if igoodi>100: continue

                    else:
                        goodi = slice(None)

                    # Interpolator
                    interpolator = self.r(xi[goodi], yi[goodi], xo, yo, listOutput = 'yes')
                    for att, val in self.ratts.items():
                        setattr(interpolator, att, val)
                    interpolator.nul = missing_value

                    # Interpolation
                    if constant:
                        sh = self._GDH.xo[yslice, xslice].shape
                        zzo = N.zeros(sh)
                        zzo += zzi[0]
                        if mmi is not None:
                            mmo = N.zeros(sh)
                            mmo += mmi[0]
                    else:

                        sh = self._GDH.xo[yslice, xslice].shape
                        zzo = interpolator.rgrd(zzi[goodi]).reshape(sh)
                        if mmi is not None:
                            mmo = interpolator.rgrd(mmi[goodi]).reshape(sh)

                    zo3d[iex, yslice, xslice] = zzo
                    if mmi is not None:
                        mo3d[iex, yslice, xslice] = mmo
                        del mmo
                    del zzo

            if hasattr(zi2d, 'mask') and zi2d[iex].mask is not N.ma.nomask:
                del zzi, xi, yi
        gc.collect()

        del zi2d

        return self._GDH.format(zo3d, mo3d, **kwargs)

    rgrd = __call__
    regrid = __call__

def cargen(xi, yi, zi, ggo, mask=None, compress=False, missing_value=None, **kwargs):
    """Interpolator from IFREMER

    :Params:

        - **xi**: Input 1D X positions.
        - **yi**: Input 1D Y positions.
        - **ggo**, optional: Output grid. Can either (xo,yo), a cdms grid or a cdms variable with a grid.
        - **mask**, optional: Mask to apply to output data [default: None]
    """
    # Helper
    GDH = _GridDataHelper_(xi, yi, ggo, mask=mask, compress=compress,
                           proj=False)
    assert GDH.grid_type == 'rect', 'cargen works only with rectangular grids'

    # Init data
    zo3d, mo3d = GDH.init_data(zi, missing_value)

    # Loop on supplementary dims
    for iex in xrange(GDH.nex):

        # Get data and mask
        get = GDH.get(iex)
        if get is None: continue
        xi, yi, zzi, mmi = get

        # Interpolate
        zo3d[iex] = _cargen_(GDH.xi, GDH.yi, zzi, GDH.xo, GDH.yo,
            GDH.missing_value).T
        if mmi is not None:
            mo3d[iex] = _cargen_(GDH.xi, GDH.yi, mmi, GDH.xo, GDH.yo, 1.).T

    # Format output
    return GDH.format(zo3d, mo3d, **kwargs)

krigdata = cargen

def krig(xi, yi, zi, ggo, mask=True, proj=True, missing_value=None, **kwargs):
    """Kriging interpolator to a grid

    :Params:

        - **xi**: Input 1D X positions.
        - **yi**: Input 1D Y positions.
        - **zi**: Input N-D with last dim as space.
        - **ggo**, optional: Output grid. Can either (xo,yo), a cdms grid or a cdms variable with a grid.
        - **mask**, optional: Mask to apply to output data [default: None]
    """
    # Helper
    assert proj is not False
    GDH = _GridDataHelper_(xi, yi, ggo, proj=proj, mask=mask, compress=True)

    # Init data
    zo3d, mo3d  = GDH.init_data(zi, missing_value)
    distfunc = 'haversine' if GDH.geo else 'simple'

    # Get data and mask
    get = GDH.get(zi2d, Ellipsis, compress, missing_value)
    if get is not None:
        gxi, gyi, zzi, mmi = get

        # Interpolate
        zo3d[:] = _krig_(gxi, gyi, zzi, GDH.xo, GDH.yo,
                         distfunc=distfunc).reshape(zo3d.shape)
        del gxi, gyi, zzi, mmi

    # Format output
    return GDH.format(zo3d, mo3d, **kwargs)



def griddata(xi, yi, zi, ggo, method='linear', cgrid=False, cache=None,
             proj=True, **kwargs):
    """Interpolation in one single shot using GridData

    :Params:

        - **xi**: 1D input x coordinates.
        - **yi**: 1D input y coordinates (same length as xi).
        - **zi**: 1D input values (same length as xi).
        - **method**, optional: Method of interpolation,
          within ('nat', 'css', 'carg', 'krig') [default: 'carg']
        - **cgrid**, optional: Output on a C-grid at U- and V-points
          deduced from ggo [default: False].
          Not available for 'carg' and 'krig' methods.

    :See also: :class:`GridData` and :func:`cargen`
    """
    if cgrid:
        return tuple([griddata(xi, yi, zi, gg, method=method, cgrid=False,
                               **kwargs) for gg in t2uvgrids(ggo)])

    # Method
    if method in ('css', 'spline', 'nat'):
        vcdwarn('Interpolation method other than nearest, linear or cubic'
                ' will not be supported in future versions')
        if method == 'nat':
            method = 'linear'
        elif method in ('css', 'spline'):
            method = 'cubic'
    assert method in _griddata_methods, ('Invalid interpolation method {}.'
                                         ' Please choose one of:').format(
                                                 method,
                                                 ' '.join(_griddata_methods))
    assert proj is not False

    # Scipy
    if method in ['nearest', 'linear', 'cubic']:

        # Helper
        GDH = _GridDataHelper_(xi, yi, ggo, compress=False, proj=proj)

        # Data
        GDH.init_data(zi, initout=False)
        if GDH.masked and not GDH.compress and method == 'cubic':
            vcwarn("Can't interpolate masked with cubic method "
                   "without compression. Switching compression on.")
            GDH.compress = True

        # Interpolator
        import scipy.interpolate as SI
        interpolator = {'nearest': SI.NearestNDInterpolator,
                        'linear': SI.LinearNDInterpolator,
                        'cubic': SI.CloughTocher2DInterpolator}[method]
        # Get stuff
        xi, yi, zi, mi = GDH.get()
        xyi = N.array([xi, yi]).T
        if method != 'nearest':
            xyi = get_tri(xyi, ttype='scipy', cache=cache)  # triangulation
        xyo = N.array([GDH.xo, GDH.yo]).T

        # Interpolate
        kw = {}
        if method != 'nearest':
            kw['fill_value'] = GDH.missing_value
        zo3d = interpolator(xyi, zi.T, **kw)(xyo).T
        if GDH.masked and not GDH.compress and method != 'nearest':
#            if method != 'nearest':
#                kw['fill_value'] = 1.
            mo3d = interpolator(xyi, mi.T, **kw)(xyo).T
        else:
            mo3d = None

        # Format output
        return GDH.format(zo3d, mo3d, **kwargs)

    if method=='krig':
        return krig(xi, yi, zi, ggo, **kwargs)

    return cargen(xi, yi, zi, ggo, **kwargs)


class _GridDataHelper_(object):

    def __init__(self, xi, yi, ggo, mask=None, compress=False, proj=None,
                 gtype=None):

        # Input grid
        assert xi.shape == yi.shape and xi.ndim == 1,   'xi and yi must be 1d arrays'
        xi = N.asarray(xi,dtype='d')
        yi = N.asarray(yi,dtype='d')
        proj = get_proj((xi, yi), proj=proj)
        self.mi = N.ones(len(xi), '?')
        if proj:
            self.xi, self.yi = proj(xi, yi)
        else:
            self.xi, self.yi = xi, yi


        # - grid and axes
        self.ongrid = isgrid(ggo)
        self.grid = get_grid(ggo, gtype=gtype)
        self.grid_type = get_grid_type(self.grid)
        self.x = self.grid.getLongitude()
        self.y = self.grid.getLatitude()
        self.xo = self.x[:]
        self.yo = self.y[:]
        if self.grid_type == "rect" and proj:
            self.xo, self.yo = N.meshgrid(self.xo, self.yo)
        if self.grid_type != 'unstruct':
            self.xo.shape = -1,
            self.yo.shape = -1,
        if proj:
            self.xo, self.yo = proj(self.xo, self.yo)
        self.axes = []
        self.nx = self.grid.shape[-1]
        self.ny = self.grid.shape[0]

        # - mask
        if mask is True:
            if hasattr(ggo, 'getMask'):
                mask = ggo.getMask()
            elif hasattr(ggo, 'mask'):
                mask = ggo.mask
        if mask is False or mask is None:
            mask = N.ma.nomask
        if mask.ndim > 2:
            mask = mask[(0,)*(mask.ndim-2)]
        self.mask = mask
        self.compress = compress
        self.inited = False

    def init_data(self, zi, missing_value=None, initout=True):
        """Init input before regridding"""

        # Convert to right dims
        self.inited = True
        self.zi = zi
        nex = zi.size/zi.shape[-1]
        if cdms2.isVariable(zi):
            zi = zi.asma()
        self.zi2d = zi.reshape(nex, zi.shape[-1]).astype('d').copy()

        # Masking
        if missing_value is None:
            if N.ma.isMA(zi):
                missing_value = zi.get_fill_value()
            else:
                missing_value = N.ma.default_fill_value(zi)
        self.missing_value = missing_value
        unmasked = (not hasattr(self.zi2d, 'mask')
                    or self.zi2d.mask is N.ma.nomask
                    or not self.zi2d.mask.any())
        self.masked = not unmasked
        self.nex = nex
        if not initout:
            return

        # Initialize output
        zo3d = N.zeros(self.zi2d.shape[:1]+self.grid.shape, zi.dtype) + missing_value
        if unmasked:
            compress = False
        else:
            compress = self.compress
        if compress or unmasked:
            mo3d = None
        else:
            mo3d = zo3d*0.
        del unmasked

        return zo3d, mo3d

    def get(self, iex=None):
        """Get a slice"""
        if not self.inited:
            raise VACUMMError('Please call .init_data() before calling .get()')

        if iex is None:
            iex = Ellipsis

        # Remove compress values or fill them
        good = self.mi # 1D
        mmi = None
        if self.masked and self.zi2d[iex].mask.any():
            if not self.compress: # No compression => interpolation of mask
                mmi = self.zi2d[iex].mask.astype('f')
                self.zi2d[iex] = self.zi2d[iex].filled(self.missing_value)
            else:
                good = good & ~self.zi2d[iex].mask
        if good.ndim==2:
            good = good.all(axis=0)

        # All data are bad
        if not good.any(): return None

        if not good.all(): # Some are good
            zzi = self.zi2d[iex][..., good]
            xi = self.xi[good]
            yi = self.yi[good]
            if mmi is not None:
                mmi = mmi[good]
        else: # There are all good
            zzi = self.zi2d[iex]
            xi = self.xi
            yi = self.yi
        if N.ma.isMA(zzi):
            zzi = zzi.filled(self.missing_value)
        del good
        return xi, yi, zzi, mmi


    def format(self, zo3d, mo3d, clipval=False, asmv=None, **kwargs):
        """Format output variable"""
        # Shape
        sho = self.zi.shape[:-1] + self.grid.shape
        if zo3d.shape != sho:
            zo = zo3d.reshape(sho)
            del zo3d
            if mo3d is not None:
                mo = mo3d.reshape(sho)
        else:
            zo = zo3d
            mo = mo3d

        # Masking
        if self.mask is not N.ma.nomask:
            mask = self.mask
            if mask.shape != zo.shape:
                mask = N.resize(mask, zo.shape)
            zo[mask] = self.missing_value
        if mo3d is not None and mo is not None:
            zo[mo>0.] = self.missing_value
        missing = closeto(zo, self.missing_value)

        # Value clipping
        if clipval:
            valmax = self.zi.max()
            valmin = self.zi.min()
            good = ~missing
            zo[good&(zo>valmax)] = valmax
            zo[good&(zo<valmin)] = valmin

        # Format
        zo = N.ma.masked_where(missing, zo, copy=0)
        if asmv is None:
            asmv = self.ongrid or self.grid_type in ('rect', 'curv')
        if not asmv:
            return zo
        zo = MV2.asarray(zo)
        set_grid(zo, self.grid)
        if cdms.isVariable(self.zi):
            cp_atts(self.zi, zo, id=True)
            for i, axis in enumerate(self.zi.getAxisList()[:-len(self.grid.shape)]):
                zo.setAxis(i, axis)
        return zo

def xy2grid(*args, **kwargs):
    """Alias for :func:`griddata`

    .. seealso::

        :class:`GridData` :func:`cargen` :func:`xy2grid`
    """
    return griddata(*args, **kwargs)

def xy2xy(xi, yi, zi, xo, yo, method='linear', proj=True, cache=None, **kwargs):
    """Interpolation between to unstructured grids using scipy

    :Params:

        - **xi/yi**: 1D input positions
        - **zi**: atleast-1D input values
        - **xo,yo**: 1D output positions
        - *proj*: convert positions to meters using mercator projection
    """
    # Helper
    assert method in ('nearest', 'linear', 'cubic')
    assert proj is not False
    GDH = _GridDataHelper_(xi, yi, (xo, yo), compress=False, proj=proj,
                           gtype='unstruct')

    # Data
    GDH.init_data(zi, initout=False)
    if GDH.masked and not GDH.compress and method == 'cubic':
        vcwarn("Can't interpolate masked with cubic method "
               "without compression. Switching compression on.")
        GDH.compress = True

    # Interpolator
    import scipy.interpolate as SI
    interpolator = {'nearest': SI.NearestNDInterpolator,
                    'linear': SI.LinearNDInterpolator,
                    'cubic': SI.CloughTocher2DInterpolator}[method]

    # Get stuff
    xi, yi, zi, mi = GDH.get()
    xyi = N.array([xi, yi]).T
    if method != 'nearest':
        xyi = get_tri(xyi, ttype='scipy', cache=cache)  # triangulation
    xyo = N.array([GDH.xo, GDH.yo]).T

    # Interpolate
    kw = {}
    if method != 'nearest':
        kw['fill_value'] = GDH.missing_value
    zo3d = interpolator(xyi, zi.T, **kw)(xyo).T
    if GDH.masked and not GDH.compress and method != 'nearest':
#        if method != 'nearest':
#            kw['fill_value'] = 1.
        mo3d = interpolator(xyi, mi.T, **kw)(xyo).T
    else:
        mo3d = None

    # Format output
    kwargs.setdefault('asmv', False)
    return GDH.format(zo3d, mo3d, **kwargs)

#    # Check input positions
#    xi = N.asarray(xi, 'd')
#    yi = N.asarray(yi, 'd')
#    xo = N.asarray(xo, 'd')
#    yo = N.asarray(yo, 'd')
#    proj = kwargs.pop('geo', proj)
#    if proj:
#        if not callable(proj):
#            proj = get_proj((xi,yi))
#        xi, yi = proj(xi, yi)
#        xo, yo = proj(xo, yo)
#
#    # Check input type
#    outtype = 0
#    if cdms2.isVariable(zi):
#        outtype = 2
#        axes = zi.getAxisList()
#        atts = get_atts(zi)
#        zi = zi.asma()
#    elif N.ma.isMA(zi):
#        if zi.mask is not N.ma.nomask and zi.mask.any():
#            outtype = 1
#
#    # Check shapes
#    zi = zi.copy()
#    si = zi.shape
#    nsi = zi.shape[-1]
#    nex = zi.size/nsi
#    zi.shape = (nex, nsi)
#    nso = len(xo)
#    zo = N.zeros((nex, nso))
#    zo[:] = N.nan
#    if outtype:
#        mo = N.zeros((nex, nso))
#        goodi = ~zi.mask
#    else:
#        mo = None
#
#    # Loop on extra dim
#    from masking import convex_hull, polygon_select
#    for iex in xrange(nex):
#
#        if outtype:
#            gi = goodi[iex]
#            mi = zi.mask[iex].astype('f')
#        else:
#            gi = slice(None)
#
#        # Build regridder only when needed
#        if iex==0 or (outtype and N.any(goodi[iex-1]!=goodi[iex])):
#
#            # Check that output points are inside convex hull
#            hull = convex_hull((xi[gi], yi[gi]), poly=True)
#            go = polygon_select(xo, yo, [hull], mask=2) ; del hull
#            if go.all():
#                del go
#                go = slice(None)
#
#            # Regridder
#            from nat import Natgrid
#            r = Natgrid(xi[gi], yi[gi], xo[go], yo[go], listOutput='yes')
#            r.igr = int(nl)
#            if outtype:
#                rm = Natgrid(xi, yi, xo[go], yo[go], listOutput='yes')
#                rm.igr = int(nl)
#
#        # Regridding
#        # - values
#        zin = zi[iex][gi]
#        if N.ma.isMA(zin): zin = zin.filled()
#        zo[iex][go] = r.rgrd(zin)
#        # - mask
#        if outtype:
#            mo[iex][go] = rm.rgrd(mi)
#    del zi
#
#    # Missing points
#    mnan = N.isnan(zo)
#    if mnan.any():
#        outtype = max(1, outtype)
#        zo = N.ma.masked_where(mnan, zo)
#        del mnan
#
#    # Return pure numeric
#    zo.shape = si[:-1]+(nso, )
#    if outtype==0: return zo
#
#    # Masking
#    if mo is not None:
#        mo.shape = zo.shape
#        zo = N.ma.masked_where(mo>.5, zo)
#        del mo, mi
#
#    # Masked arrays
#    if outtype == 1: return zo
#
#    # cdms
#    zo = MV2.asarray(zo)
#    set_atts(zo, atts)
#    for i, ax in axes[:-1]:
#        zo.setAxis(i, ax)
#    return zo



def grid2xy(vari, xo, yo, zo=None, to=None, zi=None, method='linear', outaxis=None,
            distmode='haversine'):
    """Interpolate gridded data to ramdom positions

    :Params:

        - **vari**: Input cdms variable on a grid
        - **xo**: Output longitudes
        - **yo**: Output latitudes
        - **method**, optional: Interpolation method

            - ``nearest``: Nearest neighbor
            - ``linear``: Linear interpolation

        - **zo**, optional: Output depths (negative in the ocean).
        - **to**, optional: Output times.
        - **zi**, optional: Input depths when variable in space.

        - **outaxis**, optional: Output spatial axis

            - A cdms2 axis.
            - ``None`` or ``'auto'``: Longitudes or latitudes depending
              on the range if coordinates are monotonic, else ``'dist'``.
            - ``'lon'`` or ``'x'``: Longitudes.
            - ``'lat'`` or ``'y'``: Latitudes.
            - ``'dist'`` or ``'d'``: Distance in km.

        - **distmode**, optional: Distance computation mode.
          See :func:`~vacumm.misc.grid.misc.get_distances`.

    """

    # Prefer 1D axes
    grid = get_grid(vari)
    grid = curv2rect(grid, mode=False)

    # Format

    # - input data and coordinates
    xi, yi = get_xy(grid, num=True)
    rect = xi.ndim==1
    mv = vari.getMissing()
    if mv is None:
        mv = vari.getMissing()
    order = vari.getOrder()
    vi = vari.filled(mv).astype('d')
    if 'z' not in order and zo is not None:
        zo = None
    if 't' not in order and to is not None:
        to = None
    univ = to is not None or zo is not None or True # FIXME: univ
    na = 2 # interpolation dims
    extra_axes = []
    if zo is not None:
        if zi is None:
            zi = vari.getLevel()
        zi = zi[:]
    if univ:
        if method=='nat':
            vacumm_warn('"nat" method not available with time or depth'
                        ' interpolation. Switching to linear')
            method = 'linear'

        it = 2 + ('z' in order)
        if 't' not in order and to is None: # not requested and not present
            vi = vi.reshape(vi.shape[:-it]+(1, )+vi.shape[-it:]) # insert fake dim
        elif to is not None: # requested and present
            na += 1
            taxis = vari.getTime()
            ti = taxis[:]
        else: # present but not requested
            vi = vi.reshape(vi.shape[:-it]+(1, )+vi.shape[-it:]) # insert fake dim
            extra_axes.insert(0, vari.getTime())

        if 'z' not in order and zo is None: # not requested and not present
            vi = vi.reshape(vi.shape[:-2]+(1, )+vi.shape[-2:]) # insert fake dim
        elif zo is not None: # requested and present
            na += 1
            if zi.ndim<3: # no xy assumed
                zi = zi.reshape(zi.shape+(1, 1))
            if zi.ndim < vari.ndim: # with xy: left pad with ones
                zi = zi.reshape((1, )*(vari.ndim-zi.ndim)+zi.shape)
        else: # present but not requested
            vi = vi.reshape(vi.shape[:-2]+(1, )+vi.shape[-2:]) # insert fake dim
            if extra_axes: # t as extra: ET1Z1YX -> ETZ11YX
                vi = vi.reshape(vi.shape[:-5]+vi.shape[-4:-3]+(1, 1)+vi.shape[-2:])
            else: # ETZ1YX -> EZT1YX
                vi = N.moveaxis(vi, -4, -5)
            extra_axes.append(vari.getLevel())

#        if zo is not None:
#            zi = zi.reshape((-1, )+zi.shape[-4:])
        if xi.ndim==1:
            xi = xi[None, :]
        if yi.ndim==1:
            yi = yi[:, None]


        extra_axes = ([ax for ax, o in zip(vari.getAxisList(), order) if o=='-']
                      + extra_axes)

    else:

        na = 2
        ne = vi.ndim - 2
    vi = vi.reshape((-1,)+vi.shape[-4:])
#    if zi.shape[0] != vi.shape[0]:
#        zi = N.resize(zi, vi.shape[:1] + zi.shape[1:])


    if vari.mask is not MV2.nomask:
        mi = vari.mask.astype('d').reshape(vi.shape)
    else:
        mi = None

    # - output coordinates
    isscalar = N.isscalar(xo)
    xo = N.atleast_1d(N.asarray(xo, dtype='d'))
    yo = N.atleast_1d(N.asarray(yo, dtype='d'))
    no = max(xo.size, yo.size)
    if zo is not None:
        no = max(zo.size, no)
    if to is not None:
        no = max(len(to), no)
    assert xo.ndim==1 and yo.shape==xo.shape, ('xo and yo must be scalars 1d or arrays'
                                               ' of equal lengths')
    if abs(xi.min()-xo.min())>180.:
        if xi.min()<xo.min():
            xi = xi + 360.
        else:
            xi = xi - 360.
    if zo is not None:
        zo = N.atleast_1d(zo[:])
        assert zo.shape==xo.shape,  ('zo must have the same shape as xo and yo')
    else:
        zo = N.zeros(no, 'd')
        zi = N.zeros([1]*5, 'd')
    if to is not None:
        to = A.create_time(to)
        to.toRelativeTime(taxis.units)
        to = to[:]
        assert to.shape== xo.shape,  ('"to" must have the same shape as xo and yo')
    else:
        to = N.zeros(no, 'd')
        ti = N.zeros(1, 'd')

    # Interpolate
    if method is None:
        method = 'linear'
    else:
        method = str(method)
    if univ:
            targets = [-4, -3]
            if rect:
                targets.extend([-2, -1])
            subdims = {-3:-3, -2:-2, -1:-1}
    if method == 'nearest' or (method=='nat' and mi is not None):

        if univ:
            vi, [ti, zi, yi, xi] = _monotonise_(vi, [ti, zi, yi, xi],
                targets=targets, subdims=subdims)
            vo = _nearest4dto1dxx_(xi, yi, zi, ti, vi, xo, yo, zo, to, mv)
        else:
            if rect:
                func = _nearest2dto1d_
                vi, [yi, xi] = _monotonise_(vi, [yi, xi])
            else:
                func = _nearest2dto1dc_
            vo = func(xi, yi, vi, xo, yo, mv)
            if method=='nat' and mi is not None:
                vonear = vo

    if 'linear' in method:

        if univ:
            vi, [ti, zi, yi, xi] = _monotonise_(vi, [ti, zi, yi, xi],
                targets=targets, subdims=subdims)
            vo = _linear4dto1dxx_(xi, yi, zi, ti, vi, xo, yo, zo, to, mv)
        else:
            if rect:
                func = _bilin2dto1d_
                vi, [yi, xi] = _monotonise_(vi, [yi, xi])
            else:
                func = _bilin2dto1dc_
            vi, [yi, xi] = _monotonise_(vi, [yi, xi])
            vo = func(xi, yi, vi, xo, yo, mv)


#    elif method.startswith('nat'):
#
#        # Build regridder
#        from nat import Natgrid
#        nz = vi.shape[0]
#        vo = N.zeros((nz, len(xo)))+mv
#        if xi.ndim==1:
#            xxi, yyi = xi, yi
#        else:
#            xxi, yyi = N.meshgrid(xi, yi)
#        r = Natgrid(xxi.ravel(), yyi.ravel(), xo, yo, listOutput='yes')
#        if mi is not None:
#            mo = N.zeros(len(xo))
#
#        # Z loop interpolation
#        for iz in xrange(nz):
#
#            # Base
#            vo[iz] = r.rgrd(vi[iz].ravel())
#
#            # Mask
#            if mi is not None:
#                mo[:] = r.rgrd(mi[iz].ravel())
#                vo[iz] = N.where(mo!=0., vonear[iz], vo[iz])
#
#        if mi is not None:
#            del mi, mo, vonear

    elif method != 'nearest':
        raise NotImplementedError, 'Method yet not implemented: '+method

    # Output
    vo.shape = tuple([len(ax) for ax in extra_axes]) + (vo.shape[-1], )
    varo = MV2.masked_values(vo, mv, copy=0)
    cp_atts(vari, varo, id=True)
    if outaxis=='auto':
        outaxis = None
    if outaxis is None or outaxis in ['m', 'pos', 'd', 'dist']:
        if len(xo)>1:
            #xom, yom = get_proj((xo, yo))(xo, yo)
            if outaxis is None:
                outaxis = 'lon' if N.ptp(xo)>N.ptp(yo) else 'lat'
                if ((outaxis=='lon' and (N.diff(N.sign(N.diff(xo)))!=0).any() or (N.diff(xo)==0).any()) or
                    (outaxis=='lat' and (N.diff(N.sign(N.diff(yo)))!=0).any() or (N.diff(yo)==0).any())):
                    outaxis = 'dist'
        elif isscalar:
            outaxis = None
        else:
            outaxis = 'lon'
    elif outaxis=='num' or outaxis is False:
        outaxis = varo.getAxis(-1)
    if outaxis in ['x', 'lon']:
        outaxis = A.create_lon(xo)
    elif outaxis in ['y', 'lat']:
        outaxis = A.create_lat(yo)
    elif to is not None and outaxis in ['t', 'time']:
        outaxis = A.create_time(to)
    elif zo is not None and outaxis in ['z', 'dep', 'depth']:
        outaxis = A.create_dep(zo)
    elif outaxis in ['m', 'd', 'pos', 'dist']:
        dist = N.concatenate(([0], get_distances(xo[1:], yo[1:], xo[:-1], yo[:-1],
            pairwise=True, mode=distmode))).cumsum()*0.001
        #dist = get_distances(xo[0], yo[0], xo, yo, mode=distmode)*0.001
        #dist = N.concatenate(([0],N.sqrt(N.diff(xom)**2+N.diff(yom)**2)*0.001)).cumsum()
        outaxis = cdms2.createAxis(dist, id='position')
        outaxis.long_name = 'Distance along transect'
        outaxis.units = 'km'
    if outaxis is not None and not A.isaxis(outaxis):
        outaxis = cdms2.createAxis(outaxis, id='position')
        outaxis.long_name = 'Position along transect'
    axes = extra_axes
    if outaxis is not None:
        axes.append(outaxis)
    else:
        varo = varo[...,0]
    if N.ndim(varo):
        varo.setAxisList(axes)
    return varo


def transect(var, lons, lats, depths=None, times=None, method='linear',
        subsamp=3, getcoords=False, outaxis=None, depth=None, **kwargs):
    """Make a transect in a -[T][Z]YX variable

    It calls :func:`~vacumm.misc.grid.transect_specs` to compute transect
    coordinates when not explictly specified, and :func:`grid2xy` to perform
    4D interpolations.

    :Example:

        >>> tsst = transect(sst, (1,1.6), (42.,43.), subsamp=4)
        >>> tmld = transect(mld, [1.,1.,2.], [42.,43.,43.], outaxis='dist')
        >>> tprof = transect(temp, 1., 42.)

    :Params:

        - **var**: MV2 array of at least rank 2 (YX).
        - **lons/lats**: Specification of transect, either

            - Coordinates of first and last point in degrees as
              tuples in the form ``(lon0,lon1)`` and ``(lat0,lat1)``.
              The array of coordinates is generated using :func:`transect_specs`.
            - Or explicit array of coordinates (as scalars, lists or arrays).

        - **depths**, optional: Output depths. If not a tuple, it must have
          the same size as lons and lats.
        - **times**, optional: Tuple, or time sequence or axis of the same length as
          resulting coordinates. If provided, the interpolation is first done
          in space, them onto this lagrangian time,
          and the final space-time trajectory is returned.
          If ``outaxis`` is None, ``taxis`` becomes the output axis.

        - **subsamp**, optional: Subsampling with respect to grid cell
          (only when coordinates are not explicit).
        - **method**, optional: Interpolation method (see :func:`grid2xy`).
        - **getcoords**, optional: Also get computed coordiantes.
        - **outaxis**, optional: Output spatial axis (see :func:`grid2xy`).

    :Return:

        ``tvar`` or ``tvar,tons,tlats``

    """
    # Output coordinates
    times = kwargs.get('time', times) # backward compat
    if isinstance(lons, tuple): # Find coordinates
        lon0, lon1 = lons
        lat0, lat1 = lats
        lons, lats = transect_specs(var, lon0, lat0, lon1, lat1, subsamp=subsamp)
        single = False
    else: # explicit coordinates
        single = N.isscalar(lons) and N.isscalar(lats)
        lons = N.atleast_1d(lons)
        lats = N.atleast_1d(lats)
    if depths is not None:
        if isinstance(depths, tuple):
            depths = N.linspace(depths[0], depths[1], len(lons))
        else:
            depths = N.atleast_1d(depths)
            if len(depths)!=len(lons):
                raise VACUMMError('Your depths axis must have a length of: %i (!=%i'%(
                    len(depths), len(lons)))
    if times is None and cdms2.isVariable(lons) and lons.getTime() is not None:
            times = lons.getTime()
    if times is not None:
        if isinstance(times, tuple):
            times = lindates(times[0], times[1], len(lons))
        else:
            if not A.istime(times):
                times = A.create_time(times)
            if len(times)!=len(lons):
                raise VACUMMError('Your times must have a length of: %i (!=%i'%(
                    len(times), len(lons)))

    # Interpolation
    var = grid2xy(var, lons, lats, zo=depths, to=times, zi=depth,
        method=method, outaxis=outaxis)

    # Single point
    if single: var = var[...,0]

    if not getcoords:
        return var
    coords = lons, lats
    if depths is not None:
        coords += depths,
    if times is not None:
        coords += times,
    return (var,) + coords

def refine(vari, factor, geo=True, smoothcoast=False, noaxes=False):
    """Refine a variable on a grid by a factor

    :Params:

        - **vari**: 1D or 2D variable.
        - **factor**: Refinement factor > 1
    """

    # No refinement
    if factor < 2:
        return vari

    # Initialize
    if cdms.isVariable(vari):
        op = MA
    else:
        op = N
    sho = ()
    for nxy in vari.shape:
        sho += ((nxy-1)*factor+1, )
    varo = op.zeros(sho, vari[:].dtype)
    if isinstance(vari, TransientAxis):
        varo = cdms.createAxis(varo.astype('d'))
        cp_atts(vari, varo, id=True)
        geo = False
    elif cdms.isVariable(vari):
        varo = cdms.createVariable(varo)
        cp_atts(vari, varo, id=True)
    else:
        geo = False

    # Fill
    step = 1./factor
    if vari[:].ndim == 1: # 1D linear
        varo[0::factor] = vari[:]
        for i in xrange(1, factor):
            varo[i::factor] = vari[:-1]+(vari[1:]-vari[:-1])*i*step

    else: # 2D bilinear
        #FIXME: considerer les axes 2D en input et output

        # Bilinear on pure numeric values
        if not N.ma.isMA(vari):
            xi = N.arange(vari.shape[-1])
            yi = N.arange(vari.shape[-2])
            xo = N.linspace(0., xi[-1], varo.shape[-1])
            yo = N.linspace(0., yi[-1], varo.shape[-2])
            xxo, yyo = N.meshgrid(xo, yo)
            varo[:] = _mbilin2d_(vari, xi, yi, xxo.flat, yyo.flat,
                1.e20, smoothcoast, nogeo=1-int(geo)).reshape(varo.shape)

        # Masked data
        else:

            varo.setMissing(1.e20)
            # Grids
            if cdms.isVariable(vari):
                xi = vari.getAxis(-1)
                yi = vari.getAxis(-2)
                xo = refine(xi, factor)
                yo = refine(yi, factor)
            else:
                xi = N.arange(vari.shape[-1])
                yi = N.arange(vari.shape[-2])
                xo = N.linspace(0., xi[-1], varo.shape[-1])
                yo = N.linspace(0., yi[-1], varo.shape[-2])
            xxo, yyo = N.meshgrid(xo[:], yo[:])

            # Interpolation
            varo[:] = _mbilin2d_(vari.filled(1.e20),xi[:], yi[:], xxo.flat, yyo.flat,
                1.e20, smoothcoast,len(xi),len(yi),N.size(xxo),geo).reshape(varo.shape)

            # Masking
            if vari.mask is not MV.nomask:
                fmaski = vari.mask.astype('f')
                fmasko = _mbilin2d_(fmaski,xi[:], yi[:], xxo.flat, yyo.flat,
                    1.e20, False, len(xi),len(yi),N.size(xxo),geo).reshape(varo.shape)
                varo[:] = MV.masked_where(N.greater(fmasko, .5), varo, copy=0)
            varo[:] = MV.masked_values(varo, 1.e20, copy=0)

            # Axes
            if not noaxes and cdms.isVariable(varo):
                for i in -2, -1:
#                   varo.setAxis(i, refine(vari.getAxis(i), factor))
                    varo.setAxis(i, (yo, xo)[i])
                    varo.setGrid(create_grid(xo,yo,fmasko))

    return varo
_cellave_methods = ['conservative', 'remap', 'cellave', 'conserv']
_cdat_methods = _cellave_methods+['bilinear']
_regrid2d_methods = ['nearest', 'mixt', 'interp', 'bining']+_griddata_methods+_cdat_methods

def regrid2d_method_name(method, raiseerr=True):
    """Check the method name and return its generic name"""
    if method is None or method.lower()=='auto': return 'auto'
    method = method.lower()
    if method.startswith('cons'): return 'conserv'
    if method.startswith('cell') or method.startswith('remap'): return 'cellave'
    if 'lin' in method or method=='interp': return 'bilinear'
    if method.startswith('bin'): return 'bining'
    if method.startswith('nat'): return 'nat'
    if method.startswith('near'): return 'nearest'
    if method.startswith('carg') or method.startswith('krig'): return 'carg'
    if method.startswith('dst'): return "dstwgt"
    if method.startswith('pat'): return "patch"
    if raiseerr:
        raise VACUMMError('Invalid regrid2d method. Please use for example one of these: '+
            ', '.join(_regrid2d_methods))
    else:
        return method

#: Available tools for each method
regrid2_tools = OrderedDict([
    ('cellave', dict(r2r=['regrid2', 'esmf'], c2c=['esmf'])),
    ('conserv', dict(r2r=['regrid2', 'esmf'], c2c=['esmf'])),
    ('nearest', dict(c2c=['vacumm'])),
    ('bilinear', dict(r2r=['vacumm','esmf','libcf'], r2c=['vacumm','esmf','libcf'],
        c2c=['esmf', 'libcf', 'vacumm'])),
    ('patch', dict(c2c=['esmf'])),
    ('dstwgt', dict(r2r=['vacumm'], c2c=['vacumm'])),
    ('carg', dict(r2r=['vacumm'])),
    ('nat', dict(c2r=['vacumm'])),
    ('bining', dict(c2r=['vacumm'])),
])

def regrid2d_tool_name(tool, raiseerr=True):
    """Check the tool name and return its generic name"""
    if tool is None or tool.lower()=='auto': return 'auto'
    tool = tool.lower()
    if tool.startswith('es'): return 'esmf'
    if tool.startswith('lib') or tool.startswith('g'): return 'libcf'
    if tool.startswith('cd') or tool.startswith('reg'): return 'regrid2'
    if tool.startswith('v'): return 'vacumm'
    if raiseerr:
        raise VACUMMError('Invalid regrid2d tool. Please use for example one of these: '+
            'regrid2, esmf, libcf, vacumm, auto')
    else:
        return method


def regrid2dold(vari, ggo, method='auto', mask_thres=.5, ext=False,
    bilinear_masking='dstwgt', ext_masking='poly', cdr=None, getcdr=False, usecdr=None, useoldcdr=True,
    check_mask=True, clipminmax=False, geo=None, **kwargs):
    """Regrid a variable from a regular grid to another

    If the input or output grid is curvilinear and ``method`` is set to
    ``"linear"``, ``"cellave"`` or ``"conserv"``, the :class:`CDATRegridder` is used.

    :Params:

        - **vari**: Variable cdms on regular grid
        - **ggo**: Tuple of (x,y) or a cdms grid or a cdms variable with a grid
        - **method**, optional: One of:

            - ``"auto"``: method guessed according to resolution of input and output grid (see :func:`regrid_method`)
            - ``"nearest"``: nearest neighbour
            - ``"linear"`` or ``"bilinear"``: bilinear interpolation (low res. to high res.)
            - ``"dstwgt"`` : distance weighting (low res. to high res.)
            - ``"cellave"`` : weighted regridding based on areas of cells (high res. to low res.)
            - ``"conserv"`` : same as ``cell`` but conservative (high res. to low res.)
            - ``"bining"`` : simple averaging using bining (very high res. to low res.)
            - ``"nat"`` : Natgrid interpolation (low res. to high res.) (see :class:`GridData` for more info)
            - ``"carg"`` : Interpolation with minicargen(low res. to high res.) (see :func:`cargen` for more info)

        - **cdr**, optional: :class:`CDATRegridder` instance.
        - **getcdr**, optional: Also return the computed :class:`CDATRegridder` instance.
        - **usecdr**, optional: Force the use or not of a :class:`CDATRegridder` instance,
          even for rectangular grids.
        - **useoldcdr**, optional: Force the use the old CDAT regridder (rectangular grids only).
        - **ext**, optional: Perform extrapolation when possible
        - **bilinear_masking**: the way to handle interpolation near masked values

            - ``"nearest"``: brut masking using nearest neighbor
            - ``"dstwgt"`` : distance weight data are used where interpolated mask is lower ``mask_thres``

        - **mask_thres**, optional: Threshold for masking points for some methods (~ land fraction) for
          **rectangular grids only**:

            - ``method="bilinear"`` and ``bilinear_masking="dstwght"``
            - ``method="cellave"`` or ``method="bining"``

        - **ext_masking**, optional: Manual masking method when ``ext=False`` (when needed)
            with methods ["carg",] (see :func:`~vacumm.misc.grid.masking.grid_envelop_mask`)
            if input grid is not rectangular

            - ``"poly"``: use the polygon defined by the input grid envelopp and check if output points are inside
            - ``"nearest"``: use hack with nearest 2d interpolation

        - Other keywords are passed to special interpolation functions depending on method and choices :

            - :func:`cargen` when "nat" or "carg" method is used
            - **mask_thres**, optional: Time steps when interpolated mask is greater than
              this value are masked.

    :Examples:

        >>> regrid2d(var, (lon, lat), method='bilinear', bilinear_masking='nearest')
        >>> regrid2d(var, grid, method='cellave', mask_thres=.8)
        >>> regrid2d(var, grid, method='nat', hor=.2)
    """

    # Method
    method = regrid2d_method_name(method)
    if method == 'auto':
        method = regrid2d_method_name(regrid_method(vari, ggo))

    # Check mask, grids and method
    vari = MV.asarray(vari)
    if vari.getMissing() is None:
        vari.setMissing(1.e20)
    mv = vari.getMissing()
    maski = vari.mask
    if maski is MV.nomask:
        check_mask = False
    if check_mask:
        maski = maski.astype('f')

    ggi = curv2rect(get_grid(vari), mode='none')
    nyi, nxi = ggi.shape
    nzi = vari.size/(nxi*nyi)
    ggo = get_grid(ggo)
    ggor = curv2rect(ggo, mode='none')
    nyo, nxo = ggo.shape
    loni = get_axis(ggi, -1)
    lati = get_axis(ggi, -2)
    lono = get_axis(ggor, -1)
    lato = get_axis(ggor, -2)
    xo = lono[:]
    yo = lato[:]
    xi = loni[:]
    yi = lati[:]
    curvedi = xi.ndim==2
    curvedo = xo.ndim==2
    curved = curvedi or curvedo
    curvedio = curvedi and curvedo
    xxi, xxo, yyi, yyo = None, None, None, None
    if geo is None:
        geo = A.islon(loni) or A.islat(lati) or A.islon(lono) or A.islat(lato)

    # Bounds
    if method in _cellave_methods+['nearest']:
        funci = bounds2d if curvedi else bounds1d
        loni.setBounds(funci(loni))
        lati.setBounds(funci(lati))
        funco = bounds2d if curvedo else bounds1d
        lono.setBounds(funco(lono))
        lato.setBounds(funco(lato))

    # Some checks about methods
    if curved:
        if method not in _griddata_methods+['nearest', 'bining']+_cdat_methods:
            raise NotImplementedError, 'Method not allowed with curvilinear grids: '+method
        xxi, yyi = meshgrid(loni, lati)
    maskoext = False
    if method == 'krig': method = 'carg'
    if method == 'nearest' or \
        (check_mask and method in ['bilinear', 'mixt', 'dstwgt'] and not bilinear_masking=='dstwgt') or \
            (ext is False and method == 'carg') : # Needs nearest
        xxo, yyo = meshgrid(lono, lato)
        if not isgrid(ggi, curv=True) or method=='nearest':
            xxi, yyi = meshgrid(loni, lati)
        if not ext:
            from masking import grid_envelop_mask
            maskoext = N.resize(grid_envelop_mask(ggi, ggo, poly=ext_masking=="poly"), (nzi, nyo, nxo))
    elif method == 'bining' and curvedo:
            raise NotImplementedError, 'Method not allowed with curvilinear output grid: '+method

    # 3D variable?
    if vari.ndim != 3:
        vari3d = MV.resize(vari, (nzi, nyi, nxi))
        set_grid(vari3d, ggi)
#        vari3d.getAxis(0).designateLevel() # hack against cdat bug
        maski3d = N.resize(maski, vari3d.shape)
    else:
#        if vari.getOrder()[0] not in 'tz': # TODO: OPTIMIZE REGRID2
#            vari3d = vari.clone()
#            zaxis = vari3d.getAxis(0).clone()
#            zaxis.designateLevel()
#            vari3d.setAxis(0, zaxis)
        vari3d = vari
        maski3d = maski

    # Interpolations
    if method == 'nearest':

        # Interpolation
        varo3d = _regrid2d_nearest2d_(vari3d, xxi, yyi, xxo, yyo, mv, geo, maskoext)


    elif (not curved and method in ['mixt', 'dstwgt']) or \
        (method=='bilinear' and not curved and usecdr is not True):


        # Method
        wrapper = eval('_regrid2d_%s_'%method)

        # Interpolation
        #FIXME: wrapper pour dstwgt et mask... not sure
        varo3d = wrapper(vari3d, xi, yi, xo, yo, mv, geo, ext)

        # Masking
        if check_mask:

            if bilinear_masking == 'dstwgt':
                masko = wrapper(maski3d, xi, yi, xo, yo, 1., geo, ext)
                if method != 'dstwgt':
                    varo3d_dstwgt = _regrid2d_dstwgt_(vari3d, xi, yi, xo, yo, mv, geo, ext)
                    varo3d[:] = N.where((masko<mask_thres) & (masko>0.), varo3d_dstwgt, varo3d)
                    del varo3d_dstwgt
                varo3d[masko>=mask_thres] = mv
            else:
                masko = varo3d==mv
                varo3d_nearest = _regrid2d_nearest2d_(vari3d,xxi,yyi,xxo,yyo,mv,geo,maskoext)
                varo3d = N.where(masko!=0., varo3d_nearest, varo3d)
                del varo3d_nearest
            del masko



    elif method in _cdat_methods:

        # Old regridder ?
        if useoldcdr or useoldcdr is None and not curved and method in ['cellave', 'conserv']:

            tool = 'regrid2'

        # ESMF regridder
        else:
            esmfcons = curved #and method in _cellave_methods
            if esmfcons:
                check_mask = True
                vari3d = MV2.where(maski3d, 0., vari3d)
                set_grid(vari3d, ggi)
            tool = None

        # Regridder
        if cdr is None:
            cdr = CDATRegridder(vari3d, ggor, method=method, tool=None)

        # Regrid data
        varo3d = cdr(vari3d, check_mask=check_mask)
        del maski3d

#        # Masking
#        if check_mask:# and not method.startswith('cons'):
#            good = vari3d.clone()
#            good[:] = 1-maski3d ;  del maski3d
#            goodo = cdr(good).filled(0.) ; del good
##           varo3d_nearest = _regrid2d_nearest2d_(vari3d.filled(mv),xxi,yyi,xxo,yyo,geo,maskoext,mv)
##           varo3d[:] = MV2.where(masko!=0., varo3d_nearest, varo3d)
#            varo3d[:] = MV2.masked_where((goodo<0.9999), varo3d, copy=0) ; del goodo
##            varo3d[:] = MV2.masked_where((masko>mask_thres).astype('b') & (masko>0.).astype('b'), varo3d, copy=0)
##            del varo3d_nearest



    elif method.startswith('bin'):

        # Interpolation
        varo3d = _regrid2d_bining_(vari3d, xxi, yyi, xo, yo, mv)

        # Masking
        if check_mask:
            masko = _regrid2d_bining_(maski3d, xxi, yyi, xo, yo, 1.)
            varo3d[masko>=mask_thres] = mv
            del masko

    elif method in _griddata_methods:

        # Interpolation
        varo3d = _regrid2d_griddata_(vari3d, xi, yi, xo, yo, method, ext=ext, missing_value=mv, **kwargs)

        # Masking
        if check_mask:
            masko = _regrid2d_griddata_(maski3d, xi, yi, xo, yo, method, ext=ext, missing_value=1., **kwargs)
#           varo3d_nearest = _regrid2d_nearest2d_(vari3d.filled(mv),xxi,yyi,xxo,yyo,mv,geo,False)
#           varo3d = N.where(masko!=0., varo3d_nearest, varo3d)
            varo3d[masko>=mask_thres] = mv
            if method=='carg' and not ext:
                varo3d[maskoext] = mv
#           del varo3d_nearest

    else:
        raise RuntimeError, "Well, what's this funckin' method you bastard huh: %s?"%method


    # Back to rights dims
    if vari.ndim != 3:
        varo = MV.resize(varo3d, vari.shape[:-2]+(nyo, nxo))
        del vari3d, varo3d
    else:
        varo = varo3d

    # MV2 variable
    varo[N.isnan(varo)] = mv
    varo = MV2.masked_values(varo, mv, copy=0)
#   if not isgrid(ggo, curv=True):
    set_grid(varo, ggo)
    if vari.ndim>2:
        for iaxis in range(vari.ndim-2):
            varo.setAxis(iaxis, vari.getAxis(iaxis))
    cp_atts(vari, varo, id=True)

    if getcdr: return varo, cdr
    return varo




def regrid2d(vari, ggo, method='auto', tool=None, rgdr=None, getrgdr=False,
    erode=False,  **kwargs):
    """Regrid a variable from a regular grid to another

    If the input or output grid is curvilinear and ``method`` is set to
    ``"linear"``, ``"cellave"`` or ``"conserv"``, the :class:`CDATRegridder` is used.

    :Params:

        - **vari**: Variable cdms on regular grid
        - **ggo**: Tuple of (x,y) or a cdms grid or a cdms variable with a grid
        - **method**, optional: One of:

            - ``"auto"``: method guessed between ``linear`` and ``cellave``
              according to resolution of input and output grid (see :func:`regrid_method`)
            - ``"nearest"``: nearest neighbour
            - ``"linear"`` or ``"bilinear"``: bilinear interpolation (low res. to high res.)
            - ``"dstwgt"`` : distance weighting between the four nearest grid points
              (low res. to high res.)
            - ``"patch"`` : patch recovery interpolation (low res. to high res.)
            - ``"cellave"`` : weighted regridding based on areas of cells (high res. to low res.)
            - ``"conserv"`` : same as ``cell`` but conservative (high res. to low res.)

        - **tool**, optional: Regridder. One of:

            - ``"auto"``: tool guessed depending on the method, the first available tool
              and the grids (rectangular or curvilinear).
            - ``"vacumm"``: Internal routines.
            - ``"esmf"`` and ``"libcf"``: Regridders provided by UVCDAT.
            - ``"regrid2"``: Old regridder provided by CDAT (rectangular only).

        - **rgdr**, optional: An already set up regridder instance to speed up regridding:
          :class:`CDATRegridder` instance for ``regrid2``, ``esmf`` and ``libcf`` tools,
          else a :class:`CurvedInterpolator` instance for ``vacumm`` tool with
          interpolation on curvilinear grids.
        - **getrgdr**, optional: Also return the regridder instance if it applies, or None.
        - Other keywords are passed to special interpolation functions depending on method and choices :

            - :func:`cargen` when "nat" or "carg" method is used

    :Tools/methods:

        Overview table of method availability for each tool.
        ``RECT`` means that it only works with rectangular grids.

        +----------+--------+----------+------+-------+
        | Met/Tool | Vacumm | regrid2d | ESMF | Libcf |
        +==========+========+==========+======+=======+
        | nearest  |   OK   |          |      |       |
        +----------+--------+----------+------+-------+
        | bilinear |   OK   |          |  OK  |  OK   |
        +----------+--------+----------+------+-------+
        |  dstwgt  |   OK   |          |      |       |
        +----------+--------+----------+------+-------+
        |  patch   |        |          |  OK  |       |
        +----------+--------+----------+------+-------+
        | cellave  |        |   RECT   |  OK  |       |
        +----------+--------+----------+------+-------+
        | conserv  |        |   RECT   |  OK  |       |
        +----------+--------+----------+------+-------+

    :Examples:

        >>> regrid2d(var, (lon, lat), method='linear')
        >>> regrid2d(var, grid, method='cellave')
    """

    # Check grids
    vari = MV.asarray(vari)
    mv = vari.getMissing()
    if mv is None or N.isnan(mv):
        vari.setMissing(1.e20)
    mv = vari.getMissing()
    ggi = curv2rect(get_grid(vari), mode='none')
    nyi, nxi = ggi.shape
    nzi = vari.size/(nxi*nyi)
    ggo = get_grid(ggo)
    ggor = curv2rect(ggo, mode='none')
    nyo, nxo = ggo.shape
    loni = get_axis(ggi, -1)
    lati = get_axis(ggi, -2)
    lono = get_axis(ggor, -1)
    lato = get_axis(ggor, -2)
    xo = lono[:]
    yo = lato[:]
    xi = loni[:]
    yi = lati[:]
    curvedi = xi.ndim==2
    curvedo = xo.ndim==2
    curved = curvedi or curvedo
    curvedio = curvedi and curvedo
    xxi, xxo, yyi, yyo = None, None, None, None
#    if geo is None:
#        geo = A.islon(loni) or A.islat(lati) or A.islon(lono) or A.islat(lato)
    geo = True

    # Method
    method = regrid2d_method_name(method)
    if method == 'auto':
        method = regrid2d_method_name(regrid_method(vari, ggo))
    valid_methods = ['nearest', 'bilinear', 'cellave', 'conserv', 'dstwgt', 'patch']
    if method not in valid_methods:
        raise VACUMMError('Wrong regridding method: %s. Please choose one of :'%method+
        ', '.join(valid_methods))

    # Tool verification
    tool = regrid2d_tool_name(tool)
    tools = regrid2_tools[method]
    for ttype in 'r2r', 'r2c', 'c2r', 'c2c':
        if ttype not in tools: continue
        if  (ttype=='r2r' and not curved) or \
                (ttype=='c2r' and curvedi and not curvedo) or \
                (ttype=='r2c' and curvedo and not curvedi) or \
                ttype=='c2c': # Fall back | and curved):
            if tool=='auto':
                tool = tools[ttype][0] # First available
            elif tool not in tools[ttype]:
                raise VACUMMError('Tool %s not available for these grids and this method: '%
                    (tool, method))
            break
    else:
        raise VACUMMError('No suitable tool found for these grids and this method: '+method)

    # Bounds
    if method in _cellave_methods+['nearest']:
        funci = bounds2d if curvedi else bounds1d
        loni.setBounds(funci(loni))
        lati.setBounds(funci(lati))
        funco = bounds2d if curvedo else bounds1d
        lono.setBounds(funco(lono))
        lato.setBounds(funco(lato))

    # 2D arrays
    if curvedi or method=='nearest':
        xxi, yyi = meshgrid(loni, lati)
    if curvedi or curvedo or method=='nearest':
        xxo, yyo = meshgrid(lono, lato)

    # 3D variable?
    if vari.ndim != 3:
        vari3d = MV.resize(vari, (nzi, nyi, nxi))
        set_grid(vari3d, ggi)
    else:
        vari3d = vari


    # Interpolations !
    rgdr = None
    if tool in ['regrid2', 'esmf', 'libcf']:


#        # ESMF regridder
#        if tool=='esmf' and method in ["cellave", "conserv"]:
#            check_mask = True
#            vari3d = MV2.where(maski3d, 0., vari3d)
#            set_grid(vari3d, ggi)
#        else:
#            check_mask = False

        # Regridder
        if rgdr is None:
            rgdr = CDATRegridder(vari3d, ggor, method=method, tool=tool)

        # Regrid data
        varo3d = rgdr(vari3d)#, check_mask=check_mask)
#        del maski3d

    elif method == 'nearest':

        varo3d = _regrid2d_nearest2d_(vari3d, xxi, yyi, xxo, yyo, mv, geo, False)

    elif method in ['bilinear', 'dstwgt']:

        if curvedi:
            wrapper = eval('_regrid2d_%s_c2c_'%method)
            varo3d, rgdr = wrapper(vari3d, xxi, yyi, xxo, yyo, mv, rgdr)

        elif curvedo:
            wrapper = eval('_regrid2d_%s_r2c_'%method)
            varo3d = wrapper(vari3d, xi, yi, xxo, yyo, mv)

        else:
            wrapper = eval('_regrid2d_%s_r2r_'%method)
            varo3d = wrapper(vari3d, xi, yi, xo, yo, mv, geo)


    elif method == ['mixt']:


        # Method
        wrapper = eval('_regrid2d_%s_'%method)

        # Interpolation
        varo3d = wrapper(vari3d, xi, yi, xo, yo, mv, geo, ext)

    else:
        raise RuntimeError, "Well, what's this funckin' method you bastard huh: %s?"%method


    # Back to rights dims
    if vari.ndim != 3:
        varo = MV.resize(varo3d, vari.shape[:-2]+(nyo, nxo))
        del vari3d, varo3d
    else:
        varo = varo3d

    # MV2 variable
    varo[N.isnan(varo)] = mv
    varo = MV2.masked_values(varo, mv, copy=0)
#   if not isgrid(ggo, curv=True):
    set_grid(varo, ggo)
    if vari.ndim>2:
        for iaxis in range(vari.ndim-2):
            varo.setAxis(iaxis, vari.getAxis(iaxis))
    cp_atts(vari, varo, id=True)

    if getrgdr: return varo, rgdr
    return varo



regrid2dnew = regrid2d

def _regrid2d_nearest2d_(vari3d, xxi, yyi, xxo, yyo, mv, geo, maskoext):
    """Wrapper to fortran _nearest2d_"""
    nb = min(xxi.shape)/50
    if nb == 1: nb = 0
#    nb = -1
    if N.ma.isMA(vari3d): vari3d = vari3d.filled(mv)
    varo3d = N.asarray(_nearest2d_(vari3d, xxi, yyi, xxo, yyo, nb, not geo), order='C')
    if maskoext is not False:
        varo3d[maskoext] = mv
    return varo3d

def _regrid2d_bilinear_r2r_(vari3d, xi, yi, xo, yo, mv, geo, ext=0):
    """Wrapper to fortran _bilin_"""
    if N.ma.isMA(vari3d): vari3d = vari3d.filled(mv)
    varo3d =  N.asarray(_bilin_(vari3d, xi, yi, xo, yo, mv, not geo), order='C')
    return varo3d
_regrid2d_bilinear_ = _regrid2d_bilinear_r2r_

def _regrid2d_bilinear_r2c_(vari3d, xi, yi, xxo, yyo, mv):
    """Wrapper to fortran _bilin2dto1d_"""
    if N.ma.isMA(vari3d): vari3d = vari3d.filled(mv)
    varo3d =  N.ascontiguousarray(_bilin2dto1d_(xi, yi, vari3d, xxo.ravel(), yyo.ravel(), mv))
    varo3d.shape = varo3d[:-1]+xxo.shape
    return varo3d

def _regrid2d_bilinear_c2c_(vari3d, xxi, yyi, xxo, yyo, mv, rgdr=None):
    """Wrapper to fortran _bilin2dto1dc_reduc_ through :class:`CurvedInterpolator`"""
    if rgdr is None:
        rgdr = CurvedInterpolator((xxi, yyi), (xxo, yyo), g2g=True)
    varo3d =  N.ascontiguousarray(rgdr(vari3d, method='bilinear'))
    return varo3d, rgdr

def _regrid2d_dstwgt_r2r_(vari3d, xi, yi, xo, yo, mv, geo, ext=0):
    """Wrapper to fortran _bilin_"""
    if N.ma.isMA(vari3d): vari3d = vari3d.filled(mv)
    varo3d =  N.asarray(_dstwgt_(vari3d, xi, yi, xo, yo, mv, not geo), order='C')
    return varo3d
_regrid2d_dstwgt_ = _regrid2d_dstwgt_r2r_

def _regrid2d_dstwgt_r2c_(vari3d, xi, yi, xxo, yyo, mv):
    """Wrapper to fortran _bilin2dto1d_"""
    if N.ma.isMA(vari3d): vari3d = vari3d.filled(mv)
    varo3d =  N.ascontiguousarray(_dstwgt2dto1d_(xi, yi, vari3d, xxo.ravel(), yyo.ravel(), mv))
    varo3d.shape = varo3d[:-1]+xxo.shape
    return varo3d

def _regrid2d_dstwgt_c2c_(vari3d, xxi, yyi, xxo, yyo, mv, rgdr=None):
    """Wrapper to fortran _bilin2dto1dc_reduc_ through :class:`CurvedInterpolator`"""
    if rgdr is None:
        rgdr = CurvedInterpolator((xxi, yyi), (xxo, yyo), g2g=True)
    varo3d =  N.ascontiguousarray(rgdr(vari3d, method='dstwgt'))
    return varo3d, rgdr


def _regrid2d_mixt_(vari3d, xi, yi, xo, yo, mv, geo, ext=0):
    """Wrapper to fortran _mixt2dx_"""
    if N.ma.isMA(vari3d): vari3d = vari3d.filled(mv)
    varo3d = N.asarray(_mixt2dx_(vari3d, xi, yi,  xo, yo, mv, ext).T, order='C')
    return varo3d

def _regrid2d_bining_(vari3d, xxi, yyi, xo, yo, mv):
    """Regridding using rebining"""
    xbins = meshcells(xo)
    ybins = meshcells(yo)
    varo3d = N.zeros(vari3d.shape[:-2]+(len(yo), len(xo)))
    for i, vari2d in enumerate(vari3d.filled(mv)):
        if hasattr(vari2d, 'mask') and vari2d.mask.any():
            good = ~vari2d.mask
            xc = xxi[good]
            yc = yyi[good]
            vartmp = vari2d.compressed()
        else:
            xc = xxi.ravel()
            yc = yyi.ravel()
            vartmp = vari2d.ravel()
        pos = N.asarray([yc, xc]).T
        del xc, yc
        count, edges = N.histogramdd(pos, bins=(ybins, xbins)) ; del edges
        varo3d[i], edges = N.histogramdd(pos, bins=(ybins, xbins), weights=vartmp)
        del pos, edges
        good = count!=0
        varo3d[i][good] /= count[good]
        varo3d[i][~good] = mv
        del good, count
    return varo3d


def _regrid2d_griddata_(vari3d, xi, yi, xo, yo, method, **kwargs):
    """Wrapper to griddata"""
    # Good shapes
    xxi, yyi = meshgrid(xi, yi)
    nzi, nyi, nxi = vari3d.shape
    vari2d = vari3d.reshape(nzi, nxi*nyi)
    # Restrict zone
    dxo = xo.ptp()/10.
    dyo = yo.ptp()/10.
    #if hasattr(vari2d, 'mask') and vari2d.mask is not MV2.nomask and vari2d.mask.any():
        #good = vari2d.mask
    #else:
        #good = N.zeros(vari2d.shape, '?')
    good = (xxi.ravel()>=(xo.min()-dxo)) & (xxi.ravel()<=(xo.max()+dxo)) & \
        (yyi.ravel()>=(yo.min()-dyo)) & (yyi.ravel()<=(yo.max()+dyo))
    varitmp = vari2d.compress(good.ravel(), axis=1)
    xxitmp = xxi.compress(good.ravel())
    yyitmp = yyi.compress(good.ravel())
    del vari2d
    varo2d =  griddata(xxitmp, yyitmp, varitmp, (N.asarray(xo), N.asarray(yo)), method=method, **kwargs)
    del xxitmp, yyitmp, varitmp
    varo3d = varo2d.reshape(vari3d.shape[:-2]+varo2d.shape[-2:])
    del varo2d
    return varo3d

def _regrid2d_natgridlist_(vari3d, xi, yi, xxo, yyo, ext, **kwargs):
    #FIXME: _regrid2d_natgridlist_ : missing values
    """Wrapper to Natgrid.nat with list output"""
    # Good shapes
    xxi, yyi = meshgrid(xi, yi)
    nzi, ny, nx = vari3d.shape
    nyo, nxo = xxo.shape
    vari2d = vari3d.rehape(nzi, nxi*nyi)
    # Restrict zone
    dxo = xo.ptp()/10.
    dyo = xo.ptp()/10.
    #TODO: add dxo,dxy natgridlist in regrid2d
    border = zip(xxo[0], yyo[0])+zip(xxo[1:-1, -1], yyo[1:-1, -1])+\
        zip(xxo[-1], yyo[-1])+zip(xxo[1:-1, 0], yyo[1:-1, 0])
    poly = Polygon(N.array(border))
    good = N.ones(xo.shape, '?')
    if hasattr(vari3d, 'mask') and vari3d.mask is not MV2.nomask and ~vari3d.mask.all():
        good &= ~vari2d.mask
    for i in xrange(nxi):
        for j in xrange(nyi):
            good[i, j] = Point((xxi[j, i], yyi[j, i])).within(poly)
    varitmp = vari2d.compress(good, axis=-1)
    xxitmp = xxi.compress(good).astype('d')
    yyitmp = yyi.compress(good).astype('d')
    del vari2d
    r = Natgrid(xxitmp, yyitmp, xo.astype('d').ravel(), yo.astype('d').ravel(), listOutput = 'yes')
    r.ext = int(ext)
    r.nul = -10.
    r.dup = 0
    varo2d = N.zeros((nzi, nyo, nxo), dtype=vari.dtype)
    for iz in xrange(nzi):
        varo2d[iz] = r.rgrd(varitmp[iz])
    del xxitmp, yyitmp, varitmp
    varo3d = varo2d.reshape(vari3d.shape[:-2]+(nyo, nxo))
    del varo2d
    return varo3d

def cellave2d(vari, ggo, **kwargs):
    """Shortcut to :func:`regrid2d` call with ``method="cellave"``"""
    return regrid2d(vari, ggo, method="cellave", **kwargs)
remap2d = cellave2d

def interp2d(vari, ggo, method="interp", **kwargs):
    """Shortcut to :func:`regrid2d` call with ``method="interp"`` by default"""
    return regrid2d(vari, ggo, method=method, **kwargs)

def _isaxis_(var):
    if A.isaxis(var): return True
    if isinstance(var, tuple): return False
    if isgrid(var): return False
    return True

def regrid_method(gridi, grido, ratio=.55, iaxi=-1, iaxo=-1):
    """Guess the best regridding method for passing from ``gridi`` to ``grido``

    If ``resolution(gridi) <= ratio*resolution(grido)``, ``method="cellave"``
    else  ``method="interp"``
    The :func:`~vacumm.misc.grid.misc.resol` function is used to estimate
    the resolution.
    For grids, a resolution common to X and Y axes is estimated using the following
    sequence::

        xres, yres = resol(grid)
        res = (xres**2+yres**2)**.5

    :Params:

        - **gridi**: Grid, tuple of axes or single axis or array.
        - **grido**: Grid, tuple of axes or single axis or array.
        - **ratio**, optional: Limit ratio of output grid resolution to input grid resolution.
        - **iaxi/iaxo**, optional: Dimension on which to compute the resolution with
          :func:`~vacumm.misc.grid.misc.resol` when ``gridi`` and ``grido`` are single axes
          with several dimensions.

    ..note::

        The resolution of the grids is checked in their attributes "_xres" and "_yres"
        before before trying to compute them.

    :Returns: ``'cellave'`` OR ``'linear'``
    """
    if _isaxis_(gridi) and _isaxis_(grido): # Axes

        resi = resol(gridi, axis=iaxi)
        reso = resol(grido, axis=iaxo)

    else: # Grids

        xresi, yresi = resol(gridi)
        resi = (xresi**2+yresi**2)**.5
        xreso, yreso = resol(grido)
        reso = (xreso**2+yreso**2)**.5

    # Check the ratio
    return 'cellave' if resi <= ratio*reso else 'linear'


class GriddedMerger(object):
    """Merge several gridded variables onto a grid

    :Params:

        - **grid**: Output grid
        - **id**, optional: Output id
        - **long_name**, optional: Output long name
        - **units**, optional: Output units
        - Other keywords are set a output attributes

    :Example:

        >>> gm = GriddedMerger(mygrid)
        >>> gm += var1
        >>> gm.append(var2)
        >>> gm += var3
        >>> gm -= var3
        >>> gm.insert(0, var3)
        >>> print len(gm)
        3
        >>> print gm
        ....
        >>> myvar = gm.merge(res_ratio=.4, pad=3)
    """
    def __init__(self, grid, id=None, long_name=None, units=None, **kwargs):
        self._vars = []
        self.long_name = long_name
        self.units = units
        self.id = id
        self._var = None
        self.set_grid(grid)
        for att, val in kwargs.items():
            setattr(self, '_'+att, val)

    def set_grid(self, grid):
        """Set the grid for merging"""
        self._grid = get_grid(grid)
        self._xres, self._yres = resol(self._grid)
        del self._var
        self._var = None
    def get_grid(self):
        """Get the grid for merging"""
        return self._grid
    def get_lon(self):
        """Get the longitudes of the grid"""
        return self._grid.getAxis(1)
    def get_lat(self):
        """Get the loatitudes of the grid"""
        return self._grid.getAxis(0)


    def __len__(self):
        return len(self._vars)
    def _load_(self, var, method):
        # Check grid
        grid = var.getGrid()
        assert grid is not None, 'Your variable must have a grid'
        # Check dims
        if len(self):
            firstdims = tuple([len(a) for a in var.getAxisList()[:-2]])
            assert self._firstdims == firstdims, "Incompatible number of dimension: %s. It should be: %s."%\
                (self._firstdims, firstdims)
        else:
            self._firstaxes = var.getAxisList()[:-2]
            self._firstdims = tuple([len(a) for a in self._firstaxes])
        # Method and resolution
        var._gm_method = method
        grid._xres, grid._yres = resol(grid)
        # Bounds
        for axis in var.getAxisList()[-2:]:
            axis.setBounds(bounds1d(axis))
        return var
    def add(self, *args, **kwargs):
        """Alias for :meth:`append`"""
        self.append(*args, **kwargs)
    def append(self, var, method='auto', **kwargs):
        """Append a bathymetry to the top of the merger"""
        self._vars.append(self._load_(var, method, **kwargs))
    def __iadd__(self, var):
        self.add(var)
        return self
    def remove(self, var):
        if isinstance(var, int):
            del self._vars[i]
        assert var in self._vars, 'Variable not in merger'
        self._vars.remove(var)
    def __isub__(self, var):
        self.remove(var)
        return self
    def __getitem__(self, idx):
        return self._vars[idx]
    def __setitem__(self, idx, var):
        self._vars[idx] = self._load_(var)
    def __delitem__(self, idx):
        del self._vars[idx]
    def insert(self, idx, var):
        self._vars.insert(idx, self._load_(var))

    def __str__(self):
        if not len(self): return 'No variables in the merger'
        ret = ''
        for i, var in enumerate(self._vars):
            print '\n%i) %s_n'%(i, var.id)
            for att in 'long_name', 'units', '_gm_method':
                if hasattr(var, att):
                    res += '  %s = %s'%(att, getattr(var, att))
            for att in '_gm_xres', '_gm_yres':
                if hasattr(var.getGrid(), att):
                    res += '  %s = %s'%(att, getattr(var.getGrid(), att))
        ret += '\n--\nTotal: %i variable(s)'%len(self)
        return ret


    def merge(self, res_ratio=.5, pad=5, **kwargs):
        """Merge all the variables on to a grid


        - **grid**: Out put grid or axes.
        - *res_ratio*: Resolution ratio for choosing between cell
          averaging and bilinear interpolation (see: :func:`regrid_method`).
        - *regrid_<kwparam>*: *<kwparam>* is passed to :func:`regrid2d` for interpolation.
        """
        assert len(self), 'You must add at least one variable to the merger'

        # Some useful info if the output grid
        lono = self.get_lon()
        lato = self.get_lat()
        xbo = meshcells(lono)
        ybo = meshcells(lato)
        xbmino = xbo.min()
        xbmaxo = xbo.max()
        ybmino = ybo.min()
        ybmaxo = ybo.max()
        grid = self.get_grid()
        nyo, nxo = grid.shape

        # Inits
        outmask2d = N.ones(grid.shape, '?')
        varo = MV2.zeros(self._firstdims+grid.shape)
        set_grid(varo, grid)
        total_cover = N.zeros(grid.shape)
        cover = N.zeros(grid.shape)
        weight = N.zeros(grid.shape)
        wnd = N.zeros(varo.shape)
        vo = MV2.zeros(varo.shape)
        kwregrid = kwfilter(kwargs, 'regrid')
        for var in self._vars[::-1]:

            # Axes
            loni = var.getLongitude() ; xi = loni.getValue()
            lati = var.getLatitude() ; yi = lati.getValue()
            nyi, nxi = var.getGrid().shape

            # Guess the interpolation method
            if var._gm_method != 'auto':
                method = var._gm_method
            else:
                method = regrid_method(var.getGrid(), grid, ratio=res_ratio)

            # Guess the grid bounds according to the method
            if method == 'interp':
                xmini = xi.min()
                xmaxi = xi.max()
                ymini = yi.min()
                ymaxi = yi.max()
            else:
                xmini = loni.getBounds().min()
                xmaxi = loni.getBounds().max()
                ymini = lati.getBounds().min()
                ymaxi = lati.getBounds().max()

            # Cell covering
            cover[:] = 0.
            # - limits
            if xmini>=xbmaxo or xmaxi<=xbmino or \
                ymini>=ybmaxo or ymaxi<=ybmino: continue
            ix0 = N.searchsorted(xbo, xmini, 'right')-1
            ix1 = N.searchsorted(xbo, xmaxi, 'left')-1
            iy0 = N.searchsorted(ybo, ymini, 'right')-1
            iy1 = N.searchsorted(ybo, ymaxi, 'left')-1
            # - partial+full cells
            cslice = [slice(max(iy0, 0), iy1+1), slice(max(ix0, 0), ix1+1)]
            cover[cslice[0], cslice[1]] = 1.
            # - partial cells
            if ix0>=0:
                cover[:, ix0] *= (xbo[ix0+1]-xmini)/(xbo[ix0+1]-xbo[ix0])
            if ix1<nxo:
                cover[:, ix1] *= (xmaxi-xbo[ix1])/(xbo[ix1+1]-xbo[ix1])
            if iy0>=0:
                cover[iy0, :] *= (ybo[iy0+1]-ymini)/(ybo[iy0+1]-ybo[iy0])
            if iy1<nyo:
                cover[iy1, :] *= (ymaxi-ybo[iy1])/(ybo[iy1+1]-ybo[iy1])
            # - borders
            xpad = N.ones(nyo)
            ypad = N.ones(nxo)
            ipadmax = min(ix1, nxo-1)-max(ix0, 0)
            for ipad in xrange(min(pad, ipadmax+1)):
                xcov = (ipad+.5)*xpad/(ipadmax+.5)
                cover[:, ipad] *= xcov
                cover[:, nxo-ipad-1] *= xcov
                del xcov
            jpadmax = min(iy1, nyo-1)-max(iy0, 0)
            for jpad in xrange(min(pad, jpadmax+1)):
                ycov = (jpad+.5)*ypad/(jpadmax+.5)
                cover[jpad] *= ycov
                cover[nyo-jpad-1] *= ycov
                del ycov
            del xpad, ypad

            # Regridding
            vo[:] = regrid2d(var, grid, method, **kwregrid)

            # Contribution to varo
            cover[total_cover==1] = 0.
            total_cover += cover
            wnd[:] = N.resize(cover, varo.shape)
            vo[:] *= wnd
            varo[..., cslice[0], cslice[1]] += vo[cslice[0], cslice[1]]
        varo[:] = MV2.masked_where(total_cover==0, varo/total_cover)
        del self._var, cover, total_cover, vo, wnd
        self._var = varo
        return varo

    def plot(self, **kwargs):
        """Merge and plot"""
        from ...plot import map2
        if self._var is None:
            self.merge()
        return map2(self._var)


def _shiftslicenames_(shift):
     if shift>0:
        target = firsts = 'firsts'
        out = last = 'last'
        oppos = 'first'
        neigh1 = 'last'
        neigh2 = 'lastm1'
        firsts = 'firsts'
        lasts = 'lasts'
     else:
        target = firsts = 'lasts'
        out = last = 'first'
        oppos = 'last'
        neigh1 = 'first'
        neigh2 = 'firstp1'
        firsts = 'lasts'
        lasts = 'firsts'
     return dict(target=target, out=out, oppos=oppos, last=last,
        neigh1=neigh1, neigh2=neigh2, firsts=firsts, lasts=lasts)


def shift1d(var, shift=0, bmode=None, axis=-1, copy=True, shiftaxis=True, **kwargs):
    """Interpolate data on an axis shifted by an half cell size

    :Params:

        - **var**: array.
        - **shift**, optional: Shift to operate.

            - ``0``: No shift.
            - ``<0``: Shilt toward bottom of the axis.
            - ``>0``: Shilt toward top of the axis.

        - **bmode**, optional: Boundary mode.

            - ``None``: ``"extrap"`` if axis else ``"same"``
            - ``"linear"``: Linear extrapolation.
            - ``"same"``: Constant extrapolation.
            - ``"masked"``: Mask outside data.
            - ``"cyclic"``: Cyclic extrapolation.

        - **axis**, optional: Axis on which to operate.
        - **copy**, optional: Always input copy data.
    """
    # Input data
    if copy:
        if hasattr(var, 'clone'):
            varo = var.clone()
        else:
            varo = var.copy()
    else:
        varo = var
    if shift==0: return varo
    if shift<0:
        shift = -0.5
    else:
        shift = 0.5
    bmode = kwargs.get('mode', bmode)
    if bmode is None:
        bmode = "linear" if A.isaxis(var) else "same"
    elif bmode=="nearest":
        bmode = "same"
    elif bmode=='extrap':
        bmode = 'linear'
    ax = len(var.shape)==1 and A.isaxis(var)
    if ax:
        varf = varo.getValue()
    else:
        varf = varo
    xo = None
    if cdms2.isVariable(varo):
        refvar = varo.asma().copy()
        if shiftaxis:
            xo = shift1d(var.getAxis(axis), shift, bmode='extrap')
    else:
        refvar = varf.copy()

    # Get slice specs
    ss = get_axis_slices(var, axis)
    sn = _shiftslicenames_(shift)

    # Inner data
    varf[ss[sn['target']]] = (refvar[ss['firsts']]+refvar[ss['lasts']])*0.5

    # Boundary data
    if bmode=='masked':
        varf[ss[sn['out']]] = N.ma.masked
    elif bmode=='cyclic':
        varf[ss[sn['out']]] = var[ss[sn['oppos']]]
    elif bmode=="same":
        varf[ss[sn['out']]] = refvar[ss[sn['neigh1']]]
    else:
        varf[ss[sn['out']]] = refvar[ss[sn['neigh1']]]
        varf[ss[sn['out']]] += abs(shift) * (refvar[ss[sn['neigh1']]]-refvar[ss[sn['neigh2']]])
    del refvar

    # Axes
    if ax:
        varo.assignValue(varf)
    elif xo is not None:
        varo.setAxis(axis, xo)

    return varo

def shift2d(vari, ishift=0, jshift=0, bmode=None, copy=True, **kwargs):
    """Interpolate data on an grid shifted by an half cell size in X and Y

    X and Y are supposed to be the -1 and -2 axes of var.

    :Params:

        - **var**: array.
        - **[i/j]shift**, optional: Shift to operate along X/Y.

            - ``0``: No shift.
            - ``<0``: Shilt toward bottom of the axis.
            - ``>0``: Shilt toward top of the axis.

        - **bmode**, optional: Boundary mode.

            - ``"linear"``: Linear extrapolation.
            - ``"same"``: Constant extrapolation.
            - ``"masked"``: Mask outside data.
            - ``"cyclic"``: Cyclic extrapolation.

        - **copy**, optional: Always input copy data.

    """
    if copy:
        if hasattr(vari, 'clone'):
            var = vari.clone()
        else:
            var = vari.copy()
    else:
        var = vari
    if ishift==0 and jshift==0: return var
    ny, nx = var.shape[-2:]
    ne = N.multiply.reduce(var.shape)/(nx*ny)
    sg = not A.isaxis(var) and cdms2.isVariable(var) and var.getGrid() is not None
    if sg:
        grido = shiftgrid(var.getGrid(), ishift=ishift, jshift=jshift)
    bmode = kwargs.get('mode', bmode)
    if bmode=='masked' and not N.ma.isMA(var):
        var = N.ma.asarray(var)


    # Slices
    ssx = get_axis_slices(var, -1)
    snx = _shiftslicenames_(ishift)
    ssy = get_axis_slices(var, -2)
    sny = _shiftslicenames_(jshift)

    # Shift along X only
    kwxshift = dict(axis=-1, bmode=bmode, copy=copy or jshift!=0, shiftaxis=not sg)
    if jshift==0 or ishift!=0:
        varx = shift1d(vari, ishift, **kwxshift)
        if jshift==0:
            if sg:
                set_grid(varx, grido)
            return varx

    # Shift along Y only
    kwyshift = dict(axis=-2, bmode = bmode if bmode!='cyclic' else 'extrap',
        copy=copy or ishift!=0, shiftaxis=not sg)
    if ishift==0 or jshift!=0:
        vary = shift1d(vari, jshift, **kwyshift)
        if ishift==0:
            if sg:
                set_grid(vary, grido)
            return vary

    # Inner Y
    var[:] = 0.,
    var[ssy[sny['firsts']]] += 0.5*(varx[ssy[sny['firsts']]]+varx[ssy[sny['lasts']]])

    # Top Y
    sstopyfirsts = merge_axis_slice(ssy[sny['last']], ssx[snx['firsts']])
    sstopylasts = merge_axis_slice(ssy[sny['last']], ssx[snx['lasts']])
    var[sstopyfirsts] = 0.5*(vary[sstopyfirsts]+vary[sstopylasts])

    # Shift along X and Y
    sslast = merge_axis_slice(ssy[sny['last']], ssx[snx['last']])
    sslastxy = get_axis_slices(var.ndim-1, -1)
    sslastx = sslastxy[snx['last']]
    sslasty = sslastxy[sny['last']]
    varlastx = shift1d(varx[ssx[snx['last']]], ishift, **kwxshift) #; del varx
    var[sslast] += 0.5*varlastx[sslastx] #; del varlastx
    kwyshift['axis'] = -1
    varlasty = shift1d(vary[ssy[sny['last']]], jshift, **kwyshift) #; del vary
    var[sslast] += 0.5*varlasty[sslasty] #; del varlastx

    # Grid for MV2 arrays
    if sg:
        set_grid(var, grido)

    return var


def shiftgrid(gg, ishift=0, jshift=0, bmode='linear', **kwargs):
    """Shift a grid of an half cell

    :Params:

        - **gg**: cdms2 grid.
        - **i/jshift**: Fraction cell to shift.
    """
    xx, yy = get_xy(gg)
    bmode = kwargs.get('mode', bmode)
    if len(xx.shape)==2:
        xs = shift2d(xx, ishift=ishift, jshift=jshift, bmode=bmode)
        ys = shift2d(yy, ishift=ishift, jshift=jshift, bmode=bmode)
    else:
        xs = shift1d(xx, ishift, bmode=bmode)
        ys = shift1d(yy, jshift, bmode=bmode)

    if isgrid(gg) or cdms2.isVariable(gg):
        return create_grid(xs, ys)
    return xs, ys

def _extend_init_(var, new_shape, mode):
    """Some inits for extend1d and extend2d"""
    # In and out data
    if A.isaxis(var):
        datatype = 'axis'
        varm = var.getValue()
        if len(var.shape)==1:
            datatype = 'axis1d'
            varf = N.zeros(new_shape)
        else:
            datatype = 'axis2d'
            zeros = N.ma.zeros if N.ma.isMA(var) or mode=='masked' else N.zeros
            varf = zeros(new_shape)
    elif cdms2.isVariable(var):
        varm = var.asma()
        varf = N.zeros(new_shape)
        datatype = 'mv'
    else:
        datatype = 'ma'
        zeros = N.ma.zeros if N.ma.isMA(var) or mode=='masked' else N.zeros
        varm = var
        varf = N.ma.zeros(new_shape)

    # Extrapolation mode
    if mode is None:
        mode = "extrap" if A.isaxis(var) else "same"
    elif mode=="nearest":
        mode = "same"
    elif mode=='linear': mode = 'extrap'

    return varm, varf, datatype, mode

def extend1d(var, ext=0, mode=None, axis=-1, copy=False, num=False):
    """Extrapolate data on an axis

    :Params:

        - **var**: Array.
        - **ext**, optional: Size of extension. If a tuple,
          it gives left and right extensions, else the same
          for both.

        - **mode**, optional: Interpolation mode for boundary point outside initial positions.

            - ``None``: ``"extrap"`` if axis else ``"same"``
            - ``"extrap"`` or ``"linear"``: Linear extrapolation.
            - ``"same"``: Constant extrapolation.
            - ``"cyclic"``: Cyclic extrapolation.
            - ``"masked"``: Masked.

        - **axis**, optional: Axis on which to operate.
        - **copy**, optional: Always input copy data.
    """
    # Size of extension
    if not isinstance(ext, (list, tuple)):
        ext = (ext, ext)
    if ext[0]==0 and ext[1]==0:
        if copy: return var.clone() if hasattr(var, 'clone') else var.copy()
        return var

    # Data and mode
    varo_shape = list(var.shape)
    if axis<0: axis = len(varo_shape)+axis
    ni = varo_shape[axis]
    no = ni+ext[0]+ext[1]
    varo_shape[axis] = no
    varm, varf, datatype, mode = _extend_init_(var, varo_shape, mode)

    # Get slice specs
    ss = get_axis_slices(var, axis, extinner=slice(ext[0], -ext[1] or None),
        extleft=slice(0, ext[0]), extright_for_left=slice(-ext[0], None),
        extright=slice(-ext[1], None), extleft_for_right=slice(0, ext[1]),
    )

    # Unchanged data
    varf[ss['extinner']] = varm

    # Left boundary
    if ext[0]:
        if mode=='masked':
            varf[ss['extleft']] = N.ma.masked
        elif mode=='cylic':
            varf[ss['extleft']] = varm[ss['extright_for_left']]
        else:
            if mode=='extrap':
                dv = varm[ss['first']]-varm[ss['firstp1']]
            ss['extleft'] = list(ss['extleft'])
            for i in xrange(1, ext[0]+1):
                ss['extleft'][axis] = ext[0]-i
                varf[tuple(ss['extleft'])] = varm[ss['first']]
                if mode=='extrap':
                    varf[tuple(ss['extleft'])] += i*dv

    # Right boundary
    if ext[1]:
        if mode=='masked':
            varf[ss['extright']] = N.ma.masked
        elif mode=='cylic':
            varf[ss['extright']] = varm[ss['extleft_for_right']]
        else:
            if mode=='extrap':
                dv = varm[ss['last']]-varm[ss['lastm1']]
            ss['extright'] = list(ss['extright'])
            for i in xrange(1, ext[1]+1):
                ss['extright'][axis] = ext[0]+ni+i-1
                varf[tuple(ss['extright'])] = varm[ss['last']]
                if mode=='extrap':
                    varf[tuple(ss['extright'])] += i*dv
    del varm

    # Format
    if not num:
        if datatype=='axis1d':

            varo = A.create(varf)
            cp_atts(var, varo)

        elif datatype=='axis2d':

            varo = create_axes2d(varf)
            cp_atts(var, varo)
            for iaxis in 0, 1:
                cp_atts(var.getAxis(iaxis), varo.getAxis(iaxis))

        elif datatype=='mv':

            varo = MV2.asarray(varf)
            cp_atts(var, varo)
            gridi = var.getGrid()
            for iaxis in xrange(varo.ndim):
                axi = var.getAxis(iaxis)
                if axis!=iaxis:
                    varo.setAxis(iaxis, axi)
                elif gridi is None:
                    axo = extend1d(axi, ext, mode='extrap')
                    cp_atts(axi, axo)
            if gridi is not None:
                if axis>=var.ndim-2:
                    iext = ext if axis==var.ndim-1 else 0
                    jext = ext if axis==var.ndim-2 else 0
                    grido = extendgrid(gridi, iext=iext, jext=jext)
                    set_grid(varo, extendgrid(gridi, iext=iext, jext=jext))
                else:
                    set_grid(varo, gridi)
        else:
            varo = varf
    else:
        varo = varf

    return varo

def extend2d(vari, iext=0, jext=0, mode=None, copy=False):
    """Interpolate data on an grid shifted by an half cell size in X and Y

    X and Y are supposed to be the -1 and -2 axes of var.

    :Params:

        - **var**: Array.
        - **i/jext**, optional: Size of extension along i/j.
          If a tuple,it gives left and right extensions, else the same
          for both.

        - **mode**, optional: Interpolation mode for boundary point outside initial positions.

            - ``None``: ``"extrap"`` if axis else ``"same"``
            - ``"extrap"`` or ``"linear"``: Linear extrapolation.
            - ``"same"``: Constant extrapolation.
            - ``"masked"``: Masked.

        - **copy**, optional: Always input copy data.

    """
    # Size of extensions
    if not isinstance(iext, (list, tuple)):
        iext = (iext, iext)
    if not isinstance(jext, (list, tuple)):
        jext = (jext, jext)

    # No extension
    if iext[0]==0 and iext[1]==0 and jext[0]==0 and jext[1]==0:
        if copy: return vari.clone() if hasattr(vari, 'clone') else vari.copy()
        return vari

    # Single axis
    if iext[0]==0 and iext[1]==0:
        return extend1d(vari, ext=jext, axis=-2, mode=mode, copy=copy)
    elif jext[0]==0 and jext[1]==0:
        return extend1d(vari, ext=iext, axis=-1, mode=mode, copy=copy)

    # X then Y
    varo = extend1d(vari, ext=iext, axis=-1, mode=mode, copy=False)
    varo = extend1d(varo, ext=jext, axis=-2, mode=mode, copy=False)

    # Merge with Y then X
    tmp = extend1d(vari, ext=jext, axis=-2, mode=mode, copy=False, num=True)
    varo += extend1d(tmp, ext=iext, axis=-1, mode=mode, copy=False, num=True)
    del tmp
    varo[:] *= 0.5


    return varo

def extendgrid(gg, iext=0, jext=0, mode='extrap'):
    """Extrapolate a grid

    :Params:

        - **gg**: cdms2 grid.
        - **i/jext**: Size of extrapolation along i/j.
        - **mode**, optional: Interpolation mode for boundary point outside initial positions.

            - ``"extrap"``: Linear extrapolation.
            - ``"same"``: Constant extrapolation.
            - ``"masked"``: Masked.

    """
    xx, yy = get_xy(gg)

    if len(xx.shape)==2:
        xs = extend2d(xx, iext=iext, jext=jext, mode=mode)
        ys = extend2d(yy, iext=iext, jext=jext, mode=mode)
    else:
        xs = extend1d(xx, iext, mode=mode)
        ys = extend1d(yy, jext, mode=mode)

    if isgrid(gg) or cdms2.isVariable(gg):
        return create_grid(xs, ys)
    return xs, ys


class CDATRegridder(object):
    """Regridding using CDAT regridders

    .. note:: This code is adapted from the :meth:`regrid` method of MV2 arrays.

    :Params:

        - **fromgrid**: Input grid, or variable with a grid.
        - **togrid**: Output grid.
        - **tool**, optional: One of ``"esmf"``, ``"libcf"`` and ``"regrid2"``.
        - **method**, optional: One of ``"linear"``, ``"path"``, ``"conservative"`` and ``"cellave"``.
    """

    def __init__(self, fromgrid, togrid, missing=None, order=None, mask=None, **keywords):

        if cdms2.isVariable(fromgrid) and not isgrid(fromgrid):
            fromvar = fromgrid
            fromgrid = fromvar.getGrid()
        else:
            fromvar = None

        # Verify bounds
        for g in fromgrid, togrid:
            for xy in g.getLongitude(), g.getLatitude():
                if xy.getBounds() is None:
                    func = bounds1d if len(xy.shape)==1 else bounds2d
                    xy.setBounds(func(xy))

        # Curved?
        curved = len(fromgrid.getLatitude().shape) > 1 or len(togrid.getLatitude().shape) > 1

        # Default
        regridTool = None
        regridMethod = 'linear'

#        # Some keywords for longitudes
#        for lon in fromgrid.getAxis(-1), fromgrid.getLongitude():
#            if lon.attributes.get('topology',  None) == 'circular':
#                keywords['periodicity'] = 1 # for the ESMF regridders
#                keywords['mkCyclic'] = 1    # for LibCF regridder
#                break

        # User request
        # - method
        userSpecifiesMethod = False
        for rm in 'rm', 'method', 'regridmethod', 'regrid_method', 'regridMethod':
            if keywords.has_key(rm):
                if keywords[rm] is not None:
                    regridMethod = keywords[rm]
                    userSpecifiesMethod = True
                del keywords[rm]
        # - tool
        for rt in 'rt', 'tool', 'regridtool', 'regrid_tool', 'regridTool':
            if keywords.has_key(rt):
                if keywords[rt] is not None:
                    regridTool = keywords[rt]
                del keywords[rt]

        # The method and grid determine the tool
        if not re.search('^(cons|cell)', regridMethod, re.I):
            regridTool = 'esmf'
        elif regridTool is None:
            if not curved:
                regridTool = 'regrid2'
            else:
                regridTool = 'esmf'

        # Make sure the tool can do it
        if re.search('^regrid', regridTool, re.I) and curved:
            regridTool = 'esmf'

        # Clean regrid tool names
        regridTool = regrid2d_tool_name(regridTool)
        self.tool = regridTool

        # Conservative versus cellave
        self._getdstfracs = False
        self._weidstfracs = 0
        self._mskdstfracs = False
        if re.search('esmf',regridTool, re.I) and re.search('(cons|cell)', regridMethod, re.I):
            if re.search('cell', regridMethod, re.I):
                regridMethod = 'conservative'
                self._weidstfracs = -1
            self._getdstfracs = True
            self._mskdstfracs = True
        elif re.search('^regrid',regridTool, re.I):
            if re.search('cell', regridMethod, re.I):
                regridMethod = 'conservative'
            else:
                self._weidstfracs = 1
            self._getdstfracs = True
            self._divdstfracs = True

        # Clean method names
        regridMethod = regrid2d_method_name(regridMethod, raiseerr=False)
        if regridMethod=='bilinear':
            regridMethod = 'linear'
        elif re.search('patch',regridMethod, re.I):
            regridMethod = 'patch'
        elif regridMethod=='conserv':
            regridMethod = 'conservative'
        else:
            raise VACUMMError('Wrong CDATRegridder regrid method. Please one of: '+
                'linear, conserv, patch, cellave')
        self.method = regridMethod

        # Setup regridder
        if regridTool=='regrid2': # Origignal CDMS regridder

            keywords['returnTuple'] = self._getdstfracs
            from regrid2 import Horizontal
            self._regridder = Horizontal(fromgrid, togrid)
            self._kwargs = dict(missing=missing, order=order,
                   mask=mask, **keywords)

        else: # ESMF or LIBCF regridders

            # Source mask
            srcGridMask = None
            dtype = N.float64
            if fromvar is not None:
                # set the source mask if a mask is defined with the source data
                if N.any(fromvar.mask == True):
                    srcGridMask = cdms2.avariable.getMinHorizontalMask(fromvar)
                dtype = fromvar.dtype
                if dtype.kind!='f':
                    dtype = N.float64
#            if regridMethod=='conservative':
#                # Hack for conservative
#                if srcGridMask is None:
#                    srcGridMask = N.zeros(fromgrid.shape, '?')
#                srcGridMask[:, -1] = True

            # DstAreaFractions for cell
            if self._getdstfracs:
                diag = keywords.get('diag', {})
                diag['dstAreaFractions'] = None
                keywords['diag'] = diag

            # compute the interpolation weights
            self._regridder = cdms2.mvCdmsRegrid.CdmsRegrid(fromgrid, togrid,
                            dtype = dtype,
                            regridMethod = regridMethod,
                            regridTool = regridTool,
                            srcGridMask = srcGridMask,
                            srcGridAreas = None,
                            dstGridMask = None,
                            dstGridAreas = None,
                            **keywords)
            self._kwargs = keywords

    def regrid(self, vari, weidstfracs=None, check_mask=None, csvhack=False, **keywords):
        """Regrid the variable

        :Params:

            - **vari**: Input variable.
            - **weidstfracs**, optional: Specify how to apply weights computed using
              destination fractions. Divide if <0, else multiply. It is not used with
              linear and match methods.
            - **check_mask**, optional: Mask point using masked data (the algo
              interpolate the mask and check level 0.999). MUST BE IMPROVED!
              If None, mask is checked if tool!='regrid2'
            - **csvhack**, optional: Hack to prevent a bug with conservative-like method
              of the ESMF regridder. Use it if you have strange result in the output
              most right longitude.

              .. warning:: The algo mask the must right longitude of input data before
                processing.

        """
        # Float type
        if vari.dtype.kind!='f':
            vari = vari.astype('d')

        # Hack for conservative method bug
        if csvhack and self.method=='conservative':
            mask = N.zeros(vari.shape, '?')
            mask[..., -1] = True
            vari = MV2.masked_where(mask, vari)
            del mask

        # Regrid
        kwargs = self._kwargs.copy()
        kwargs.update(keywords)
        res = self._regridder(vari, **kwargs)

        # Outputs
        if not isinstance(res, tuple):
            varo = res
            if self._getdstfracs:
                wo = N.resize(kwargs['diag']['dstAreaFractions'], varo.shape)
        else: # regrid2 old regridder
            varo, wo = res

        # Working with destination fractions
        if self._getdstfracs:
            mask = wo==0.

            # Masking
            if self._mskdstfracs:
                varo[:] = MV2.masked_where(mask, varo, copy=0)

            # Scaling
            if weidstfracs is None:
                weidstfracs = self._weidstfracs
            if weidstfracs:
                if weidstfracs==1:
                    varo[:] *= wo
                else:
                    wo[mask] = 1.
                    varo[:] /= wo
            del wo, mask

        # Check mask
        if check_mask is None:
            check_mask = self.tool != 'regrid2'
        if check_mask:
            good = vari.clone()
            good[:] = 1-N.ma.getmaskarray(vari)
            goodo = self.regrid(good, weidstfracs=weidstfracs, check_mask=False, csvhack=csvhack, **keywords)
            del good
            if isinstance(goodo, tuple): goodo = goodo[0]
            goodo = goodo.filled(0.)
            varo[:] = MV2.masked_where((goodo<0.9999), varo, copy=0) ; del goodo

        # Finalize
        cp_atts(vari, varo)
        return varo

    __call__ = regrid



class CurvedInterpolator(object):
    """Interpolator from a curved grid to random points or another grid

    :Params:

        - **fromgrid**: Input grid, or variable with a grid.
        - **topts**: Output coordinates or grid. It a tuple of coordinates, it is
          assumed that it refers to random points, and NOT a grid.
        - **g2g**, optional: Force the interpretation of ``topts`` as a grid or axes
          of a grid.

    :Examples:

        >>> interpolator = CurvedInterpolator(ssti.getGrid(), (lono, lato)) # lono/lato 1D
        >>> interpolator = CurvedInterpolator(ssti.getGrid(), grido) # grid
        >>> ssto = interpolator(ssti)
    """
    valid_methods = ['bilinear', 'nearest', 'dstwgt']

    def __init__(self, fromgrid, topts, g2g=False):

        # Input coordinates
        xxi, yyi = get_xy(fromgrid, mesh=True, num=True)
        self._shapei = xxi.shape

        # Output coordinates
        if cdms2.isVariable(topts): topts = topts.getGrid()
        if g2g:
            topts = get_grid(topts)
        if isgrid(topts):
            xxo, yyo = meshgrid(topts.getLongitude(), topts.getLatitude())
            xo = xxo.ravel()
            yo = yyo.ravel()
            self._grido = topts
            self._shapeo = topts.shape
        else:
            xo, yo = topts
            if len(xo)!=len(yo):
                raise VACUMMError('Output axes must have the same length: %i!=%i'%
                    (len(xo), len(yo)))
            xo = N.asarray(xo)
            yo = N.asarray(yo)
            self._grido = None
            self._shapeo = xo.shape

        # Find relative positions
        self._p, self._q = _curv2rel_(xxi, yyi, xo, yo)

    def interp(self, vari, method='bilinear'):
        """Interpolate

        :Params:

            - **vari**: Variable to interpolate.
            - **method**, optional: Interpolation method. One of : %s.
        """%', '.join(self.valid_methods)
        # Method and function
        method = regrid2d_method_name(method)
        if method not in self.valid_methods:
            raise VACUMMError('Invalid interpolation method. Choose of: '+
                ', '.join(self.valid_methods))
        if method=='bilinear': method = 'bilin'
        func = eval('_%s2dto1dc_reduc_'%method)

        # Interpolate
        mv = vari.get_fill_value()
        zzi = vari.asma() if cdms2.isVariable(vari) else N.ma.asarray(vari)
        zzi = zzi.reshape((-1, )+self._shapei).filled(mv)
        zo = func(self._p, self._q, zzi, mv).reshape(vari.shape[:-2]+self._shapeo)
        zo = N.ascontiguousarray(zo)
        varo = N.ma.masked_values(zo, mv)

        # Finalize the cdms variable
        if cdms2.isVariable(vari):
            varo = MV2.asarray(varo)
            cp_atts(vari, varo)
            for i, ax in enumerate(vari.getAxisList()[:-2]):
                varo.setAxis(i, ax)
            if self._grido:
                set_grid(varo, self._grido)

        return varo

    regrid = __call__ = interp

def _monotonise_(vari, axes, targets=None, back=False, subdims=None):
    if targets is None:
        targets = range(-len(axes), 0)
    if not isinstance(targets, list):
        targets = [targets]
    if not targets:
        if back:
            return vari
        return vari, axes
    n = vari.ndim
    targets = [(t if t < 0 else t-n) for t in targets]
    axes = [ax[:] for ax in axes]
    if not back:
        oaxes = list(axes)
    varo = vari
    for tg in targets:
        ax = axes[tg]

#        # 2D case
#        #FIXME: add case when ax is 2d but not curvilinear like variable depths
#        ash = ax.shape
#        if len(ash)==2:
#            if min(ash)>1:
#                continue
#            ax = ax.ravel()

        if ax.size>1:
            subdim = (subdims[tg] if isinstance(subdims, dict)
                and tg in subdims else -1)
            if (N.ma.diff(ax, axis=subdim)<0).any():
                sel = varo.ndim*[slice(None)]
                sel[tg] = slice(None, None, -1)
                varo = varo[tuple(sel)]
                if not back:
                    sel = ax.ndim*[slice(None)]
                    sel[subdim] = slice(None, None, -1)
                    oaxes[tg] = ax[tuple(sel)]
#USE SUBDIM

#        ax.shape = ash

    if back:
        return varo
    return varo, oaxes




######################################################################
######################################################################
