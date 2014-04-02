# -*- coding: utf8 -*-
"""Regridding utilities

.. seealso::

    Tutorials: :ref:`user.tut.misc.grid.regridding`
"""
# Copyright or © or Copr. Actimar (contributor(s) : Stephane Raynaud) (2010)
# 
# raynaud@actimar.fr
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
# 
import gc, os, subprocess
import re
import warnings

import numpy as N, cdms2,  MV2,  regrid2
from cdms2.axis import TransientAxis
#MA = N.oldnumeric.ma
MA = N.ma
import genutil
import tempfile, os, shutil
from copy import deepcopy
from _geoslib import Point, Polygon
MV=MV2
cdms=cdms2
from vacumm import VACUMMError, VACUMMWarning
from kriging import krig as _krig_

# Python functions
__all__ = ['fill1d', 'regular', 'regular_fill1d', 'cellave1d', 'spline_interp1d', 
    'refine', 'GridData', 'griddata', 'cargen', 'fill2d', 'regrid2d', 
    'SCRIP', 'scrip', 'regrid1d', 'interp1d', 'nearest1d', 'cubic1d',
    'xy2grid', 'grid2xy', 'fill1d', 'GriddedMerger', 'regrid_method', 
    'cellave2d', 'interp2d', 'xy2xy', 'shift1d', 'shift2d', 
    'shiftgrid',  'transect', 'CDATRegridder', 'extend1d', 'extend2d', 
    'extendgrid', 'regrid2d_method_name', 'fill1d2', 'krig']
__all__.sort()

# Fortran functions
_interp_funcs = ['interp1d', 'interp1dx', 'interp1dxx', 
    'remap1d', 'remap1dx', 'remap1dxx', 'nearest2d', 'bilin', 'dstwgt',
    'mbilin2d', 'mixt2dx', 'cargen', 'bilin2dto1d', 
    'nearest2dto1d', 'nearest2dto1dc', 'bilin2dto1dc']
    
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
    #cmd = ["f2py", "-c","-m","_interp_"]
    cmd = ["make", "_interp_.so"]
    #if None in [os.environ.get('F77'), os.environ.get('F90')]:
    #    cmd.append("--fcompiler=gnu95")
    #cmd.append("interp.f90")
    out = subprocess.Popen(cmd, cwd=os.path.dirname(__file__), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if out[1]!='':
        raise ImportError("Can't build _interp_ for importation:\n%s"%('\n'.join(out)))
    exec import_interp


# Interpolation methods
_griddata_methods = ['nat', 'natgrid', 'css', 'splines', 'carg', 'krig']
_cellave_methods = ['conservative', 'remap', 'cellave', 'conserv']
_cdat_methods = _cellave_methods+['bilinear', 'patch']
_regrid2d_methods = ['nearest', 'mixt', 'interp', 'bining']+_griddata_methods+_cdat_methods
_interp1d_methods = ['nearest', 'linear', 'cubic', 'hermit']
_regrid1d_methods = _interp1d_methods+_cellave_methods

def regrid1d_method_name(method, raiseerr=True):
    """Check the method name and return its generic name"""
    if method is None or method.lower()=='auto': return 'auto'
    method = method.lower()
    if method.startswith('cons'): return 'conserv'
    if method.startswith('cell') or method.startswith('remap'): return 'cellave'
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


def regrid1d(vari, axo, method='auto', axis=None, xmap=None, xmapper=None, mask_thres=.5, 
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
        if method in _interp1d_methods:
            method = _interp1d_methods.index(method)
        else:
            method = -2 # Wrong method
    else:
        method = int(method)
    assert method >= -1 and method < 4, 'Wrong method'
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
        

def _regrid1dnew_(vari, axo, method='auto', axis=None, axi=None, iaxo=None, iaxi=None, 
    xmap=None, xmapper=None, mask_thres=.5, extrap=0):
    """Interpolation along one axis
    
    :Params:
    
        - **vari**: Input cdms array.
        - **axo**: Output cdms axis or array. It can be of any dimensions.
        - **method**:
        
            - ``"nearest"``: Nearest neighbor
            - ``"linear"``: Linear interpolation
            - ``"cubic"``: Cubic interpolation
            - ``"cellave"``: Cell averaging
            - ``"conserv"``: Conservative cel averaging (like ``cellave`` but with integral preserved)
            
        - **axis**, optional: Dimension (int) on which the interpolation is performe. 
          If not specified, it is guessed from the 
          input and output axis types, or set to ``0``.
        - **axi**, optional: Input axis. It defaults to the axis-th axis of ``vari``.
          Like ``axo``, it can be of any dimensions.
        - **iaxo**, optional: Dimension of ``axo`` on which the interpolation is performed
          when ``axo`` has more than one dimension.
        - **iaxi**, optional: Same as ``iaxo`` but for ``axi`.
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
    if missing_value is None:
        vari.setMissing(1.e20)
    missing_value = vari.getMissing()
    
    # Output axis
    axon = axo[:] 
    if N.ma.isMA(axon): axon = axon.filled(missing_value)

    # Working data axis
    if axis is None:# Guess it
    
        if A.isaxis(axo): # From axis type
            axis = order.find(A.axis_type(axo))
            
        elif axon.ndim==1:  # From axis length
            axis = _subshape_(vari, axo)
            if axis is None: 
                raise VACUMMError('Please, specify the "axis" parameter for interpolation')
            
        # Not found, so 0 without verification now (later)
        if axis == -1 or axis is None: axis = 0
        
    elif axis == -1:      # Last axis
        axis = vari.ndim-1
    elif axis < 0 or axis >= vari.ndim:
        raise VACUMMError('Wrong "axis" parameter for interpolation')
    if axis in [vari.ndim-2, vari.ndim-1]:
        grid = None
        
    # Input axis
    # - get it
    if axi is None: axi = vari.getAxis(axis)
    # - special case for time axis: must have same units !
    if (A.axis_type(axo), order[axis]) == ('t', 't') and \
        hasattr(axi, 'units') and hasattr(axo, 'units') and \
        not are_same_units(axi.units, axo.units):
        axi = ch_units(axi, axo.units, copy=True)
    # - numeric version
    axin = axi[:]

    # Subaxes
    nxo = axon.ndim
    if nxo>1:
        if iaxo is None:
            iaxo = _getiax_(vari, axo, axis)
        if iaxo is None: 
            raise VACUMMError("Please, specifiy the 'iaxo' parameter")
    else:
        iaxo = 0
    nxi = axin.ndim
    if nxi>1:
        if iaxi is None:
            iaxi = _getiax_(vari, axi, axis)
        if iaxi is None: 
            raise VACUMMError("Please, specifiy the 'iaxi' parameter")
    else:
        iaxi = 0
            
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
    if axind.ndim>3:
        axind = axind.reshape((-1, axin.shape[iaxi]))
    axond, bakmapo = _toright_(axon, iaxo)
    if axond.ndim>3:
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
        if method in _interp1d_methods:
            method = _interp1d_methods.index(method)
        else:
            method = -2 # Wrong method
    else:
        method = int(method)
    assert method >= -1 and method < 4, 'Wrong method'
    if method == 2: method = 3 # Hermit is always better
    
    # Routine and arguments
    if nxi==1 and nxo==1: # 1D->1D
        interp_func = _interp1d_
        remap_func = _remap1d_
    elif nxo==1: # 1D->ND
        interp_func = _interp1dx_
        remap_func = _remap1dx_        
    else: # ND->ND
        interp_func = _interp1dxx_
        remap_func = _remap1dxx_
    if method == -1: # Cellave
        regrid_func = remap_func
        args = [conserv]
    else: # Interp
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
    
        
    # First guess
    print "regrid var"
    varo2d = regrid_func(vari2d, axind, axond, missing_value, *args)

    # Mask
    maski2d = vari2d==missing_value
    if method and N.any(maski2d):
        
        # Float mask
        maski2df = maski2d.astype('f')
        masko2d = regrid_func(maski2df, axind, axond, missing_value, *args)
        masko2d[masko2d==missing_value] = 1.
        
        # Masking
        if method == -1:
            
            # Cell case: threshold
            varo2d[:] = N.where(masko2d>mask_thres, missing_value, varo2d)
    
        else:
            
            # Lower method
            varolow = interp_func(vari2d, axind, axond, missing_value, 
                min(abs(method), 2)-1, extrap=extrap)
        
            # Cubic case: lower order again
            if method >= 2:
                masko2dc = interp_func(maski2df, axind, axond, missing_value, 1)
                masko2dc[masko2dc==missing_value] = 1.
                varolowc = interp_func(vari2d, axind, axond, missing_value, 0)
                varolow[:] = N.where(masko2dc!=0., varolowc, varolow)
                del masko2dc, varolowc
            
            # Select between nearest and linear, or linear and cubic
            varo2d[:] = N.where(masko2d!=0., varolow, varo2d)
            del varolow
            
        del maski2df, masko2d
    
    # Reshape back
    varos = varo2d.reshape(varis.shape[:-1]+axond.shape[-1:]) ; del varo2d
    varon = varos.transpose(*bakmapv) ; del varos

    # Convert to cdms
    varo = MV2.masked_values(varon, missing_value)
    cp_atts(vari, varo, id=True)
    varo.setMissing(missing_value)
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
    gc.collect()
    return varo


def nearest1d(vari, axo, **kwargs):
    """Interpolation along an axes
    
    - **vari**: Input cdms array
    - **axo**: Output cdms axis
    - *axis*: Axis on wich to operate
    - Other keywords are passed to :func:`regrid1d`
    
    .. note::
    
        This is an wrapper to :func:`regrid1d` using ``nearest`` as a default method.
        See its help for more information.
    """
    return regrid1d(vari, axo, 'nearest', **kwargs)

def interp1d(vari, axo, method='linear', **kwargs):
    """Linear or cubic interpolation along an axes
    
    
    - **vari**: Input cdms array
    - **axo**: Output cdms axis
    - *axis*: Axis on wich to operate 
    - Other keywords are passed to :func:`regrid1d`
    
    .. note::
    
        This is an wrapper to :func:`regrid1d` using ``linear`` as a default method.
        See its help for more information.
    """
    return regrid1d(vari, axo, method, **kwargs)
    
def cubic1d(vari, axo, **kwargs):
    """Cubic interpolation along an axes
    
    
    - **vari**: Input cdms array
    - **axo**: Output cdms axis
    - *axis*: Axis on wich to operate 
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
    cxi = MA.masked_array(xi,mask=mask).compressed()
    cvi = vi.compressed()
    cii = MA.masked_array(N.arange(len(xi)),mask=mask).compressed()
    nc = len(cii)
    cjj = N.arange(nc)

    # Identify gaps
    gaps = cii[1:]-cii[:-1]
    igaps = MA.masked_where(gaps == 1,cjj[:-1]).compressed()

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
            gapstart = ii[dd==1]+1
            gapend = ii[dd==-1]
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

        # Helper
        self._GDH = _GridDataHelper_(xi, yi, ggo, geo=geo, compress=compress)
        
        # Init interpolator
        if method in ['css', 'splines']:
            assert not self._GDH.curv, 'css method does not work with curvilinear output grids'
            from css.css import Cssgrid
            self.r = Cssgrid
            self.sub = None
            self.method = 'css'
        elif method in ['nat', None, 'natgrid']:
            from nat import Natgrid
            self.r = Natgrid
            if self._GDH.curv:
                self.sub = None
            else:
                self.sub = sub
            self.method = 'nat'
        self.margin = margin
        
        # Attributes of interpolator
        if self._GDH.curv:
            kwargs['igr'] = 0
        else:
            kwargs.setdefault('igr', nl)
        kwargs['ext'] = ext
        kwargs.setdefault('hor', -1)
        kwargs.setdefault('ver', -1)
        self.ratts = kwargs
        
    def __call__(self, zi, missing_value=None, **kwargs):
        """Interpolate zi on output grid
        
        - **zi**: At least a 1D array.
        """

        # Init
        zi2d, zo3d, mo3d, nex, compress, missing_value = self._GDH.init(zi, missing_value)
        
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
            get = self._GDH.get(zi2d, iex, compress)
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
                    if self._GDH.curv:
                        xo = self._GDH.x[yslice, xslice].ravel().astype('d')
                        yo = self._GDH.y[yslice, xslice].ravel().astype('d')
                    else:
                        xo = self._GDH.x[xslice].astype('d')
                        yo = self._GDH.y[yslice].astype('d')
                        
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
                    if self._GDH.curv: # Curvilinear
                        interpolator = self.r(xi[goodi], yi[goodi], xo, yo, listOutput = 'yes')
                    else: # Rectangular
                        interpolator = self.r(xi[goodi], yi[goodi], xo, yo)
                    for att, val in self.ratts.items():
                        setattr(interpolator, att, val)
                    interpolator.nul = missing_value

                    # Interpolation
                    if constant:
                        if self._GDH.curv:
                            sh = self._GDH.x[yslice, xslice].shape
                        else:
                            sh = len(self._GDH.y[yslice]), len(self._GDH.x[xslice])
                        zzo = N.zeros(sh)
                        zzo += zzi[0]
                        if mmi is not None:
                            mmo = N.zeros(sh)
                            mmo += mmi[0]
                    elif self._GDH.curv: # Curvilinear
                    
                        sh = self._GDH.x[yslice, xslice].shape
                        zzo = interpolator.rgrd(zzi[goodi]).reshape(sh)
                        if mmi is not None:
                            mmo = interpolator.rgrd(mmi[goodi]).reshape(sh)
                            
                    else: # Rectangular
                        try:
                            zzo = interpolator.rgrd(zzi[goodi])
                            if mmi is not None:
                                mmo = interpolator.rgrd(mmi[goodi])
                        except:
                            continue
                        del interpolator
                        if self.method == 'nat':
                            zzo = zzo.T
                            if mmi is not None:
                                mmo = mmo.T
                    zo3d[iex, yslice, xslice] = zzo
                    if mmi is not None:
                        mo3d[iex, yslice, xslice] = mmo
                        del mmo
                    del zzo

            if hasattr(zi2d, 'mask') and zi2d[iex].mask is not MA.nomask: 
                del zzi, xi, yi
        gc.collect()
        
        del zi2d
        
        return self._GDH.format(zi, zo3d, mo3d, missing_value, **kwargs)
        
    rgrd = __call__
    regrid = __call__

def cargen(xi, yi, zi, ggo, mask=None, geo=None, compress=False, missing_value=None, **kwargs):
    """Interpolator from IFREMER
    
    :Params:
    
        - **xi**: Input 1D X positions.
        - **yi**: Input 1D Y positions.
        - **ggo**, optional: Output grid. Can either (xo,yo), a cdms grid or a cdms variable with a grid.
        - **mask**, optional: Mask to apply to output data [default: None]
    """
    # Helper
    GDH = _GridDataHelper_(xi, yi, ggo, geo=geo, mask=mask, compress=compress)
    assert not GDH.curv, 'cargen does not work with output curvilinear grids'
    
    # Init data
    zi2d, zo3d, mo3d, nex, compress, missing_value = GDH.init(zi, missing_value)
    
    # Loop on supplementary dims
    for iex in xrange(nex):
        
        # Get data and mask
        get = GDH.get(zi2d, iex, compress)
        if get is None: continue
        xi, yi, zzi, mmi = get
            
        # Interpolate
        zo3d[iex] = _cargen_(xi, yi, zzi, GDH.x[:], GDH.y[:], 1.e20).transpose()
        if mmi is not None:
            mo3d[iex] = _cargen_(xi, yi, mmi, GDH.x[:], GDH.y[:], 1.).transpose()
    
    # Format output
    return GDH.format(zi, zo3d, mo3d, missing_value, **kwargs)
    
krigdata = cargen

def krig(xi, yi, zi, ggo, mask=None, geo=None, compress=False, missing_value=None, **kwargs):
    """Kriging interpolator
    
    :Params:
    
        - **xi**: Input 1D X positions.
        - **yi**: Input 1D Y positions.
        - **ggo**, optional: Output grid. Can either (xo,yo), a cdms grid or a cdms variable with a grid.
        - **mask**, optional: Mask to apply to output data [default: None]
    """
    # Helper
    GDH = _GridDataHelper_(xi, yi, ggo, geo=geo, mask=mask, compress=compress)
    assert not GDH.curv, 'cargen does not work with output curvilinear grids'
    
    # Init data
    zi2d, zo3d, mo3d, nex, compress, missing_value = GDH.init(zi, missing_value)
    xo, yo = N.meshgrid(GDH.x[:], GDH.y[:])
    xo.shape = -1
    yo.shape = -1
    
    # Loop on supplementary dims
    for iex in xrange(nex):
        
        # Get data and mask
        get = GDH.get(zi2d, iex, compress)
        if get is None: continue
        xi, yi, zzi, mmi = get
            
        # Interpolate
        zo3d[iex] = _krig_(xi, yi, zzi, xo, yo).reshape(zo3d.shape[-2:])
        if mmi is not None:
            mo3d[iex] = _krig_(xi, yi, mmi, xo, yo).reshape(zo3d.shape[-2:])
    
    # Format output
    return GDH.format(zi, zo3d, mo3d, missing_value, **kwargs)
    


def griddata(xi, yi, zi, ggo, method='carg', cgrid=False, **kwargs):
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
        ggu, ggv = t2uvgrids(ggo)
        
    if method in ('nat', 'css', 'spline'):
        kwout = {}
        for att in 'outtype', 'missing_value', 'clipval':
            if kwargs.has_key(att):
                kwout[att] = kwargs.pop(att)
        if cgrid:
            zou = GridData(xi, yi, ggu, method=method, **kwargs)(zi, **kwout)
            zov = GridData(xi, yi, ggv, method=method, **kwargs)(zi, **kwout)
            return zou, zov
        return GridData(xi, yi, ggo, method=method, **kwargs)(zi, **kwout)
        
    if method=='krig':
        return krig(xi, yi, zi, ggo, **kwargs)

    if cgrid:
        zou = cargen(xi, yi, zi, ggu,**kwargs)
        zov = cargen(xi, yi, zi, ggv,**kwargs)
        return zou, zov
    return cargen(xi, yi, zi, ggo, **kwargs)


class _GridDataHelper_(object):
    def __init__(self, xi, yi, ggo, geo=None, mask=None, compress=False, outtype=-1):
        # Input grid
        assert xi.shape == yi.shape and xi.ndim == 1,   'xi and yi must be 1d arrays'
        xi = N.asarray(xi,dtype='d')
        yi = N.asarray(yi,dtype='d')
        self.mi = N.ones(len(xi), '?')
        #self.mi = N.zeros(len(xi), '?')
        #xy = []
        #for i, (x, y) in enumerate(zip(xi, yi)): # Check unicity
            #if (x, y) not in xy:
                #xy.append((x, y))
                #self.mi[i] = True
        
        # - mask
        if mask is True:
            if hasattr(ggo, getMask):
                mask = ggo.getMask()
            elif hasattr(ggo, mask):
                mask = ggo.mask
        if mask in [False, None]:
            mask = MA.nomask
            
        # - axes
        self.x, self.y = get_xy(ggo, m=False)
        if geo: # Force geographic grid
            if self.x[:].ndim == 2 or self.y[:].ndim == 2:
                self.x, self.y = axes2d(self.x, self.y)
            else:
                self.x = A.create_lon(self.x)
                self.y = A.create_lat(self.y)
        if A.isaxis(self.x) or A.isaxis(self.y): # Real axes
            ggo = get_grid(ggo) 
            self.x = get_axis(ggo, -1)
            self.y = get_axis(ggo, -2)
            self.outtype = 2
        elif outtype == -1:
            if mask is not MA.nomask:
                self.outtype = 1
            else:
                self.outtype = -1
        else:
            self.outtype = outtype

        # - grid object
        if isgrid(ggo):
            self.grid = ggo
        elif A.islon(self.x) and A.islat(self.y):
            self.grid = create_grid(self.x, self.y, mask=mask)
            self.geo = True
        else:
            self.grid = None
        self.gridshape = (self.y[:].shape[0], self.x[:].shape[-1])
        self.curv = self.x[:].ndim == 2 and self.y[:].ndim == 2 
            
        self.geo = geo
        self.mask = mask
        self.compress = compress
        self.xi = xi
        self.yi = yi
        self.axes = []
        self.nx = self.x.shape[-1]
        self.ny = self.y.shape[0]
    
    def init(self, zi, missing_value):
        """Init input before regridding"""
        # Convert to right dims
        nex = zi.size/zi.shape[-1]
        if cdms2.isVariable(zi): zi = zi.asma()
        zi2d = zi.reshape(nex, zi.shape[-1]).astype('d')
        
        # Initialize output
        if missing_value is None:
            try:
                missing_value = zi.getMissing()
            except:
                missing_value = 1.e20
        zo3d = N.zeros(zi2d.shape[:1]+self.gridshape, zi.dtype)+missing_value
        unmasked = not hasattr(zi2d, 'mask') or zi2d.mask is MA.nomask or not zi2d.mask.any()
        if unmasked:
            compress = False
        else:
            compress = self.compress
        if compress or unmasked:
            mo3d = None
        else:
            mo3d = zo3d*0.
        del unmasked
        return zi2d, zo3d, mo3d, nex, compress, missing_value

    def get(self, zi2d, iex, compress):
        """Get a slice"""
        # Remove compress values or fill them
        good = self.mi
        mmi = None
        if hasattr(zi2d, 'mask') and zi2d[iex].mask is not MA.nomask and zi2d[iex].mask.any():
            if not compress: # No compression => interpolation of mask
                mmi = zi2d[iex].mask.astype('f')
                zi2d[iex] = zi2d[iex].filled(1.e20)
            else:
                good = good & ~zi2d[iex].mask
        
        # All data are bad
        if not good.any(): return None
        
        if not good.all(): # Some are good
            zzi = zi2d[iex].compress(good)
            xi = self.xi.compress(good)
            yi = self.yi.compress(good)
            if mmi is not None:
                mmi = mmi.compress(good)
        else: # There are all good
            zzi = zi2d[iex]
            xi = self.xi
            yi = self.yi
        del good
        return xi, yi, zzi, mmi

        
    def format(self, zi, zo3d, mo3d, missing_value, clipval=False, **kwargs):
        """Format output variable"""
        # Shape
        sho = zi.shape[:-1]+zo3d.shape[-2:]
        if zo3d.shape != sho:
            zo = zo3d.reshape(sho)
            del zo3d
            if mo3d is not None:
                mo = mo3d.reshape(sho)
        else:
            zo = zo3d
            mo = mo3d
            
        
        # Masking
        if self.mask is not MA.nomask:
            
            mask = self.mask
            if mask.shape != zo.shape:
                mask = N.resize(mask, zo.shape)
            zo[mask] = missing_value
        if mo3d is not None and mo is not None:
            zo[mo!=0.] = missing_value
        missing = closeto(zo, missing_value)
        if missing.any() and self.outtype==-1: 
            self.outtype = 1
        if self.outtype==-1:
            self.outtype = 0
    
        # Value clipping
        if clipval:
            valmax = zi.max()
            valmin = zi.min()
            good = ~missing
            zo[good&(zo>valmax)] = valmax
            zo[good&(zo<valmin)] = valmin

        # Pure numeric
        if self.outtype == 0: return zo
        
        # Masking
        zo = MA.masked_values(zo, missing_value, copy=0)
        if self.outtype == 1: return zo
        
        # Gridding
        zo = cdms.createVariable(zo)
        if self.grid is not None:
            set_grid(zo, self.grid)
        else:
            zo.setAxis(-1, self.x)
            zo.setAxis(-2, self.y)
        if cdms.isVariable(zi):
            cp_atts(zi, zo, id=True)
            for i, axis in enumerate(zi.getAxisList()[:-2]):
                zo.setAxis(i, axis)
        return zo

def xy2grid(*args, **kwargs):
    """Alias for :func:`griddata`
    
    .. seealso::
    
        :class:`GridData` :func:`cargen` :func:`xy2grid`
    """
    return griddata(*args, **kwargs)
    
def xy2xy(xi, yi, zi, xo, yo, nl=False, proj=True, **kwargs):
    """Interpolation between to unstructured grids using :mod:`Natgrid`
    
    :Params:
    
        - **xi/yi**: 1D input positions
        - **zi**: atleast-1D input values
        - **xo,yo**: 1D output positions
        - *nl*: Non linear interpolation using natural neighbours
        - *proj*: convert positions to meters using mercator projection
    """
    # Check input positions
    xi = N.asarray(xi, 'd')
    yi = N.asarray(yi, 'd')
    xo = N.asarray(xo, 'd')
    yo = N.asarray(yo, 'd')
    proj = kwargs.pop('geo', proj)
    if proj:
        if not callable(proj):
            proj = get_proj((xi,yi))
        xi, yi = proj(xi, yi)
        xo, yo = proj(xo, yo)
        
    # Check input type
    outtype = 0
    if cdms2.isVariable(zi):
        outtype = 2
        axes = zi.getAxisList()
        atts = get_atts(zi)
        zi = zi.asma()
    elif N.ma.isMA(zi):
        if zi.mask is not N.ma.nomask and zi.mask.any():
            outtype = 1
        
    # Check shapes
    zi = zi.copy()
    si = zi.shape
    nsi = zi.shape[-1]
    nex = zi.size/nsi
    zi.shape = (nex, nsi)
    nso = len(xo)
    zo = N.zeros((nex, nso))
    zo[:] = N.nan
    if outtype:
        mo = N.zeros((nex, nso))
        goodi = ~zi.mask
    else:
        mo = None
    
    # Loop on extra dim
    for iex in xrange(nex):
        
        if outtype:
            gi = goodi[iex]
            mi = zi.mask[iex].astype('f')
        else:
            gi = slice(None)
        
        # Build regridder only when needed
        if iex==0 or (outtype and N.any(goodi[iex-1]!=goodi[iex])):
            
            # Check that output points are inside convex hull
            hull = masking.convex_hull((xi[gi], yi[gi]), poly=True)
            go = masking.polygon_select(xo, yo, [hull], mask=True) ; del hull
            if go.all():
                del go
                go = slice(None)
            
            # Regridder
            from nat import Natgrid
            r = Natgrid(xi[gi], yi[gi], xo[go], yo[go], listOutput='yes')
            r.igr = int(nl)
            if outtype:
                rm = Natgrid(xi, yi, xo[go], yo[go], listOutput='yes')
                rm.igr = int(nl)
        
        # Regridding
        # - values
        zin = zi[iex][gi]
        if N.ma.isMA(zin): zin = zin.filled()
        zo[iex][go] = r.rgrd(zin)
        # - mask
        if outtype:
            mo[iex][go] = rm.rgrd(mi)
    del zi  
    
    # Missing points
    mnan = N.isnan(zo)
    if mnan.any():
        outtype = max(1, outtype)
        zo = N.ma.masked_where(mnan, zo)
        del mnan
    
    # Return pure numeric
    zo.shape = si[:-1]+(nso, )
    if outtype==0: return zo

    # Masking
    if mo is not None:
        mo.shape = zo.shape
        zo = N.ma.masked_where(mo>.5, zo)
        del mo, mi
    
    # Masked arrays
    if outtype == 1: return zo
        
    # cdms
    zo = MV2.asarray(zo)
    set_atts(zo, atts)
    for i, ax in axes[:-1]:
        zo.setAxis(i, ax)
    return zo
        
    

def grid2xy(vari, xo, yo, method='bilinear', outaxis=None):
    """Interpolate gridded data to ramdom positions
    
    :Params:
    
        - **vari**: Input cdms variable on a grid
        - **xo**: Output longitudes
        - **yo**: Output latitudes
        - **method**, optional: Interpolation method
        
            - ``nearest``: Nearest neighbor
            - ``bilinear``: Linear interpolation
            - ``nat``: Natgrid interpolation
        
        - **outaxis**, optional: Output spatial axis
        
            - A cdms2 axis.
            - ``None`` or ``'auto'``: Longitudes or latitudes depending 
              on the range if coordinates are monotonic, else ``'dist'``.
            - ``'lon'`` or ``'x'``: Longitudes.
            - ``'lat'`` or ``'y'``: Latitudes.
            - ``'dist'`` or ``'d'``: Distance in km.
        
    """
    # Prefer 1D axes
    grid = get_grid(vari)
    grid = curv2rect(grid)

    # Format
    # - input
    xi, yi = get_xy(grid, num=True)
    rect = xi.ndim==1
    mv = vari.getMissing()
    if mv is None:
        mv = 1.e20
        vari.setMissing(mv)
    zi = vari.filled().astype('d').reshape((-1,)+vari.shape[-2:])
    nz = zi.shape[0]
    if vari.mask is not MV2.nomask:
        mi = vari.mask.astype('d').reshape(zi.shape)
    else:
        mi = None
    # - output
    isscalar = N.isscalar(xo)
    xo = N.atleast_1d(N.asarray(xo, dtype='d'))
    yo = N.atleast_1d(N.asarray(yo, dtype='d'))
    
    # Interpolate
    if method == 'nearest' or mi is not None:
        
        func = _nearest2dto1d_ if rect else _nearest2dto1dc_
        zo = func(xi, yi, zi, xo, yo, mv)
            
        if mi is not None: zonear = zo
        
    if method == 'bilinear':
        
        # Base
        func = _bilin2dto1d_ if rect else _bilin2dto1dc_
        zo = func(xi, yi, zi, xo, yo, mv)
                    
        # Mask
        if mi is not None:
            mo = func(xi, yi, mi, xo, yo, mv)
            zo[:] = N.where(mo!=0., zonear, zo)
            del mi, mo
        
    elif method.startswith('nat'):
        
        # Build regridder
        from nat import Natgrid
        zo = N.zeros((nz, len(xo)))+mv
        if xi.ndim==1:
            xxi, yyi = xi, yi
        else:
            xxi, yyi = N.meshgrid(xi, yi)
        r = Natgrid(xxi.ravel(), yyi.ravel(), xo, yo, listOutput='yes')
        if mi is not None:
            mo = N.zeros(len(xo))
        
        # Z loop interpolation
        for iz in xrange(nz):
            
            # Base
            zo[iz] = r.rgrd(zi[iz].ravel())
            
            # Mask
            if mi is not None:
                mo[:] = r.rgrd(mi[iz].ravel())
                zo[iz] = N.where(mo!=0., zonear[iz], zo[iz])
        
        if mi is not None: del mi, mo
        
    elif method != 'nearest':
        raise NotImplementedError, 'Method yet not implemented: '+method
        
    # Output
    zo.shape = vari.shape[:-2]+(zo.shape[-1], )
    varo = MV2.masked_values(zo, mv, copy=0)
    cp_atts(vari, varo, id=True)
    if outaxis=='auto': outaxis = None
    if outaxis is None or outaxis in ['m', 'pos', 'd', 'dist']:
        if len(xo)>1:
            xom, yom = get_proj((xo, yo))(xo, yo)
            if outaxis is None:
                outaxis = 'lon' if xom.ptp()>yom.ptp() else 'lat'
                if (outaxis=='lon' and (N.diff(N.sign(N.diff(xom)))!=0).any() or (N.diff(xom)==0).any()) or \
                    (outaxis=='lat' and (N.diff(N.sign(N.diff(yom)))!=0).any() or (N.diff(yom)==0).any()):
                    outaxis = 'dist'
        elif isscalar:
            outaxis = None
        else:
            outaxis = 'lon'
    if outaxis in ['x', 'lon']:
        outaxis = A.create_lon(xo)
    elif outaxis in ['y', 'lat']:
        outaxis = A.create_lat(yo)
    elif outaxis in ['m', 'd', 'pos', 'dist']:
        dist = N.concatenate(([0],N.sqrt(N.diff(xom)**2+N.diff(yom)**2)*0.001)).cumsum()
        outaxis = cdms2.createAxis(dist, id='position')
        outaxis.long_name = 'Distance along transect'
        outaxis.units = 'km'
    if outaxis is not None and not A.isaxis(outaxis):
        outaxis = cdms2.createAxis(outaxis, id='position')
        outaxis.long_name = 'Position along transect'
    axes = vari.getAxisList()[:-2]
    if outaxis is not None:
        axes.append(outaxis)
    else:
        varo = varo[...,0]
    varo.setAxisList(axes)
    return varo


def transect(var, lons, lats, times=None, method='bilinear', subsamp=3, 
    getcoords=False, outaxis=None, **kwargs):
    """Make a transect in a -YX variable
    
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
        
        - **times**, optional: Time sequence or axis of the same length as 
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
    time = kwargs.get('time', times) # backward compat
    if isinstance(lons, tuple): # Find coordinates
        lon0, lon1 = lons
        lat0, lat1 = lats
        lons, lats = transect_specs(var, lon0, lat0, lon1, lat1, subsamp=subsamp)
        single = False
    elif times is None and cdms2.isVariable(lons) and lons.getTime() is not None:
        times = lons.getTime()
    else: # explicit coordinates
        single = N.isscalar(lons) and N.isscalar(lats)
        lons = N.atleast_1d(lons)
        lats = N.atleast_1d(lats)
    
    
    # Spatial interpolation
    if outaxis=='time':
        outaxis = None 
        ko = False
    else:
        ko = outaxis is not None
    var = grid2xy(var, lons, lats, outaxis=outaxis)
    
    # Interpolate to a lagrangian time axis
    if times is not None and var.getTime() is not None:
        
        if ko: outaxis = var.getAxis(-1)
        
        # Time axis
        if not A.istime(times):
            times = A.create_time(times)
        
        # Interpolate
        if len(times)!=len(lons):
            raise VACUMMError('Your time axis must have a length of: %i (!=%i'%(
                len(times), len(lons)))
        var_square = regrid1d(var, times) # (nl,nz,nl)
        
        # Init out
        iaxis = var.getOrder().index('t')
        sel = [slice(None)]*var_square.ndim
        if outaxis is None:
            sel[-1] = 0
            oaxis = iaxis
        else:
            sel[iaxis] = 0
            oaxis = -1
        
        # Select the diagnonal the ~square matrix
        var = var_square[tuple(sel)].clone() # (nz,nl)
        varsm = var_square.asma()
        ii = N.arange(len(times))
        sel = [slice(None)]*varsm.ndim
        sel[iaxis] = sel[-1] = ii
        varsm = varsm[tuple(sel)]
        if oaxis!=0:
            varsm = N.rollaxis(varsm, 0, oaxis)
        var[:] = varsm
        
        
    # Single point
    if single: var = var[...,0]
    
    if not getcoords: return var
    return var, lons, lats

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
        if not MA.isMA(vari):
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
            varo[:] = _mbilin2d_(xi[:], yi[:], vari.filled(1.e20), xxo.flat, yyo.flat, geo, 
                1.e20, smoothcoast).reshape(varo.shape)
                
            # Masking
            if vari.mask is not MV.nomask:
                fmaski = vari.mask.astype('f')
                fmasko = _mbilin2d_(xi[:], yi[:], fmaski, xxo.flat, yyo.flat, geo, 
                    1.e20, False).reshape(varo.shape)
                varo[:] = MV.masked_where(N.greater(fmasko, .5), varo, copy=0)
            varo[:] = MV.masked_values(varo, 1.e20, copy=0)

            # Axes
            if not noaxes and cdms.isVariable(varo):
                for i in -2, -1:
#                   varo.setAxis(i, refine(vari.getAxis(i), factor))
                    varo.setAxis(i, (yo, xo)[i])

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

def regrid2d(vari, ggo, method='auto', mask_thres=.5, ext=False, 
    bilinear_masking='dstwgt', ext_masking='poly', cdr=None, getcdr=False, usecdr=None, useoldcdr=True, 
    mixt_fill=True,  check_mask=True, clipminmax=False, geo=None, **kwargs):
    """Regrid a variable from a regular grid to another
    
    If the input or output grid is curvilinear and ``method`` is set to
    ``"linear"``, ``"cellave"`` or ``"conserv"``, the :class:`CDATRegridder` is used.
    
    :Params:
    
        - **vari**: Variable cdms on regular grid
        - **ggo**: Tuple of (x,y) or a cdms grid or a cdms variable with a grid
        - **method**, optiona: One of:
        
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
          **rectangular grids only*:
            
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
            maskoext = N.resize(masking.grid_envelop_mask(ggi, ggo, poly=ext_masking=="poly"), (nzi, nyo, nxo))
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
        varo3d = _regrid2d_nearest2d_(vari3d.filled(mv), xxi, yyi, xxo, yyo, geo, maskoext, mv)
        
    
    elif (not curved and method in ['mixt', 'dstwgt']) or \
        (method=='bilinear' and not curved and usecdr is not True):
        
        
        # Method
        wrapper = eval('_regrid2d_%s_'%method)
        
        # Interpolation
        #FIXME: wrapper pour dstwgt et mask... not sure
        varo3d = wrapper(vari3d.filled(mv), xi, yi, xo, yo, geo, mv, ext)
        
        # Masking
        if check_mask:
            
            if bilinear_masking == 'dstwgt':
                masko = wrapper(maski3d, xi, yi, xo, yo, geo, 1., ext)
                if method != 'dstwgt':
                    varo3d_dstwgt = _regrid2d_dstwgt_(vari3d.filled(mv), xi, yi, xo, yo, geo, mv, ext)
                    varo3d[:] = N.where((masko<mask_thres) & (masko>0.), varo3d_dstwgt, varo3d)
                    del varo3d_dstwgt
                varo3d[masko>=mask_thres] = mv
            else:
                masko = varo3d==mv
                varo3d_nearest = _regrid2d_nearest2d_(vari3d.filled(mv),xxi,yyi,xxo,yyo,geo,maskoext,mv)
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
#           varo3d_nearest = _regrid2d_nearest2d_(vari3d.filled(mv),xxi,yyi,xxo,yyo,geo,False,mv)
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


def _regrid2d_nearest2d_(vari3d, xxi, yyi, xxo, yyo, geo, maskoext, mv):
    """Wrapper to fortran _nearest2d_"""
    varitmp = N.asarray(vari3d, order='F')
    xxitmp = N.asarray(xxi, order='F')
    yyitmp = N.asarray(yyi, order='F')
    xxotmp = N.asarray(xxo, order='F')
    yyotmp = N.asarray(yyo, order='F')
    nb = min(xxi.shape)/50 
    if nb == 1: nb = 0
#    nb = -1
    varo3d = N.asarray(_nearest2d_(varitmp, xxitmp, yyitmp, xxotmp, yyotmp, nb,  nogeo=not geo), order='C')
    del varitmp, xxitmp, yyitmp, xxotmp, yyotmp
    if maskoext is not False:
        varo3d[maskoext] = mv
    return varo3d

def _regrid2d_bilinear_(vari3d, xi, yi, xo, yo, geo, mv, ext):
    """Wrapper to fortran _bilin_"""
    varitmp = N.asarray(vari3d, order='F')
    varo3d =  N.asarray(_bilin_(varitmp, xi, yi, xo, yo, mv, nogeo=not geo), order='C')
    del varitmp
    return varo3d
    
def _regrid2d_dstwgt_(vari3d, xi, yi, xo, yo, geo, mv, ext):
    """Wrapper to fortran _dstwgt_"""
    varitmp = N.asarray(vari3d, order='F')
    varo3d =  N.asarray(_dstwgt_(varitmp, xi, yi, xo, yo, mv, nogeo=not geo), order='C')
    del varitmp
    return varo3d

def _regrid2d_mixt_(vari3d, xi, yi, xo, yo, geo, mv, ext):
    """Wrapper to fortran _mixt2dx_"""
    varitmp = N.asarray(vari3d.T, order='F')
    varo3d = N.asarray(_mixt2dx_(varitmp, xi, yi,  xo, yo, mv, ext).T, order='C')
    del varitmp
    return varo3d

def _regrid2d_bining_(vari3d, xxi, yyi, xo, yo, mv):
    """Regridding using rebining"""
    xbins = meshcells(xo)
    ybins = meshcells(yo)
    varo3d = N.zeros(vari3d.shape[:-2]+(len(yo), len(xo)))
    for i, vari2d in enumerate(vari3d):
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
   

def shift1d(var, shift=0, mode=None, axis=-1, copy=True, shiftaxis=True):
    """Interpolate data on an axis shifted by an half cell size
    
    :Params:
    
        - **var**: array.
        - **shift**, optional: Shift to operate.
        
            - ``0``: No shift.
            - ``<0``: Shilt toward bottom of the axis.
            - ``>0``: Shilt toward top of the axis.
            
        - **mode**, optional: Interpolation mode for boundary point outside initial positions.
        
            - ``None``: ``"extrap"`` if axis else ``"same"``
            - ``"extrap"``: Linear extrapolation.
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
    if mode is None:
        mode = "extrap" if A.isaxis(var) else "same"
    elif mode=="nearest": 
        mode = "same"
    ax = len(var.shape)==1 and A.isaxis(var)
    if ax:
        varf = varo.getValue()
    else:
        varf = varo
    xo = None
    if cdms2.isVariable(varo):
        refvar = varo.asma().copy()
        if shiftaxis:
            xo = shift1d(var.getAxis(axis), shift, mode='extrap')
    else:
        refvar = varf.copy()
        
    # Get slice specs
    ss = get_axis_slices(var, axis)
    sn = _shiftslicenames_(shift)
        
    # Inner data
    varf[ss[sn['target']]] = (refvar[ss['firsts']]+refvar[ss['lasts']])*0.5
    
    # Boundary data
    if mode=='masked':
        varf[ss[sn['out']]] = N.ma.masked
    elif mode=='cyclic':
        varf[ss[sn['out']]] = var[ss[sn['oppos']]]
    elif mode=="same":
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
        
def shift2d(vari, ishift=0, jshift=0, mode=None, copy=True):
    """Interpolate data on an grid shifted by an half cell size in X and Y
    
    X and Y are supposed to be the -1 and -2 axes of var.
    
    :Params:
    
        - **var**: array.
        - **[i/j]shift**, optional: Shift to operate along X/Y.
        
            - ``0``: No shift.
            - ``<0``: Shilt toward bottom of the axis.
            - ``>0``: Shilt toward top of the axis.
            
        - **mode**, optional: Interpolation mode for boundary point outside initial positions.
        
            - ``"extrap"``: Linear extrapolation.
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
    if mode=='masked' and not N.ma.isMA(var):
        var = N.ma.asarray(var)

    
    # Slices
    ssx = get_axis_slices(var, -1)
    snx = _shiftslicenames_(ishift)
    ssy = get_axis_slices(var, -2)
    sny = _shiftslicenames_(jshift)
        
    # Shift along X
    kwxshift = dict(axis=-1, mode=mode, copy=copy or jshift!=0, shiftaxis=not sg)
    if jshift==0 or ishift!=0:
        varx = shift1d(vari, ishift, **kwxshift)
        if jshift==0: 
            if sg:
                set_grid(varx, grido)
            return varx
            
    # Shift along Y
    kwyshift = dict(axis=-2, mode = mode if mode!='cyclic' else 'extrap', 
        copy=copy or ishift!=0, shiftaxis=not sg)
    if ishift==0 or jshift!=0:
        vary = shift1d(vari, jshift, **kwyshift)
        if ishift==0:
            if sg:
                set_grid(vary, grido)
            return vary
            
    # Inner Y
    var[:] = 0.
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
            
        
def shiftgrid(gg, ishift=0, jshift=0, mode='extrap'):
    """Shift a grid of an half cell
    
    :Params:
    
        - **gg**: cdms2 grid.
        - **i/jshift**: Fraction cell to shift.
    """
    xx, yy = get_xy(gg)
    
    if len(xx.shape)==2:
        xs = shift2d(xx, ishift=ishift, jshift=jshift, mode=mode)
        ys = shift2d(yy, ishift=ishift, jshift=jshift, mode=mode)
    else:
        xs = shift1d(xx, ishift, mode=mode)
        ys = shift1d(yy, jshift, mode=mode)
    
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
    ss = get_axis_slices(var, axis, extinner=slice(ext[0], no-ext[1]), 
        extleft=slice(0, ext[0]), extright=slice(no-ext[1], None))
    
    # Unchanged data
    varf[ss['extinner']] = varm
        
    # Left boundary
    if ext[0]:
        if mode=='masked':
            varf[ss['extleft']] = N.ma.masked
        else:
            if mode=='extrap':
                dv = varm[ss['first']]-varm[ss['firstp1']]
            for i in xrange(1, ext[0]+1):
                ss['extleft'][axis] = ext[0]-i
                varf[ss['extleft']] = varm[ss['first']]
                if mode=='extrap':
                    varf[ss['extleft']] += i*dv

    # Right boundary
    if ext[1]:
        if mode=='masked':
            varf[ss['extright']] = N.ma.masked
        else:
            if mode=='extrap':
                dv = varm[ss['last']]-varm[ss['lastm1']]
            for i in xrange(1, ext[1]+1):
                ss['extright'][axis] = ext[0]+ni+i-1
                varf[ss['extright']] = varm[ss['last']]
                if mode=='extrap':
                    varf[ss['extright']] += i*dv
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
        regridTool = 'esmf'   
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
        userSpecifiesTool = False
        for rt in 'rt', 'tool', 'regridtool', 'regrid_tool', 'regridTool':
            if keywords.has_key(rt):
                if keywords[rt] is not None:
                    regridTool = keywords[rt]
                    userSpecifiesTool = True
                del keywords[rt]

        # The method determines the tool
        if re.search('(cons|cell|patch)', regridMethod, re.I):
            if not curved and re.search('cons', regridMethod, re.I):
                regridTool = 'regrid2'
            else:
                regridTool = 'esmf'
            
        # Make sure the tool can do it
        if re.search('^regrid', regridTool, re.I) and curved:
            regridTool = 'esmf'

        # Clean regrid tool names
        if re.search('(libcf|gs|gridspec)',regridTool, re.I):
            regridTool = 'libcf'
        elif re.search('^regrid',regridTool, re.I):
            regridTool = 'regrid2'
        else:
            regridTool = 'esmf'
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
            
    def regrid(self, vari, weidstfracs=None, check_mask=True, csvhack=False, **keywords):
        """Regrid the variable
        
        :Params:
        
            - **vari**: Input variable.
            - **weidstfracs**, optional: Specify how to apply weights computed using 
              destination fractions. Divide if <0, else multiply. It is not used with
              linear and match methods.
            - **check_mask**, optional: Mask point using masked data (the algo 
              interpolate the mask and check level 0.999). MUST BE IMPROVED!
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


######################################################################
######################################################################
from misc import axis1d_from_bounds, get_xy, isgrid, t2uvgrids, get_grid, \
    set_grid, bounds1d, bounds2d, get_axis, \
    meshgrid, create_grid, resol, meshcells, curv2rect, merge_axis_slice, \
    get_axis_slices, get_axis, transect_specs, create_axes2d
from .. import axes as A
#from ...misc.axml.base import Base as XMLBase
from ...misc.misc import cp_atts, intersect, kwfilter, get_atts, set_atts, closeto
from ...misc.atime import are_same_units, ch_units
import masking
from basemap import get_proj
#from kriging import krig as _krig2d

