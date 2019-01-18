# -*- coding: utf8 -*-
"""
Grid utilities.

It deals with bounds, areas, interpolations...

See: :ref:`user.tut.misc.grid`.
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
from __future__ import absolute_import
import operator
from gc import collect
import six
from six.moves import range

import numpy as N
import MV2
import cdms2
from numpy import ma as MA
from cdms2.coord import FileAxis2D
from cdms2.hgrid import AbstractCurveGrid, TransientCurveGrid
from cdms2.grid import AbstractGrid


from vacumm import VACUMMError, vacumm_warn, VACUMMWarning
from .misc import numod, cp_atts
from .units import deg2m
from .axes import (islat, islon, isaxis, create_lon, create_lat, isdep,
                   create_axes2d)
from .geo import get_distances, haversine  # compat


cdms = cdms2
MV = MV2

__all__ = ['isoslice', 'isgrid', 'get_resolution',
           'get_closest', 'get_closest_depth', 'get_geo_area',
           'bounds2d', 'bounds1d', 'get_grid', 'get_grid_axes', "gridsel",
           'grid2d', 'var2d', 'axis1d_from_bounds', 'get_xy',
           'deg2xy', 'set_grid', 'check_xy_shape',
           'meshgrid', 'meshbounds', 'meshweights',
           't2uvgrids', 'meshcells', 'cells2grid',
           'curv_grid', 'bounds2mesh', 'resol',
           'create_grid', 'rotate_grid', 'isregular',
           'monotonic', 'transect_specs',
           'depth2dz', 'isdepthup', 'makedepthup',
           'makealtitudeup', 'dz2depth', 'get_axis_slices',
           'xextend', 'xshift', 'curv2rect',  'isrect',
           'create_grid2d', 'create_var2d',
           'merge_axis_slice', 'merge_axis_slices',
           'get_zdim', 'coord2slice', 'mask2ind',
           'varsel',
           'clone_grid', 'are_same_grids', 'get_axis']


def isoslice(var, prop, isoval=0, axis=0, masking=True):
    """
    result = isoslice(variable,property[,isoval=0])

    result is a a projection of variable at property == isoval in the first
    nonsingleton dimension.  In the case when there is more than one zero
    crossing, the results are averaged.

    Example
    -------
    Assume two three dimensional variable, s (salt), z (depth), and
    u (velicity), all on the same 3D grid.  x and y are the horizontal
    positions, with the same horizontal dimensions as z (the 3D depth
    field).  Here, assume the vertical dimension, that will be projected,
    is the first

    >>> s_at_m5  = isoslice(s,z,-5);        # s at z == -5
    >>> h_at_s30 = isoslice(z,s,30);       # z at s == 30
    >>> u_at_s30 = isoslice(u,s,30);       # u at s == 30

    Get from OCTANT library : https://github.com/hetland/octant
    """
    if (len(var.squeeze().shape) <= 2):
        raise ValueError('variable must have at least two dimensions')
    if not prop.shape == var.shape:
        raise ValueError('dimension of var and prop must be identical')
    var = var.swapaxes(0, axis)
    prop = prop.swapaxes(0, axis)
    prop = prop-isoval
    sz = N.shape(var)
    var = var.reshape(sz[0], -1)
    prop = prop.reshape(sz[0], -1)
    # find zero-crossings (zc == 1)
    zc = N.where((prop[:-1, :]*prop[1:, :]) < 0.0, 1.0, 0.0)
    varl = var[:-1, :]*zc
    varh = var[1:, :]*zc
    propl = prop[:-1, :]*zc
    proph = prop[1:, :]*zc
    result = varl - propl*(varh-varl)/(proph-propl)
    result = N.where(zc == 1., result, 0.)
    szc = zc.sum(axis=0)
    szc = N.where(szc == 0., 1, szc)
    result = result.sum(axis=0)/szc
    if masking:
        result = N.ma.masked_where(zc.sum(axis=0) == 0, result)
        if all(result.mask):
            raise VACUMMWarning('property==%f out of range (%f, %f)' %
                                (isoval, (prop+isoval).min(), (prop+isoval).max()))
    result = result.reshape(sz[1:])
    return(result)


def isgrid(gg, curv=None):
    """Check if gg is a grid

    Parameters
    ----------
    gg: A cdms grid.
    curv: optional
        If True, restrict to curvilinear grids.
    """
    if curv:
        return isinstance(gg, AbstractCurveGrid)
    if curv is False:
        return isinstance(gg, AbstractGrid) and not isinstance(gg, AbstractCurveGrid)
    return isinstance(gg, AbstractGrid)


def get_resolution(mygrid, lon_range=None, lat_range=None):
    """Get the mean resolution of a grid

    .. warning:: Deprecated: please use :meth:`resol` instead.

    """

    # Get data
    if isinstance(mygrid, tuple):
        if len(mygrid) == 2:
            lon, lat = mygrid
            mask = None
        else:
            lon, lat, mask = mygrid
    else:
        lon = mygrid.getLongitude().getValue()
        lat = mygrid.getLatitude().getValue()
        mask = mygrid.getMask()

    # Work on 2D arrays
    if lon.shape != lat.shape:
        lon, lat = N.meshgrid(lon, lat)

    # Initial mask
    if mask is MA.nomask:
        lon[:] = N.where(mask, 999., lon)
        lat[:] = N.where(mask, 999., lat)

    # Geographic selection mask
    if lon_range is not None:
        mask_lon = N.logical_and(N.greater_equal(
            lon, lon_range[0]), N.less_equal(lon, lon_range[1]))
        if mask is MA.nomask:
            mask = mask_lon
        else:
            mask = N.logical_and(mask, mask_lon)
    if lat_range is not None:
        mask_lat = N.logical_and(N.greater_equal(
            lat, lat_range[0]), N.less_equal(lat, lat_range[1]))
        if mask is MA.nomask:
            mask = mask_lat
        else:
            mask = N.logical_and(mask, mask_lon)
        lon[:] = N.where(mask, 999., lon)
        lat[:] = N.where(mask, 999., lat)

    # Select valid data
    if mask is MA.nomask:
        lon = N.compress((1-mask).ravel(), lon.ravel())
        lat = N.compress((1-mask).ravel(), lat.ravel())
    else:
        lon = lon.ravel()
        lat = lat.ravel()

    # Compute areas
    xx = deg2m(lon, lat=lat)
    yy = deg2m(lat)
    dx = xx[1:]-xx[:-1]
    dy = yy[1:]-yy[:-1]
    ds = dx*dy
    del xx, yy, dx, dy
    collect()

    # Return mean resolution
    return N.sqrt(ds).mean()


def get_closest(xx, yy, xp, yp, distmode=None, mask=None, gridded=True, **kwargs):
    """Find the closest unmasked point on a grid and return indices

    Parameters
    ----------
    xx:
        1D or 2D X axis, or random positions.
    yy:
        1D or 2D Y axis, or random positions.
    xp:
        X position of the point (float)
    yp:
        Y //
    geo: optional
        If True, force to treat the grid as geographical,
        thus convert coordinates to meters using
          :func:`~vacumm.misc.phys.units.deg2m`.
    gridded: optional
        Treat input as gridded points,
        otherwise treat them as random points (flatten ``xx`` and ``yy``).

    Returns
    -------
    tuple, int

        - ``(i,j)`` 2-element tuple of indices along y and x for gridded data,
        - OR ``i`` if not gridded.
    """
    # Longitude/Latitude grid?
    proj = kwargs.pop('proj', None)
    proj = kwargs.pop('geo', proj)
    if proj:
        distmode = True
    elif distmode is None:
        distmode = 'haversine' if islon(xx) or islat(yy) else 'direct'

    # Gridded
    xx = N.asarray(xx[:])
    yy = N.asarray(yy[:])
    if gridded:
        xx, yy = meshgrid(xx, yy, copy=0)
    else:
        xx = xx.ravel()
        yy = yy.ravel()

    # Compute the distance to all grid points
    if mask is not None:
        xx = MA.array(xx, mask=mask, copy=False)
        yy = MA.array(yy, mask=mask, copy=False)

    # Distances
    dd = get_distances(xp, yp, xx, yy, mode=distmode)

    # Find the smallest distance
    imin = N.ma.argmin(dd)
    if gridded:
        return N.unravel_index(imin, dd.shape)[::-1]
    return imin


def get_closest_depth(zz, zp, mask=None, **kwargs):
    """Find the closest unmasked point on a depth vector and return indices

    Parameters
    ----------
    zz:
        1D Z axis, or random positions.
    zp:
        Z position of the point (float)

    Returns
    -------
    tuple, int

        - ``(i,)`` 1-element tuple of indices along z for gridded data,
        - OR ``i`` if not gridded.
    """

    # Gridded
    zz = N.asarray(zz[:])

    # Compute the distance to all grid points
    if mask is not None:
        zz = MA.array(zz, mask=mask, copy=False)
    dd = N.abs(zz - zp)

    # Find the smallest distance
    imin = N.ma.argmin(dd)
    return imin


def get_geo_area(grid, mask=None):
    """Compute cell areas on the a regular geographical grid.

    Parameters
    ----------
    grid: cdms2 grid
    mask: array
        Force the use of this mask

    Returns
    -------
    array
        2D array of areas in m^2

    .. TODO:: get_geo_area: treat 2D grid + use standard projection
    """

    if mask is None and grid.mask is not False:
        mask = grid.getMask()

    latBounds, lonBounds = grid.getBounds()
    lat = grid.getLatitude()

    latWeights = deg2m(latBounds[:, 1] - latBounds[:, 0], lat=lat)
    lonWeights = deg2m(lonBounds[:, 1] - lonBounds[:, 0])

    area = N.outerproduct(latWeights, lonWeights)

    if mask.any():
        area = MA.array(area, mask=mask)

    return area


def isregular(axis, tol=.05, iaxis=None, dx=None):
    """Check is an 1S or 2D axis is regular

    Parameters
    ----------
    axis:
        A :mod:`cdms2` axis or :mod:`numpy` array
    tol: optional
        Relative tolerance
    iaxis: optional
        On which direction to operate for 2D axis

    Example
    -------
    >>> isregular(lon1d, dx=2., tol=.01)
    >>> isregular(lon2d, iaxis=1)
    """
    cdaxis = isaxis(axis)
    if cdaxis:
        xx = axis.getValue()
    else:
        xx = axis[:]
    if iaxis is None:
        if xx.ndim == 1 or not cdaxis:
            iaxis = 0
        else:
            iaxis = int(islon(axis))
    thisdx = N.diff(xx, axis=iaxis)
    dxx = N.abs(N.diff(thisdx, axis=iaxis))
    if dx is None:
        dx = N.median(thisdx)
    return not ((dxx/dx) % 1. > tol).any()


def bounds1d(xx, cellwidth=None):
    """Compute bounds on a linear sequence or axis.
    It is based on :func:`cdms2.genGenericBounds` of CDAT.

    Example
    -------
    >>> bounds1d(xx).shape
    360, 2

    Returns
    -------
    array(nx,2)
    """

    try:
        xx = xx.getValues().astype('d')
    except:
        if not isinstance(xx, N.ndarray):
            xx = N.asarray(xx[:], 'd')
    M = numod(xx)

    if len(xx) > 1:
        leftPoint = M.array([1.5*xx[0]-0.5*xx[1]])
        midArray = (xx[0:-1]+xx[1:])/2.0
        rightPoint = M.array([1.5*xx[-1]-0.5*xx[-2]])
        bnds = N.concatenate((leftPoint, midArray, rightPoint))
    else:
        assert cellwidth is not None
        delta = cellwidth/2.0
        bnds = M.array([xx[0]-delta, xx[0]+delta])

    # Transform to (n,2) array
    retbnds = M.zeros((len(xx), 2), 'd')
    retbnds[:, 0] = bnds[:-1]
    retbnds[:, 1] = bnds[1:]

    return retbnds


def bounds2d(*xyaxes):
    """Compute bounds on a rectangular grid.

    Parameters
    ----------
    xyaxes:
        2D arrays(ny,nx) (or 2x1D arrays)

    Example
    -------
    >>> xx_bounds,yy_bounds = bounds2d(xx,yy)

    Returns
    -------
    xx(ny,nx,4)
    ...
    """

    results = []
    half = N.array(0.5, 'd')
    quart = N.array(0.25, 'd')

    if len(xyaxes) == 2:
        xyaxes = meshgrid(*xyaxes)

    for xy in xyaxes:

        if hasattr(xy, 'getValue'):
            xy = xy.getValue()
        assert xy.ndim == 2, 'You must pass 2d arrays as argument'
        ny, nx = xy.shape

        M = numod(xy)

        bounds = M.zeros((ny, nx, 4), 'd')

        inside = quart * (xy[0:-1, 0:-1] + xy[0:-1, 1:] +
                          xy[1:, 0:-1] + xy[1:, 1:])
        left = half * (xy[0:-1, 0] - half * (xy[0:-1, 1] - xy[0:-1, 0]) +
                       xy[1:, 0] - half * (xy[1:, 1] - xy[1:, 0]))
        right = half * (xy[0:-1, -1] + half * (xy[0:-1, -1] - xy[0:-1, -2]) +
                        xy[1:, -1] + half * (xy[1:, -1] - xy[1:, -2]))
        bot = half * (xy[0, 0:-1] - half * (xy[1, 0:-1] - xy[0, 0:-1]) +
                      xy[0, 1:] - half * (xy[1, 1:] - xy[0, 1:]))
        top = half * (xy[-1, 0:-1] + half * (xy[-1, 0:-1] - xy[-2, 0:-1]) +
                      xy[-1, 1:] + half * (xy[-1, 1:] - xy[-2, 1:]))
        left_bot = xy[0, 0] + half * (xy[0, 0] - xy[1, 1])
        right_bot = xy[0, -1] + half * (xy[0, -1] - xy[1, -2])
        right_top = xy[-1, -1] + half * (xy[-1, -1] - xy[-2, -2])
        left_top = xy[-1, 0] + half * (xy[-1, 0] - xy[-2, 1])

        bounds[1:, 1:, 0] = inside
        bounds[1:, 0:-1, 1] = inside
        bounds[0:-1, 0:-1, 2] = inside
        bounds[0:-1, 1:, 3] = inside

        bounds[1:, 0, 0] = left
        bounds[0:-1, 0, 3] = left
        bounds[1:, -1, 1] = right
        bounds[0:-1, -1, 2] = right
        bounds[0, 1:, 0] = bot
        bounds[0, 0:-1, 1] = bot
        bounds[-1, 0:-1, 2] = top
        bounds[-1, 1:, 3] = top

        bounds[0, 0, 0] = left_bot
        bounds[0, -1, 1] = right_bot
        bounds[-1, -1, 2] = right_top
        bounds[-1, 0, 3] = left_top

        results.append(bounds)

    if len(results) == 1:
        return results[0]
    else:
        return results


def meshgrid(x, y, copy=1):
    """Convert pure numeric x/y axes to 2d arrays of same shape.
    It works like :func:`numpy.meshgrid` but is more flexible.

    Parameters
    ----------

    x:
        1D or 2D array
    x:
        2D or 2D array

    Return
    ------
    array(ny,nx)
        x positons
    array(ny,nx)
        y positions
    """

    # Be sure to have numpy arrays and make a copy
    x = _cdat2num_(x)
    Mx = numod(x)
    My = numod(y)
    if N.ma.isMA(x):
        if copy:
            x = x.copy()
        xmask = N.ma.getmaskarray(x)
    else:
        x = N.array(x[:], copy=copy)
        xmask = None
    y = _cdat2num_(y)
    if N.ma.isMA(y):
        if copy:
            y = y.copy()
        ymask = N.ma.getmaskarray(y)
    else:
        y = N.array(y[:], copy=copy)
        ymask = None

    # All is ok
    if x.ndim == 2 and x.shape == y.shape:
        return x, y

    # From regular grid
    if x.ndim == 1 and y.ndim == 1:
        x, y = N.meshgrid(x, y)
        if xmask is not None or ymask is not None:
            if xmask is None:
                xmask = N.zeros(len(x), '?')
            elif ymask is None:
                ymask = N.zeros(len(y), '?')
            xmask, ymask = N.meshgrid(xmask, ymask)
            x = N.ma.masked_where(xmask, x, copy=False)
            y = N.ma.masked_where(ymask, y, copy=False)
        return x, y

    # From x 1d
    if x.ndim == 1 and y.ndim == 2 and y.shape[1] == x.shape[0]:
        xx = Mx.resize(x, y.shape)
        return xx, y

    # From y 1d
    if y.ndim == 1 and x.ndim == 2 and x.shape[0] == y.shape[0]:
        yy = My.resize(y, x.shape[::-1]).transpose()
        return x, yy

    raise VACUMMWarning("Unable to safely convert to 2D axes")


def bounds2mesh(xb, yb=None):
    """Convert 2D 4-corners cell bounds arrays (results from :func:`bounds2d`) to 2D arrays

    Parameters
    ----------

    xb:
        array(ny,nx,4) or array(nx,2)
    yb:
        None or array(ny,nx,4) or array(ny,2)

    Return
    ------
    array(ny+1,nx+1)
        x mesh bounds
    array(ny+1,nx+1)
        y mesh bounds
    """
    res = ()
    for xyb in xb, yb:
        if not isinstance(xyb, N.ndarray):
            xyb = N.asarray(xyb)
        M = numod(xyb)

        # 1D
        if xyb.ndim == 2:
            xy = M.zeros((len(xyb)+1, ), dtype=xb.dtype)
            xy[:-1] = xyb[:, 0]
            xy[-1] = xyb[-1, 1]

        # 1D
        else:
            ny, nx, nb = xyb.shape
            # - center + lower + left
            xy = M.zeros((ny+1, nx+1), dtype=xyb.dtype)
            xy[:-1, :-1] = xyb[:, :, 0]
            # - right
            xy[:-1, -1] = xyb[:, -1, 1]
            # - top
            xy[-1, :-1] = xyb[-1, :, 3]
            # - top-left
            xy[-1, -1] = xyb[-1, -1, 2]

        if yb is None:
            return xy
        res += xy,
    return res


def meshcells(x, y=None, xwidth=None, ywidth=None):
    """Return a 1D or 2D array corresponding the cell corners

    Parameters
    ----------
    x:
        1D or 2D array
    y:
        1D or 2D array

    Returns
    -------
    array(nx+1)[,array(ny+1,nx+1)]
        x[,y] mesh cells arrays

    TODO: Add full support to meshcells for xwidth and ywidth params
    """
    if not isinstance(x, N.ndarray):
        x = N.asarray(x)
    if y is None and len(x.shape) == 1:
        return bounds2mesh(bounds1d(x, xwidth))
    if y is not None:
        if not isinstance(y, N.ndarray):
            y = N.asarray(y)
        xy = meshgrid(x, y)
        return bounds2mesh(*bounds2d(*xy))
    else:
        return bounds2mesh(bounds2d(x))


def meshweights(x, y=None, distmode='haversine', axis=-1, proj=None):
    """Return a 1D or 2D array corresponding the cell weights

    Parameters
    ----------
    x:
        1D or 2D array
    y:
        1D or 2D array
    proj: optional
        Geographic projection before computing weights
        (for 2D case only with both x and y)

        - ``True``: use default mercator projection (see :func:`~vacumm.misc.grid.basemap.merc`),
        - else, directly apply projection.

    Returns
    -------
    array(nx) OR array(ny,nx)

    """
    if not isinstance(x, N.ndarray):
        x = N.asarray(x)

    # Without y
    if y is None:

        if axis < 0:
            axis += x.ndim
        dx = x.copy()
        M = numod(x)
        dd = M.diff(x, axis=axis)
        sl = [slice(None)]*x.ndim

        def ss(*sss):
            s = list(sl)
            s[axis] = slice(*sss) if len(sss) != 1 else sss[0]
            return tuple(s)
        dx[ss(1, -1)] = .5*(dd[ss(None, -1)]+dd[ss(1, None)])
        dx[ss(0)] = dd[ss(0)]
        dx[ss(-1)] = dd[ss(-1)]
        del dd
        return dx

    # With y
    if not isinstance(y, N.ndarray):
        y = N.asarray(y)

    # 1D
    if x.ndim == 1:
        return meshweights(x, axis=axis), meshweights(y, axis=axis)

    # 2D
    if distmode is None:
        distmode = 'haversine' if islon(x) or islat(y) else 'direct'
    xxb, yyb = meshcells(x, y)
    dx = get_distances(xxb[:, :-1], yyb[:, :-1], xxb[:, 1:], yyb[:, 1:],
                       distmode=distmode, pairwise=True)
    dx = 0.5*(dx[1:]+dx[:-1])
    dy = get_distances(xxb[:-1, :], yyb[:-1, :], xxb[1:, :], yyb[1:, :],
                       distmode=distmode, pairwise=True)
    dy = 0.5*(dy[:, 1:]+dy[:, :-1])
    return dx*dy


def cells2grid(xxb, yyb):
    """Convert a grid of cells (grid of corners, like results from :func:`meshcells`) to a grid (grid of centers)

    Parameters
    ----------
    xxb:
        X of corners (ny+1,nx+1)
    yyb:
        Y of corners (ny+1,nx+1)

    Returns
    -------
    array(ny,nx), array(ny,nx)
        x and y positions
    """
    M = numod(xxb)
    ny, nx = xxb.shape
    xxyy = [M.zeros((ny+1, nx+1), dtype=xxb.dtype),
            M.zeros((ny+1, nx+1), dtype=yyb.dtype)]
    bb = (xxb, yyb)
    ij2ic = M.reshape([0, 1, 3, 2], (2, 2))
    for ib in 0, 1:
        # Center
        xxyy[ib][1:-1, 1:-1] = 0.25 * (bb[ib][1:, 1:, 0] + bb[ib][:-1, 1:, 1] +
                                       bb[ib][:-1, :-1, 2] + bb[ib][1:, :-1, 3])
        # Sides
        for i in 0, -1:
            xxyy[ib][i, 1:-1] = 0.5 * (bb[ib][i, 1:, ij2ic[i, 0]] +
                                       bb[ib][i, :-1, ij2ic[i, 1]])  # Bottom/top
            xxyy[ib][1:-1, i] = 0.5 * (bb[ib][1:, i, ij2ic[0, i]] +
                                       bb[ib][:-1, i, ij2ic[1, i]])  # Left/right
        # Corners
        for j in 0, -1:
            for i in 0, -1:
                xxyy[ib][j, i] = bb[ib][j, i, ij2ic[j, i]]
    return xxyy[0], xxyy[1]


def meshbounds(*args, **kwargs):
    """A shortcut to :func:`meshcells`"""
    return meshcells(*args, **kwargs)


def coord2slice(gg, lon=None, lat=None, mode='slice', mask=None, assubmask=True,
                squeeze=False, **kwargs):
    """Convert from coordinates to indices on grid or axes

    Parameters
    ----------
    gg:
        A grid, a cmds variable with a grid,
        a tuple of axes or an axis.
    lon/lat:
        Interval of coordinates or slice. If single coordinates, the
        interval becomes for instance ``(lon,lon,'ccb')``.
    mode: optional
        Output mode. You can pass only the first letter.

        - ``"slice"``: Return a :class:`slice` object.
        - ``"indices"``: Return a 3-element tuple as in a :class:`slice` context
          ``(start,stop,step)``.
        - ``"array"``: Return an array of valid indices of shape ``(2,nvalid)``,
          where ``out[0]`` gives the X indices and ``out[1]`` gives the Y indices.

    mask: optional
        Also use this mask to get indices
        (for 2D axes only).
    squeeze: optional
        If ``asslice`` is False ad
    asslice: optional
        DEPRECATED. Use ``mode='slice'`` instead.

    Return
    ------
    tuple(2,nvalid) OR (islice, jslice, mask) OR slice
        In array mode, it always returns an array of indices of shape ``(2,nvalid)``,
        possibly empty if not intersection if found.

        Here ``i/jslice`` is a slice or a tuple of indices,
        depending on asslice. It can also be ``None`` if
        not intersection are found.

        If ``gg`` is grid, a tuple of axes or a 2D axis:
        ``islice, jslice, mask``.
        In case of a tuple of 1D axes, ``mask`` is ``None`` because not relevant.

        If ``gg`` is a 1D axis: ``ijslice``, just as the
        :meth:`mapIntervaExt`` method of 1D axes.

    Examples
    --------
    >>> ijk = coord2slice(lon1d, lon=(lon0,lon1), mode='i')
    >>> xijk, yijk, mask = coord2slice((lon1d,lat1d),
        lon=(lon0,lon1), lat=slice(0,3), mode='i')
    >>> islice, jslice, mask = coord2slice(gridcurv, lat=(lat0,lat0,'ccb'), mode='s')
    >>> xijk, yijk, mask = coord2slice(lon2d, lon=(lon0,lon1), mode='i')
    >>> ij = coord2slice(grid, lon, lat, mode='a')

    """
    # Output mode
    if 'asslice' in kwargs:
        mode = 'slice' if kwargs['asslice'] else 'indices'
    if mode is None:
        mode = 'slice'
    else:
        mode = str(mode).lower()
        if mode[:1] not in 'ias':
            raise VACUMMError(
                'Wrong mode: %s. Choose one of: slice, indices, array' % mode)
    if mode[0] == 'a':
        assubmask = True

    # Format selector specs
    if lon == ':' or lon == slice(None):
        lon = None
    if lon == ':' or lon == slice(None):
        lon = None
    elif isinstance(lon, tuple):
        if len(lon) == 2:
            lon += 'cc',
        if len(lon) == 3:
            lon += None,
    elif isinstance(lon, (int, float)):
        lon = (lon, lon, 'cob')
    if lat == ':' or lat == slice(None):
        lat = None
    elif isinstance(lat, tuple):
        if len(lat) == 2:
            lat += 'cc',
        if len(lat) == 3:
            lat += None,
    elif isinstance(lat, (int, float)):
        lat = (lat, lat, 'cob')

    # CDMS variable
    if cdms2.isVariable(gg) and not isaxis(gg):

        gg = get_grid(gg)

    # Grid
    if isgrid(gg):

        # Get axes
        xx, yy = get_xy(gg)

        # 2D
        if len(xx.shape) == 2:

            # As axes
            if not isaxis(xx) or not isaxis(yy):
                xx, yy = create_axes2d(xx, yy)

            # Get indices
            xii, xjj, mask = coord2slice(
                xx, lon=lon, mode='indices', mask=mask, assubmask=False)
            if xii is None:
                return None, None, mask
            yii, yjj, mask = coord2slice(
                yy, lat=lat, mode='indices', mask=mask, assubmask=False)
            if yii is None:
                return None, None, mask

            # Merge
            imin = max(xii[0], yii[0])
            imax = min(xii[1], yii[1])
            jmin = max(xjj[0], yjj[0])
            jmax = min(xjj[1], yjj[1])
            xijk = (imin, imax, 1)  # max(xii[2],yii[2]) instead of 1?
            yijk = (jmin, jmax, 1)
            sxijk = slice(*xijk)
            syijk = slice(*yijk)

            # Submask and slices
            if mode[0] == 's':
                xijk = sxijk
                yijk = syijk
            if assubmask:
                mask = mask[syijk, sxijk]
            if mode[0] == 'a':  # array of valid indices
                iijj = N.indices(xx.shape)[:, syijk, sxijk]
                iijj = iijj.reshape(2, -1)
                return iijj[:, ~mask.ravel()]
            return xijk, yijk, mask

        # 1D

        # - as axes
        if not isaxis(xx):
            xx = create_lon(xx)
        if not isaxis(yy):
            yy = create_lon(yy)

        # - indices
        cmode = mode if mode[0] != 'a' else 'slice'
        xijk = coord2slice(xx, lon=lon, mode=cmode)
        if xijk is None:
            if mode[0] == 'a':
                return N.zeros((2, 0))
            return None, None, None
        yijk = coord2slice(yy, lat=lat, mode=cmode)
        if yijk is None:
            if mode[0] == 'a':
                return N.zeros((2, 0))
            return None, None, None

        # - mask
        mask = None
#        mask = N.ones((len(yy), len(xx)))
#        mask[slice(*yijk), slice(*xijk)] = False

        # - return
        if mode == 'a':
            iijj = N.indices((len(yy), len(xx)))[:, yijk, xijk]
            iijj.shape = 2, -1
            return iijj
        return xijk, yijk, mask

    # Single Axis

    # - selector
    if not isaxis(gg):
        raise VACUMMError(
            "If gg is not a grid or an tuple, it must be an axis")
    axis = int(gg.isLongitude())
    sel, osel = (lat, lon)[::1-2*axis]
    if axis:
        minval, maxval = -720., 720
        modeval = 'mask'
    else:
        minval, maxval = -90, 90
        modeval = 'clip'

    # - 2D
    if len(gg.shape) == 2:

        # Selection with provided coordinates
        if isinstance(sel, slice):
            ijk = sel.indices(gg.shape[axis])
            oijk = (0, gg.shape[1-axis], 1)
            mask = N.ones(gg.shape, '?')
#            mask[(ijk, None)[::1-2*axis]] = False
            mask[(sel, slice(None))[::1-2*axis]] = False
        else:
            xijk, yijk, mask = _coord2ind2d_(gg, sel, mask=mask,
                                             minval=minval, maxval=maxval, modeval=modeval)
            if xijk is None:
                return None, None, mask
            ijk, oijk = (yijk, xijk)[::1-2*axis]

        # Selection along other coordinates with slice
        if isinstance(osel, slice):
            i, j, k = osel.indices(gg.shape[1-axis])
            oijk = (max(i, oijk[0]), min(j, oijk[1]), 1)
            bad = N.ones(gg.shape, '?')
            bad[(slice(None), slice(*oijk))[::1-2*axis]] = False
            mask |= bad
            del bad

        # Submask and slices
        xijk, yijk = (oijk, ijk)[::1-2*axis]
        if mode[0] == 's':
            xijk = slice(*xijk)
            yijk = slice(*yijk)
            if assubmask:
                mask = mask[yijk, xijk]
        elif assubmask:
            mask = mask[slice(*yijk), slice(*xijk)]
        if mode[0] == 'a':  # array of valid indices
            iijj = N.indices(xx.shape)[:, yijk, xijk]
            iijj.shape = 2, -1
            return iijj[:, ~mask.ravel()]
        return xijk, yijk, mask

    # - 1D
    if isinstance(sel, slice):
        if mode[0] == 's':
            return sel
        ijk = sel.indices(len(gg))
        if mode[0] == 'i':
            return ijk
    else:
        ijk = gg.mapIntervalExt(sel)
    if ijk is None:
        if mode[0] == 'a':
            return N.array((0,))
        return
    if mode[0] == 's':
        ijk = slice(*ijk)
    elif mode[0] == 'a':
        return N.arange(len(gg))[ijk]
    return ijk


def mask2ind(mask, extended=False):
    """Convert a mask to min and max valid indices

    Parameters
    ----------
    mask:
        N-dimensional array of logical where
        valid points are at False.
    extended: optional
        Extend indices. Accepts
        ``True``, ``False``, ``int``, ``(intmin, intmax)``.

    Return
    ------
    array(ndim,2)
        array of (min,max) for each dimension
    """
    inds = N.indices(mask.shape)
    ndim = inds.shape[0]
    inds.shape = ndim, -1
    mask = mask.ravel()
    inds = inds[:, ~mask]
    mm = N.empty((ndim, 2), 'l')
    if extended == 0:
        extended = False
    if extended is True:
        extended = 1
    elif extended is False:
        extended = 0
    if isinstance(extended, int):
        extended = (extended, extended)
    mm[:, 0] = inds.min(axis=1) - extended[0]
    mm[:, 1] = inds.max(axis=1) + extended[1]
    del inds
    mm[mm < 0] = 0
    return mm


def _mask4tomask_(mask4):
    """Convert mask at corners to mask at center of cell using the AND operator"""
    mask = mask4[:-1, :-1].copy()
    mask &= mask4[:-1, 1:]
    mask &= mask4[1:, :-1]
    mask &= mask4[1:, 1:]
    return mask


def _coord2ind2d_(cc2d, sel, mask=None, maskonly=False,
                  minval=None, maxval=None, modeval='mask'):
    """Convert lon/lat 2d coord to index respectively with lon/lat selection

    Parameters
    ----------
    cc2d: array(ny,nx):
        Coordinates.
    sel:
        Coordinates interval.

    Returns
    -------
    (None, None, mask) OR ((imin,imax,istep),(jmin,jmax,jstep),mask(ny,nx))

    """

    #  No selection
    if sel is None:
        nj, ni = cc2d.shape
        if mask is None:
            return [(0, nj, 1), (0, ni, 1), N.zeros((nj, ni), '?')]
        [(jmin, jmax), (imin, imax)] = mask2ind(mask)
        return [(imin, imax+1, 1), (jmin, jmax+1, 1), mask]

    # Selection specs
    if len(sel) == 2:
        sel = sel+('cc', )
    b = sel[2]
    e = b[2:]
    if isinstance(cc2d, FileAxis2D):
        cc2d = cc2d.clone()
    if cdms2.isVariable(cc2d):
        cc2d = cc2d.asma()

    # Limits
    if minval is not None or maxval is not None:
        if minval is None:
            minval = -N.inf
        if maxval is None:
            maxval = N.inf
        if modeval == 'mask':
            cc2d = N.ma.masked_outside(cc2d, minval, maxval, False)
        else:
            cc2d = N.ma.clip(cc2d, minval, maxval)

    # According to bounds of cells
    if 'b' in e:

        # Bounds become coordinates
        cc2db = meshcells(cc2d)

        # Remove 'b' from interval specifications
#        selb = sel[:2]+(b[:2],)   # no 'b' #FIXME: what is this selb for?

        # Select min first
        selmin = (sel[0], N.inf)+(b[:2],)
        maskmin = _coord2ind2d_(cc2db, selmin, mask=None, maskonly=True)
        if maskmin is None:
            return mask if maskonly else (None, None, mask)
        xymask = _mask4tomask_(maskmin)
#        del maskmin

        # Select max
        selmax = (-N.inf, sel[1])+(b[:2],)
        maskmax = _coord2ind2d_(cc2db, selmax, mask=None, maskonly=True)
        if maskmax is None:
            return mask if maskonly else (None, None, mask)
        xymask |= _mask4tomask_(maskmax)
#        del maskmax

        # Merge
        if mask is not None:
            xymask |= mask
        if maskonly or xymask.all():
            return xymask if maskonly else (None, None, xymask)
        [(jmin, jmax), (imin, imax)] = mask2ind(xymask, 'n' in e)
        return [(imin, imax+1, 1), (jmin, jmax+1, 1), xymask]

    # Operators
    Less = N.ma.less if b[0] == 'c' else N.ma.less_equal
    Greater = N.ma.greater if b[1] == 'c' else N.ma.greater_equal

    # Mask
    mask = False if mask is None else mask.copy()
    mask |= Less(cc2d, sel[0]).filled(True)
    mask |= Greater(cc2d, sel[1]).filled(True)
    if maskonly or mask.all():
        return mask if maskonly else (None, None, mask)

    # Indices
    [(jmin, jmax), (imin, imax)] = mask2ind(mask, 'e' in e)
    return [(imin, imax+1, 1), (jmin, jmax+1, 1), mask]


def get_grid(gg, geo=True, intercept=False, strict=False):
    """Get a cdms grid from gg

    Parameters
    ----------

    g:
        A cdms variable with a grid OR cdms grid OR a tuple like (xx,yy) where xx and yy are numpy arrays or cdms axes
    geo: optional
        Output grid must have longitude and latitude
    intercept: optional
        Raise an error in case of problem

    Examples
    --------
    >>> grid = get_grid((lon, lat))
    >>> grid = get_grid(var)
    >>> grid = get_grid(grid) # does nothing

    Return
    ------
    cdms2 grid

    """
    # From a grid
    if isgrid(gg):
        gg.getAxis(-1).designateLongitude()
        gg.getAxis(-2).designateLatitude()
        return gg

    # From a variable
    if cdms.isVariable(gg):
        if gg.getGrid() is not None or strict:
            return gg.getGrid()
        if len(gg.shape) < 2:
            return
        xx = gg.getLongitude()
        yy = gg.getLatitude()
        if xx is None:
            xx = gg.getAxis(-1)
        if yy is None:
            yy = gg.getAxis(-2)
        gg = (xx, yy)

    # From a tuple of coordinates
    if isinstance(gg, tuple):

        xx, yy = gg

        # Check axis is compatible
        if isaxis(xx):
            if getattr(xx, 'axis', 'x').upper() != 'X':
                if not intercept:
                    return
                raise VACUMMError(
                    "Can't make a grid with %s axis as longitude" % xx.axis)
            xx.designateLongitude()
        if isaxis(yy):
            if getattr(yy, 'axis', 'y').upper() != 'Y':
                if not intercept:
                    return
                raise VACUMMError(
                    "Can't make a grid with %s axis as latitude" % yy.axis)
            yy.designateLatitude()

        # Create grid
        if (isinstance(xx,  tuple) or xx[:].ndim == 1) and \
                (isinstance(yy, tuple) or yy[:].ndim == 1):  # Rectangular

            return create_grid(xx, yy)

        elif xx[:].ndim == 2 or yy[:].ndim == 2:  # Curvilinear

            return curv_grid(xx, yy)

        else:  # Generic

            return cdms.createGenericGrid(yy[:].ravel(), xx[:].ravel())

    if intercept:
        raise VACUMMError('No way to guess the grid')


def set_grid(var, gg, axes=True, intercept=False):
    """Set a grid to a variable

    Parameters
    ----------
    var:
        A cdms variable with at least 2 dimensions
    gg:
        A cdms grid or tuple or axes (see :func:`get_grid`)
    axes: optional
        Set grid raw axes to the variable (gg.getAxis(i))
        which are different from real axes with curvilinear grids

    Return
    ------
    array
        The input variable
    """
    if not cdms2.isVariable(var):
        var = MV.asarray(var)
    gg = get_grid(gg, intercept=intercept)
    if gg is not None:
        for i in -2, -1:
            axis = gg.getAxis(i)
            setattr(axis, 'axis', ['Y', 'X'][i])
            var.setAxis(i, axis)
    var.setGrid(gg)
    return var


def get_grid_axes(gg, raw=False):
    """Get the (lat,lon) axes from a grid

    Parameters
    ----------
    gg:
        A cdms grid or a tuple or axes (see :func:`get_grid`)
    raw: optional
        If True, return raw axes which are different from real
        axes with curvilinear grids

    Return
    ------
    axis, axis
       loon and lat axis
    """
    gg = get_grid(gg)
    if raw:
        return gg.getAxis(0), gg.getAxis(1)
    return gg.getLatitude(), gg.getLongitude()


def curv_grid(xaxis, yaxis, xatts=None, yatts=None, id=None, mask=None,
              **kwargs):
    """Create a curvilinear 2D grid from 1D or 2D axes

    Parameters
    ----------
    xaxis:
        Numeric or formatted X axis
    yaxis:
        Numeric or formatted Y axis
    xatts: optional
        Attributes of X axis
    yatts: optional
        Attributes of Y axis
    id: optional
        Id of the grid
    **kwargs
        Other keywords are passed to :func:`create_axes2d`

    Example
    -------
    >>> curvgrid = curv_grid(lon2d, lat2d)

    Return
    ------
    cdms2 curved grid

    See also
    --------
    :func:`create_grid2d`

    """
    xaxis2d, yaxis2d = create_axes2d(xaxis, yaxis, xatts=xatts,
                                     yatts=yatts, numeric=False, **kwargs)
    if id is None:
        id = 'grid2d'
    return TransientCurveGrid(yaxis2d, xaxis2d, id=id, tempmask=mask)


def create_grid2d(*args, **kwargs):
    """Alias for :func:`curv_grid`"""
    return curv_grid(*args, **kwargs)


def grid2d(*args, **kwargs):
    """Alias for :func:`curv_grid`"""
    return create_grid2d(*args, **kwargs)


def create_var2d(var, xaxis=None, yaxis=None,
                 xatts=None, yatts=None, gid=None, copy=1,
                 lonid=None, latid=None, iid=None, jid=None, **kwargs):
    """Create 2D cdms variable with on a proper 2d curvilinear grid

    You should generally call this routine when you want to attach 2D axes
    to a variable. This may happen with netcdf that doesn't follow CF
    conventions, as in the example bellow.

    Parameters
    ----------
    var:
        Numeric or formatted X axis
    xaxis: optional
        Numeric or formatted X axis. Mandatory if var is not a cdms variable!
    yaxis: optional
        Numeric or formatted Y axis. Mandatory if var is not a cdms variable!
    atts: optional
        Attributes of the variable .
    xatts: optional
        Attributes of X axis.
    yatts: optional
        Attributes of Y axis.
    gid: optional
        Id of the grid.
        - All other keywords are passed to :func:`cdms.createVariable`

    Return
    ------
    :class:`MV2.array`
        A :mod:`MV2` array on a curvilinear grid.

    Example
    -------
    >>> f = cdms2.open('myfile.nc')
    >>> lon2d = f('X')
    >>> lat2d = f('Y')
    >>> sst = f('sst')
    >>> f.close()
    >>> sstc = create_var2d(sst, lon2d, lat2d)


    See also
    --------
    :func:`curv_grid` :func:`get_axis`

    """

    assert len(var.shape) >= 2, 'Input variable must have a at least 2d shape'

    if cdms.isVariable(var):
        if yaxis is None:
            yaxis = get_axis(var, 0)
        if xaxis is None:
            xaxis = get_axis(var, 1)
    var = MV2.array(var, copy=copy, **kwargs)

    mygrid2d = curv_grid(xaxis, yaxis, xatts=xatts, yatts=yatts, id=gid,
                         iid=iid, jid=jid, lonid=lonid)  # ,mask=MV.getmask(var))

    set_grid(var, mygrid2d)

    return var


def var2d(*args, **kwargs):
    """Alias for :func:`create_var2d`"""
    return create_var2d(*args, **kwargs)


def create_grid(lon, lat, mask=None, lonatts={}, latatts={}, curv=None, **kwargs):
    """Create a cdms rectangular or curvilinear grid from axes

    Parameters
    ----------

    lon:
        Array or axis of longitudes or any argument passed to :func:`~vacumm.misc.axes.create_lon`.
    lat:
        Array or axis of latitudes or any argument passed to :func:`~vacumm.misc.axes.create_lat`.
    mask: optional
        Grid mask.
    (lon/lat)atts: optional
        Attributes to set for axes.
    curv: optional

            - ``None``: If one axis is 2D, the grid will be curvilinear.
            - ``True|False``: Force the grid to be or not to be curvilinear.

    Return
    ------
    A :mod:`cdms2` grid object.

    Example
    -------

        >>> create_grid([1.,3.,5.,7], numpy.arange(45., 60., .5))
        >>> create_grid((.1, 8., 1.), (45., 60., .5))

    See also
    --------
    :func:`~vacumm.misc.axes.create_lon` :func:`~vacumm.misc.axes.create_lat`
               :func:`get_grid` :func:`set_grid`
    """

    # Arrays
    if isinstance(lon, list):
        lon = N.asarray(lon)
    if isinstance(lat, list):
        lat = N.asarray(lat)

    # Rectangular or curvilinear?
    if curv is None and not isinstance(lon, tuple) and not isinstance(lat, tuple):
        curv = lon[:].ndim != 1 or lat[:].ndim != 1

    if not curv:  # Rectangular
        # Get axes
        if not islon(lon):
            lon = create_lon(lon, **lonatts)
        if not islat(lat):
            lat = create_lat(lat, **latatts)
        # Create grid
        return cdms2.createRectGrid(lat, lon, mask=mask, **kwargs)

    else:  # Curvilinear
        return create_grid2d(lon, lat, xatts=lonatts, yatts=latatts, mask=mask, **kwargs)


def get_axis(gg, iaxis=0, geo=True, strict=False):
    """A robust way to get an axis from a variable or a grid

    This is a generic way to get a 1D and 2D axes: the only
    way to get longitude and latitude axes of a curvilinear grid are the
    :meth:`getLongitude` and :meth:`getLatitude` methods
    (you can't use :meth:`getAxis` with such grid).

    Parameters
    ----------
    gg:
        The variable or a grid (see :func:`get_grid`)
    iaxis: optional
        The axis number [default: 0]

    Examples
    --------
    >>> lon = get_axis(grid, -1)
    >>> lat = get_axis(var, 0)

    Returns
    -------
    cdms2 axis
    """

    grid = get_grid(gg, intercept=False, geo=geo, strict=strict)
    if grid is not None:
        getfrom = grid
    else:
        getfrom = gg
    ndim = len(getfrom.shape)
    if iaxis < 0:
        iaxis = ndim+iaxis
    if iaxis <= 2 and isgrid(grid, curv=True):
        if iaxis == 1:
            return getfrom.getLongitude()
        else:
            return getfrom.getLatitude()
    return getfrom.getAxis(iaxis)


def clone_grid(grid):
    """Clone a grid"""
    if grid is None:
        return
    id = grid.id

    lon = grid.getLongitude().clone()
    lat = grid.getLatitude().clone()
    if cdms2.isVariable(lon):
        iaxis = lon.getAxis(1).clone()
        jaxis = lon.getAxis(0).clone()
        lon.setAxisList([jaxis, iaxis])
        lat.setAxisList([jaxis, iaxis])

    grid = create_grid(lon, lat, curv=len(lon.shape) == 2)
    grid.id = id
    return grid


def gridsel(gg, lon=None, lat=None):
    """Extract a subregion from an axis or a grid

    Lat and lon generic selection are used for
    spatial selections.

    .. note:: Properly works on curved grids thanks to :func:`coord2slice`.

    Parameters
    ----------

    gg:
        1D or 2D cdms2 axis, or a tuple of them , or a cdms2 grid.
    lon/lat: optional
        A slice, or a tuple of coordinates, or ':'.

    Return
    ------

        An extraction in the same format as input format

    Examples
    --------

        >>> gg = gridsel(gg, lon=(-6,4), lat=slice(0,34))
        >>> lon = gridsel(lon, lon=(-5,8))
        >>> lon,lat = gridsel((lon,lat), lat=':', lon=(0,2,'ccb'))
    """

    # Format selector specs
    if lon == ':' or lon == slice(None):
        lon = None
    elif isinstance(lon, tuple):
        if len(lon) == 2:
            lon += 'cc',
        if len(lon) == 3:
            lon += None,
    if lat == ':' or lat == slice(None):
        lat = None
    elif isinstance(lat, tuple):
        if len(lat) == 2:
            lat += 'cc',
        if len(lat) == 3:
            lat += None,
    if lon is None and lat is None:
        gg

    # Format input to always treat a grid
    if isinstance(gg, tuple):  # (lon, lat) -> grid
        outtype = 'tuple'
        gg = get_grid(gg)
    elif isgrid(gg):  # grid
        outtype = 'grid'
    elif isaxis(gg):  # axis
        #        if len(gg.shape)==2 and (
        #            (lon is not None and not isinstance(lon, slice)) or
        #            (lat is not None and not isinstance(lat, slice))):
        #            raise VACUMMError("Can't apply geographical selection to a lonely 2D axis")
        if islon(gg, ro=True):
            outtype = 'lon'
        elif islat(gg, ro=True):
            outtype = 'lat'
        else:
            raise VACUMMError("Axis must be longitude or latitude")
    else:
        raise TypeError('Please provide a grid,  a tuple of axes or an axis')

    # Selection
    if len((gg if outtype[0] not in ['g', 't'] else gg.getLongitude()).shape) == 2:  # Curvilinear

        # Convert to slice
        islice, jslice, mask = coord2slice(gg, lon=lon, lat=lat)

        # Grid
        if outtype[0] in ['g', 't']:

            islice = islice.indices(gg.shape[1])
            jslice = jslice.indices(gg.shape[0])
            gg = gg.subSlice(jslice, islice)
            if mask.any():
                gg.setMask(mask)

        else:  # Axis

            gg = gg[jslice, islice]

#        # Selection using slices
#        if isinstance(lon, slice):
#            if outtype=='lon':
#                gg = gg[:, lon]
#            else:
#                gg = create_grid(gg.getLongitude()[:, lon], gg.getLatitude()[:, lon])
#            lon = None
#        if isinstance(lat, slice):
#            if outtype=='lat':
#                gg = gg[lat, :]
#            else:
#                gg = create_grid(gg.getLongitude()[lat, :], gg.getLatitude()[lat, :])
#            lat = None
#
#        # Selection using coordinates
#        if lon is not None or lat is not None:
#            mask,select = gg.intersect([lon, lat, None])
#            gg = gg.subSlice(**select)
#            gg.setMask(mask)

    else:  # Rectangular

        if outtype in ['grid', 'tuple']:
            axlon = gg.getLongitude()
            axlat = gg.getLatitude()
        elif outtype == 'lon' and isinstance(lon, tuple):
            ijk = gg.mapIntervalExt(lon)
            if ijk is None:
                raise VACUMMError('Empty selection on axis: '+str(lon))
            lon = slice(*ijk)
        elif outtype == 'lat' and isinstance(lat, tuple):
            ijk = gg.mapIntervaExt(lat)
            if ijk is None:
                raise VACUMMError('Empty selection on axis: '+str(lon))
            lat = slice(*ijk)

        # Selection using slices
        if isinstance(lon, slice):
            if outtype == 'lon':
                ii = lon.indices(len(gg))
                gg = gg.subaxis(*ii)
            elif outtype == 'grid':
                ii = lon.indices(len(axlon))
                axlon = axlon.subaxis(*ii)
                gg = create_grid(axlon, axlat)
            lon = None
        if isinstance(lat, slice):
            if outtype == 'lat':
                ii = lat.indices(len(gg))
                gg = gg.subaxis(*ii)
            elif outtype == 'grid':
                ii = lat.indices(len(axlat))
                axlat = axlat.subaxis(*ii)
                gg = create_grid(axlon, axlat)
            lat = None

        # Selection using coordinates
        if lon is not None or lat is not None:
            if lat is None:
                lat = (-90, 90)
            if lon is None:
                lon = (-180, 180)
            gg = gg.subGridRegion(lat, lon)

    # Out
    if outtype in ['lon', 'lat', 'grid']:
        return gg
    return (gg.getLongitude(), gg.getLatitude())


select_region = gridsel


def varsel(var, lon=None, lat=None, **kwargs):
    """Extract a subregion of a variable

    Lat and lon generic selection are used for
    spatial selections.

    If `var` has no grid or a rectangular grid, it simply
    equivalent to::

        var(lon=lon, lat=lat, **kwargs)

    .. note:: Properly works on curved grids thanks to :func:`coord2slice`.

    Parameters
    ----------

    var:
        MV2 array.
    lon/lat:
        Coordinates intervals or slices.
        - Extra keywords are passed to the variable.

    """
    grid = var.getGrid()

    # Usual way
    if not isgrid(grid, curv=True):
        return var(lon=lon, lat=lat, **kwargs)

    # Curved grid
    islice, jslice, mask = coord2slice(grid, lon=lon, lat=lat)
    xid = grid.getAxis(1).id
    yid = grid.getAxis(0).id
    var = var(**{xid: islice, yid: jslice})
    if mask.any():
        var[:] = MV2.masked_where(mask, var, copy=0)
    if kwargs:
        var = var(**kwargs)
    return var


def axis1d_from_bounds(axis1d, atts=None, numeric=False):
    """Create a numeric of formatted 1D axis from bounds

    axis1d:
        Input 1d axis from which we get bounds

    numeric:
        Return a simple numeric array
    atts:
        Attributes for outputs axis

    Example
    -------

    >>> xx = axis1d_from_bounds(xxb)
    """

    # Get bounds
    mybounds1d = axis1d.getBounds()
    if mybounds1d is None:
        mybounds1d = bounds1d(axis1d)
        axis1d.setBounds(mybounds1d)

    # Numerical values
    numax = N.zeros(len(axis1d)+1, 'd')
    numax[0] = mybounds1d[0, 0]
    numax[-1] = mybounds1d[-1, 1]
    numax[1:-1] = 0.5*(mybounds1d[1:, 0]+mybounds1d[:-1, 1])
    if numeric:
        return numax

    # Formatted axis
    ax = cdms.createAxis(numax, savespace=1)
    if atts is not None:
        for att, val in atts.items():
            setattr(ax, att, val)
    return ax


def get_xy(gg, proj=False, mesh=None, num=False, checklims=True, **kwargs):
    """Get axes from gg

    Parameters
    ----------

    gg:
        (xx,yy) or grid (1D or 2D), cdms variable (see :func:`get_grid`),
        or a dict of lon/lat/x/y, or a (2,npts) numpy array.
    proj: optional
        If True or basemap instance, convert to meters.
        If None, check if lon and lat axes to force conversion. [default: False]
    mesh: optional
        Get axes as 2D arrays.
    num: optional
        Get axes as numpy arrays.
    checklims: optional
        For 2D axes only, mask longitudes outside (-720,720),
        and clip latitude outside (-90,90).
    Example
    -------

        >>> lon, lat = get_xy((xaxis,yaxis))
        >>> lon, lat = get_xy(grid)
        >>> x, y = get_xy(dict(lon=(-10,0), lat=(42,43,'cc')), proj=True)
    """
    # Get axes
    if hasattr(gg, 'xy') and callable(gg.xy):
        xx, yy = gg.xy()
    elif isinstance(gg, (tuple, list)):
        if len(gg) == 4:
            xmin, ymin, xmax, ymax = gg
            xx = N.array([xmin, xmax])
            yy = N.array([ymin, ymax])
        else:
            xx, yy = gg
            if isinstance(xx, tuple):  # (xmin,xmax,'ccb')
                xx = N.array(xx[:2])
                yy = N.array(yy[:2])
    elif isinstance(gg, N.ndarray) and not cdms2.isVariable(gg) and \
            gg.ndim == 2 and gg.shape[0] == 2:
        xx, yy = gg
    elif isinstance(gg, dict):
        xy = []
        for keys in [['x', 'lon', 'longitude'], ['y', 'lat', 'latitude']]:
            for key in keys:
                if key in gg:
                    pp = gg[key]
                    if isinstance(pp, tuple):
                        pp = pp[:2]
                    xy.append(N.asarray(pp))
                    break
            else:
                raise KeyError('Dictionay has no %s key' % ' or '.join(keys))
        xx, yy = xy

    else:
        if cdms.isVariable(gg) and gg.getGrid() is not None:
            gg = gg.getGrid()
        try:
            xx = gg.getLongitude()
            if xx is None:
                xx = gg.getAxis(-1)
            yy = gg.getLatitude()
            if yy is None:
                yy = gg.getAxis(-2)
        except:
            gg = cdms.createVariable(gg)
            assert gg.ndim > 1,  'Input must be at least a 2D array in this case'
            xx = gg.getAxis(-1)
            yy = gg.getAxis(-2)

    # Pure numeric?
    xxn = _cdat2num_(xx)
    yyn = _cdat2num_(yy)
    if num:
        xx = xxn
        yy = yyn

    # Check limits of 2D axes
    if checklims:
        if xxn[:].ndim == 2:
            if isaxis(xx):
                xx = xx.clone()
                xx[:] = N.ma.masked_outside(xx[:], -720., 720., copy=False)
            else:
                xx = N.ma.masked_outside(xx[:], -720., 720., copy=False)
        if yyn[:].ndim == 2:
            if isaxis(yy):
                yy = yy.clone()
                yy[:] = N.ma.clip(yy[:], -720., 720.)
            else:
                yy = N.ma.clip(yy[:], -90., 90.)

    # Convert to meters?
    proj = kwargs.pop('m', proj)
    if proj is None:  # and islon(xx) and islat(yy):
        proj = True
    if proj is not False:
        xx, yy = deg2xy(xx, yy, proj=proj, mesh=mesh, **kwargs)

    # Mesh ?
    if mesh is True:  # and (xx[:].ndim == 1 or  yy[:].ndim == 1):
        return meshgrid(xx[:], yy[:])
    elif mesh is False:
        if xxn[:].ndim == 2:
            xx = xxn[0]
        if yyn[:].ndim == 2:
            yy = xxn[:, 0]

    return xx, yy


def _cdat2num_(xy):
    """Convert CDAT axis or variable to numeric (masked) form"""
    if hasattr(xy, 'asma'):
        xy = xy.asma()
    elif hasattr(xy, 'getValue'):
        xy = xy.getValue()
    if not isinstance(xy, N.ndarray):
        xy = N.asarray(xy)
    return xy


def deg2xy(lon, lat, proj=None, inverse=False, mesh=None, **kwargs):
    """Convert from degrees to map (m) coordinates using map projection, and reverse

    Parameters
    ----------

    lon:
        Longitudes in degrees
    lat:
        Latitude in degrees
    proj: optional
        Proj object for projection. If False, returns (lon,lat).
        If None, a new instance using :func:`~vacumm.misc.basemap.get_proj` is created,
        where proj is passed as a parameter.
    inverse: optional
        Inverse transform (from meters to degrees)

    Example
    -------

        >>> x, y = deg2xy(lon, lat)
        >>> x2d, y2d = deg2xy(lon, lat, mesh=True)
    """

    # Check shapes
    lon, lat = check_xy_shape(lon, lat, mesh=mesh)

    # Nothing to do
    proj = kwargs.pop('m', proj)
    if proj is False:
        return lon, lat

    # Check if we are sure to not have degrees
    xmin = lon.min()
    ymin = lat.min()
    xmax = lon.max()
    ymax = lat.max()
    if not inverse and xmin >= -0.1 and ymin >= -0.1 and (xmax > 720. or ymax > 90.):
        return lon, lat

    # Need a projection instance
    if not callable(proj):
        from .basemap import get_proj
        proj = get_proj((lon, lat), proj=proj)

    # Transform
    return proj(lon, lat, inverse=inverse)


def _dist2x2d_(xx, yy, mode):
    kw = dict(pairwise=True, mode=mode)
    return (get_distances(xx[:, :-1], yy[:, :-1], xx[:, 1:], yy[:, 1:], **kw),
            get_distances(xx[:-1], yy[:-1], xx[1:], yy[1:],  **kw))


# def _dist1x1d_(xxyy, xy, axis, mode):
#    kw = dict(pairwise=True, mode=mode)
#    dslices = get_axis_slices(xxyy, axis)
#    if axis==-1:
#        return

def resol(axy, mode='median',  axis=None, meters=False, cache=True, lat=None,
          checklims=True, **kwargs):
    """Get the resolution of an axis or a grid

    Parameters
    ----------
    axy:
        It can be either

        - a 1D axis or array
        - a grid of 1D or 2D axes or tuple of (lon,lat)

    mode: optional

        - ``"raw"``: Get the local resolution between "grid points".
        - ``"averaged"``: Return an averaged resolution (do not use if the grid highly anisotropic!!).
        - ``"local"``: The local resolution is averaged and extrapolated to grid points.
        - A string: An :mod:`numpy` attribute is used to compute the resolution.
          For instance ``"median"`` mode implies the use of :func:`numpy.median`.
        - A callable function: Directly used to compute the resolution.

    axis: optional
        Direction on which to compute resolution if a single
        2D axis is passed.
    meters: optional
        Get the resolution in meters instead of degrees.
    axis: optional
        Direction on which to compute resolution if a single
          2D axis is passed.
    lat: optional
        Latitude to use for projection of 1D zonal axis.

        .. warning:: If not provided but needed, it defaults to 45. and
            a warning is emitted.

    .. warning::

        If you work on a pair of 2D axes, resolution is computed along
        grid axes, and NOT along X and Y.
        In this case, X and Y must have consistent units,
        and resolution along X and Y are defined by:

        .. math::

            dx_{i,j} = \sqrt{(x_{i+1,j}-x_{i,j})^2+ (y_{i+1,j}-y_{i,j})^2}

            dy_{i,j} = \sqrt{(x_{i,j+1}-x_{i,j})^2+ (y_{i,j+1}-y_{i,j})^2}

    Examples
    --------

        >>> dx = resol(lon)
        >>> dx,dy = resol(grid)
        >>> dx2d = resol(x2d, mode='loc')
        >>> dx2d, dy2d = resol((x1d, y2d))
    """

    # Get what to inspect
    if cdms2.isVariable(axy):
        if not isaxis(axy):
            axy = axy.getGrid()
            assert axy is not None, 'You must pass a variable with a grid'
        else:
            axy = axy.asma()

    # Meters?
    proj = kwargs.get('proj', False)
    if proj:
        meters = True
    if meters is None:
        meters = False

    # Check chache (for cdms grid and axes)
    if kwargs.get('averaged', False):  # compat
        mode = 'averaged'
    if cache and not mode.startswith('loc'):
        pres = 'res'+'m' if meters else ''
        if hasattr(axy, '_'+pres):
            return getattr(axy, '_'+pres)
        if hasattr(axy, '_x'+pres) and hasattr(axy, '_y'+pres):
            return getattr(axy, '_x'+pres), getattr(axy, '_y'+pres)
        if int(cache) > 1:
            return

    # Numerical values
    single = None
    if isinstance(axy, tuple) and len(axy) == 1:
        axy = axy[:]
        single = False
    if not isgrid(axy) and not isinstance(axy, tuple):  # single axis
        if single is None:
            single = True

        atype = 'x' if islon(axy) else 'y'
        if isaxis(axy):
            axy = axy.getValue()
        xy = axy,

    else:  # grid or pair of axes > convert to 2D arrays
        single = False

        xy = get_xy(axy, num=True, proj=False, mesh=True, checklims=checklims)

    # Local distances
    res = ()
    distmode = 'haversine' if meters else 'simple'
    kwdist = dict(mode=distmode)
    warnlatmsg = ("Your axis resolution along X is computed with a default "
                  "latitude of 45. You should also provide a valid Y axis or "
                  "at least set this latitude properly with the lat parameter.")
    if len(xy) == 2 and (xy[0][:].ndim == 2 or xy[1][:].ndim == 2):  # 2x2D

        xy = meshgrid(*xy)
        res = _dist2x2d_(*xy,  **kwdist)

    else:  # Single 1D or 2D

        kwdist.update(pairwise=True)

        if N.ndim(xy[0][:]) == 2:  # 2D

            if axis is None:
                if atype == 'y':
                    axis = -2
                else:
                    axis = -1
            if axis < 0:
                axis += xy[0].ndim

            ds = get_axis_slices(2, axis)
            if atype == 'x':
                if lat is None:
                    lat = 45.
                    vacumm_warn(warnlatmsg)
                res = get_distances(axy[ds['firsts']],
                                    lat, axy[ds['lasts']], lat, **kwdist),
            else:
                res = get_distances(
                    0., axy[ds['firsts']], 0., axy[ds['lasts']], **kwdist),

        else:

            if atype == 'x':
                if lat is None:
                    lat = 45.
                    vacumm_warn(warnlatmsg)
                res = get_distances(axy[:-1], lat, axy[1:], lat, **kwdist),
            else:
                res = get_distances(0., axy[:-1], 0., axy[1:], **kwdist),

    # Averages
    if mode.startswith('loc'):
        if xy[0].ndim == 2:  # 2D
            eres = ()
            if len(xy) == 1:
                kk = [axis]
            else:
                kk = [1, 0]
            for ik, k in enumerate(kk):

                # Init output
                eshape = list(res[ik].shape)
                eshape[k] += 1
                eres += N.ma.resize(res[ik], eshape),

                # Slices specs
                sl = get_axis_slices(res[ik].shape, k)
                sle = get_axis_slices(eshape, k)

                # Core
                eres[ik][sle['mid']] = .5 * \
                    (res[ik][sl['firsts']]+res[ik][sl['lasts']])

                # Limits
                eres[ik][sle['first']] = 1.5 * \
                    res[ik][sl['first']]-0.5*res[ik][sl['firstp1']]
                eres[ik][sle['last']] = 1.5 * \
                    res[ik][sl['last']]-0.5*res[ik][sl['lastm1']]

            res = eres

        else:  # 1D

            eres = tuple([N.ma.resize(r, (r.shape[0]+1, )) for r in res])
            for ik in range(len(xy)):
                eres[ik][1:-1] = .5*(res[ik][:-1]+res[ik][1:])
                eres[ik][0] = 1.5*res[ik][0]-0.5*res[ik][1]
                eres[ik][-1] = 1.5*res[ik][-1]-0.5*res[ik][-2]
            res = eres

    elif mode != 'raw':

        if not isinstance(mode, tuple):
            mode = mode,
        if len(mode) < len(res):
            mode *= 2
        eres = []
        for ik, dcell in enumerate(res):
            m = mode[ik]
            if m.startswith('ave'):
                m = 'mean'
            if isinstance(m, str):
                func = getattr(N, m)
            else:
                func = m
            if N.ma.isMA(dcell):
                dcell = dcell.compressed()
            eres.append(func(dcell))
        res = tuple(eres)

    # Caching
    if cache and (isinstance(mode, tuple) or not mode.startswith('loc')):
        if isgrid(axy):
            setattr(axy, '_x'+pres, res[0])
            setattr(axy, '_y'+pres, res[1])
            setattr(axy, '_mres', mode)
        elif isaxis(axy):
            setattr(axy, '_'+pres, res)
            setattr(axy, '_mres', mode)

    if single:
        return res[0]
    return res


def check_xy_shape(xx, yy, mesh=None):
    """Check that xx and yy have the same shape

    xx:
        X positions (in meters or degrees).
    yy:
        Y positions (in meters or degrees).
    mesh:
        Return a 2D axes if True, or if None and xx and yy have not the same shape [default: None]

    Example
    -------

    >>> xx2d, yy2d = check_xy_shape(xx, yy, mesh=True)
    """
    # Numpy arrays
    xx = N.asarray(xx[:])
    yy = N.asarray(yy[:])

    # Check if we want 2D axes
    if xx.shape != yy.shape:
        mesh = True
    if mesh is None:
        mesh = xx.shape == yy.shape

    # Go
    if not mesh:
        return xx, yy
    if xx.ndim == 2 and yy.ndim == 2 and xx.shape == yy.shape:
        return xx, yy
    if xx.ndim != 1 or yy.ndim != 1:
        raise VACUMMError('xx and yy must be 1D for meshgrid transform.')
    return N.meshgrid(xx, yy)


def t2uvgrids(gg, getu=True, getv=True, mask=None):
    """Convert a (C) grid at T-points to a grid at U- and/or V- points

    Parameters
    ----------

    gg:
        A (x,y) or a cdms grid or a cdms variable with a grid (see :func:`get_grid`)
    getu:
        Get the grid at U-points [default: True]
    getv:
        Get the grid at V-points [default: True]
    mask:
        If False, do not try to guess masks ; if None, try to get them from
        original grid by conversion to U and V points with t2uvmask() [default: None]

    Return
    ------
    ``ugrid,vgrid`` OR ``ugrid`` OR ``vgrid`` depending on getu and getv

    Example
    -------

        >>> gridu, gridv = t2uvgrids(gridt)
    """
    if not getu and not getv:
        return
    # Get rights axes
    if cdms.isVariable(gg):
        gg = gg.getGrid()
    if isgrid(gg) and mask is not False:
        mask = gg.getMask()
    if mask is False:
        mask is None
    xx, yy = get_xy(gg)
    if not islon(xx):
        xx = create_lon(xx)
    if not islat(yy):
        yy = create_lon(yy)
    # U-points
    if mask is not None:
        from .masking import t2uvmasks
    if getu:
        xu = xx.clone()
        xu[:-1] = .5*(xx[:-1]+xx[1:])
        xu[-1] = xx[-1]+.5*(xx[-1]-xx[-2])
        xu.long_name += ' at U-points'
        yu = yy.clone()
        gu = cdms.createRectGrid(yu, xu)
        if mask is not None:
            gu.setMask(t2uvmasks(mask, getv=False))
        if not getv:
            return gu
    # V-points
    yv = yy.clone()
    yv[:-1] = .5*(yy[:-1]+yy[1:])
    yv[-1] = yy[-1]+.5*(yy[-1]-yy[-2])
    yv.long_name += ' at V-points'
    xv = xx.clone()
    gv = cdms.createRectGrid(yv, xv)
    if mask is not None:
        gv.setMask(t2uvmasks(mask, getu=False))
    if not getu:
        return gv
    return gu, gv


def rotate_grid(ggi, angle, pivot='center', **kwargs):
    """Rotate a grid

    Parameters
    ----------
    ggi:
        Cdms grid ou (lon,lat) or variable (see :func:`get_grid`).
    angle:
        Angle in degrees.
    pivot: optional
        it can be either

            - A tuple of (lon, lat)
            - A string that specify the vertical and horizontal position.
              Use the following keys : ``top``, ``bottom``,
              ``left``, ``right``, ``center``.

        - Other keywords are passed to :func:`~vacumm.misc.grid.curv_grid`

    Returns
    -------
    A curvilinear cdms grid.

    Example
    -------
    >>> mygrid = rotate_grid((lon,lat), 60., 'top right')
    >>> mygrid2 = rotate_grid(mygrid, 60., (-5.,45))
    """
    # Get 2D axes
    xxi, yyi = meshgrid(*get_xy(ggi))

    # Find pivot
    if isinstance(pivot, tuple):
        xpivot, ypivot = pivot
    else:
        nx, ny = xxi.shape
        pivot = pivot.lower()
        if pivot == 'center':
            pivot = 'center center'
        sxpivot, sypivot = pivot.split()
        ipivot = [0, int(nx/2), -1][['left', 'center', 'right'].index(sxpivot)]
        jpivot = [0, int(ny/2), -1][['bottom', 'center', 'top'].index(sypivot)]
        xpivot = xxi[jpivot, ipivot]
        ypivot = yyi[jpivot, ipivot]

    # Rotate
    angle = angle*N.pi/180.
    xxo = xpivot + (xxi-xpivot)*N.cos(angle) - (yyi-ypivot)*N.sin(angle)
    yyo = ypivot + (xxi-xpivot)*N.sin(angle) + (yyi-ypivot)*N.cos(angle)

    return curv_grid(xxo, yyo, **kwargs)


def monotonic(xx, xref=None):
    """Make monotonic an array of longitudes

    Example
    -------

    >>> xxm = monotonic(xx, xref=120.)
    >>> print (N.diff(xxm)<0).any()
    False
    """

    # Longitude reference
    modulo = xx.modulo if hasattr(xx, 'modulo') else 360.
    if xref is not None:
        xx[:] = (xx[:]-xref) % modulo+xref

    # Monotonic modulo
    steps = N.diff(xx, axis=-1)
    bigsteps = N.abs(steps) > modulo/2.
    if not N.any(bigsteps):
        return xx
    fixsteps = N.where(bigsteps, -modulo*N.sign(steps), 0.)
    fixsteps.cumsum(axis=-1, out=fixsteps)
    xx[..., 1:] += fixsteps
    return xx


def xextend(vari, n, mode=None):
    """Extend a variable along x axis

    .. note:: x is the last axis and must be 1D

    .. todo::

        Make possible to use 2D axes with :func:`xextend`

    Parameters
    ----------

    vari:
        :mod:`cdms2` variable
    n:
        size of the extention

            - if an integer: ``nleft = nright = n``
            - else: nleft, nright = n

    mode: optional
        mode of extension

            - ``None``: choose ``"cyclic"`` is x axis
              has a "modulo" attribute
            - ``"cyclic"``: assume x axis is cyclic
            - ``"constant"``: extrapolate with first or last values
            - else, extend with masked values
              and linearly extrapolate positions

    Example
    -------

        >>> varo = xextend(vari, (5, 7), mode='cyclic')
    """
    # Get extension sizes and mode
    xi = vari.getAxis(-1)
    if mode is None and hasattr(xi, 'modulo'):
        mode = 'cyclic'
    if not isinstance(n, int):
        nleft, nright = n
    else:
        nleft = nright = n

    # Intialize
    nxi = vari.shape[-1]
    nxo = nleft+nxi+nright
    varo = MV2.zeros(vari.shape[:-1]+(nxo, ), dtype=vari.dtype)
    varo[:] = MV2.masked
    varo[..., nleft:nleft+nxi] = vari
    cp_atts(vari, varo)
    axes = vari.getAxisList()

    # Extend
    if mode == 'cyclic':
        if not hasattr(xi, 'modulo'):
            modulo = 360.
        else:
            modulo = xi.modulo
        varo[..., :nleft] = vari[..., nxi-nleft:]
        varo[..., nxo-nright:] = vari[..., :nright]
        xleft = xi[nxi-nleft:]-modulo
        xright = xi[:nright]+modulo
    else:
        dx0 = xi[1]-xi[0]
        dx1 = xi[-1]-xi[-2]
        xleft = xi[0]-N.arange(1, nleft+1)*dx0
        xright = xi[-1]+N.arange(1, nright+1)*dx1
        if mode == 'constant':
            for ileft in range(nleft):
                varo[..., ileft] = vari[..., 0]
            for iright in range(nright):
                varo[..., nxo-iright-1] = vari[..., -1]

    # Insert axes
    xo = cdms2.createAxis(N.concatenate(
        (xleft, xi[:], xright)).astype(xi[:].dtype))
    cp_atts(xi, xo)
    axes[-1] = xo
    varo.setAxisList(axes)
    return varo


def xshift(vari, n, mode=None, monot=True):
    """Shift a cyclic array along x axis

    .. note:: x is the last axis and must be 1D

    .. todo::

        Make possible to use 2D axes with :func:`xshift`

    Parameters
    ----------

    vari:
        :mod:`cdms2` variable
    n:
        size of the shift

            Positive means east bound is moved toward east.
            If ``mode is None``:

            - if an integer, use it as grid points
            - else (a float), use at as degrees.

    Example
    -------

        >>> varo = xshift(vari, 50)
    """
    # Get shift value
    xi = vari.getAxis(-1).getValue()
#    if mode=='left':
#
#    if mode is None:
#        if not isinstance(n, int):
#            n = float(n)

    # Intialize
    varo = vari.clone()
    xo = xi.copy()
    if n == 0:
        return varo
    if n > 0:
        varo[..., :-n] = vari[..., n:]
        varo[..., -n:] = vari[..., :n]
        xo[:-n] = xi[n:]
        xo[:-n] = xi[n:]
    else:
        varo[..., n:] = vari[..., :-n]
        varo[..., :n] = vari[..., -n:]
        xo[n:] = xi[:-n]
        xo[n:] = xi[:-n]
    if xi[-1] > xi[0]:
        xi = monotonic(xi)
    varo.getAxis(-1).assignValue(xo)
    return varo


def curv2rect(gg, mode="warn", tol=1.e-2, f=None):
    """Convert a curvilinear grid to a rectangular grid

    .. warning::

        This not an interpolation.
        Longitudes and latitudes are converted from 2D arrays
        to 1D arrays using axis averages.

    Parameters
    ----------

  gg:
      Grid, tuple of axes or MV2 variable.
  mode: optional
      what ot do when it does not seems to be rectangular

              - ``"warn"``: simply show a warning,
              - ``"raise"``: raise a :class:`VACUMMError`,
              - ``"none"`` or ``False``: don't do it.
              - else just convert it at your own risks.

    tol: optional
        Tolerance (passed to :func:`isrect`).
    f: optional
        File object or name (passed to :func:`isrect`).

    Example
    -------

        >>> rgrid = curv2rect(cgrid)
        >>> rgrid = curv2rect((lon2d,lat2d))

    Parameters
    ----------

    gg:
        A variable or a curvilinear grid (see :func:`get_grid`)
    """

    # MV2 variable
    if cdms2.isVariable(gg):
        grid = get_grid(gg)
        if not isgrid(grid, curv=True):
            return gg
        rgrid = curv2rect(grid, mode=mode, tol=tol, f=f)
        if rgrid is not None and grid is not rgrid:
            set_grid(gg, rgrid)
        return gg

    # Is it rectangular?
    gg = get_grid(gg)
    if gg is None:
        return
    if not isrect(gg, tol=tol, f=f):
        if mode == "warn":
            vacumm_warn("Grid seems not trully rectangular. Not converted.")
            return gg
        elif mode == "raise":
            raise VACUMMError('Cannot convert to rectangular grid')
        if mode is False or mode == "none":
            return gg

    # Convert it
    xx, yy = get_xy(gg, mesh=True, num=True)
    x = xx.mean(axis=0)
    y = yy.mean(axis=1)
    grido = create_grid(x, y)
    cp_atts(gg.getLongitude(), grido.getLongitude())
    cp_atts(gg.getLatitude(), grido.getLatitude())
    return grido


def isrect(gg, tol=1.e-2, mode="real", f=None, nocache=False):
    """Check wether a grid is trully rectangular

    Parameters
    ----------

    gg:
        cdms2 grid, tuple of cdms axes or variabe with a grid.
    tol: optional
        Tolerance for coordinates deformation.
    mode: optional
        Simple or real check.

            - ``"simple"``: Simply check if axes are 1D.
            - ``"real"``: Also check that axes along one direction do
              vary along the other direction with tolerance tol.

    f: optional
        netcdf file object.

    """

    # Check that it's a grid
    gg = get_grid(gg, intercept=False)

    # Already a rectangular
    is1d = isgrid(gg, curv=False)
    if not mode.startswith('r'):
        return is1d
    if is1d:
        return True

    # File grid
    if not nocache and f is not None:
        from .io import ncget_fgrid
        fgrid = ncget_fgrid(f, gg)
    else:
        fgrid = None

    # From cache
    if not nocache:
        for grid in fgrid, gg:
            if hasattr(grid, '_vacumm_isrect'):
                return grid._vacumm_isrect

    # Check if this grid is really rectangular
    # - resolution from cache
    for grid in fgrid, gg:
        if grid is None:
            continue
        res = resol(grid, cache=2)
        if res is not None:
            break
    else:
        res = resol(gg)
        if fgrid is not None:
            fgrid._res = res
    dx, dy = res
    # - relative resolution variations
    xx, yy = get_xy(gg, num=True)
    xx = N.ma.masked_outside(xx, -720., 720., copy=False)
    yy = N.ma.clip(yy, -90, 90.)
    Dx = N.ma.median(xx.ptp(axis=0))
    Dy = N.ma.median(yy.ptp(axis=-1))
    res = Dx < tol*dx and Dy < tol*dy

    # Cache result
    for grid in gg, fgrid:
        if grid is not None:
            grid._vacumm_isrect = res

    return res


def transect_specs(gg, lon0, lat0, lon1, lat1, subsamp=3, getxy=False,
                   m=None, getmap=False, **kwargs):
    """Get specs for a transect given a grid and starting and ending points

    Parameters
    ----------

    lon/lat0/1:
        Coordinates of first and last point in degrees.
    subsamp: optional
        Subsampling with respect to grid cell.
    getxy: optional
        Also return projected coordinates.
    getproj: optional
        Also return projection map.

    Return
    ------
    lons, lats[, xx, yy][, proj]

    See also
    --------
    :func:`resol`, :meth:`mpl_toolkits.basemap.Basemap.gcpoints`.
    """

    # Bounds and resolution
    from .basemap import get_map  # FIXME: haversine and geo for transect_specs
    dx, dy = resol(gg, proj='merc')
    if m is None:
        m = get_map(gg, resol=None, proj='merc')
    x0, y0 = m(lon0, lat0)
    x1, y1 = m(lon1, lat1)
    Dx = N.abs(x1-x0)
    Dy = N.abs(y1-y0)

    # Number of points
    if dy*Dx > dx*Dy:
        npts = Dx/dx
    else:
        npts = Dy/dy
    npts *= subsamp
    npts = int(N.ceil(npts))

    # Coordinates
    xx, yy = m.gcpoints(lon0, lat0, lon1, lat1, npts)  # FIXME: module geo
    lons, lats = m(xx, yy, inverse=True)
    ret = N.array(lons), N.array(lats)
    if getxy:
        ret += N.array(xx), N.array(yy)
    if 'getproj' in kwargs:
        getmap = kwargs['getproj']
    if getmap:
        ret += m
    return ret


def depth2dz(depth, axis=None, mode=None, masked=True):
    """Convert from depth to layer thickness

    Parameters
    ----------

    depth:
        1D or ND array of depth (variable or axis).
    axis: optional
        Axis of vertical dimension.
        It is guessed of not provided.
    mode: optional
        Mode for defining the reference layer.
        Possible values: ``'first'``, ``'last'``, ``center`` or ``None``.
        It is either the first or the last layer. The layer is
        computed normally like other layer by difference,
        whereas its opposite (hidden) layer is either masked or has its
        thickness set the the thickness of its adjacent layer.
        If ``mode`` is set to ``None``, it set to ``'last'``
        if depths are positive up.
    masked: optional
        If True, hidden layer is masked.

    """
    # Init
    if cdms2.isVariable(depth):
        dz = depth.clone()
        dz.long_name = "Layer thickness"
        dz.id = 'dz'
        dz.units = 'm'
    elif isaxis(depth):
        dz = N.empty(depth.shape)
    else:
        dz = depth.copy()
    if masked and not N.ma.isMA(dz):
        dz = N.ma.asarray(dz)
        dz[:] = N.ma.masked

    # Guess axis to work on
    if axis is None:
        axis = _getzdim_(depth)
    slices = get_axis_slices(depth, axis)

    # Priority mode
    if isinstance(mode, six.string_types):
        mode = mode.lower()
    if mode not in ['last', 'first', None, 'center']:
        mode = None
    if mode is None:
        mode = ['first', 'last'][int(isdepthup(depth, axis=axis))]

    # Compute
    if depth.shape[axis] == 1:
        dz[:] = 0
        return dz
    if mode == 'last':

        dz[slices['lasts']] = depth[slices['lasts']]-depth[slices['firsts']]
        if not masked:
            dz[slices['first']] = dz[slices['firstp1']]

    elif mode == 'first':

        dz[slices['firsts']] = depth[slices['firsts']]-depth[slices['lasts']]
        if not masked:
            dz[slices['last']] = dz[slices['lastm1']]

    else:  # center

        dz[slices['mid']] = 0.5 * \
            (depth[slices['lastsp1']]-depth[slices['firstsm1']])
        if not masked:
            dz[slices['first']] = dz[slices['firstp1']]
            dz[slices['last']] = dz[slices['lastm1']]

    dz[:] = abs(dz)
    return dz


def get_axis_slices(ndim, axis, **kwargs):
    """Get standard slices for an axis of a ndim array

    Parameters
    ----------

    ndim:
        The number of dimensions. It can also be
        a tuple (like an array shape) or an array.
    axis:
        Index of the axis.

    Return
    ------
    A dictionary of tuples of slices. All tuples have a
        length of ndim, and can be used has a slice for the array
        (see example).

        - "all": Select everything.
        - "first"/"last": First and last.
        - "firstp1": Second element.
        - "firstp2": Third element.
        - "lastm1": Element before the last one.
        - "lastm2": Second element before the last one.
        - "firsts": All but the last.
        - "lasts": All but the first.
        - "firstsm1": All but the last two.
        - "lastsp1": All but the first two.
        - "mid": All but the first and last.

    Example
    -------

        >>> var = N.arange(3*4).reshape(3, 4)
        >>> ss = get_axis_slices(var, axis=0)
        >>> var_north = var[ss['last']]
        >>> var_south = var[ss['first']]
   """
    if hasattr(ndim, 'shape'):
        ndim = ndim.shape
    if hasattr(ndim, '__len__'):
        ndim = len(ndim)
    sel = [slice(None)]*ndim
    selmid = list(sel)
    selmid[axis] = slice(1, -1)
    selmid = tuple(selmid)
    sellasts = list(sel)
    sellasts[axis] = slice(1, None)
    sellasts = tuple(sellasts)
    selfirsts = list(sel)
    selfirsts[axis] = slice(0, -1)
    selfirsts = tuple(selfirsts)
    sellastsp1 = list(sel)
    sellastsp1[axis] = slice(2, None)
    sellastsp1 = tuple(sellastsp1)
    selfirstsm1 = list(sel)
    selfirstsm1[axis] = slice(0, -2)
    selfirstsm1 = tuple(selfirstsm1)
    sellast = list(sel)
    sellast[axis] = -1
    sellast = tuple(sellast)
    selfirst = list(sel)
    selfirst[axis] = 0
    selfirst = tuple(selfirst)
    sellastm1 = list(sel)
    sellastm1[axis] = -2
    sellastm1 = tuple(sellastm1)
    sellastm2 = list(sel)
    sellastm2[axis] = -3
    sellastm2 = tuple(sellastm2)
    selfirstp1 = list(sel)
    selfirstp1[axis] = 1
    selfirstp1 = tuple(selfirstp1)
    selfirstp2 = list(sel)
    selfirstp2[axis] = 2
    selfirstp2 = tuple(selfirstp2)
    if kwargs:
        for key, val in kwargs.items():
            if isinstance(val, (list, tuple)):
                val = slice(*val)
            ksel = list(sel)
            ksel[axis] = val
            kwargs[key] = ksel
    return dict(all=sel, mid=selmid, lasts=sellasts, firsts=selfirsts,
                lastsp1=sellastsp1, firstsm1=selfirstsm1,
                last=sellast, first=selfirst, lastm1=sellastm1, firstp1=selfirstp1,
                lastm2=sellastm2, firstp2=selfirstp2, **kwargs)


def merge_axis_slices(slices1, slices2):
    """Merge standard tuples of slices stored in dictionaries created with :func:`get_axis_slices`:

    Parameters
    ----------
    slice1/2:
        Dictionaries of tuples of slices.

    Example
    -------

        >>> slicesx = get_axis_slices(3, -1)  # for X
        >>> slicesy = get_axis_slices(3, -2)  # for Y
        >>> slices12 = merge_axis_slices(slicesx, slicesy) # for X and Y
        """
    slicesm = {}
    for key in slices1:
        slicesm[key] = merge_axis_slice(slices1[key], slices2[key])
    return slicesm


def merge_axis_slice(sel1, sel2):
    """Merge a single tuple of slices

    Theses slices may have been created with :func:`get_axis_slices`

    Parameters
    ----------

    sel1/2:
        Tuples of slices.

    Example
    -------

        >>> sel1 = (slice(None), -2)
        >>> sel2 = (slice(1,4), slice(None)
        >>> merge_axis_slice(sel1, sel2)
        (slice(1, 4, None), -2)

    """
    ndim1 = len(sel1)
    ndim2 = len(sel2)
    sel1 = list(sel1)
    sel2 = list(sel2)
    selm = []
    if ndim1 < ndim2:
        sel1 = [slice(None)]*(ndim2-ndim1)+sel1
    elif ndim1 > ndim2:
        sel2 = [slice(None)]*(ndim1-ndim2)+sel2
    for axis in range(len(sel1)):
        if sel1[axis] is None or sel1[axis] == slice(None) or \
                sel1[axis] is Ellipsis or sel1[axis] == ':':
            selm.append(sel2[axis])
        else:
            selm.append(sel1[axis])
    return tuple(selm)


def get_zdim(var, axis=None, default=None, strict=True):
    """Get the index of the Z dimension"""
    if isinstance(axis, int):
        return axis
    ndim = len(var.shape)
    if ndim == 1 and (not strict or isdep(var)):  # a single axis
        return 0
    if isaxis(axis) and cdms2.isVariable(var):  # find axis index in var axes
        i = var.getAxisList().index(axis)
        if i != -1:
            axis = i
    # find the same axis length
    if axis is not None and not strict and hasattr(axis, '__len__') and len(axis) != 0:
        good = [i_ for i_ in var.shape if i_ == len(axis)]
        if len(good) == 1:
            return good[0]
        axis = None
    if axis is None:
        if cdms2.isVariable(var):  # search for z in order
            axis = var.getOrder().find('z')
            if axis == -1:
                axis = None
        if not strict and axis is None and ndim > 2:
            axis = ndim-3
            if cdms2.isVariable(var) and var.getOrder()[axis] != '-':
                axis = None
    if axis is None:
        axis = default
    return axis


_getzdim_ = get_zdim


def isdepthup(depth, axis=None, ro=True):
    """Guess if depth is positive up

    It first check "positive" attribute, then guess from values:
    more positive values means positive down.

    ..warning:: Bad values must be masked

    Parameters
    ----------

    depth:
        Depth arrays or axis.
    axis: optional
        Z axis index.
    ro: optional
        Read-only? If False and if depth
        is an axis, its 'positive' attribute is marked
        'up' or 'down' upon the results.

    Return
    ------
    ``None`` if not Z axis found, else ``True`` or ``False``
    """
    # From attribute
    positive = getattr(depth, 'positive', None)
    if isinstance(positive, six.string_types):
        if positive.lower() in ['down', 'up']:
            return positive == 'up'
        positive = None

    # Slices
    axis = get_zdim(depth, axis=axis, default=None)
    if axis is None:
        return
    slices = get_axis_slices(depth, axis)

    # From values
    firstval = depth[slices['first']].max()
    lastval = depth[slices['last']].max()
    isup = N.abs(firstval) > N.abs(lastval)

    # Set attribute
    if not ro and isaxis(depth):
        depth.positive = 'up' if isup else 'down'
    return isup


def _makelevup_(negative, vv, levels=None, axis=None, default=None, ro=False, strict=True):
    # A list of variables
    single = not isinstance(vv, (list, tuple))
    if single:
        vv = [vv]

    # Get depth
    getlevels = isaxis(levels) or isinstance(levels, N.ndarray)
    if levels is None:
        if isaxis(vv[0]):
            levels = levels  # ?
        else:
            levels = vv[0].getLevel()
    if axis is None and isaxis(levels):
        axis = levels
    if isinstance(levels, bool):
        isup = levels
    else:

        axis = get_zdim(vv[0], axis, default=default, strict=strict)

        # No levels found
        if axis is None:
            vv = vv[0] if single else vv
            if getlevels:
                return vv, levels
            return vv

        # Positive up?
        isup = isdepthup(levels, axis=axis, ro=ro)

        # Revert levels
        if not isup and getlevels:
            if negative:
                levels[:] *= -1
            if isaxis(levels):
                levels.assignValue(levels[::-1])
                levels.positive = 'up'
            else:
                levels[:] = levels.take(
                    list(range(levels.shape[axis]-1, -1, -1)), axis=axis)

    # Revert variables or levels
    if not isup:
        for ivar, var in enumerate(vv):
            if isaxis(var):
                var.assignValue(var[::-1] * (-1 if negative else 1))
                var.positive = 'up'
            else:
                axis = get_zdim(var, axis=axis)
                if axis is None:
                    continue
                var[:] = var.take(
                    list(range(var.shape[axis]-1, -1, -1)), axis=axis)
                depaxis = var.getAxis(axis)
                depaxis.assignValue(depaxis[::-1] * (-1 if negative else 1))

    vv = vv[0] if single else vv
    if getlevels:
        return vv, levels
    return vv


def makedepthup(vv, depth=None, axis=None, default=None, ro=False, strict=True):
    """Make depth and variables positive up

    Parameters
    ----------

    vv:
        A single variable or a list of them.
    depth: optional
        Explicit depths to not guess
        it with :meth:`getLevel`. If True or False, simply revert along Z
        dimension.
    axis: optional
        Z dimension (else guessed with :func:`get_zdim`).

    """
    return _makelevup_(True, vv, levels=depth, axis=axis, default=default,
                       ro=ro, strict=strict)


def makealtitudeup(vv, altitude=None, axis=None, default=None, ro=False, strict=True):
    """Make altitude and variables positive up

    Parameters
    ----------

    vv:
        A single variable or a list of them.
    altitude: optional
        Explicit altitudes to not guess
        it with :meth:`getLevel`. If True or False, simply revert along Z
        dimension.
    axis: optional
        Z dimension (else guessed with :func:`get_zdim`).

    """
    return _makelevup_(False, vv, levels=altitude, axis=axis, default=default,
                       ro=ro, strict=strict)


def dz2depth(dz, ref=None, refloc=None, copyaxes=True, mode='edge',
             zerolid=False):
    """Conversion from layer thickness to depths

    Parameters
    ----------

    dz:
        Layer thickness from bottom to top.
    ref:
        Reference depth for integration, which depends
        on ``refloc``:

            - ``"top"`` or ``"eta"`` or ``"ssh"``: Sea surface height.
            - ``"bottom"`` or ``"depth"`` or ``"bathy"``: Bottom depth.
            - Else, estimated by checking the standard_name of the
              variable, if some values are negatives ('top')
              or if maximal value is greater then 15. ('bottom').

    refloc: optional
        ``"top"`` | ``"eta"``,
        ``"bottom"`` | ``"depth"``,  or ``None``.
    mode: optional

            - ``"edge"`` or ``"edge+"``: Compute depths at layer edges (interfaces).
              Adding a + includes the bottom layer and add a vertical level.
            - ``"middle"``: Compute depths at the middle of layers

        - **zerolid**, optional: Force the surface to be at a zero depth

    """

    # Init depths
    if mode is None:
        mode = 'edge'
    mode = str(mode)
    ext = '+' in mode
    mode = mode[:3]
    assert mode in ('edg', 'mid'), ("Invalid mode: Please choose one of: "
                                    "edge or middle")
    ext = ext and mode == 'edg'
    withtime = dz.getTime() is not None
    nt = dz.shape[0] if withtime else 1
    nz = dz.shape[int(withtime)]
    if ext:
        nz += 1
    shape = (nt, nz) + dz.shape[int(withtime)+1:]
    depths = MV2.zeros(shape, dz.dtype)
    depths.long_name = 'Depths'
    depths.units = 'm'
    depths.id = 'depths'

    # Format dz
    dzm = dz.asma()
    if not withtime:
        dzm = dzm.reshape((1, )+dzm.shape)
    # TODO: dz2depth: dzshift
    # TODO: dz2depth: positive up/down

    # Guess reference
    if ref is None:
        ref = 0.
    if isinstance(refloc, six.string_types):
        refloc = refloc.lower()
    if refloc in ["eta", "ssh"]:
        refloc = 'top'
    elif refloc in ["depth", "bathy"]:
        refloc = 'bottom'
    if refloc not in ['top', 'bottom']:
        refloc = None
    refsn = getattr(ref, 'standard_name', None)
    ref = ref.asma() if cdms2.isVariable(ref) else N.ma.asarray(ref)
    if refloc is None:

        # From standard_name (using cf)
        if refsn is not None:
            try:
                from vacumm.data.cf import VAR_SPECS
                if refsn in VAR_SPECS['ssh']['standard_names']:
                    refloc = 'top'
                elif True in [(refsn in VAR_SPECS[vn]['standard_names'])
                              for vn in ['bathy', 'bathy_u', 'bathy_v']]:
                    refloc = 'bottom'
            except:
                pass

        # From values
        if refloc is None:
            if (ref <= 0).any():
                refloc = 'top'
            else:
                refloc = 'bottom' if ref.max() > 15 else 'top'

    # Integrate
    depths[:, int(ext):] = dzm.cumsum(axis=1)

    # Add reference
    if refloc == 'top':  # zero at top
        for it in range(nt):
            depths[it] -= depths[it, -1]
        if zerolid:
            ref = 0
    depths[:] += ref

    # Substract surface for zerolid
    if zerolid and refloc != 'top':
        for dep in depths:
            dep[:] -= dep[-1]

    # Middle of layers
    if mode == 'mid':
        if refloc == 'top':
            depths[:] += dzm * 0.5
        else:
            depths[:] -= dzm * 0.5
    del dzm

    # Format axes
    if not withtime:
        depths = depths[0]
    if copyaxes:
        axes = dz.getAxisList()
        if ext:
            iaxis = int(withtime)
            from .regridding import extend1d
            zaxis = extend1d(dz.getAxis(iaxis), (1, 0), mode='linear')
            axes[iaxis] = zaxis
        depths.setAxisList(axes)
        grid = dz.getGrid()
        if grid is not None:
            depths.setGrid(grid)

    return depths


def are_same_grids(grid0, grid1):
    """Check if two grids are similars"""
    if grid0.shape != grid1.shape:
        return False
    lon0 = grid0.getLongitude()[:]
    lat0 = grid0.getLatitude()[:]
    lon1 = grid1.getLongitude()[:]
    lat1 = grid1.getLatitude()[:]
    if (lon0.shape != lon1.shape) or (lat0.shape != lat0.shape):
        return False
    for a0, a1 in [(lon0, lon1), (lat0, lat1)]:
        if not N.ma.allclose(a0, a1):
            return False
    return True
