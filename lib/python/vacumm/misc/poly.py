# -*- coding: utf8 -*-
"""Low level polygon utilities

"""
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
import numpy as N
from matplotlib.pyplot import gca
from _geoslib import Point, Polygon, LineString

from vacumm import VACUMMError
from .misc import nduniq

GEOS_POINT = GEOS_POINTS = 0
GEOS_LINESTRING = GEOS_LINES = GEOS_LINE = GEOS_POLYLINE = GEOS_POLYLINES = 1
GEOS_POLYGON = GEOS_POLYGONS = GEOS_POLYS = GEOS_POLY = 2

#: GEOS shape classes
SHAPES = (Point, LineString, Polygon)

#:GEOS shape names
SHAPE_NAMES = ['point', 'linestring', 'polygon']

def get_geos_type(gtype, mode='int'):
    """Get the shape type from gtype

    Params
    ------
    gtype: int, string, :class:`Point`, :class:`Polygon`, :class:`LineString`
        Input type to inspect. Possible types are those of mode keyword.
    mode: string, optional
        Output mode:

        - int: 0 = points, 1 = lines, 2 = polygons
        - geos: A :class:`Point`, :class:`Polygon` or :class:`LineString` class.
        - string: one of :data:`SHAPE_NAMES`
    """
    if isinstance(gtype, list):
        gtype = gtype[0]

    # Get it as int
    if isinstance(gtype, int):
        assert gtype in (GEOS_POINT, GEOS_LINE, GEOS_POLYGON)
    elif isinstance(gtype, string):
        att = 'GEOS_'+gtype.upper()
        from . import poly as vcpl
        assert hasattr(vcpl, att)
        gtype = getattr(vcpl, att)
    elif gtype in SHAPES:
        gtype = SHAPES.index(gtype)
    elif isinstance(gtype, SHAPES):
        gtype = SHAPES.index(gtype.__class__)
    else:
        raise VACUMMError('Invalid geos shape type')

    # Return it
    if mode=='int':
        return gtype
    if mode=='geos':
        return SHAPES[gtype]
    if mode=='string':
        return SHAPE_NAMES[gtype]
    raise VACUMMError('Invalid output mode')

def is_point(gobj):
    return get_geos_type(gobj) == GEOS_POINT

def is_linestring(gobj):
    return get_geos_type(gobj) == GEOS_LINESTRING

def is_polygon(gobj):
    return get_geos_type(gobj) == GEOS_POLYGON

def proj_shape(shape, proj):
    """Project a Point, a Polygon or a LineString shape using a geographical projection

    :Params:

        - **shape**: A Point, Polygon, or a LineString instance with coordinates in degrees.
        - **proj**: A projection function of instance.
          See :func:`~vacumm.misc.grid.basemap.get_proj`.

    :Return: A similar instance with its coordinates converted to meters
    """
    if not callable(proj): return shape

    # Point
    if isinstance(shape, Point):
        return Point(*proj(shape.boundary))

    # LineString and Polygon
    return shape.__class__(proj(*shape.boundary.T).T)

def clip_shape(shape, clip=None):
    """Clip a :class:`Point`, a :class:`Polygon` or a :class:`LineString`
    shape using clipping polygon

    :Params:

        - **shape**: A valid :class:`Point`, a :class:`Polygon` or a
          :class:`LineString` instance.
        - **clip**, optional: A clipping polygon.

    :Return: A possible empty list of intersection shapes.
    """
    if clip is None: return [shape]
    clip = create_polygon(clip)
    shapes = []
    if shape.within(clip):
        return [shape]
    if shape.intersects(clip):
        try:
            return shape.intersection(clip)
        except:
            pass
    return shapes

def clip_shapes(shapes, clip=None):
    """Same as :func:`clip_shape` but applied to a list of shapes"""
    clip = create_polygon(clip)
    shapes = []
    for shape in shapes:
        shapes.extend(clip_shape(shape, clip))
    return shapes

def create_polygon(data, proj=False, mode='poly'):
    """Create a simple :class:`Polygon` instance using data

    :Param:

        - **data**: Can be either a ``(N,2)`` or ``(2,N)`` array,
          a ``(xmin,ymin,xmax,ymax)`` list, tuple or array, or a Polygon.
        - **proj**, optional: Geographical projection function.
        - **mode**, optional: Output mode. ``"line"`` = :class:`LineString`,
          ``"verts"`` = numpy vertices, else :class:`Polygon`.

    """
    if isinstance(data, Polygon):
        if callable(proj):
            data = data.boundary
        else:
            return data

    # Convert to numeric
    data = N.asarray(data, 'float64')

    # xmin,ymin,xmax,ymax form
    if data.ndim == 1:
        assert len(data) == 4, '1D form must have 4 elements (xmin,ymin,xmax,ymax), not %i'%len(poly)
        xmin, ymin, xmax, ymax = data
        data =  N.asarray([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])

    # Check order
    if data.shape[0] == 2:
        data = data.T

    # Projection
    if callable(proj):
        xdata, ydata = proj(*data.T)
        data = N.asarray([xdata, ydata]).T

    # Create Polygon or LineString
    mode = str(mode)
    if mode.startswith('v'): return data
    shaper = LineString if mode.startswith('l') else Polygon
    return shaper(data)

def plot_polygon(poly, ax=None, **kwargs):
    """Simply plot a polygon on the current plot using :func:`matplotlib.pyplot.plot`

    :Params:

        - **poly**: A valid :class:`Polygon` instance.
        - Extra keyword are passed to plot function.

    """
    xx = poly.boundary[:, 0]
    yy = poly.boundary[:, 1]
    xx = N.concatenate((xx, xx[:1]))
    yy = N.concatenate((yy, yy[:1]))
    if ax is None:
        ax = gca()
    return ax.plot(xx, yy, **kwargs)


def polygons(polys, proj=None, clip=None, shapetype=2, **kwargs):
    """Return a list of Polygon instances

    :Params:

        - **polys**: A tuple or list of polygon proxies (see examples).
        - **shapetype**: 1 = Polygons, 2=Polylines(=LineStrings) [default: 2]
        - **proj**: Geographical projection to convert positions. No projection
          by default.
        - **clip**, optional: Clip all polygons with this one.

    :Example:

        >>> import numpy as N
        >>> X = [0,1,1,0]
        >>> Y = N.array([0,0,1,1])
        >>> polygons( [(X,Y)]] )                                # from like grid bounds
        >>> polygons( [zip(X,Y), [X,Y], N.array([X,Y])] )       # cloud
        >>> polygons( [polygons(([X,Y],), polygons(([X,Y],)])   # from polygins
        >>> polygons( [[min(X),min(Y),max(X),max(Y)],] )        # from bounds
    """

    # Input
#    if hasattr(polys, 'get_shapes'):
#        polys = polys.get_shapes()
#    elif hasattr(polys, 'getLongitude') and hasattr(polys, 'getLatitude'):
#        polys = (polys.getLongitude()[:], polys.getLatitude()[:])
#        if N.ndim(polys[0])==1:
#            polys = N.meshgrid(*polys)
    if isinstance(polys, (Polygon, N.ndarray, tuple)) or hasattr(polys, 'get_shapes'):
        polys = [polys]

    if kwargs.has_key('m'): proj = kwargs['m']
    kwclip = vcm.kwfilter(kwargs, 'clip_')

    # Numeric or geos polygons
    out_polys = []
#    gmt_polys = []
    if shapetype == 1:
        shaper = LineString
    else:
        shaper = Polygon
    if clip is not None:
        clipl = clip # degrees
        kwclip.setdefault('proj', proj)
        clip = create_polygon(clip, **kwclip)
        del kwclip['proj']
        kwclip['mode'] = 'verts'
    else:
        clipl = None

    # Loop on polygon data
#    from grid import isgrid, isrect, curv2rect, get_grid
#    from axes import islon, islat
    for poly in polys:

#        # Grid (cdms var, grid or tuple of axes) -> get envelop
#        if (cdms2.isVariable(poly) or isgrid(poly) or
#                (isinstance(poly, tuple) and len(poly)==2 and
#                islon(poly[1]) and islat(poly[0]))):
#            grid = get_grid(poly)
#            if grid is None: continue
#            poly = grid_envelop(grid)

        # It's already a polygon
        if isinstance(poly, Polygon):
            pp = [poly]

#        # Polygon from GMT
#        if isinstance(poly, (str, Basemap)):
##            poly = GSHHS_BM(poly)
##           gmt_polys.append(poly)
#            out_polys.extend(GSHHS_BM(poly, clip=clipl, proj=proj).get_shapes())
#            continue

        # Shapes instance
        if hasattr(poly, 'get_shapes'):

            # Good polygon?
            if poly.get_type() == 0 or poly.get_type() != shapetype:
                continue

            # Clip
            poly = poly.clip(clip)

            # Append to list
            out_polys.extend(poly.get_shapes())
            continue

        # Make sure to have a polygon with the right projection
        poly = create_polygon(poly, proj=proj,
            mode='line' if shaper is LineString else 'poly')

        # Clip
        if clip is not None:
            out_polys.extend(clip_shape(poly, clip))
        else:
            out_polys.append(poly)

    return out_polys


def sort_shapes(shapes, reverse=True):
    """Sort shapes according to their surface or length

    Params
    ------
    reverse: optional, boolean
        Reverse the sorted order with smallest first.
    """
    gtype = get_geos_type(shapes)
    if is_polygon(gtype):

        shapes.shapes.sort(cmp=lambda p0, p1: cmp(p0.area(), p1.area()),
            reverse=reverse)
        sorted = 1-2*int(reverse)

    elif is_line(gtype):

        shapes.sort(cmp=lambda p0,
            p1: cmp(N.diff(p0.boundary).sum(), N.diff(p1.boundary).sum()),
            reverse=reverse)
        sorted = 1-2*int(reverse)

    else:

        sorted = 0

def convex_hull(xy, poly=False, method='delaunay'):
    """Get the envelop of cloud of points

    :Params:

        - **xy**: (x,y) or array of size (2,nxy) (see :func:`~vacumm.misc.grid.misc.get_xy`).
        - *poly*:

            - ``True``: Return as Polygon instance.
            - ``False``: Return two 1D arrays ``xpts,ypts``

            - *method*:

                - ``"angles"``: Recursive scan of angles between points.
                - ``"delaunay"``: Use Delaunlay triangulation.

    """

    # Coordinates
    if cdms2.isGrid(xy):
        xx = xy.getLongitude().getValue()
        yy = xy.getLatitude().getValue()
        if xx.ndim==1:
            xx, yy = N.meshgrid(xx, yy)
    else:
        xx, yy = N.asarray(xy)
    if xx.ndim>1:
        xx = xx.ravel()
        yy = yy.ravel()

    # Various methods
    if method.startswith('delau'):

        # Delaunay
        from matplotlib.delaunay import Triangulation
        xy = nduniq(N.asarray([xx, yy]).T)
        hull = Triangulation(xy[:, 0], xy[:, 1]).hull
        xe = xy[:, 0][hull]
        ye = xy[:, 1][hull]

    elif method == 'quarters':

        np = len(xx)

        # Most south point
        ip = N.argmin(yy)
        xe = [xx[ip]]
        ye = [yy[ip]]

        # SW
        while True:
            good = xx>xe[-1]
            if not good.any(): break
            ip = N.argmin(yy[good])
            xe.append(xx[ip])
            ye.append(yy[ip])

        # NW
        while True:
            good = yyx>ye[-1]
            if not good.any(): break
            ip = N.argmax(xx[good])
            xe.append(xx[ip])
            ye.append(yy[ip])

        pass
        #TODO: finish convex_hull with quaters

    elif method=='angles':

        # Angles
        np = len(xx)
        xx0 = N.resize(xx, (np, np))
        xx1 = N.resize(xx, (np, np)).transpose()
        dx = xx1-xx0 ; del xx0, xx1
        yy0 = N.resize(yy, (np, np))
        yy1 = N.resize(yy, (np, np)).transpose()
        dy = yy1-yy0 ; del yy0, yy1
        angles = N.arctan2(dx, dy) ; del dx, dy
        idx = N.arange(np)
        angles[idx, idx] = 10.

        # Most south point
        ip0 = ip = N.argmin(yy)
        xe = [xx[ip]]
        ye = [yy[ip]]

        # Recursive search
        ic = 0
        while True:
            ip = N.argmin(angles[ip])
            if ip == ip0: break
            xe.append(xx[ip])
            ye.append(yy[ip])
            ic += 1
            if ic > np: break
        xe = N.asarray(xe)
        ye = N.asarray(ye)

    else:
        raise NotImplementedError

    # Polygon or positions?
    if poly:
        return create_polygons((xe, ye))
    return xe, ye

