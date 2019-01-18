# -*- coding: utf8 -*-
"""Low level polygon utilities

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
import six
import numpy as N
from matplotlib.pyplot import gca
from matplotlib.collections import PolyCollection, LineCollection
from _geoslib import Point, Polygon, LineString

import cdms2

from vacumm import VACUMMError
from .misc import nduniq, kwfilter, ArgList

__all__ = ['get_shape_type', 'is_point', 'is_linestring', 'is_polygon',
           'proj_shape', 'clip_shapes', 'clip_shapes', 'create_polygon',
           'as_polygons', 'polygons', 'sort_shapes', 'filter_shapes',
           'convex_hull', 'get_shapes_bounds', 'get_convex_hull',
           'undersamp_shapes', 'get_shape_name', 'transform_shape_coords',
           'is_shape_type', 'plot_polygon', 'create_shape',
           'plot_shapes', 'plot_shapes_on_map', 'create_shapes',
           'SHAPE_POINT', 'SHAPE_POINTS', 'SHAPE_LINESTRING', 'SHAPE_LINES',
           'SHAPE_POLYLINE', 'SHAPE_POLYLINES',
           'SHAPE_POLYGON', 'SHAPE_POLYGONS', 'SHAPE_POLYS', 'SHAPE_POLY',
           'SHAPES', 'SHAPE_NAMES']

SHAPE_POINT = SHAPE_POINTS = 0
SHAPE_LINESTRING = SHAPE_LINES = SHAPE_LINE = SHAPE_POLYLINE = \
    SHAPE_POLYLINES = 1
SHAPE_POLYGON = SHAPE_POLYGONS = SHAPE_POLYS = SHAPE_POLY = 2

#: GEOS shape classes
SHAPES = (Point, LineString, Polygon)

#: GEOS shape names
SHAPE_NAMES = ['point', 'linestring', 'polygon']

SHAPE_NAME_ALIASES = {'Point': 'point', 'MultiPoint': 'point',
                      'LineString': 'linestring',
                      'MultiLineString': 'linestring',
                      'Polygon': 'polygon', 'MultiPolygon': 'polygon',
                      }


def get_shape_type(shape_type, mode='int'):
    """Get the shape type from shape_type

    Params
    ------
    shape_type: int, string, :class:`Point`, :class:`Polygon`, :class:`LineString`
        Input type to inspect. Possible types are those of mode keyword.
    mode: string, optional
        Output mode:

        - int: 0 = points, 1 = lines, 2 = polygons
        - geos: A :class:`Point`, :class:`Polygon`
          or :class:`LineString` class.
        - string: one of :data:`SHAPE_NAMES`
    """
    if isinstance(shape_type, list):
        shape_type = shape_type[0]

    # Get it as int
    if isinstance(shape_type, int):
        assert shape_type in (SHAPE_POINT, SHAPE_LINE, SHAPE_POLYGON)
    elif isinstance(shape_type, six.string_types):
        if shape_type in SHAPE_NAME_ALIASES:
            shape_type = SHAPE_NAME_ALIASES[shape_type]
        if shape_type in SHAPE_NAMES:
            shape_type = SHAPE_NAMES.index(shape_type)
        else:
            raise VACUMMError('Invalid shape type: {}'.format(shape_type))
    elif shape_type in SHAPES:
        shape_type = SHAPES.index(shape_type)
    elif isinstance(shape_type, SHAPES):
        shape_type = SHAPES.index(shape_type.__class__)
    else:
        raise VACUMMError('Invalid shape type: {}'.format(shape_type))

    # Return it
    if mode == 'int':
        return shape_type
    if mode == 'geos':
        return SHAPES[shape_type]
    if mode == 'string':
        return SHAPE_NAMES[shape_type]
    raise VACUMMError('Invalid output mode')


def get_shape_name(shape_type):
    """Get the shape name"""
    return get_shape_type(shape_type, mode='string')


def is_point(gobj):
    return get_shape_type(gobj) == SHAPE_POINT


def is_linestring(gobj):
    return get_shape_type(gobj) == SHAPE_LINESTRING


def is_polygon(gobj):
    return get_shape_type(gobj) == SHAPE_POLYGON


def is_shape_type(gobj, gtype):
    return get_shape_type(gobj) == get_shape_type(gtype)


def create_shape(data, shape_type, transform=None):
    """Create a Point, a Polygon or a LineString

    Parameters
    ----------
    data: array_like or tuple
        For a point, a (x, y) tuple or a Point.
        For lines and polygons, a ``(N,2)`` or ``(2,N)`` array,
        a ``(xmin,ymin,xmax,ymax)`` list, tuple or array, or a shape.
    mode: optional
        Output mode. ``"line"`` = :class:`LineString`,
        `"verts"`` = numpy vertices, else :class:`Polygon`.

    """
    if shape_type is None:
        if not isinstance(data, SHAPES):
            raise VACUMMError('Cannot guess shape type to create shape')
        shape_type = data
    shape_type = get_shape_type(shape_type, mode='string')

    # Point
    if shape_type.startswith('point'):
        if isinstance(data, Point):
            return transform_shape_coords(data, transform)
        elif isinstance(data, SHAPES):
            raise VACUMMError('Cannot create a Point with a Polygon '
                              'or a LineString')
        if six.callable(transform):
            data = transform(*data)
        return Point(data)

    # Direct
    if isinstance(data, LineString):
        if shape_type == 'linestring':
            return data
        data = data.boundary
    elif isinstance(data, Polygon):
        if shape_type == 'polygon':
            return data
        data = data.boundary

    # From basemap instance
    if hasattr(data, 'llcrnrlon'):
        data = [data.llcrnrlon, data.llcrnrlat,
                data.urcrnrlon, data.urcrnrlat]

    # From coords
    data = N.asarray(data, 'float64')

    # From bounds (xmin, ymin, xmax, ymax)
    if data.ndim == 1:
        assert len(data) == 4, ('1D form must have 4 elements '
                  '(xmin,ymin,xmax,ymax), not %i') % len(data)
        xmin, ymin, xmax, ymax = data
        data = N.asarray([[xmin, ymin], [xmax, ymin],
                          [xmax, ymax], [xmin, ymax]])

    # Check order
    if data.shape[0] == 2 and data.shape[1] != 2:
        data = data.T

    # Transform
    if six.callable(transform):
        data = transform(*data.T).T

    # Create shape
    return SHAPES[SHAPE_NAMES.index(shape_type)](data)


def create_polygon(data, mode='poly', transform=None):
    """Create a simple :class:`Polygon` instance using data

    .. note:: This function is deprecated.

    Parameters
    ----------
    data:
        Can be either a ``(N,2)`` or ``(2,N)`` array,
        a ``(xmin,ymin,xmax,ymax)`` list, tuple or array, a Polygon,
        or a :class:`~mpl_toolkits.basemap.Basemap` instance.
    proj: optional
        Geographical projection function.
    mode: optional
        Output mode. ``"line"`` = :class:`LineString`,
        ``"verts"`` = numpy vertices, else :class:`Polygon`.

    """
    # Create the shape
    shape_type = 'linestring' if mode.startswith('l') else 'polygon'
    shape = create_shape(data, shape_type, transform=transform)

    # Output
    if mode.startswith('v'):
        return shape.boundary
    return shape


def create_shapes(data, shape_type=None, transform=None, clip=None, samp=None):
    """Create a list of shapes
    """
    # Get a list of shape compatible data
    if shape_type is not None:
        shape_type = get_shape_type(shape_type, mode='string')
    if hasattr(data, 'get_shapes'):  # Shapes object
        data = data.get_shapes()
    else:  # Other
        if (isinstance(data, SHAPES + (tuple, )) or
                (isinstance(data, N.ndarray) and
                 (shape_type != 'point' or data.ndim==1)) and
                (isinstance(data, list) and N.isscalar(data[0]))):
            data = [data]

    # Create shapes
    shapes = [create_shape(dat, shape_type, transform=transform)
        for dat in data]

    # Undersample
    if samp is not None and samp > 1:
        shapes = undersamp_shapes(shapes, samp)

    # Clip
    if clip is not None:
        shapes = clip_shapes(shapes, clip)

    return shapes


def as_polygons(polys, clip=None, shapetype=2, **kwargs):
    """Return a list of Polygon instances

    .. note:: This function is deprecated. Please use :func:`create_shapes`
        instead, ``shape_type='polygons``.

    Parameters
    ----------
    polys:
        A tuple or list of polygon proxies (see examples).
    shapetype:
        2=Polygons, 1=Polylines(=LineStrings) [default: 2]
    proj:
        Geographical projection to convert positions. No projection
        by default.
    clip: optional
        Clip all polygons with this one.

    Example
    -------
    >>> import numpy as N
    >>> X = [0,1,1,0]
    >>> Y = N.array([0,0,1,1])
    >>> as_polygons( [(X,Y)]] )                                # from like grid bounds
    >>> as_polygons( [zip(X,Y), [X,Y], N.array([X,Y])] )       # cloud
    >>> as_polygons( [polygons(([X,Y],), polygons(([X,Y],)])   # from polygins
    >>> as_polygons( [[min(X),min(Y),max(X),max(Y)],] )        # from bounds
    """
    return create_shapes(polys, shape_type=shapetype, clip=clip, **kwargs)

#    # Input
##    if hasattr(polys, 'get_shapes'):
##        polys = polys.get_shapes()
##    elif hasattr(polys, 'getLongitude') and hasattr(polys, 'getLatitude'):
##        polys = (polys.getLongitude()[:], polys.getLatitude()[:])
##        if N.ndim(polys[0])==1:
##            polys = N.meshgrid(*polys)
#    if isinstance(polys, (Polygon, N.ndarray, tuple)) or hasattr(polys, 'get_shapes'):
#        polys = [polys]
#
#    if 'm' in kwargs:
#        proj = kwargs['m']
#    kwclip = kwfilter(kwargs, 'clip_')
#
#    # Numeric or geos polygons
#    out_polys = []
##    gmt_polys = []
#    if shapetype == 1:
#        shaper = LineString
#    else:
#        shaper = Polygon
#    if clip is not None:
#        #        clipl = clip # degrees (backup)
#        kwclip.setdefault('proj', proj)
#        clip = create_polygon(clip, **kwclip)
#        del kwclip['proj']
#        kwclip['mode'] = 'verts'
##    else:
##        clipl = None
#
#    # Loop on polygon data
##    from grid import isgrid, isrect, curv2rect, get_grid
##    from axes import islon, islat
#    for poly in polys:
#
#        #        # Grid (cdms var, grid or tuple of axes) -> get envelop
#        #        if (cdms2.isVariable(poly) or isgrid(poly) or
#        #                (isinstance(poly, tuple) and len(poly)==2 and
#        #                islon(poly[1]) and islat(poly[0]))):
#        #            grid = get_grid(poly)
#        #            if grid is None: continue
#        #            poly = grid_envelop(grid)
#
#        # It's already a polygon
#        if isinstance(poly, Polygon):
#            poly = [poly]
#
##        # Polygon from GMT
##        if isinstance(poly, (str, Basemap)):
###            poly = GSHHS_BM(poly)
## gmt_polys.append(poly)
##            out_polys.extend(GSHHS_BM(poly, clip=clipl, proj=proj).get_shapes())
##            continue
#
#        # Shapes instance
#        if hasattr(poly, 'get_shapes'):
#
#            # Good polygon?
#            if poly.get_type() == 0 or poly.get_type() != shapetype:
#                continue
#
#            # Clip
#            poly = poly.clip(clip)
#
#            # Append to list
#            out_polys.extend(poly.get_shapes())
#            continue
#
#        # Make sure to have a polygon with the right projection
#        poly = create_polygon(poly, proj=proj,
#                              mode='line' if shaper is LineString else 'poly')
#
#        # Clip
#        if clip is not None:
#            out_polys.extend(clip_shape(poly, clip))
#        else:
#            out_polys.append(poly)
#
#    return out_polys


polygons = as_polygons


def plot_shapes(shapes, ax=None, transform=None, fill=None, points=None,
                color=None, facecolor=None, edgecolor=None, linewidth=None,
                alpha=1, s=None,
                **kwargs):
    """Basic plot a single or a list of Polygons, LineStrings or Points

    Points are plotted with :func:` ~matplotlib.pyplot.scatter`,
    linestrings are plotted with :func:` ~matplotlib.pyplot.plot`,
    and polygons are plotted with :func:` ~matplotlib.pyplot.fill`

    Parameters
    ----------
    shapes: list
        List of shapes
    transform: callable
        Function to transform coordinates
    fill: bool or None
        Fill polygons or even lines?

    """
    # Inits
    if ax is None:
        ax = gca()
    if not isinstance(shapes, (list, tuple)):
        shapes = [shapes]
    shape_type = get_shape_type(shapes[0], 'string')
    kwpoints = kwfilter(kwargs, 'points_')
    kwlines = kwfilter(kwargs, 'lines_')
    kwfill = kwfilter(kwargs, 'fill_')
    facecolor = kwargs.pop('fillcolor', facecolor) or facecolor  # compat
    if facecolor is not None:
        kwfill.setdefault('facecolor', facecolor)
        kwpoints.setdefault('c', facecolor)
    if edgecolor is not None:
        for kw in kwfill, kwpoints, kwlines:
            kw.setdefault('edgecolor', edgecolor)
    if color is not None:
        kwpoints.setdefault('c', color)
        kwpoints.setdefault('edgecolors', color)
        kwlines.setdefault('color', color)
        kwfill.setdefault('edgecolor', color)
        kwfill.setdefault('facecolor', color)
    if s is not None:
        kwpoints.setdefault('s', s)
    if linewidth is not None:
        kwpoints.setdefault('linewidths', linewidth)
        kwlines.setdefault('linewidth', linewidth)
        kwfill.setdefault('linewidth', linewidth)
    kwlines.setdefault('linestyle', 'solid')
    kwlines.setdefault('alpha', alpha)
    kwpoints.setdefault('alpha', alpha)
    kwfill.setdefault('alpha', alpha)
    for kw in kwlines, kwpoints, kwfill:
        kw.update(kwargs)

    # Coordinates
    if shape_type == 'point':
        coords = N.asarray([p.boundary for p in shapes]).T
        if transform:
            coords = transform(*coords)
    else:
        coords = []
        for shape in shapes:
            x = shape.boundary[:, 0]
            y = shape.boundary[:, 1]
            if transform:
                x, y = transform(x, y)
            coords.append(N.transpose([x, y]))
#            coords.append(x)
#            coords.append(y)

    # Lines and polygons
    if shape_type != 'point':
        if shape_type == 'linestring' or fill is False:
            collec = LineCollection
            kw = kwlines
        else:
            collec = PolyCollection
            kw = kwfill
        ax.add_collection(collec(coords, **kw))
        ax.autoscale_view()
#        func(*coords, **kw)
#        func(coords[0], coords[1])

    # Points
    if shape_type == 'point' or points:
        if shape_type != 'point':
            coords = (N.concatenate(coords[::2]), N.concatenate(coords[1::2]))
        ax.scatter(coords[0], coords[1], **kwpoints)


def plot_shapes_on_map(shapes, m, **kwargs):
    """Plot shapes on an existing Basemap instance"""
    if six.callable(m):
        if hasattr(m, 'map'):
            m = m.map
        kwargs.update(ax=m.ax, transform=m)
    return plot_shapes(shapes, **kwargs)


def plot_polygon(poly, ax=None, **kwargs):
    """Simply plot a polygon on the current plot using
    :func:`matplotlib.pyplot.plot`

    .. note:: Please use the :func:`plot_shapes` for a more generic
        plotting function.

    Parameters
    ----------
    poly: :class:`Polygon`
    **kwargs:
        Extra keyword are passed to plot function.

    """
    xx = poly.boundary[:, 0]
    yy = poly.boundary[:, 1]
    xx = N.concatenate((xx, xx[:1]))
    yy = N.concatenate((yy, yy[:1]))
    if ax is None:
        ax = gca()
    return ax.plot(xx, yy, **kwargs)


def sort_shapes(shapes, key=None, reverse=True):
    """Sort shapes according to their surface or length

    Params
    ------
    key: None or callable
        If None, sorting is performed with key :func:`get_shape_type`
    reverse: optional, boolean
        Reverse the sorting order with smallest first.
    """
    gtype = get_shape_type(shapes)

    if is_point(gtype):
        return 0

    shapes.sort(key=get_shape_size, reverse=reverse)

    return 1


def get_shape_size(shape):
    """Get the size of shape, depending on its type

    - A Point: 0.
    - A Linestring: total length
    - A Polygon: area

    Parameter
    ---------
    shape: a Geos shape

    Return
    ------
    float
        The "size"
    """
    gtype = get_shape_type(shape)

    if is_point(gtype):
        return 0

    if is_polygon(gtype):
        return shape.area()

    return N.diff(shape.boundary).sum()


def get_convex_hull(xy, poly=False, method='delaunay'):
    """Get the envelop of cloud of points

    Parameters

        xy**: (x,y) or array of size (2,nxy) (see :func:`~vacumm.misc.grid.misc.get_xy`).
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
        if xx.ndim == 1:
            xx, yy = N.meshgrid(xx, yy)
    else:
        xx, yy = N.asarray(xy)
    if xx.ndim > 1:
        xx = xx.ravel()
        yy = yy.ravel()

    # Various methods
    if method.startswith('delau'):

        # Delaunay
        from scipy.spatial import ConvexHull
        xy = N.unique(N.vstack((xx, yy)).T, axis=0)
        hull = ConvexHull(xy)
        xe = xy[hull.vertices, 0]
        ye = xy[hull.vertices, 1]
        del xy

    elif method == 'quarters':

        np = len(xx)

        # Most south point
        ip = N.argmin(yy)
        xe = [xx[ip]]
        ye = [yy[ip]]

        # SW
        while True:
            good = xx > xe[-1]
            if not good.any():
                break
            ip = N.argmin(yy[good])
            xe.append(xx[ip])
            ye.append(yy[ip])

        # NW
        while True:
            good = yy > ye[-1]
            if not good.any():
                break
            ip = N.argmax(xx[good])
            xe.append(xx[ip])
            ye.append(yy[ip])

        pass
        # TODO: finish convex_hull with quaters

    elif method == 'angles':

        # Angles
        np = len(xx)
        xx0 = N.resize(xx, (np, np))
        xx1 = N.resize(xx, (np, np)).transpose()
        dx = xx1-xx0
        del xx0, xx1
        yy0 = N.resize(yy, (np, np))
        yy1 = N.resize(yy, (np, np)).transpose()
        dy = yy1-yy0
        del yy0, yy1
        angles = N.arctan2(dx, dy)
        del dx, dy
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
            if ip == ip0:
                break
            xe.append(xx[ip])
            ye.append(yy[ip])
            ic += 1
            if ic > np:
                break
        xe = N.asarray(xe)
        ye = N.asarray(ye)

    else:
        raise NotImplementedError

    # Polygon or positions?
    if poly:
        return create_polygon((xe, ye))
    return xe, ye


convex_hull = get_convex_hull


def get_shapes_bounds(shapes):
    """Get the coordinates bounds of shapes


    Parameters
    ----------
    shapes: list of shapes or shape
        A single shape or a list of themes.
        A shape must have the :attr:`boundary` attribute.

    Return
    ------
    tuple
        ``xmin,ymin,xmax,ymax``
    """
    al = ArgList(shapes)
    ss = al.get()
    xmin = min([shape.boundary[..., 0].min() for shape in ss])
    xmax = max([shape.boundary[..., 0].max() for shape in ss])
    ymin = min([shape.boundary[..., 1].min() for shape in ss])
    ymax = max([shape.boundary[..., 1].max() for shape in ss])
    return xmin, ymin, xmax, ymax


def clip_shape(shape, clip=None):
    """Clip a :class:`Point`, a :class:`Polygon` or a :class:`LineString`
    shape using clipping polygon

    Parameters
    ----------
    shape: :class:`Point`, a :class:`Polygon`, :class:`LineString`
    clip: optional
        A clipping polygon or a ``[xmin, ymin, xmax, ymax]``.

    Return
    ------
    A possible empty list of intersection shapes.
    """
    if clip is None:
        return [shape]
    clip = create_polygon(clip)
    shapes = []
    if shape.within(clip):
        return [shape]
    if shape.intersects(clip):
        try:
            return shape.intersection(clip)
        except Exception as e:
            print('error !', e.message)
            pass
    return shapes


def clip_shapes(shapes, clip=None):
    """Same as :func:`clip_shape` but applied to a list of shapes"""
    if clip is None:
        return list(shapes)
    clip = create_polygon(clip)
    oshapes = []
    for shape in shapes:
        oshapes.extend(clip_shape(shape, clip))
    return oshapes


def transform_shape_coords(shape, transform=None, **kwargs):
    """Transform shape coordinates

    A typical transformation is for converting from degrees to meters.

    Parameters
    ----------
    shape: Point, LineString or Polygon
    transform: callable
        Must be comptatible with: ``x, y = transform(x, y, **kwargs)``

    Return
    ------
    Point, LineString or Polygon
    """
    if not six.callable(transform):
        return shape

    # Point
    if isinstance(shape, Point):
        return Point(*transform(shape.boundary))

    # LineString and Polygon
    return shape.__class__(N.transpose(transform(*shape.boundary.T)))


proj_shape = transform_shape_coords


def undersamp_shapes(shapes, samp):
    """Undersample shapes


    Parameters
    ----------
    shapes: list of geos shapes or shape
        A single shape or a list of themes.
        A shape must have the :attr:`boundary` attribute.
    samp: int

    Return
    ------
    geos shape or list of them
    """
    al = ArgList(shapes)
    ss = al.get()

    if not shapes:
        return al.put([])
    if not samp:
        return al.put(list(ss))

    gtype = get_shape_type(shapes, mod='geos')

    outss = []
    for shape in ss:

        if samp > 1 and shape.shape[0] > (2*samp+1):
            shape = gtype(shape.coord[::samp])

        outss.append(shape)

    return al.put(outss)


def filter_shapes(shapes, checker, transform=None):
    """Select shapes that satisfy the checker function

    For each selected shape: ``checker(shape) is True``


    Example
    -------
    >>> checker = lambda shape: shape.area() > 10e3 * 103
    >>> shapes = filter_shapes(shapes, checker, transform=proj)
    """
    al = ArgList(shapes)
    ss = al.get()

    if not shapes:
        return al.put([])

    outss = []
    for shape in ss:
        if transform:
            shape = transform_shape_coords(shape, transform)
        if checker(shape):
            outss.append(shape)

    return al.put(outss)
