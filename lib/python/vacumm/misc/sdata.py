# -*- coding: utf8 -*-
"""High level spatial data tools"""
# Copyright or © or Copr. Actimar/IFREMER (2016-2016)
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

from .misc import kwfilter
from .poly import clip_shape
from .io import read_shapefile
from .basemap import get_proj
from .grid import resol
from .regridding import xy2xy, regrid2d, regrid_method
from .color import get_cmap, land
from .core_plot import Map
from .plot import map2, _colorbar_

class Shapes(object):
    """A class to read shapefiles and return GEOS objects
    Inspired from basemap.readshapefile

    Here are the conversion rules from shapefile to GEOS objects :

        - Points and multipoints are interpreted as :class:`Points`.
        - Polylines are interpreted as :class:`LineString`.
        - Polygons are interpreted as :class:`Polygons`.

    :Params:

        - **input**: Refers to a shapefile or is a shapes isntance ;
          if a shapefile, it assumes that <input>.shp contains points,
          multipoints, lines or polygons, and that <input>.dbf contains their attributes.
        - **proj*, optional: A projection function to convert coordinates. It must accept
          the "inverse" keyword.
        - **m*, optional: A Basemap instance for converting for plotting.
        - *inverse*, optional: Inverset the conversion with proj .
        - **clip*, optional: If in the form ``(xmin,ymin,xmax,ymax)``,
          clips to this box ; if a polygon like argument,
          it clips to this polygon
          (see :func:`~vacumm.misc.grid.masking.polygons()` for arguments).
          If simply ``True`` and *m* is present, it clips to the bounds of *m*.
        - **min_area*, optional: Minimal area to keep a polygon
        - **samp**, optional: An integer to undersample coordinates of polygons and lines.
        - **shapetype**, optional:

            - If 0, it must only deal with points ;
            - if 1, only polylines ;
            - if 2, only polygons (conversion 1<->2 is automatic).
    """
    POINTS = POINT = 0
    LINES = LINE = 1
    POLYGONS = POLYS = POLY = POLYGON = 2
    INPUT_POINTS = 1
    INPUT_MULTIPOINTS = 8
    INPUT_POLYLINES = 3
    INPUT_POLYGONS = 5
    def __init__(self, input, m=None, proj=False, inverse=False, clip=True,
            shapetype=None, min_area=None, sort=True, reverse=True, samp=1,
            clip_proj=True):

        # Load data
        sdata = read_shapefile(input, m=m, proj=proj, inverse=inverse, clip=clip,
            shapetype=shapetypre, min_area=min_area, sort=sort, samp=samp,
            clip_proj=clip_proj)

        # Set attributes
        for key in ['shapes', 'shaper', 'proj', 'sorted', 'm', 'm_projsync']:
            setattr('_'+key, sdata[key])
        for key in ['xmin', 'xmax', 'ymin', 'ymax', 'xpmin',
                'xpmax', 'ypmin', 'ypmax']:
            setattr(key, sdata[key])


    def clip(self, zone, copy=True, sort=True, reverse=True, **kwargs):
        """Clip to zone

        :Params:

            - **zone**: ``[xmin, ymin, xmax, ymax]``
            - **copy**, optional: If ``True``, make a copy of current instance,
              else simply rehandled the list of shapes.
            - If ``copy==True``, Other parameters are passed to the initialization
              of the new instance.
        """
        if not copy:
            zone = self._clip_zone_(zone)
            if zone is None: return self
            if self._type > 0:
                newshapes = []
                for shape in self:
                    if zone.intersects(shape):
                        intersections = shape.intersection(zone)
                        newshapes.extend(filter(lambda s: isinstance(s,  self._shaper), intersections))
                self._shapes = newshapes
                if sort: self.sort(reverse=reverse)
            return self
        return self.__class__(self, clip=zone, **kwargs)

#   def join_lines(self, poly=True):
#       # Array of end points
#       nsh = len(self)
#       ends = N.zeros((2*nsh, 3))
#       for i in xrange(nsh):
#           ends[i, 0:1] = self._shapes[i][0]
#           ends[i+1, 0:1] = self._shapes[i][-1]
#           ends[i:i+2, 2] = float(i)
#       # Distances
#       xx = ends[:, 0].reshape((nsh*2, nsh*2))
#       yy = ends[:, 1].reshape((nsh*2, nsh*2))
#       dst = (xx-xx.transpose())**2+(yy-yy.transpose())**2
#       new_shapes = []


    def __len__(self):
        return len(self._shapes)

    def __getitem__(self, key):
        return self._shapes[key]

#    def project(self, proj=True, inverse=False):
#        """Project shapes using proj
#
#        - **proj**: A Basemap instance or pure projection instance (from Basemap).
#          If True, a mercator projection is used.
#        - *inverse*: Inverse projection [default: False]
#        """
#        if proj is True:
#            proj=merc(lon=(self.xmin, self.xmax), lat=(self.ymin, self.ymax))
#        for i, shape in enumerate(self):
#            bb = shape.boundary
#            if self._type == 0:
#                bb = proj(b[0], b[1], inverse=inverse)
#            else:
#                bb[:, 0], bb[:, 1] = proj(bb[:, 0], bb[:, 1], inverse=inverse)
#            self._shapes[i] = self._shaper(bb)
#        if isinstance(proj, Basemap): proj = proj.projtran
#        self._proj.append((proj, inverse))

    def sort(self, reverse=True):
        """Sort shapes according to their surface or length

        - *reverse*: If True, greater polygons are first [default: True]
        """
        self._sorted = sort_shapes(self._shapes, reverse=reverse)

    @property
    def sorted(self):
        return self._sorted

    def get_type(self):
        """Return the type of shapes

        - :attr:`POINTS` = Points,
        - :attr:`LINES` = LineStrings = PolyLines
        - :attr:`POLYGONS` = Polygons
        """
        return self._type
    type = property(fget=get_type)

    def is_type(self, type):
        """Check type

        :Example:

            >>> self.is_type(self.POLYS)
        """
        return self.get_type()==type

    def _get_proj_(self, proj=None):
        """Get a valid projection function to operate on shapes coordinates"""
        if proj is None: return
        if proj is True or isinstance(proj, basestring):

            if hasattr(self, '_proj') and callable(self._proj): # already projected
                return

            if proj is True and callable(self._m): # from map

                return self._m

            # using grid.basemap.get_proj
            if self.xmin>self.xmax:
                gg = None
            else:
                gg = ([self.xmin,self.xmax],
                    N.clip([self.ymin,self.ymax], -89.99, 89.99))
            kw = dict(proj=proj) if isinstance(proj, basestring) else {}
            return get_proj(gg, **kw)

        if callable(self._proj): # already projected

            if proj is False: # no projection -> project back
                return lambda x, y, inverse=False: self._proj(x, y, inverse=not inverse)

            if proj is not self._proj: #re-projection
                old_proj = proj
                def proj(x, y, inverse=False):
                    if inverse:
                        return self._proj(*(old_proj(x, y, True)+(False, )))
                    return old_proj(*self._proj(x, y, True))

            else: # no need to project
                proj = None

        return proj



    def get_shapes(self, key=None, proj=None):
        """Get the list of geos objects (polygons, etc)

        :Param:

            - **key**: A slice selector applied to the list.
            - **proj**: ``True``, or a callable to project or re-project coordinates.
        """
        # Objects to work on
        if key is None:
            shapes = self._shapes
        else:
            shapes = self._shapes[key]
        single = not isinstance(shapes, list)
        if single: shapes = [shapes]

        # Projection
        proj = self._get_proj_(proj)

        # Loop on shapes
        polys = []
        for poly in shapes:
            if proj:
                poly = proj_shape(poly, proj)
            polys.append(poly)

        if single: return polys[0]
        return polys

    def get_data(self, key=None, proj=None):
        """Get the numeric version of the list of geos objects (polygons, etc)

        :Param:

            - **key**: A slice selector applied to the list.
            - **proj**: ``True``, or a callable to project or re-project coordinates.
        """
        if not len(self): return []
        # Projection
        proj = self._get_proj_(proj)

        # Shapes to work on
        if key is None:
            shapes = self._shapes
        else:
            shapes = self._shapes[key]
        single = not isinstance(shapes, list)
        if single: shapes = [shapes]

        # Loop on shapes
        data = []
        for poly in shapes:
            xy = poly.boundary
            if proj:
                if self.is_type(self.POINTS): # Points
                    xy = N.array(proj(*xy))
                else: # Lines
                    xx, yy = proj(*xy.T)
                    xy = N.asarray([xx, yy]).T
                    del xx, yy
            data.append(xy)

        if single: return data[0]
        return data

    def get_points(self, key=None, split=True, proj=None):
        """Get all the points from all the shapes as a tuple (x,y)"""
        if not len(self): return N.array([[],[]])

        # Projection
        proj = self._get_proj_(proj)

        # Shapes to work on
        if key is None:
            shapes = self._shapes
        else:
            shapes = self._shapes[key]
        single = not isinstance(shapes, list)
        if single: shapes = [shapes]

        # Loop in shapes
        xx, yy = [], []
        for poly in shapes:
            xy = poly.boundary
            if self.is_type(self.POINTS):
                xx.append(xy[0])
                yy.append(xy[1])
            else:
                xx.extend(xy[:, 0].tolist())
                yy.extend(xy[:, 1].tolist())
        if split:
            if proj: return proj(N.asarray(xx), N.asarray(yy))
            return N.asarray(xx), N.asarray(yy)
        if proj: return N.asarray(proj(xx, yy))
        return N.asarray([xx, yy])

    def get_xy(self, key=None, proj=None):
        """Shortcut to ``get_points(split=false)``"""
        return self.get_points(split=False, key=key, proj=proj)
    xy = property(get_xy, doc="XY coordinates as a (2,npts) array")


    def resol(self, deg=True):
        """Compute the mean "resolution" of the shapes based on the first shape

        - **deg**:

            - if ``False``: return a resolution in meters has a the median distance between points
            - if ``True``: return the median distance between points as a resolution in degrees ``(xres,yres)``

        """
        if not len(self): return 0,0
        x, y = self.get_xy(key=0)
        if deg and callable(self._proj): # m->deg
            dx, dy = resol((x, y), proj=False)
            x0 = x.mean()
            y0 = y.mean()
            x1, y1 = self._proj(x0+dx, y0+dx, inverse=True)
            return x1-x0, y1-y0
        elif not deg and not callable(self._proj):
            return resol((x, y), proj=True)
        return resol((x, y), proj=False)


    def get_map(self):
        """Return the associated basemap instance if set"""
        if hasattr(self._m, 'map'): return self._m.map
        return self._m

    def plot(self, select=None, ax=None, fill=None, points=False, lines=True,
            fillcolor=None, color='k', s=None, linewidth=None, m=None, show=True,
            alpha=1, autoscale=True, title=None, **kwargs):
        """Plot shapes

        :Params:

            - **select**, optional: argument for selecting shapes in the list [defaul: None].
            - **fill**, optional: Force filling (True/False), else guessed from shpe type, ie filling for polygons only [default: None]
            - **ax**, optional: Axes instance.
            - **m**, optional: :class:`~vacumm.misc.core_plot.Map` instance
              (created with :func:`~vacumm.misc.plot.map2`) or a :class:`~mpl_toolkits.basemap.Basemap` instance.
            - **points**, optional: Plots shapes as points.
            - **lines**, optional: Plot shapes as lines (if a of type :attr:`POINTS`).
            - **fill_<params>**, optional: ``<param>`` is passed to
              :class:`~matplotlib.collections.PolyCollection`.
            - **lines_<params>**, optional: ``<param>`` is passed to
              :class:`~matplotlib.collections.LineCollection` or to :class:`~matplotlib.collections.PolyCollection`.
            - **points_<params>**, optional: ``<param>`` is passed to
              :class:`~matplotlib.pyplot.scatter`.
            - **m_<params>**, optional: ``<param>`` is passed to
              :class:`~vacumm.misc.plot.map2` if ``m is True``.
            - **autoscale**, optional: Autoscale axis limits?
        """

        # Keywords
        kwpoints = kwfilter(kwargs, 'points')
        kwlines = kwfilter(kwargs, 'lines')
        kwpoints.setdefault('c', color)
        if s is not None: kwpoints.setdefault('s', s)
        if linewidth is not None:
            kwlines.setdefault('linewidth', linewidth)
        kwlines.setdefault('linestyles', 'solid')
        kwlines.setdefault('color', color)
        kwlines.setdefault('alpha', alpha)
        kwpoints.setdefault('alpha', alpha)
        kwfill = kwfilter(kwargs, 'fill')
        kwfill.update(kwlines)
        if fill is None and self.is_type(self.POLYGONS): fill = True


        # Map
        if m is True or m=='auto':
            if m!='auto' and getattr(self, '_m', None):
                m = self._m
            else:
                m = Map.get_current(axes=ax) or True
        if m is True:
            kwmap = kwfilter(kwargs,'m_')
            if not len(self):
                warn('No shape found, thus nothing to plot')
            else:
                if not kwmap.has_key('lon'):
                    dx = (self.xmax-self.xmin)*.05
                    kwmap['lon'] = (self.xmin-dx, self.xmax+dx)
                if not kwmap.has_key('lat'):
                    dy = (self.ymax-self.ymin)*.05
                    kwmap['lat'] = (self.ymin-dy, self.ymax+dy)
            kwmap.setdefault('res', None)
            kwmap.setdefault('proj', 'merc')
            kwmap.update(show=False, axes=ax, title=title)
            m = map2(**kwmap)
            ax = m.axes
        isbm = isinstance(m, Basemap)

        # Plot on what?
        if ax is None:
            ax = getattr(m, 'axes',  None) or P.gca()

        # Plot what?
        if points:
            xx, yy = self.get_points(split=True, proj=m)
        if not self.is_type(self.POINTS):
            data = self.get_data(select, proj=m)

        # Polygons or lines
        oo = []
        if not self.is_type(self.POINTS):
            if (fill is None and self.is_type(self.POLYGONS)) or fill is True: # Polygons
                if fillcolor is None: fillcolor=land
                for kv in dict(facecolor=fillcolor).items():
                    kwfill.setdefault(*kv)
                cc = PolyCollection(data, **kwfill)
            else:
                cc = LineCollection(data, **kwlines)
            ax.add_collection(cc)
            oo.append(cc)
            if isbm:
                m.set_axes_limits(ax=ax)

        # Points
        if points:
            cc = ax.scatter(xx, yy, **kwpoints)
            oo.append(cc)
            if isbm:
                m.set_axes_limits(ax=ax)

        # Special properties
        for key in ['label']:
            if key in kwargs and hasattr(cc, 'set_'+key):
                getattr(cc, 'set_'+key)(kwargs[key])

        # Finalize
        if title:
            ax.set_title(title)
        if autoscale:
            ax.autoscale_view()

        if show: P.show()

        return oo


class GSHHSBM(Shapes):
    """Shoreline from USGS using Basemap

    Initialized with a valid Basemap instance with resolution not equal to None,
    or thanks to arguments passed to :func:`create_mapplot.map`

    - *input*: Basemap or Shapes instance [default: None]

    """
    def __init__(self, input=None, clip=None, sort=True, reverse=True, proj=False, **kwargs):

        # From another Shapes instance
        if isinstance(input, Shapes):
            if clip is None: # m is not None and
                clip = [input.xmin, input.ymin, input.xmax, input.ymax]
            input = input._m

        # Get the map without projection
        if not isinstance(input, Basemap):

            # Base to create the map
            kwmap = kwargs.copy()
            if isinstance(input, str):
                kwmap['res'] = input
            elif isinstance(input, dict):
                kwmap.update(input)

            # Clipping zone
            if clip is not None:

                # Vertices
                clip = create_polygon(clip, mode='verts')


                # Map extension from clip bounds
                kwmap.setdefault('lon_min', clip[:, 0].min())
                kwmap.setdefault('lon_max', clip[:, 0].max())
                kwmap.setdefault('lat_min', clip[:, 1].min())
                kwmap.setdefault('lat_max', clip[:, 1].max())

            # Default resolution is 'i' if nothing to estimate it
            if not 'res' in kwmap and not 'resolution' in kwmap and \
                (   (not 'lon' in kwmap and
                        (not 'lon_min' in kwmap or not 'lon_max' in kwmap)) or
                    (not 'lat' in kwmap and
                        (not 'lat_min' in kwmap or not 'lat_max' in kwmap))):
                kwmap['res'] = 'i'

            # Check lats
            if 'lat_min' in kwmap: kwmap['lat_min'] = max(kwmap['lat_min'], -89.99)
            if 'lat_max' in kwmap: kwmap['lat_max'] = min(kwmap['lat_max'], 89.99)

            # Build the map
            m = create_map(**kwmap)
            self.res = m.resolution

        else:

            clip = False
            m = input

        # Get unprojected polygons vertices
        self.res = m.resolution
        assert m.resolution is not None, 'Your map needs its resolution to be set'
        all_verts = []
        for i, verts in enumerate(m.coastpolygons):
            if m.coastpolygontypes[i] in [2,4]: continue # Skip lakes
            if m.projection!='cyl': # Project back
                verts = m.projtran(verts[0], verts[1], inverse=True)
            all_verts.append(N.asarray(verts).T)

        # Final initialization
        Shapes.__init__(self, all_verts, m=m, clip=clip, sort=sort, proj=proj,
            shapetype=Shapes.POLYGON, **kwargs)

GSHHS_BM = GSHHSBM


class XYZ(object):
    """Class to manipulate xyz data (randomly spaced)

    - **xyz**: It can be either

        - a .xyz ascii file, or a netcdf/grd file with variables ``x``, ``y`` and ``z``,
        - a (x,y,z) tuple,
        - a (3, npts) array,
        - another XYZ instance.

    - *long_name*: Long name
    - *units* Units
    - *tranform*: It can be either

        - a factor applied to z at initialisation
        - a fonction that takes z as the only argument to filter its data.

    - *exc*: Polygons to exclude data (see :meth:`exclude`).
             Several polygons must be passed as a tuple (poly1, poly2, ...).
    - *sel*: Polygons to select data (see :meth:`select`).
             Several polygons must be passed as a tuple (poly1, poly2, ...).
    - load_<keywords>: keywords are passed to :func:`numpy.loadtxt`
    - rsamp_<keywords>: keywords are passed to :func:`rsamp`
    - Other keywords are set as atrributes.


    Slicing:

    - len(obj): number of xyz data
    - obj[1:3]: [(x0,y0,z0),(x1,y1,z1)]

    Operations :

    >>> xyz += 2
    >>> xyz3 = xyz1 + xyz2/2. # concatenation
    """

    def __init__(self, xyz, m=None,  units=None, long_name=None, transp=True,
        trans=False, magnet=0, rsamp=0, id=None, **kwargs):
        # Load data
        self._selections = []
        self._exclusions = []
        if isinstance(xyz, XYZ):
            x = xyz._x.copy()
            y = xyz._y.copy()
            z = xyz._z.copy()
            if units is None: units = xyz.units
            if long_name is None: long_name = xyz.long_name
#           if m is None: m = XYZ._m
            self._selections = copy.deepcopy(xyz._selections)
            self._exclusions = copy.deepcopy(xyz._exclusions)
        elif hasattr(xyz, 'xyz'):
            xyz = xyz.xyz
        elif isinstance(xyz, (tuple, N.ndarray)):
            # Direct
            x, y, z = xyz
        elif os.path.exists(xyz):
            # Read from file
            if xyz.endswith('.nc') or xyz.endswith('.grd'):
                # Netcdf or grid file
                f = cdms2.open(xyz)
                x = f('x').filled()
                y = f('y').filled()
                z = f('z').filled()
                f.close()
            else:
                # Ascii file
                data = N.loadtxt(xyz, **kwfilter(kwargs, 'load'))
                x = data[:, 0]
                y = data[:, 1]
                z = data[:, 2]
        else:
            raise TypeError, 'xyz must be either a .xyz file or a tuple of (x, y, z) values'
        # Float64 are needed for polygons
        self._x = N.asarray(x, 'float64')
        self._y = N.asarray(y, 'float64')
        self._z = N.asarray(z)
        if trans is not False:
            if operator.isNumberType(trans):
                self._z *= trans
            elif callable(trans):
                self._z[:] = trans(self._z)
        self._m = None
        self.units = units
        self.long_name = long_name
        self.set_transp(transp)
        self.set_magnet(magnet)
        if rsamp is False: rsamp = 0
        self._rsamp = rsamp
        for att, val in kwargs.items(): setattr(self, att, val)
        self._mask = None
        self._xres_man = self._yres_man = 0
        self._proj = self._proj_(True)
        self.id = None

        # Now update arrays
        self._update_(**kwargs)

    def __str__(self):
        s = 'XYZ:'
        for att in 'long_name', 'units':
            if getattr(self, att) is not None:
                s += '\n %s: %s'%(att, getattr(self, att))
        s += '\n total npts: %i'%self._x.shape[0]
        s += '\n full extension: xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g'%\
            (self._x.min(), self._x.max(), self._y.min(), self._y.max(), self._z.min(), self._z.max())
        s += '\n %i selection polygon(s)'%len(self._selections)
        s += '\n %i exclusion polygon(s)'%len(self._exclusions)
        if len(self._selections) or len(self._exclusions):
            s += '\n filtered extension: xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g'%\
                (self._xmin, self._xmax, self._ymin, self._ymax, self._zmin, self._zmax)
        return s

    def copy(self):
        """Deep copy"""
        self._m = None
        return copy.deepcopy(self)

    def __iadd__(self, other):
        if isinstance(other, (float, int)):
            self._z += other
        else:
            # Consolidate to remove exceptions and selections
            self.consolidate()
            # Load
            other = self._load_(other)
            # Transparency
            if not other.get_transp():
                self.exclude(other)
                self.consolidate()
            # Concatenate arrays
            for l in 'xyz':
                setattr(self, '_%s'%l, N.concatenate((getattr(self, l), getattr(other, l))))
            # Now update everything
            self._update_()
        return self
    def __isub__(self, other):
        assert operator.isNumberType(other), 'Only scalars are allowed for such arithmetic operation'
        self._z -= other
        return self
    def __imul__(self, other):
        assert operator.isNumberType(other), 'Only scalars are allowed for such arithmetic operation'
        self._z *= other
        return self
    def __idiv__(self, other):
        assert operator.isNumberType(other), 'Only scalars are allowed for such arithmetic operation'
        self._z /= other
        return self

    def __add__(self, other):
        newxyz = self.copy()
#       if not isinstance(other, (float, int)):
#           newxyz = self.copy().consolidate()
#           if hasattr(other, '_transp') and not other._transp:
#               # Handle transparency
#               newxyz.exclude(other)
#               newxyz.consolidate()
#           else: # Simply load with no opacity
#               other = self._load_(other)
#       else:
#           newxyz = self.copy()
        newxyz += other
        return newxyz
    def __sub__(self, other):
        newxyz = self._class__(self)
        newxyz -= other
        return newxyz
    def __mul__(self, other):
        newxyz = self._class__(self)
        newxyz *= other
        return newxyz
    def __div__(self, other):
        newxyz = self._class__(self)
        newxyz /= other
        return newxyz

    def _load_(self, d):
        # Load for adding
        # - already a XYZ
        if isinstance(d, self.__class__): return d
        # - XYZ via .xyz() (like shapes/shorelines or mergers)
        if hasattr(d, 'xyz'):
            return d.xyz
        # - load it (from file or array)
        return self.__class__(d)

    def _update_(self,**kwargs):

        # Mask (1==good)
        del self._mask
        npt = len(self._x)
        self._mask = N.ones(npt, '?')
        ii = N.arange(npt, dtype='i')

        # Radius underspampling
        if self._rsamp is False: self._rsamp = 0
        if not hasattr(self, '_rsamp_old'):
            self._rsamp_old = self._rsamp
            self._rsamp_mask = None
        if self._rsamp==0 and self._rsamp_mask is not None:
            del self._rsamp_mask
        elif (self._rsamp!=0 and self._rsamp_mask is None) or self._rsamp != self._rsamp_old:
            del self._rsamp_mask
            kwrsamp=kwfilter(kwargs, 'rsamp')
            self._rsamp_mask = rsamp(self._x, self._y, abs(self._rsamp), proj=self._rsamp<0, getmask=True,**kwrsamp)
        if self._rsamp_mask is not None:
            self._mask &= self._rsamp_mask
        self._rsamp = self._rsamp_old

        # Selections
        if len(self._selections):
            smask = N.zeros(npt, '?')
            rectmask = N.zeros(npt, '?')
            for pl in self._selections:

                # Mask is good only within N/S/W/E limits of the polygon
                rectmask[:] = \
                    (self._x>=pl.boundary[:, 0].min()) & (self._x<=pl.boundary[:, 0].max()) & \
                    (self._y>=pl.boundary[:, 1].min()) & (self._y<=pl.boundary[:, 1].max())

                # Mask is good only within the polygon
                for i in ii[rectmask]:
                    smask[i] |= Point((self._x[i], self._y[i])).within(pl)
            del rectmask
            self._mask &= smask
            del smask

        # Exclusions
        for exc in self._exclusions:

            # Check if inclusions and magnet zone
            if isinstance(exc, tuple):
                exc, incs, magnet = exc
            else:
                incs = []
                magnet = None

            # We check only good points within the N/S/W/E limits of the polygon
            good = self._mask & \
                ((self._x>=exc.boundary[:, 0].min()) & (self._x<=exc.boundary[:, 0].max()) & \
                 (self._y>=exc.boundary[:, 1].min()) & (self._y<=exc.boundary[:, 1].max()))

            # Mask is bad within the polygon
            for i in ii[good]:
                out = not Point((self._x[i], self._y[i])).within(exc)

                # Check inclusions and magnet zone
                if not out:

                    # Inclusions = exclusions of exclusions
                    for inc in incs:
                        out = Point((self._x[i], self._y[i])).within(inc)
                        if out: break

                    # Magnet zone
                    if magnet != None:
                        radius, xmgt, ymgt = magnet
                        dst = N.sqrt((xmgt-self._x[i])**2+(ymgt-self._y[j])**2)
                        out = dst.min() > radius
                        del dst

                # Apply
                self._mask[i] = out
            del good
        del ii

        # Limits
        x = self.x
        y = self.y
        z = self.z
        if len(x):
            self._xmin = x.min()
            self._xmax = x.max()
            self._ymin = y.min()
            self._ymax = y.max()
            self._zmin = z.min()
            self._zmax = z.max()
            self._xmean = x.mean()
            self._ymean = y.mean()
        else:
            warn('No more points after exclusion')
            self._xmin = self._xmax = self._ymin = self._ymax = self._zmin = self._zmax = self._xmean = self._ymean = None

        # Resolution
        if len(x):
            self._xres_mauto, self._yres_mauto = self.resol(deg=False)
            self._xres_dauto, self._yres_dauto = self.resol(deg=True)
        else:
            self._xres_mauto = self._yres_mauto = self._xres_dauto = self._yres_dauto = 0
#       del self._m
#       self._m = None

    def consolidate(self):
        """Apply radius undersampling and all exclusions and selections to data and reset them"""
        x = self.x
        y = self.y
        z = self.z
        self._x = x
        self._y = y
        self._z = z
        self.reset_rsamp()
        self.reset_selections()
        self.reset_exclusions()
        return self

    def mask(self):
        """Get the current mask due to exclusion and selection polygons

        .. seealso::

            :meth:`exclude` :meth:`select`
        """
        return ~self._mask

    def _clip_mask_(self, zone, inverse):
        zone = polygons([zone])[0]
        good = self._mask.copy()
        good = good & \
              (self._x>=zone.boundary[:, 0].min()) & (self._x<=zone.boundary[:, 0].max()) & \
              (self._y>=zone.boundary[:, 1].min()) & (self._y<=zone.boundary[:, 1].max())
        ii = N.arange(len(good), dtype='i')
        for i in ii.compress(good):
            good[i] = Point((self._x[i], self._y[i])).within(zone)
        del ii
        if inverse:
            good = ~good
        return good

    def clip(self, zone=None, margin=None, inverse=False, mask=False, id=None, **kwargs):
        """Geographical selection of part of the data

        - *zone*: (xmin,ymin,xmax,ymax) or a float/int a complex polygon (see :func:`~vacumm.misc.grid.masking.polygons`).
        - *margin*: Margin around ``zone`` relative to the resolution (see :meth:`resol`)
        - *inverse*: Inverse the selection.
        - *mask*: ``zone`` must be interpreted as a mask
        """
        if zone is None or zone == (None, )*4:
            return self
        if margin is not None and isinstance(zone, tuple) and len(zone)==4:
            zone = list(zone)
            if zone[0] is None: zone[0] = self.xmin
            if zone[2] is None: zone[2] = self.xmax
            if zone[1] is None: zone[1] = self.ymin
            if zone[3] is None: zone[3] = self.ymax
            xres, yres = self.get_res()
            xmin = zone[0]-margin*xres
            xmax = zone[2]+margin*xres
            ymin = zone[1]-margin*yres
            ymax = zone[3]+margin*yres
            zone = (xmin, ymin, xmax, ymax)
        if mask is False:
            mask = self._clip_mask_(zone, inverse)
        x = self._x[mask]
        y = self._y[mask]
        z = self._z[mask]
        kwargs.setdefault('units', self.units)
        kwargs.setdefault('long_name', self.long_name)
        kwargs.setdefault('id', id)
        xyz = self.__class__((x, y, z), **kwargs)
#       xyz.include(*self._inclusions)
        xyz.exclude(*self._exclusions)
        return xyz

    def zone(self, poly=False, mask=True):
        """Get xmin,ymin,xmax,ymax

        - *poly*: if True, return zone as a Polygon instance"""
        zone = self.get_xmin(mask), self.get_ymin(mask), self.get_xmax(mask), self.get_ymax(mask)
        if poly: return polygons([zone])[0]
        return zone

#   def get_map(self):
#       """Return the map instance or None"""
#       return self._m


    def set_transp(self, transp):
        """Set :attr:`transp`

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        self._transp = transp
    def get_transp(self):
        """Get :attr:`transp`

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        return self._transp
    transp = property(get_transp, set_transp, doc="Transparency boolean attribute")

    def set_magnet(self, magnet):
        """Set the magnet integer attribute.
        If set to ``0``, no magnet effect.

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        if magnet is False: magnet = 0
        self._magnet = magnet
    def get_magnet(self):
        """Get the magnet integer attribute

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        return self._magnet
    magnet = property(get_magnet, set_magnet, doc="Magnet integer attribute")

    def set_rsamp(self, rsamp):
        """Set the radius sampling :attr:`rsamp`
        If set to ``0``, no sampling."""
        if rsamp is False: rsamp = 0
        self._rsamp = rsamp
        if not hasattr(self, '_rsamp_old') or rsamp != self._rsamp_old:
            self._update_()
    def get_rsamp(self):
        """Get the radius sampling :attr:`rsamp`"""
        return self._rsamp
    def reset_rsamp(self):
        """Reset :attr:`rsamp` without affecting data """
        self._rsamp = self._rsamp_old = 0
        if hasattr(self, '_rsamp_mask'): del self._rsamp_mask
        self._rsamp_mask = None
    del_rsamp = reset_rsamp
    rsmap = property(get_rsamp, set_rsamp, del_rsamp, "Radius of unsersampling")

    def tocfg(self, cfg, section, param=None):
        """Dump one or all parameters as options to a cfg section

        - **cfg**: ConfigParser object
        - **section**: Section of cfg
        - *param*: A single or a list of parameter names
        """
        # List of params
        allowed = ['xmin', 'xmax', 'ymin', 'ymax','zmin','zmax', 'long_name', 'units', 'transp',
                   'xres','yres','exclusions', 'selections']
        if param is None:
            param = allowed
        elif isinstance(param, str):
            param = [param]
        # Get string
        for param in param:
          if param.endswith('res'):
             val = getattr(self, '_'+param+'_mauto')
          elif param.endswith('sions'):
             # selections, exclusions
             val = []
             for var in getattr(self, '_'+param):
                 if isinstance(var, tuple):
                    # exclusions and inclusions
                    exc, incs = var
                    if len(incs) == 0:
                       # no inclusions
                       val.append(exc.get_coords().tolist())
                    else:
                       polys = []
                       for inc in incs:
                          polys.append(inc.get_coords().tolist())
                          val.append((exc, polys))
                 else:
                     # normal case
                     val.append(var.get_coords().tolist())
          elif param in ['units', 'long_name']:
              val = getattr(self, param)
          else:
              val = getattr(self, '_'+param)
          # Check section
          if not cfg.has_section(section):
              cfg.add_section(section)
          # Dump
          cfg.set(section, param, str(val))


    def get_xmin(self, mask=True):
        if mask is True or mask == 'masked': return self._xmin
        return self.get_x(mask).min()
    xmin = property(get_xmin, doc="X min")
    def get_xmax(self, mask=True):
        if mask is True or mask == 'masked': return self._xmax
        return self.get_x(mask).max()
    xmax = property(get_xmax, doc="X max")
    def get_ymin(self, mask=True):
        if mask is True or mask == 'masked': return self._ymin
        return self.get_y(mask).min()
    ymin = property(get_ymin, doc="Y min")
    def get_ymax(self, mask=True):
        if mask is True or mask == 'masked': return self._ymax
        return self.get_y(mask).max()
    ymax = property(get_ymax, doc="Y max")
    def get_zmin(self, mask=True):
        if mask is True or mask == 'masked': return self._zmin
        return self.get_z(mask).min()
    zmin = property(get_zmin, doc="Z min")
    def get_zmax(self, mask=True):
        if mask is True or mask == 'masked': return self._zmax
        return self.get_z(mask).max()
    zmax = property(get_zmax, doc="Z max")

    def exclude(self, *zones):
        """Add one or more zones where data are not used.

        A zone can be :

        - an argument to :func:`~vacumm.misc.grid.masking.polygons` to get a :class:`_geoslib.Polygon` instance,
        - another :class:XYZ` instance from which the convex hull (see :meth:`hull`) is used as a delimiting area

        :Usage:

        >>> xyz.exclude([[-8,43],[-5.5,43],[-6,45.]],[[-10,45],[-7,47],[-10,49.]])
        >>> xyz.exclude(polygon1,polygon2)
        >>> xyz.exclude(xyz1,[-5,42,-3,48.])

        .. seealso::

            :meth:`select` :meth:`exclusions`
        """
        zones = list(zones)
        for i, zone in enumerate(zones):
            if isinstance(zone, XYZ):
                zone = zone.shadows()
            if isinstance(zone, tuple): # inclusions
                exc = polygons([zone[0]])[0]
                if len(zone) == 1 or len(zone[1])==0:
                    zone = exc
                else:
                    zone = exc, polygons(zone[1]), zone[2]
            else:
                zone = polygons([zone])[0]
            zones[i] = zone
        self._exclusions.extend(zones)
        self._update_()
    def select(self, *zones):
        """Add one or more zone (polygons) where only these data are used

        A zone is an argument to :func:`~vacumm.misc.grid.masking.polygons` to get a :class:`_geoslib.Polygon` instance.

        :Usage:

        >>> xyz.select([[-8,43],[-5.5,43],[-6,45.]],[[-10,45],[-7,47],[-10,49.]])
        >>> xyz.select(polygon1,polygon2)

        .. seealso::

            :meth:`exclude`  :meth:`selections`
        """
        self._selections.extend(polygons(list(zones)))
        self._update_()
    def reset_selections(self):
        """Remove all selections"""
        del self._selections
        self._selections = []
        self._update_()
    def reset_exclusions(self):
        """Remove all exclusions"""
        del self._exclusions
        self._exclusions = []
        self._update_()
    def selections(self):
        """Get all selection polygons as a tuple"""
        return self._selections
    def exclusions(self):
        """Get all exclusion polygons as a tuple"""
        return self._exclusions


    def _filter_(self, var, mask):
        """
        - var can be a list
        - mask can be
            - 'valid', True
            - False, None
            - 'revert', 'masked'
            - mask array
        """
        if mask == 'valid': mask = True
        assert self._mask is not None
        if mask is False :
            return var
        if isinstance(mask, str) and (mask.startswith('rever') or mask.startswith('inver') or
            mask.startswith('mask')):
            mask = ~self._mask
        else:
            mask = self._mask
        if isinstance(var, list):
            return [v.compress(mask) for v in var]
#       if isinstance(var, list):
#           return [var[i] for i, m in enumerate(self._mask) if m]
        return var.compress(mask)
    def get_x(self, mask=True):
        """Get valid X positions"""
        return self._filter_(self._x, mask)
    x = property(get_x, doc='Valid X positions')
    def get_y(self, mask=True):
        """Get valid Y positions"""
        return self._filter_(self._y, mask)
    y = property(get_y, doc='Valid Y positions')
    def get_z(self, mask=True):
        """Get valid Z values"""
        return self._filter_(self._z, mask)
    z = property(get_z, doc='Valid Z values')
    def get_xy(self, mask=True):
        """Return coordinates as a (2, npts) array :attr:`xy`

        - ``xy()[0]``: X
        - ``xy()[1]``: Y
        """
        return N.asarray([self.get_x(mask), self.get_y(mask)])
    xy = property(get_xy, doc='Coordinates as a (2, npts) array')
    def get_xyz(self, mask=True, split=False):
        """Return coordinates and data as a (3, npts) array :attr:`xyz`

        - ``xy()[0]``: X
        - ``xy()[1]``: Y
        - ``xy()[2]``: Z
        """
        if split: return self.get_x(mask), self.get_y(mask), self.get_z(mask)
        return N.asarray([self.get_x(mask), self.get_y(mask), self.get_z(mask)])
    xyz = property(get_xyz, doc='Coordinates and data as a (3, npts) array')

    def astuple(self, mask=True):
        """Shortcut to ``xyz(split=True)`` (see :meth:`xyz`)"""
        return self.get_xyz(mask=mask, split=True)


    def __len__(self):
        return self._mask.sum()
    def __getitem__(self, key):
        x = self.x[key]
        y = self.y[key]
        z = self.z[key]
        if isinstance(key, int):
            return x, y, z
        return zip(x, y, z)

    def interp(self, xyo, xyz=False, **kwargs):
        """Interpolate to (xo,yo) positions using :class:`nat.Natgrid`

        :Params:

            - **xo**: Output X
            - **yo**: Output Y
            - *xyz*: If True, return a :class:`XYZ` instance instead of a :mod:`numpy` array
            - **interp_<param>**, optional: ``<param>`` is passed to the
              :func:`~vacumm.misc.grid.regridding.xy2xy` interpolation routine.
            - Other params are passed to XYZ initialization for the output dataset.

        :Returns: An XYZ instance
        """
        # FIXME: Natgrid still required by this module ??
        # from nat import Natgrid
        if isinstance(xyo, (tuple, N.ndarray)):
            xo, yo = xyo
        elif hasattr(xyo, 'x'):
            xo = xyo.x
            yo = xyo.y
        elif hasattr(xyo, 'xy'):
            xo, yo = xyo.xy
        elif hasattr(xyo, 'xyz'):
            xo, yo, tmp = xyo.xyz
        else:
            raise TypeError, 'Wrong input type'
        kwinterp = kwfilter(kwargs, 'interp_')
        zo = vcgr.xy2xy(self.x, self.y, self.z, xo, yo, **kwinterp)
        if not xyz: return zo
        kwargs.setdefault('units', self.units)
        kwargs.setdefault('long_name', self.long_name)
        return self.__class__((xo, yo, zo), **kwargs)

    def hull(self, out='xy', mask=True):
        """Return the convex hull

        :Returns: Depends on ``out``

        - ``"xy"``: (xhull, yhull)
        - ``"ind"``: indices of points
        - ``"poly"``: :class:`_geoslib.Polygon` instance
        """
        return vacumm.misc.grid.masking.convex_hull((self.get_x(mask), self.get_y(mask)), poly=out=='poly')

    def shadows(self):
        """Get the polygons defining the 'shadow' of this dataset.

        It consists of a tuple of two elements:

            - the convex hull as a polygon,
            - a list of exclusion polygons that intersect the convex hull.

        Therefore, a point in the shadow must be inside the convex hull polygon,
        and outside the exclusion polygons.

        :Returns: (hull_poly, [exclusion_poly1,...])
        """
        hull_poly = self.hull(out='poly')
        exc_polys = []
        for exc in self.exclusions():
            if isinstance(exc, tuple): exc = exc[0]
            if hull_poly.intersects(exc):
                exc_polys.append(exc)
            #FIXME: inclusions
        magnet = None
        if self.get_magnet():
            if self.get_magnet() < 0:
                radius = -self.get_magnet()
            else:
                xres, yres = self.get_res(deg=True)
                radius = N.sqrt((xres**2+yres**2)/2.)*self.get_magnet()
            magnet = radius, self.x, self.y
        return hull_poly, exc_polys, magnet


    def contains(self, x, y):
        """Check if one or several points are within a the convex hull

        - **x,y**: X,Y positions as floats or lists or an :mod:`numpy` arrays.

        .. seealso:

            :meth:`hull`
        """
        hull = self.hull(out='poly')
        try:
            len(x)
            x = N.asarray(x)
            y = N.asarray(y)
            good = (x>=self.xmin) & (x<=self.xmax) & \
                  (y>=self.ymin) & (y<=self.ymax)
            xg = x.compress(good)
            yg = y.compress(good)
            out = []
            for xp, yp in zip(xg, yg):
                out.append(Point(xp, yp).within(hull))
            if isinstance(x, list): return out
            return N.array(out)
        except:
            if (x<self.xmin)|(x>self.xmax)|(y<self.ymin)|(y>self.ymax):
                return False
            return Point(x, y).within(hull)


    def resol(self, convex_hull_method='delaunay', exc=[], deg=False):
        """Return the mean resolution.

        Algorithm: Median distances between facets of triangles

        :Returns: (xres,yres)

        """
        # Coordinates
        x = self.x
        y = self.y
        if not deg:
            x, y = self._proj_(True)(x, y)

#       if method.startswith('tri'): # Using the median length of triangles
        from scipy.spatial import Delaunay

        coord=N.vstack((x,y)).T
        t=Delaunay(coord)
        distances = []
        for p1,p2,p3 in t.vertices:
           distances.append(N.sqrt((x[p1]-x[p2])**2+(y[p1]-y[p2])**2))
           distances.append(N.sqrt((x[p1]-x[p3])**2+(y[p1]-y[p3])**2))
           distances.append(N.sqrt((x[p2]-x[p3])**2+(y[p2]-y[p3])**2))

        xres = yres = N.median(distances)
        del t, distances

#       else: # Using the convex hull
#
#           # Convex hull polygon
#           hull = convex_hull((x, y), poly=True, method=convex_hull_method)
#
#           # Area
#           # - max = convex hull
#           area = hull.area()
#           # - substract intersections with exclusions
#           #FIXME: xyz.resol: non overlapping exclusions+inclusions
#           for e in exc+self.exclusions():
#               if isinstance(e, tuple):
#                   e, incs = e
#               else:
#                   incs = []
#               if hull.intersects(e):
#                   for i in hull.intersection(e):
#                       area -= i.area()
#
#           # Area and mean resolution
#           xres = yres = N.sqrt(area/len(x))

        return xres, yres

    def set_res(self, xres, yres=None):
        """Set the resolution of the dataset

        If ``yres`` is not, it is set to ``xres``.
        When a value is **negative**, it is supposed to be in **meters** (not in degrees)"""
        if yres is None: yres = xres
        self._xres_man = xres
        self._yres_man = yres

    def get_res(self, deg=False, auto=None):
        """Get the mean X and Y resolutions in meters or degrees"""
        # Get manual or auto resolutions
        if self._xres_man and auto is not True:
            xres = self._xres_man
        elif deg:
            xres = self._xres_dauto
        else:
            xres = self._xres_mauto
        if self._yres_man and auto is not True:
            yres = self._yres_man
        elif deg:
            yres = self._yres_dauto
        else:
            yres = self._yres_mauto
        # Conversions
        return self._get_res_(xres, yres, deg)

    def _get_res_(self, xres, yres, degres):
        if degres: # need degrees
            if xres<0: xres = -vcpu.m2deg(xres, self._ymean)
            if yres<0: yres = -vcpu.m2deg(yres)
            return xres, yres
        else: # need meters
            if xres>0: xres = vcpu.deg2m(xres, self._ymean)
            if yres>0: yres = vcpu.deg2m(yres)
            return xres, yres


    def _proj_(self, proj):
        """Get a proper projection or None"""
        if proj is None: proj = False
        if not proj: return
        if not callable(proj):
            proj = get_proj((self._x, self._y))
        return proj

    def get_grid(self, res=None, xmin=None, xmax=None, ymin=None, ymax=None, relres=.5, degres=False, id='xyz_grid'):
        """Generate a rectangular grid based on x/y positions and resolution

        - *res*: Resolution. It can be:

                - a float where then ``xres=yres=res``
                - a tuple as ``(xres,yres)``
                - else it is guessed using :meth:`get_res` (and maybe :meth:`resol`)` and multiplied by ``relres``

        - *relres*: Relative resolution factor applied to ``res`` when resolution is guessed (``res=None``)
        - *degres*: When ``res`` is explicitly given, it interpreted as degrees is ``degres`` is True.
        - *xmin,xmax,ymin,ymax*: Bounds of the grid. If not specified, bounds of the dataset are used (see :meth:`xmin`, etc).

        .. note::

            Resolutions are adjusted when they are not mutiple of grid extensions (slightly decreased).
            Therefore, extensions of the grid are always preserved.

        .. seealso::

            :meth:`resol`, :meth:`togrid`
        """
        # Resolution
#       proj = self._proj_(proj)
        if res is None: # Auto
            # Meters
#           dx, dy = tuple([r*relres for r in self.get_res(deg=False, auto=True)])
            dx, dy = tuple([r*relres for r in self.get_res(deg=True, auto=True)])
#           # To degrees
#           dx, dy = self._get_res_(-dx, -dy, True)
        else: # Manual
            if isinstance(res, tuple):
                xres, yres = res
            else: # Same for X and Y
                xres = yres = res
            if not degres: # Given in meters
                xres = -xres
                yres = -yres
            # Now in degrees
            dx, dy = self._get_res_(xres, yres, True)

        # Bounds
        if xmin is None: xmin = self.xmin
        if xmax is None: xmax = self.xmax
        if ymin is None: ymin = self.ymin
        if ymax is None: ymax = self.ymax

        # Rectifications
        xratio = (xmax-xmin)/dx
        yratio = (ymax-ymin)/dy
        dx += dx*(xratio%1)/int(xratio)
        dy += dy*(yratio%1)/int(yratio)

        # Grid
        grid = vacumm.misc.grid.create_grid(N.arange(xmin, xmax+dx/2., dx), N.arange(ymin, ymax+dy/2., dy))
        grid.id = id
        return grid
    grid = property(get_grid, doc="Rectangular grid based on x/y positions and resolution")

    def togrid(self, grid=None, mask=False, cgrid=False,  **kwargs):
        """Interpolate to a regular grid

        - **grid**: The output grid. It can be either:

            - a (x,y) tuple or a grid or a :mod:`MV2` variable with a grid,
            - ``None``, thus guessed using :meth:`grid`

        - *mask*: It can be either:

            - ``None``, ``False`` or ``MV2.nomask``: no masking
            - an array: this mask array is directly applied
            - a :class:`Shapes` instance (or :class:`~vacumm.bathy.shorelines.ShoreLine`)
              or a single char GSHHS resolution (and optionally 's' for Histolitt)
            - a callable fonction so that ``mask = thisfunc(mask, **kwmask)``
            - a float: data with this value are masked

        - *mask_<param>*: <param> is passed to :func:`~vacumm.misc.grid.masking.polygon_mask`
          for evaluation of mask thanks to the polygons.
        - *grid_<param>*: <param> is passed to :func:`grid`.
        - *cgrid*: If ``True``, returns bathy at U- and V-points, else at T-points
        - Other keyparam are passed to :func:`~vacumm.misc.grid.regridding.griddata`
          for regridding.

        Return: ``(Zx,Zy)`` OR ``Z`` depending on cgrid.

        .. seealso:

            :func:`~vacumm.misc.grid.masking.polygon_mask`
            :class:`~vacumm.misc.grid.basemap.GSHHS_BM`
            :class:`Shapes`  :class:`~vacumm.bathy.shorelines.ShoreLine`
        """
        # Grid
        kwgrid = kwfilter(kwargs, 'grid_')
        if grid is None:
            grid = self.get_grid(**kwgrid)

        # Interpolation
        kwmask = kwfilter(kwargs, 'mask_')
        kwargs.setdefault('method', 'carg')
        kwargs.setdefault('ext', True)
        var = vcgr.griddata(self.x, self.y, self.z, grid, mask=None, cgrid=cgrid, **kwargs)
        if not cgrid: var = var,

        # Mask from polygons
        if isinstance(mask, (Shapes, str)):

            # Clip for polygons
            clip = kwmask.pop('clip', .1)
            if isinstance(clip, float):
                xx, yy = get_xy(grid)
                xmin = xx[:].min() ; xmax = xx[:].max()
                ymin = yy[:].min() ; ymax = yy[:].max()
                dx = (xmax-xmin)*clip
                dy = (ymax-ymin)*clip
                clip = [xmin-dx, ymin-dy, xmax+dx, ymax+dy]

            # Get polygons
            if isinstance(mask, Shapes): # Direct
                shapes = mask.clip(clip)
            else: # GSHHS
                shapes = GSHHS_BM(mask, clip=clip)

            # Create mask
            mask = vacumm.misc.grid.masking.polygon_mask(grid, shapes.get_shapes(), **kwmask)

        # Masking
        if mask is not None and mask is not False and mask is not MV2.nomask:
            for vv in var:
                if hasattr(mask, 'ndim'): # array
                    vv[:] = MV.masked_where(mask, vv, copy=0)
                elif callable(mask): # function
                    vv[:] = mask(vv, **kwmask)
                else: # value
                    vv[:] = MV.masked_values(vv, mask, copy=0)

        # Attributes
        for vv in var:
            if self.units is not None:
                vv.units = self.units
            if self.long_name is not None:
                vv.long_name = self.long_name
        if cgrid: return var
        return var[0]

    def toxy(self, xo, yo, mask=None, outtype='tuple'):
        """Interpolate on random points using :func:`~vacumm.misc.grid.regridding.xy2xy`

        - **xo,yo**: Output positions
        - *mask*: It can be either:

            - ``None``, ``False`` or ``MV2.nomask``: no masking
            - a :class:`Shapes` instance (or :class:`~vacumm.bathy.shorelines.ShoreLine`) or a single char GSHHS resolution (and optionally 's' for Histolitt)
        - *outtype*: Define output type

            - ``"tuple"``: as a tuple (x, y, z)
            - ``"xyz"``: as xyz block
            - ``"XYZ"``: as an :class:`XYZ` (or subclass) instance
        """

        # Interpolate
        zo = vcgr.xy2xy(self.x, self.y, self.z, xo, yo, **kwargs)


        # Mask from polygons
        if isinstance(mask, (Shapes, str)):

            # Clip for polygons
            clip = kwmask.pop('clip', .1)
            if isinstance(clip, float):
                dx = (xo.max()-xo.min())*clip
                dy = (yo.max()-yo.min())*clip
                clip = [xo.min()-dx, yo.min()-dy, xo.max()+dx, yo.max()+dy]

            # Get polygons
            if isinstance(mask, Shapes): # Direct
                shapes = mask.clip(clip)
            else: # GSHHS
                shapes = GSHHS_BM(mask, clip=clip)

            # Maskit
            xo, yo, zo = polygon_select(xo, yo, shapes.get_shapes(), zz=zo)

        # Out
        if outype=='xyz':
            return N.asarray([xo, yo, zo])
        if outtype=='XYZ':
            return self.__class__((xo, yo, zo))
        if callable(outtype):
            return outtype([xo, yo, zo])
        return xo, yo, zo



    def plot(self, size=5., color=None, alpha=1., masked_alpha=.3,
        masked_size=None, linewidth=0., show=True, savefig=None,
        savefigs=None, m=None, colorbar=True, title=None,
        units=None,  cmap=None, mode='valid', zorder=100,
        masked_zorder=50, margin=2,
        xmin=None, xmax=None, ymin=None, ymax=None,
        xres=None, yres=None, **kwargs):
        """Scatter plot of bathymetry points

        :Params:

            - **mode**, optional: 'valid', 'masked' or 'both'.
            - **size**, optional: Size of markers.
            - **color**, optional: Color of markers.
            - **alpha**, optional: Alpha transparency of markers.
            - **zorder**, optional: zorder of markers.
            - **m**, optional: Use this :class:`~mpl_toolkits.basemap.Basemap` instance
              to plot the points.
            - **masked_size**, optional: Size of masked markers.
            - **masked_alpha**, optional: Alpha transparency of masked markers.
            - **masked_zorder**, optional: zorder of masked markers.
        """

        # Params
        kwm = kwfilter(kwargs, 'map')
        kwhull = kwfilter(kwargs, 'hull')
        kwcb = kwfilter(kwargs, 'colorbar')
        kwplot = dict(linewidth=linewidth, cmap=get_cmap(cmap))
        kwplot.update(kwargs)
        pts = None

        # Limits
        if m is True:
            if xmin is None: xmin = self.xmin
            if xmax is None: xmax = self.xmax
            if ymin is None: ymin = self.ymin
            if ymax is None: ymax = self.ymax
        if margin != 0:
            if xres is None or yres is None:
                xxres, yyres = self.resol(deg=True)
            if xres is None: xres = xxres
            if yres is None: yres = yyres
            if xmin is not None: xmin -= margin*xres
            if xmax is not None: xmax += margin*xres
            if ymin is not None: ymin -= margin*yres
            if ymax is not None: ymax += margin*yres

        # Map
        if m is True:
            kwm.setdefault('lon', (xmin, xmax))
            kwm.setdefault('lat', (ymin, ymax))
            kwm['projection'] = 'merc'
            m = map2(show=False, savefig=None, **kwm)
        if m:
            G = m.map if hasattr(m, 'map') else m
            self._m = G
        else:
            G = P

        # Max values
        if mode == 'both':
            zmin = self._z.min()
            if not kwplot.has_key('vmin'): kwplot['vmin'] = zmin
            zmax = self._z.max()
            if not kwplot.has_key('vmax'): kwplot['vmax'] = zmax

        # Valid points
        if mode != 'masked' or mode == 'both':
            # Points
            pts = self._plot_('valid', G, m, size, color, alpha=alpha, zorder=zorder,
                label=self.long_name, **kwplot)
#           # Convex hull
#           if hull:
#               xhull, yhull = self.hull(out='xy')
#               if callable(m):
#                   xhull, yhull = m(xhull, yhull)
#               xhull = xhull.tolist()+[xhull[0]]
#               yhull = yhull.tolist()+[yhull[0]]
#               kwhull.setdefault('alpha', masked_alpha)
#               kwhull.setdefault('zorder', zorder)
#               if color is None:
#                   hullcolor = '#888888'
#               else:
#                   hullcolor = color
#               kwhull.setdefault('color', hullcolor)
#               G.plot(xhull, yhull, '-', **kwhull)
        # Masked points
        if mode == 'masked' or mode == 'both':
            if mode == 'masked':
                masked_alpha = alpha
                masked_size = size
                masked_zorder = zorder
            elif masked_size is None: masked_size = size
            p = self._plot_('masked', G, m, masked_size, color, alpha=masked_alpha, zorder=masked_zorder, **kwplot)
            if pts is None: pts = p

        # Limits
        if m:
            G.set_axes_limits(P.gca())
        else:
            if xmin is not None: P.xlim(xmin=xmin)
            if xmax is not None: P.xlim(xmax=xmax)
            if ymin is not None: P.ylim(ymin=ymin)
            if ymax is not None: P.ylim(ymax=ymax)

        # Decorations
        if title is None and self.long_name is not None:
            title = self.long_name
        if title is not None:
            if self.long_name is None:
                self.long_name = title
            P.title(title)
#           for pt in pts:
#               if pt is not None:
#                   pt.set_label(title)
#                   break
        if units is None and self.units is not None:
            units = self.units
        if units is not None:
            if self.units is None:
                self.units = units
        if colorbar:
            _colorbar_(pts, units=units, **kwcb)
        if savefig:
            P.savefig(savefig)
        elif savefigs:
            savefigs(savefigs)
        if show:
            P.show()
        return pts

    def _plot_(self, mask, G, proj, size, color, **kwargs):
        """Generic plot"""
        x, y = self.get_x(mask), self.get_y(mask)
        if not len(x): return None
        if callable(proj): x, y,  = proj(x, y)
        if color is None: color = self.get_z(mask)
        kwargs.setdefault('linewidth', 0)
        if kwargs.get('marker', None) is None:
            kwargs['marker'] = 'o'
        return  G.scatter(x, y, s=size, c=color, **kwargs)

    def save(self, xyzfile, **kwargs):
        """Save to a file

        - **xyzfile**: Output file name

            - write a netcdf file if it ends with ".nc" or ".grd"
            - write a sinux file if it ends with ".snx"
            - else write an ascii file with 3 columns

        - Other keywords are passed to :func:`numpy.savetxt` for ascii saving
        """
        if xyzfile.endswith('.nc') or xyzfile.endswith('.grd'):
            x = cdms2.createVariable(self.x, id='x')
            y = cdms2.createVariable(self.y, id='y')
            z = cdms2.createVariable(self.z, id='z')
            i = cdms2.createAxis(range(len(x)), id='i')
            f = cdms2.open(xyzfile, 'w')
            for var in x, y, z:
                var.setAxis(0, i)
                f.write(var)
            f.close()
        elif xyzfile.endswith('.snx'):
            write_snx(self.xyz.T, type='point', **kwargs)
        else:
            N.savetxt(xyzfile, self.xyz.T, **kwargs)


class XYZMerger(object):
    """Mix different bathymetries"""
    def __init__(self, *datasets, **kwargs):
        self._xmin = kwargs.pop('xmin', None)
        self._xmax = kwargs.pop('xmax', None)
        self._ymin = kwargs.pop('ymin', None)
        self._ymax = kwargs.pop('ymax', None)
        self._datasets = list(datasets)
        self._XYZ = XYZ
        self.long_name = kwargs.get('long_name', None)
        self.units = kwargs.get('units', None)

    def __str__(self):
        s = 'Merger:'
        for att in 'long_name', 'units':
            if getattr(self, att) is not None:
                s += '\n %s: %s'%(att, getattr(self, att))
        n = len(self)
        if not n:
            return s+'\n no datasets in the merger'
        xyz = self.xyz
        s += '\n extension: xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g'%\
            (xyz.xmin, xyz.xmax, xyz.ymin, xyz.ymax, xyz.zmin, xyz.zmax)
        if n==1:
            s += '\n there is 1 dataset in the merger:'
        else:
            s += '\n there are %i datasets in the merger:'%n
        sep = '\n'+'-'*3
        s += sep
        s += sep.join(['\n'+str(d) for d in self])
        s += sep
        return s


#   def __copy__(self):
#       return self.__class__(xmin=self._xmin, xmax=self._xmax, ymin=self._ymin, ymax=self._ymax, *self._datasets)
    def copy(self):
        return copy.deepcopy()
#       return self.__copy__()

    def _load_(self, d):
        # XYZMerger instance
        if isinstance(d, self.__class__):
        # xyzmerger._load_: pas clair
#           if d is self: return
#           return d.copy()
            return [self._XYZ(dd) for dd in d._datasets]
        # XYZ instance
        if isinstance(d, self._XYZ):
            return d
        # Get XYZ from b
        if hasattr(d, 'xyz'):
            return d.xyz
        # New XYZ instance
        return self._XYZ(d)

    def tolist(self):
        """Return the merger as a list of datasets"""
        return self._datasets

    def ids(self):
        return [d.id for d in self._datasets]

    def append(self, d):
        """Append a dataset to the merger"""
        self += d
    def remove(self, d):
        """Remove a dataset from the merger"""
        self -= d
    def clean(self):
        """Remove all current dataset"""
        del self._datasets
        self._datasets = []

    def __iadd__(self, datasets):
        if not isinstance(datasets, list): datasets = [datasets]
        for d in datasets:
            d = self._load_(d)
            if d is not None and d not in self._datasets:
                self._datasets.append(d)
        return self
    def __isub__(self, d):
        if d not in self._datasets:
            warn('Dataset not in merger')
        self._datasets.remove(d)
        return self

    def __getitem__(self, key):
        return self._datasets[self._key_(key)]
    def __delitem__(self, key):
        del self._datasets[self._key_(key)]
    def __setitem__(self, key, b):
        self._datasets[self._key_(key)] = self._load_(b)
    def __len__(self):
        return len(self._datasets)
    def _key_(self, key):
        if isinstance(key, str):
            key = self.ids().find(key)
        return key

    def get_xyz(self, mask=True, **kwargs):
        """Merge current dataset"""
        assert len(self), 'You must add at least one dataset to the merger'
        xyz = self._XYZ(self._datasets[0], **kwargs)
        if mask is False:
            xyz.reset_exclusions()
            xyz.reset_selections()
        if mask:
            pass
        for d in self[1:]:
            if mask is False:
                d = copy.deepcopy(d)
                d.reset_exclusions()
                d.reset_selections()
            xyz += d
            if mask:
                pass
        xyz.long_name = self.long_name
        xyz.units = self.units
        return xyz
    xyz = property(get_xyz, doc='Coordinates and data as a (3, npts) array')

    def merge(self, **kwargs):
        """Shortcut to :meth:`xyz`"""
        return self.get_xyz(**kwargs)

    def togrid(self, *args, **kwargs):
        """Interpolate merged bathymetries to a grid"""
        return self.xyz.togrid(*args, **kwargs)

    def plot(self, color=None, marker=None, mode='cluster', title='XYZ merger',
        show=True, colorbar=True, savefig=None, savefigs=None, legend=True,
        xmin=None, xmax=None, ymin=None, ymax=None, margin=5, xres=None, yres=None, **kwargs):
        """

        - *alpha*: Alpha transparency:

            - applied to **all** points if ``mode="cluster"``
            - applied to **hidden** points if ``mode="data"``

        - *mode*: Display mode:

            - ``"cluster"``: Points from different datasets have different colors and markers,
                and hidden points are transparent.
            - ``"data"``: Points have the same marker, colors depends on Z value and hidden
                points are masked.

        - *marker*: Define a single or several markers to be used.
        - *legend*: Show a legend if ``mode="cluster"``.
        - *title*: Title of the plot.
        - *m*: :class:`~mpl_toolkits.basemap.Basemap` instance.
        - *m_margin*: Margin for ``m``, relative to the mean resolution (see :meth:`XYZ.resol`)
        - *m_<keywords>*: Keywords are passed to :func:`~vacumm.misc.plot.map`.
        - Extra keywords are passed to :meth:`XYZ.plot`.
        """


        # Limits
        if None in [xmin, xmax, ymin, ymax] or margin != 0.:
            xyz = self.get_xyz(mask=False)
        if xmin is None: xmin = xyz.get_xmin(False)
        if xmax is None: xmax = xyz.get_xmax(False)
        if ymin is None: ymin = xyz.get_ymin(False)
        if ymax is None: ymax = xyz.get_ymax(False)
        if margin != 0. and None in [xres, yres]:
            xxres, yyres = xyz.resol(deg=True)
            if xres is None: xres = xxres
            if yres is None: yres = yyres

        # Arguments to plots
        kwleg = kwfilter(kwargs, 'legend')
        kwsf = kwfilter(kwargs, 'savefig')
        kwsfs = kwfilter(kwargs, 'savefigs')
        kwplot = dict(show=False, colorbar=False, title=False,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            margin=margin, xres=xres, yres=yres)
        kwplot.update(kwargs)

        # Select mode
        pp = []
        if mode == 'cluster':

            # Colors and symbols
            if color is None:
                color = simple_colors
            color = broadcast(color, len(self))
            if marker is None:
                marker = Markers
            marker = broadcast(marker, len(self))

            # Copy all datasets
            datasets = [copy.deepcopy(d) for d in self._datasets]

            # Loop on datasets
            m = kwplot.pop('m', None)
            for j, i in enumerate(xrange(len(datasets)-1, -1, -1)):

                # Plot current dataset
                pp.append(datasets[i].plot(mode='both', color=color[i], marker=marker[i],
                    zorder=200+i, masked_zorder=50+i, m=m, **kwplot))
                m = datasets[i]._m

                # Mask next datasets with current one if not transparent
                if i != 0 and not datasets[i].get_transp():
                    for k in xrange(i-1, -1, -1):
                        datasets[k].exclude(datasets[i])

            # Legend
            if legend:
                hh = ()
                ll = ()
                for d, p in zip(datasets[::-1], pp):
                    if d.long_name is not None:
                        hh += p,
                        ll += d.long_name,
                legend_alpha = kwleg.pop('alpha', .3)
                kwleg.setdefault('loc', 'best')
                kwleg.setdefault('shadow', False)
                leg = P.legend(hh, ll, **kwleg)
                leg.legendPatch.set_alpha(legend_alpha)
                leg.set_zorder(1000)

        else:

            # Physical case
            pp.append(self.xyz.plot(color=color, marker=marker, **kwplot))
            if colorbar:
                _colorbar_(pp[0])

        # Plot attributes
        if title is None and self.long_name is not None:
            title = self.long_name
        if title is not False:
            P.title(title)
        if savefig:
            P.savefig(savefig, **kwsf)
        if savefigs:
            Savefigs(savefigs, **kwsfs)
        if show:
            P.show()
        return pp


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
        if self._var is None:
            self.merge()
        return map2(self._var)

