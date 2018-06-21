# -*- coding: utf8 -*-
"""Utilities derived from mpl_toolkits.basemap"""

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
__all__  = ['gshhs_reslist', 'gshhs_autores', 'cached_map', 'cache_map', 'get_map',
'GSHHS_BM', 'merc', 'clean_cache', 'reset_cache', 'get_map_dir', 'get_proj',
'create_map', 'RSHPERE_WGS84', 'GSHHS_RESLIST']
__all__.sort()

import os
import numpy as N
from mpl_toolkits.basemap import Basemap, __version__ as basemap_version
from mpl_toolkits.basemap.proj import Proj
from matplotlib import get_configdir
import cPickle, stat
#FIXME:imports
from ...misc.io import Shapes
from misc import get_xy
from ...misc.phys.constants import R as rsphere_mean
from ...misc.phys.units import deg2m, m2deg
from ...misc.misc import kwfilter, dict_check_defaults,  zoombox
from vacumm.config import get_config_value
from ...__init__ import vacumm_warn

#: Earth radius of wgs84 ellipsoid
RSHPERE_WGS84 = (6378137.0,6356752.3141)
rshpere_wgs84 = RSHPERE_WGS84

#: GSHHS shorelines letters
GSHHS_RESLIST = ['f', 'h', 'i', 'l', 'c']
gshhs_reslist = GSHHS_RESLIST

def gshhs_autores(lon_min, lon_max, lat_min, lat_max, asindex=False, shift=None):
    """Guess best resolution from lon/lat bounds"""
    testresol=((lon_max-lon_min)+(lat_max-lat_min))/2.0
    ires = N.array([-1.,1. ,5.,15.,50.]).searchsorted(testresol)-1
    if isinstance(shift, int):
        ires += shift
        ires = N.clip(ires, 0, len(GSHHS_RESLIST)-1)
    if asindex:
        return ires
    return GSHHS_RESLIST[ires]

# Cached maps
def cached_map(m=None, mapdir=None, verbose=False, **kwargs):
    """Check if we have a cached map

    - *m*: A Basemap instance [Default: None]
    - *mapdir*: Where are stored the cached maps. If ``None``,
      :func:`matplotlib.get_configdir` is used as a parent directory,
      which is the matplotlib configuration directory
      (:file:`~/.matplotlib` undex linux), and
      :file:`basemap/cached_maps` as the subdirectory.

    :Example:

    >>> m = cached_map(lon_min=-5, lon_max=6, lat_min=40, lat_max=50, projection='lcc', resolution='f')
    >>> m = cached_map(m) # Does only caching of map
    """

    # We already have one in live memory
    if isinstance(m, Basemap):
        # Save it
        cache_map(m, mapdir=mapdir)
        # Get it
        return m
    # Guess
    file = _cached_map_file_(mapdir=mapdir, **kwargs)
    if file is None: return None
    if verbose: print 'Checking', file, os.path.exists(file)
    if not os.path.exists(file): return None
    if verbose: print 'Loadind cached map from '+os.path.basename(file)
    try:
        f = open(file)
        m = cPickle.load(f)
        f.close()
        return m
    except:
        vacumm_warn('Error while loading cached basemap instance from dir: '+
            os.path.dirname(file))
        os.remove(file)
        return None

def cache_map(m, mapdir=None):
    """Cache a map if still not cached"""
    if m is None or m.resolution is None: return
    file = _cached_map_file_(m, mapdir=mapdir)
    if file is None: return
    if not os.path.exists(file):

        # Dump
        try:
            f = open(file, 'wb')
            m.ax = None
            cPickle.dump(m, f)
            f.close()
        except:
            vacumm_warn('Error while trying to cache basemap instance into: '+
                os.path.dirname(file))
            return

        # Access to all if not in user directory
        if not file.startswith(os.path.expanduser("~")):
            os.chmod(file, stat.S_IROTH+stat.S_IWOTH+stat.S_IWGRP+stat.S_IRGRP+stat.S_IWUSR+stat.S_IRUSR)

        # Clean
        clean_cache()


def clean_cache(mapdir=None, maxsize=None):
    """Clean cache directory by checking its size

    :Params:

        - **mapdir**, optional: Directory where maps are cached
        - **maxsize**, optional: Maximal size of directory in bytes.
          Default value from :confopt:`[vacumm.misc.grid.basemap]max_cache_size`
          configuration value.
    """
    from ...misc.misc import dirsize
    mapdir = get_map_dir(mapdir)
    if mapdir is None:
        mapdir = os.path.join(get_configdir(), 'basemap', 'cached_maps')
    cache_size = dirsize(mapdir)
    if maxsize is None:
        maxsize = eval(get_config_value('vacumm.misc.grid.basemap', 'max_cache_size'))
    if cache_size>maxsize:
        files = [os.path.join(mapdir, ff) for ff in os.listdir(mapdir)]
        files.sort(cmp=lambda f1, f2: cmp(os.stat(f1)[8], os.stat(f2)[8]))
        for ff in files:
            cache_size -= os.path.getsize(ff)
            try:
                os.remove(ff)
            except:
                vacumm_warn('Error while trying to clean basemap cache in: '+
                    os.path.dirname(ff))
                return
            if cache_size<=maxsize: break

def reset_cache(mapdir=None):
    """Remove all cached maps"""
    mapdir = get_map_dir(mapdir)
    for file in [os.path.join(mapdir, ff) for ff in os.listdir(mapdir)]:
        os.remove(file)

def get_map_dir(mapdir=None):
    """Get the directory where cqched maps are stored"""
    if mapdir is None:
        mapdir = os.path.join(get_configdir(), 'basemap', 'cached_maps')
    return mapdir

def _cached_map_file_(m=None, mapdir=None, **kwargs):
    mapdir = get_map_dir(mapdir)
    if not os.path.exists(mapdir):
        os.makedirs(mapdir)
    if m is None:
        if kwargs.has_key('resolution') and kwargs['resolution'] is None:
            return None
        res = kwargs['resolution']
        kwargs['resolution'] = None
        m = Basemap(**kwargs)
    elif m.resolution is None:
        return None
    else:
        res = m.resolution
    srs = m.srs.replace(' ', '')+'+res='+res
    szone = '+%.5f+%.5f+%.5f+%.5f' % (m.llcrnrlon, m.llcrnrlat, m.urcrnrlon, m.urcrnrlat)
#    bversion = '.'.join(basemap_version.split('.')[:2])
    return os.path.join(mapdir, 'basemap-%s.%s.%s.pyk' % (basemap_version, srs, szone))


def create_map(lon_min=-180., lon_max=180., lat_min=-90., lat_max=90.,
        projection='cyl', resolution='auto', epsg=None,
        lon_center=None, lat_center=None, lat_ts=None, zoom=None, ax=None,
        overlay=False, fullscreen=False, nocache=False, cache_dir=None, **kwargs):
    """Generic creation of a :class:`Basemap` instance with caching

    .. todo:: Merge :func:`get_map` with :func:`create_map`
    """
    kwmap = kwfilter(kwargs, 'basemap', defaults={'area_thresh':0.})
    kwmap.update(kwfilter(kwargs, 'map_'))

    # Map arguments
    kwargs.setdefault('area_thresh', 0.)
    kwargs.setdefault('rsphere', RSHPERE_WGS84) # WGS-84
    if kwargs['rsphere'] in [None, False, True]: del kwargs['rsphere']
    projection = kwargs.pop('proj', projection)
    if lon_center is None: lon_center = .5*(lon_min+lon_max)
    if lat_center is None: lat_center = .5*(lat_min+lat_max)
    if lat_ts is None: lat_ts = lat_center
    if lon_max-lon_min<1.e-5:
        lon_min = lon_center-1
        lon_max = lon_center+1
    if lat_max-lat_min<1.e-5:
        lat_min = N.clip(lat_center-1, 0, 90)
        lat_max = N.clip(lat_center+1, 0, 90)
    if isinstance(zoom, (int, float)):
        lon_min, lat_min, lon_max, lat_max = zoombox(
            [lon_min, lat_min, lon_max, lat_max], zoom)

    # Special overlay case
    if overlay:
        projection = 'merc'
        resolution = None
        lat_center = 0
        lon_center = 0
    elif projection == None:
        projection = 'cyl'

    # Guess resolution
    res = kwargs.pop('res', resolution)
    if res is True:
        res = 'auto'
    elif res is False or res=='None':
        res = None
    elif isinstance(res, int):
        if res < 0:
            res= 'auto'
        else:
            res = GSHHS_RESLIST[4-res]
    if res == 'auto':
        res = gshhs_autores(lon_min, lon_max, lat_min, lat_max)
    if res in GSHHS_RESLIST:
        kwargs.setdefault('resolution', res)
    else:
        kwargs['resolution'] = None

    # Basemap args
    if isinstance(projection, str) and projection.lower() == 'rgf93' :
        # RGF93
        kwargs.update(lon_0=3, lat_0=46.5, lat_1=44, lat_2=49,
            rsphere=RSHPERE_WGS84, projection='lcc')
    else:
        # standard
        kwargs.setdefault('lon_0', lon_center)
        kwargs.setdefault('lat_0', N.clip(lat_center, -90, 90))
        kwargs.setdefault('lat_1', kwargs['lat_0'])
        kwargs.setdefault('lat_2', kwargs['lat_1'])
        kwargs.setdefault('lat_ts', N.clip(lat_center, -90, 90))
        kwargs.setdefault('projection', projection)
    kwargs.setdefault('llcrnrlon', lon_min)
    kwargs.setdefault('urcrnrlon', lon_max)
    kwargs.setdefault('llcrnrlat', N.clip(lat_min, -90, 90))
    kwargs.setdefault('urcrnrlat', N.clip(lat_max, -90, 90))
    kwargs['epsg'] = epsg

    # Check the cache
    kwcache = kwargs.copy()
    if cache_dir is not None:
        kwcache['mapdir'] = cache_dir
    if not nocache:
        mymap = cached_map(**kwcache)
    else:
        mymap = None

    # Create the map object
    if mymap is None:
        mymap = Basemap(ax=ax, **kwargs)

        # Cache it?
        if int(nocache)<=1:
            if cache_dir is not None:
                kwcache = {'mapdir':cache_dir}
            else:
                kwcache = {}
            cache_map(mymap, **kwcache)
    elif ax is not None:
        mymap.ax = ax
    mymap.res = res
    return mymap


def get_map(gg=None, proj=None, res=None, auto=False, **kwargs):
    """Quickly create :class:`Basemap` instance

    :Params:

        - **gg**, optional: cdms grid or variable, or (xx,yy).
        - **res**, optional: Resolution.
        - **proj**, optional: Projection [default: None->'merc']
        - **auto**, optional: If True, get geo specs according to grid. If False, whole earth.
          If None, auto = res is None.

    .. todo:: Merge with :func:`create_map`
    """
    from vacumm.misc.grid import misc
    if proj is None:
        proj = 'merc'
    if auto is None:
        auto = res is not None
    if gg is None: auto = False
    kwmap = dict(resolution=res, projection=proj)
    if auto:
        xx, yy = misc.get_xy(gg, proj=False)
        lat_center = yy.mean()
        lon_center = xx.mean()
        kwmap.update(
            llcrnrlon = xx.min(),
            urcrnrlon = xx.max(),
            llcrnrlat = yy.min(),
            urcrnrlat = yy.max())
    else:
        lat_center = 0.
        lon_center = 0.
    return Basemap(lat_ts=lat_center, lat_0=lat_center, lon_0=lon_center,  **kwmap)

class GSHHS_BM(Shapes):
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


#    def __init__(self, input=None, clip=None, sort=True, reverse=True, proj=False, **kwargs):
#        # Clipping argument
#        if clip is not None:
#            clip = create_polygon(clip)
#
#        from_map = not isinstance(input, Shapes)
#        if not from_map:
#            # Already a Shapes instance
#            self._m = input._m
#            self._proj = input.get_proj(proj)
#            polys = input.get_shapes(proj=proj)
#
#        else:
#            # Get the map
#            if isinstance(input, Basemap):
#                assert input.resolution is not None, 'Your map needs its resolution to be set'
#                m = input
#            else:
#                if isinstance(input, str):
#                    kwargs['res'] = input
#                elif isinstance(input, dict):
#                    kwargs.update(input)
#
#                # Map extension from clip bounds
#                if clip is not None:
#                    bb = clip.boundary
#                    kwargs.setdefault('lon_min', bb[:, 0].min())
#                    kwargs.setdefault('lon_max', bb[:, 0].max())
#                    kwargs.setdefault('lat_min', bb[:, 1].min())
#                    kwargs.setdefault('lat_max', bb[:, 1].max())
#
#                # Default resolution is 'i' if nothing to estimate it
#                if not kwargs.has_key('res') and not kwargs.has_key('resolution') and \
#                    (   (not kwargs.has_key('lon') and
#                            (not kwargs.has_key('lon_min') or not kwargs.has_key('lon_max'))) or
#                        (not kwargs.has_key('lat') and
#                            (not kwargs.has_key('lat_min') or not kwargs.has_key('lat_max')))):
#                    kwargs['res'] = 'i'
#
#                # Check lats
#                if kwargs.has_key('lat_min'): kwargs['lat_min'] = max(kwargs['lat_min'], -90)
#                if kwargs.has_key('lat_max'): kwargs['lat_max'] = min(kwargs['lat_max'], 90)
#
#                # Build the map
#                m = create_map(**kwargs)
#
#
#            polys = m.coastpolygons
#            self._m = m
#            self._proj = proj
#
#        # Convert to GEOS polygons and clip
#        self._shapes = []
#        for i, pp in enumerate(polys):
#
#            # Get the polygon with good projection
#            if not from_map:
#                poly = pp
#            else:
#                if m.coastpolygontypes[i] in [2,4]: continue # Skip lakes
#                if callable(proj) and m.projection!='cyl': # Project back for reprojection
#                    pp = m.projtran(pp[0], pp[1], inverse=True)
#                poly = create_polygon(pp, proj=proj)
#
#            # Clip it
#            self._shapes.extend(clip_shape(poly, clip))
#
#        # Save some info
#        self._info = []
#        self._type = 2
#        self._shaper = Polygon
#        if self._shapes:
#            xy = N.concatenate([s.boundary for s in self._shapes])
#            self.xmin = xy[:, 0].min()
#            self.xmax = xy[:, 0].max()
#            self.ymin = xy[:, 1].min()
#            self.ymax = xy[:, 1].max()
#            del xy
#        else:
#            xmin = N.inf
#            xmax = -N.inf
#            ymin = N.inf
#            ymax = -N.inf
#
#        # Sort polygons?
#        if sort:
#            self.sort(reverse=reverse)

def merc(lon=None, lat=None, **kwargs):
    """Mercator map

    - Extra keywords are passed to :class:`mpl_toolkits.basemap.Basemap`
    - *lon*: Longitudes to define ``llcrnrlon`` and ``urcrnrlon``
    - *lat*: Latitudes to define  ``lat_ts``, ``llcrnrlat`` and ``urcrnrlat``
    """
    kwargs.setdefault('resolution', None)
    if lon is not None:
        lon = N.asarray(lon)
        kwargs.setdefault('llcrnrlon', lon.min())
        kwargs.setdefault('urcrnrlon', lon.max())
    if lat is not None:
        lat = N.asarray(lat)
        kwargs.setdefault('llcrnrlat', lat.min())
        kwargs.setdefault('urcrnrlat', lat.max())
        lat_ts = N.median(lat)
    else:
        lat_ts = 0.
    kwargs.setdefault('lat_ts',lat_ts)
    return Basemap(projection='merc', **kwargs)

#proj = merc()

def basic_proj(xx, yy, inverse=False):
    """A basic projection using :func:`vacumm.misc.phys.units.deg2m`
    and :func:`vacumm.misc.phys.units.deg2m`
    """
    if inverse:
        yy = m2deg(yy)
        return m2deg(xx, yy), yy
    return deg2m(xx, yy), deg2m(yy)

def get_proj(gg=None, proj=None, **kwargs):
    """Setup a default projection using x,y coordinates and
    :class:`~mpl_toolkits.basemap.proj.Proj` or :func:`basic_proj`

    Projection is set by default to "basic".

    :Params:

        - **gg**, optional: Grid or coordinates (see :func:`~vacumm.misc.grid.misc.get_xy`).
          If not provided, lon bounds are set to (-180,180) and lat bounds to (-89.99,89.99).
        - Other keywords are passed to :class:`~mpl_toolkits.basemap.proj.Proj`. One of
          them is the projection type, which defaults to configuration option
          :confopt:`[vacumm.misc.grid.basemap] proj`.

    :Return:
        A :class:`mpl_toolkits.basemap.proj.Proj` instance
        or :func:`basic_proj`

    :Examples:

        >>> proj = get_proj(sst.getGrid(), proj='laea')
        >>> x, y = proj(lon, lat)

        >>> proj = get_proj((lon, lat))
        >>> xx, yy = N.meshgrid(lon, lat)
        >>> xx, yy = proj(xx, yy)
        >>> print proj(xx, yy, inverse=True)

        >>> proj = get_proj(R=6000000.)
    """
    if proj is False:
        return False
    if proj is None or proj is True:
        proj = 'basic'
    if proj == 'basic':
        return basic_proj
    if callable(proj): return proj
    if gg is not None:
        x,y = get_xy(gg, num=True)
        xmin, ymin, xmax, ymax = x.min(), y.min(), x.max(), y.max()
    else:
        xmin, ymin, xmax, ymax = -180, -90, 180, 90
        y = [0]
    projparams = kwargs.copy()
    ymax = min(ymax, 89.99)
    ymin = max(ymin, -89.99)
    if not isinstance(proj, str):
        proj = get_config_value('vacumm.misc.grid.basemap', 'proj')
    dict_check_defaults(projparams, R=rsphere_mean, units='m',
        proj=proj,
        lat_ts = N.median(y) if len(y)>10 else N.mean(y),
        lon_0 = N.median(x) if len(x)>10 else N.mean(x))
    dict_check_defaults(projparams, lat_0=projparams['lat_ts'])
    return Proj(projparams, xmin, ymin, xmax, ymax)


from masking import polygons, create_polygon, proj_shape, clip_shape
