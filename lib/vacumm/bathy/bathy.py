# -*- coding: utf8 -*-
"""
Bathmetry tools
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
from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
import six
import os.path
from io import StringIO
from collections import OrderedDict
from operator import itemgetter

import cdms2, MV2
from matplotlib.ticker import Formatter
import matplotlib.pyplot as P
from configobj import ConfigObj
import numpy as N

from vacumm import VACUMMError, vcwarn, VACUMM_CFG
from vacumm.config import check_data_file
from vacumm.misc.config import ConfigManager
from vacumm.misc.misc import kwfilter, auto_scale
from vacumm.misc.color import cmap_custom, auto_cmap_topo
from vacumm.misc.cf import format_var
from vacumm.misc.axes import create_axes2d
from vacumm.misc.grid import (get_xy, resol, create_grid, meshgrid,
                              meshbounds, get_axis)
from vacumm.misc.masking import polygon_mask
from vacumm.misc.regridding import refine, regrid2d, GriddedMerger
from vacumm.misc.plot import map2
from vacumm.misc.filters import norm_atan, deriv2d, shapiro2d
from vacumm.misc.sdata import XYZ, XYZMerger, Shapes

from .shorelines import (get_best_shoreline, get_shoreline, ShoreLine,
                         get_best_res)


__all__ = ['print_bathy_list', 'Bathy', 'plot_bathy', 'fmt_bathy',
           'bathy_list', 'XYZBathy', 'XYZBathyBankClient', 'XYZBathyBank',
           'XYZBathyMerger', 'GriddedBathy', 'GriddedBathyMerger',
           'NcGriddedBathy', 'NcGriddedBathyError', 'format_bathy',
           ]


def bathy_list(name=None, cfgfile=None, minstatus=0):
    """Get the list of bathymetries as dict

    :Params:

        - **name**, optional: Request this bathymetry.
        - **cfgfile**, optional: Load this bathymetry config file.
        - **minstatus**, optional: Minimal status of availability
          of the bathymetry (see :func:`vacumm.config.check_data_file`):

            - ``0``: Bathy file not on disk and not downloadable
              (however, user can specify the path manually).
            - ``1``: File is not on disk, but downloadable.
            - ``2``: File is on disk.

    :Return:

        - If ``name`` is not provided, a dict of subdicts like this one::

            {'etopo2':{'file':'/home/me/bathy.nc', 'var':'btdata', ...}, ...}

        - Else, the subdict of the requested bathy::

            {'file':'/home/me/bathy.nc', 'var':'btdata', ...}

    """

    # Read config file
#    sections, cfg = get_config_sections(cfg=cfgfile, parent_section=__name__,
#        parent_option='cfgfile_gridded', getcfg=True)

#    if not cfgfile: cfgfile = get_cfgfiles_gridded() # File name
#    dvars = get_dir_dict(__names__) # Various directory names
#    bathys = {}
#    f = ConfigParser()
#    f.read(cfgfile)


    # Loop on section
    bathies = {}
    for section in VACUMM_CFG['gridded'].sections:
        status = _check_bathy_file_(section, avail=True)
#        if not os.path.exists(f.get(section, 'file', vars=dvars)): continue
        if status<=minstatus:
            continue
        bathies[section] = {'status':status}
        bathies[section].update(VACUMM_CFG['gridded'][section])

    # Found any bathy?
    if name is not None:
        if not len(bathies): raise NcGriddedBathyError('No bathy are available')
        if not name in bathies:
            raise NcGriddedBathyError('Wrong bathymetry name: "%s". Please choose one of these: %s.'%
                (name, ', '.join(list(bathies.keys()))))
        return bathies[name]
    return bathies

def _check_bathy_file_(bathyname, **kwargs):
    """Make sure to have a bathy file (see :func:`~vacumm.config.check_data_file`)"""
    return check_data_file(VACUMM_CFG['gridded'][bathyname], **kwargs)


#def get_cfgfiles_gridded():
#    """Name of the configuration file containing the list of netcdf gridded bathy files"""
#    return get_config_value(__name__, 'cfgfile_gridded')

def print_bathy_list(cfgfile=None):
    """Print available bathymetries"""
#    if not cfgfile: cfgfile = get_cfgfiles_gridded()
    print('List of bathymetries:')
    for bathyname,vv in bathy_list(cfgfile=cfgfile).items():
        print('\n[%s]' % bathyname)
        print('  status: '+['Not available (need user input)',
            'Downloadable', 'On disk'][vv['status']])
        for key in 'file','file_url','file_license','var','lon','lat':
            if key in vv:
                print('  - %s: %s' % (key,vv[key]))
list_bathy = print_bathy_list

class Bathy(object):
    """Get a portion of the gridded bathymetry

    .. warning::

        This class is **DEPRECATED**.
        Please use :class:`NcGriddedBathy` instead.

    :Params:

        - *lon*: Longitude range [default: None]
        - *lat*: Latitude range [default: None]
        - *sampling*: Sampling in both directions [default: 1]
        - *fill*: If True, fill positive values with 0. [default: True]
        - *mask_function*: To mask points. Should be 'greater' or 'less'.
          If empty or None, no mask is applied [default: 'greater']
        - *mask_value*: Value used with mask_function [default: 0.]
        - *name*: Name of bathymetry within is list given by print_bathy_list()

    :Example:


    >>> from vacumm import Bathy
    >>> gg = Bathy((-9,10),(40,50),'etopo2') # object creation
    >>> gg.show() # quick plot
    >>> gg.zoom((-5,0),(42,50)) # Select another zoom (lon/lat)
    >>> mybathyvar = gg() # Get the current zoom
    >>> mybathyvar = gg((-9,-8),(45,50)) # Get another zoom
    >>> gg.write_ascii('/home/yoman/mybathy.dat') # dump to file
    >>> gg.write_netcdf('/home/yoman/mybathy.nc') # same for a netcdf file

    """


    def __init__(self, lon=None, lat=None, name='etopo2', **kwargs):
        # Choose the bathymetry
        bathy_name = kwargs.get('bathy_name',name)
        mybathy = bathy_list(bathy_name)
        self.file_name = mybathy['file']
        self._lon_name = mybathy['lon']
        self._lat_name = mybathy['lat']
        self._var_name = mybathy['var']
        self.bathy_name = bathy_name

        # Default values
        self._lon_range = None
        self._lat_range = None
        self._fill = True
        self._sampling =1

        # Zoom
        cdms2.axis.longitude_aliases.append(self._lon_name)
        cdms2.axis.latitude_aliases.append(self._lat_name)
        self.zoom(lon,lat,**kwargs)

    def __call__(self):

        self.zoom(self._lon_range,self._lat_range,self._sampling,self._fill)
        return self.bathy


    def zoom(self,lon=None,lat=None,sampling=None,fill=None,mask_function='greater',mask_value=0.):
        """Select a zoom and set the current variable

        :Params:

        - *lon*: Longitude range [default: None]
        - *lat*: Latitude range [default: None]
        - *sampling*: Sampling in both directions [default: 1]
        - *fill*: If True, fill positive values with 0. [default: True]
        - *mask_function*: To mask points. Should be 'greater' or 'less'.
          If empty or None, no mask is applied [default: 'greater']
        - *mask_value*: Value used with mask_function [default: 0.]
        """

        if lon is None: lon = self._lon_range
        if lat is None: lat = self._lat_range
        if lon == self._lon_range and lat == self._lat_range: return
        if sampling is not None: self._sampling = sampling
        if fill is not None: self._fill = fill

        self._lon_range = lon
        self._lat_range = lat
#        kwargs = {self._lon_name:self._lon_range,self._lat_name:self._lat_range}

        f = cdms2.open(self.file_name)
        ilonmin,ilonmax,ilonsamp = f[self._var_name].getLongitude().mapIntervalExt(self._lon_range)
        ilatin,ilatmax,ilatsamp = f[self._var_name].getLatitude().mapIntervalExt(self._lat_range)
        ilonsamp = ilatsamp = self._sampling
        self.bathy = f[self._var_name][ilatin:ilatmax:ilatsamp,ilonmin:ilonmax:ilonsamp]
        f.close()
        self.bathy.setMissing(mask_value)
        if mask_function not in ['',None]:
            self.bathy[:] = getattr(MV,'masked_'+mask_function)(self.bathy[:],mask_value)
        if self._fill:
            self.bathy[:] = self.bathy.filled()
#           self.bathy.unmask()
        self.bathy.id = 'bathy' ; self.bathy.name = self.bathy.id
        self.bathy.long_name = self.bathy_name[0].upper()+self.bathy_name[1:]+' bathymetry'
        if 'units' not in self.bathy.attributes:
            self.bathy.units = 'm'
        lon_axis = self.bathy.getAxis(1)
        lon_axis.id = 'lon'
        lon_axis.long_name = 'Longitude'
        lon_axis.units = 'degrees_east'
        lat_axis = self.bathy.getAxis(0)
        lat_axis.id = 'lat'
        lat_axis.long_name = 'Latitude'
        lat_axis.units = 'degrees_north'
        self.bathy.setAxisList([lat_axis,lon_axis])

        self._lon2d,self._lat2d = N.meshgrid(lon_axis.getValue(),lat_axis.getValue())
        self._lon2dp,self._lat2dp = create_axes2d(lon_axis.getValue(),lat_axis.getValue(),bounds=True,numeric=True)
        self._xxs = self._yys = None


    def lon(self):
        """ Longitude axis """
        return self.bathy.getLongitude()

    def lat(self):
        """ Latitude axis """
        return self.bathy.getLatitude()


    def show(self,lon=None,lat=None,**kwargs):
        """Plot the current bathymetry

        - *lon*: Longitude range.
        - *lat*: Latitude range.

        .. seealso::

            :func:`plot_bathy`
        """

        # Zoom if requested
        self.zoom(lon,lat)

        return plot_bathy(self, **kwargs)


    plot=show


    def write_netcdf(self,file):
        """Write the bathymetry as netcdf file

        - **file**: Output netcdf file name
        """
        assert file.endswith('.nc'),'Your netcdf file name must end with .nc'
        f = cdms2.open(file,'w')
        f.write(self.bathy)
        f.from_bathy_file = self.file_name
        f.close()


    def write_ascii(self,file,fmt='%10f %10f %10f'):
        """Write the bathymetry as an ascii file in the format : x y z

        - **file**: Output file name

        - *fmt*: Format of lines [default: '%10f %10f %10f']
        """

        bathy1d = self.bathy.ravel()
        f = open(file,'w')
        for i in range(len(bathy1d)):
            f.write((fmt+'\n') % (self._lon2d.flat[i],self._lat2d.flat[i],bathy1d[i]))
        f.close()
        del bathy1d
        print('Bathy wrote to ascii file',file)


def plot_bathy(bathy, shadow=True, contour=True, shadow_stretch=1.,
               shadow_shapiro=False, show=True, shadow_alpha=1.,
               shadow_black=.3, white_deep=False, nmax_levels=30, m=None, alpha=1.,
               zmin=None, zmax=None, **kwargs):
    """Plot a bathymetry

    - *lon*: Longitude range.
    - *lat*: Latitude range.
    - *show*:Display the figure [default: True]
    - *pcolor*: Use pcolor instead of contour [default: False]
    - *contour*: Add line contours [default: True]
    - *shadow*:Plot south-west shadows instead of filled contours.
    - *nmax_levels*: Max number of levels for contours [default: 30]
    - *white_deep*: Deep contours are white [default: False]
    - All other keyword are passed to :func:`~vacumm.misc.plot.map2`
    """

    # Input
    bb = bathy
    if isinstance(bathy, GriddedBathy):
        bathy = bathy.bathy()
        if shadow:
            xxs = getattr(bb, '_xxs', None)
            yys = getattr(bb, '_yys', None)
            if xxs is None:
                lon2d = bb._lon2d
                lat2d = bb._lat2d
    elif shadow:
        xxs = yys = None
        lon2d,lat2d = meshgrid(get_axis(bathy, -1).getValue(),
                               get_axis(bathy, -2).getValue())

    # Masking
    if 'maxdep' in kwargs:
        zmin = -kwargs['maxdep']
    if 'maxalt' in kwargs:
        zmax = kwargs['maxalt']
    if zmin is not None:
        bathy[:] = MV2.masked_less(bathy, zmin)
    if zmax is not None:
        bathy[:] = MV2.masked_greater(bathy, zmax)

    # Default arguments for map
    if hasattr(bathy, 'long_name'):
        kwargs.setdefault('title', bathy.long_name)
    if 'cmap' not in kwargs:
        vmin = bathy.min()
        vmax = bathy.max()
#        print 'cmap topo', vmin, vmax
        kwargs['cmap'] = auto_cmap_topo((kwargs.get('vmin', vmin),
              kwargs.get('vmax', vmax)))
        kwargs.update(vmin=kwargs['cmap']._vacumm_topomin,
                      vmax=kwargs['cmap']._vacumm_topomax,
                      levels_keepminmax=True)
#       kwargs.setdefault('ticklabel_size','smaller')
    kwargs.setdefault('clabel_fontsize', 8)
    kwargs.setdefault('clabel_alpha',.7*alpha)
    kwargs.setdefault('clabel_glow_alpha', kwargs['clabel_alpha'])
    kwargs.setdefault('fill', 'contourf')
    kwargs['nmax_levels'] = kwargs.get('nmax', nmax_levels)
    kwargs['show'] = False
    kwargs['contour'] = contour
    if shadow:
        kwargs.setdefault('alpha',.5*alpha)
    kwargs.setdefault('projection', 'merc')
    kwargs.setdefault('fmt', BathyFormatter())
    kwargs.setdefault('colorbar_format', BathyFormatter())
    kwargs.setdefault('units', False)
    kwargs.setdefault('levels_mode', 'normal')
    kwargs.setdefault('bgcolor', '0.8')
    kwargs.setdefault('contour_linestyle', '-')
    savefig = kwargs.pop('savefig', None)
    kwsavefig = kwfilter(kwargs, 'savefig_')


    # White contour when dark
    if contour and white_deep:
        if hasattr(kwargs['cmap'], '_vacumm_topomin'):
            data = (kwargs['cmap']._vacumm_topomin,
                    kwargs['cmap']._vacumm_topomax)
        else:
            data = bathy
        levels = auto_scale(data, nmax=nmax)
        colors = []
        nlevel = len(levels)
        for i in range(nlevel):
            if i < nlevel/2:
                colors.append('w')
            else:
                colors.append('k')
        kwargs['levels'] = levels
        kwargs.setdefault('contour_colors',tuple(colors))

    # Call to map
    m = map2(bathy, m=m, **kwargs)

    # Add shadow
    if shadow:

        # Filter
        data = MV2.array(bathy,'f',fill_value=0.)
        if shadow_shapiro:
            data = shapiro2d(data,fast=True).shape

        # Gradient
        grd = deriv2d(data,direction=45.,fast=True,fill_value=0.).filled(0.)
        grdn = refine(grd, 3)
        grdn = norm_atan(grdn,stretch=shadow_stretch).clip(0,1.) ; del grd

        # Grid
#           im = m.map.imshow(grdn,cmap=P.get_cmap('gist_yarg'),alpha=1) # gist_yarg , YlGnBu
        if xxs is None or yys is None:
            xx, yy = m(lon2d,lat2d)
            xxr = refine(xx, 3)
            yyr = refine(yy, 3)
            xxs, yys = meshbounds(xxr, yyr)
            if isinstance(bb, GriddedBathy):
                bb._xxs = xxs
                bb._yys = yys
            del xx, yy, xxr, yyr

        # Cmap
        cmap = cmap_custom(( ((1, )*3, 0), ((shadow_black, )*3, 1) ))

        # Plot
        pp = m.map.pcolormesh(xxs, yys, grdn, cmap=cmap)#P.get_cmap('gist_yarg'))
        pp.set_zorder(.9)
        pp.set_linewidth(0)
        pp.set_alpha(shadow_alpha*N.clip(alpha*2, 0, 1))
        del grdn

    # Show it?
    if savefig:
        m.savefig(savefig, **kwsavefig)
    if show:
        P.show()
    return m

bathy=Bathy


class BathyFormatter(Formatter):
    """Label formatter for contour"""

    def __init__(self, fmt='%i m', signed=True):
        self.signed = signed
        self.fmt = fmt

    def __mod__(self, val):
        s = self.fmt % abs(val)
        if self.signed and val < 0:
            s = '-'+s
        return s

    def __call__(self, val, i=None):
        return self % val



class XYZBathy(XYZ):
    """For random bathymetries, like mnt files

        - **xyz**: A (x,y,z) tuple or a xyz file.
        - *long_name*: Long name (like a title) for this bathymetry.

    .. note::

        Bathymetry is supposed to be positive over water.

    :Usage:

        >>> xzy1 = XZY('bathy1.xzy', long_name='My bathy')           # Load
        >>> xyz1.select([-6,45,-4,48])
        >>> xyz1 += XZY('bathy2.xzy')
        >>> print xyz1.x
        >>> bathy = xyz1.togrid(mygrid)

    .. seealso::

        :class:`~vacumm.misc.io.XYZ` :class:`~vacumm.bathy.shoreline` :class:`~vacumm.tide.station_info.MeanSeaLevel`
    """
    def __init__(self, xyz, check_sign=False, *args, **kwargs):
        # Initialize
        kwargs.setdefault('units', 'm')
        kwargs.setdefault('long_name', 'Bathymetry')
        XYZ.__init__(self, xyz, *args, **kwargs)

        # Negative bathy
        if check_sign:
            zmax = self.zmax
            zmin = self.zmin
            if zmax > 0. and (zmin > 0. or abs(zmin) < zmax):
                self._z[:] *= -1.

    def togrid(self, grid=None, mask=None, cgrid=False,  **kwargs):
        """Interpolate data to a grid"""

        # Grid
        kwgrid = kwfilter(kwargs, 'grid_')
        if grid is None:
            grid = self.get_grid(**kwgrid)


        # Shoreline for masking
        if (mask is not None and mask is not False and
                mask is not MV2.nomask):

            # Clip for polygons
            kwmask = kwfilter(kwargs, 'mask')
            clip = kwmask.pop('clip', .1)
            if isinstance(clip, float):
                xx, yy = get_xy(grid)
                xmin = xx[:].min() ; xmax = xx[:].max()
                ymin = yy[:].min() ; ymax = yy[:].max()
                dx = (xmax-xmin)*clip
                dy = (ymax-ymin)*clip
                clip = [xmin-dx, ymin-dy, xmax+dx, ymax+dy]

            # Auto resolution
            if mask is True:
                mask = 'auto'
            if mask == 'auto' and False:
                if grid is None:
                    grid = self.get_grid(**kwfilter(kwargs, 'grid'))
                mask = get_best_shoreline(grid)

            # Direct handling
            if isinstance(mask, (str, Shapes)):
                # Get mask
                mask = get_shoreline(mask, clip=clip)

        # Normal case
        return XYZ.togrid(self, grid=grid, mask=mask, cgrid=cgrid,  **kwargs)
    togrid.__doc__ = XYZ.togrid.__doc__

    def bathy(self, *args, **kwargs):
        """Alias for :meth:`togrid`"""
        return self.togrid(*args, **kwargs)

    def plot(self, m=True, **kwargs):
        kwargs.update(m=m)
        return XYZ.plot(self, **kwargs)
    plot.__doc__ = XYZ.plot.__doc__+"""

    .. seealso::

        For value of default parameters: :meth:`~vacumm.misc.io.XYZ.plot`
    """


XYZ_BATHY_BANK_CFG_SPECS = """
[__many__]
xyzfile=path(default=None,expand=True)
xmin=float(default=None)
xmax=float(default=None)
ymin=float(default=None)
ymax=float(default=None)
transp=boolean(default=True)
long_name=string(default=None)
zorder=integer(default=0)
"""

XYZ_BATHY_RANK_CFGM = ConfigManager(XYZ_BATHY_BANK_CFG_SPECS.split('\n'))


class XYZBathyBankClient(object):
    """An single element (bathy infos) of the :class:`XYZBathyBank`

    - **bank**: A bank instance.
    - **id**:Id of this bathy.
    - Other keywords are set as attributes (like long_name, etc).

    :Usage:

    Let ``xyzbc`` being an element of the bank.

    >>> print xyzbc
    >>> xyzbc.long_name = 'New long name'
    >>> xyzbc.xmax = -2.
    >>> xyzbc.transp = False
    >>> xyzbc.id = 'iroise'             # Rename in the bank
    >>> xyzbc.xyzfile = 'newfile.mnt'   # xmin, xmax, ymin, ymax updated
    >>> xyz = xyzbc.load()
    """

    def __init__(self, bank, id, **kwargs):
        self.bank = bank
        if id not in self.bank.cfg:
            self.bank.cfg[id] = {}
        self._id = id
        self._xyz = None
        if kwargs:
            self.update(**kwargs)
        self._validate_()
        pass

    @property
    def cfg(self):
        return self.bank.cfg[self.id]

    def _validate_(self):
        self.bank.validate_cfg()
        for key, val in self.cfg.items():
            if val is None or val == '':
                self.cfg[key] = self[key]



    def __getitem__(self, key):

        if key == 'id':
            return self.id
        value = self.cfg[key]
        if key == 'long_name':
            return value or self.id.capitalize()
        if key in ('xmin', 'xmax', 'ymin', 'ymax'):
            if value is None:
                xyz = self.xyz
                if xyz is not None:
                    value = self.cfg[key] = getattr(xyz, key)
        return value

    def __setitem__(self, key, value):
#        assert key in self.cfg
        if key == 'id':
            self.bank.rename(self._id, value)
        else:
            self.cfg[key] = value
            self._validate_()

    def get_id(self):
        return self._id

    def set_id(self, value):
        self['id'] = value

    id = property(get_id, set_id, doc='Bathymetry id')

    def get_xmin(self):
        return self['xmin']

    def set_xmin(self, value):
        self['xmin'] = value

    xmin = property(get_xmin, set_xmin, doc='X min')

    def get_xmax(self):
        return self['xmax']

    def set_xmax(self, value):
        self['xmax'] = value

    xmax = property(get_xmax, set_xmax, doc='X max')

    def get_ymin(self):
        return self['ymin']

    def set_ymin(self, value):
        self['ymin'] = value

    ymin = property(get_ymin, set_ymin, doc='Y min')

    def get_ymax(self):
        return self['ymax']

    def set_ymax(self, value):
        self['ymax'] = value

    ymax = property(get_ymax, set_ymax, doc='Y max')

    def get_long_name(self):
        return self['long_name']

    def set_long_name(self, value):
        self['long_name'] = value

    long_name = property(get_long_name, set_long_name, doc='Long name')

    def get_transp(self):
        return self['transp']

    def set_transp(self, value):
        self['transp'] = value

    transp = property(get_transp, set_transp, doc='Transparency')

    def get_zorder(self):
        return self['zorder']

    def set_zorder(self, value):
        self['zorder'] = value

    zorder = property(get_zorder, set_zorder, doc='Z order')

    def update(self, **kwargs):
        for key, value in kwargs.items():
            self[key] = value

    def get_xyz(self):
        """Load the file as a :class:`XYZBathy` instance"""
        if hasattr(self, '_xyz') and self._xyz is not None:
            return self._xyz
        if not self['xyzfile'] or not os.path.exists(self['xyzfile']):
            vcwarn('Not a valid xyzfile name: ' + str(self['xyzfile']))
            return None
        self._xyz = XYZBathy(self['xyzfile'], long_name=self['long_name'],
                             transp=self['transp'], id=self.id)
        return self._xyz

    xyz = property(get_xyz, doc='XYZBathy instance')

    def __str__(self):
        sio = StringIO()
        cfg = ConfigObj(self.cfg.dict())
        cfg.write(sio)
        lines = sio.getvalue()
        sio.close()
        return lines

    def save(self):
        """Save the bank"""
        self.bank.save()

    def load(self, zone=None, margin=None):
        """Load the bathymetry as a :class:`XYZBathy` instance

        .. seealso::

            :meth:`~vacumm.misc.io.XYZ.clip`
        """
        xyz = self.xyz
        if xyz is not None and zone is not None:
            return xyz.clip(zone=zone, margin=margin)
        return xyz

    def plot(self, id, **kwargs):
        """Load and plot this bathy"""
        return self.load().plot(**kwargs)


class XYZBathyBank(object):
    """Bank of XYZ bathymetry infos

    .. seealso::

        :class:`XYZBathyBankClient`

    :Usage:

    >>> xyzb = XYZBathyBank()

    """

    def __init__(self, cfg=None):

        # Config
        self.load_cfg(cfg)
        self.cfgfile = self.cfg.filename

        # Populate
        self._clients = OrderedDict()
        for id in self.cfg.sections:
            self._clients[id] = XYZBathyBankClient(self, id)
        self._order = None

    def load_cfg(self, cfg):
        """Load a configuration and validate it"""
        cfg_update = XYZ_BATHY_RANK_CFGM.load(cfg, validate='check')
        if hasattr(self, 'cfg'):
            self.cfg.merge(cfg_update)
        else:
            self.cfg = cfg_update
        return self.cfg

    def validate_cfg(self):
        """Revalidate the configuration"""
        cfg = self.cfg
        del self.cfg
        return self.load_cfg(cfg)

    @property
    def ids(self):
        """Return the list of bathy ids"""
        return self.cfg.sections

    def __str__(self):
        return '\n'.join(ConfigObj(self.cfg.dict()).write())

    def __len__(self):
        return len(self.cfg.sections)

    def check_id(self, id):
        assert id in self._clients, ('No such XYZBathyBankClient '
                                     'id: {}'.format(id))

    def __getitem__(self, id):
        if isinstance(id, XYZBathyBankClient):
            id = id._id
        self.check_id(id)
        return self._clients[id]

    def __setitem__(self, id, xyzfile):
        self.cfg[id] = {'xyzfile': xyzfile}
        self.validate_cfg()
        self._clients[id] = XYZBathyBankClient(self, id)
        self._update_order_()

    def __delitem__(self, id):
        if isinstance(id, XYZBathyBankClient):
            id = id._id
        self.check_id(id)
        del self.cfg[id], self._clients[id]
        self._update_order_()

    def __contains__(self, id):
        return id in self._clients

    def add(self, xyzfile, id=None, force=False, **kwargs):
        """Add a bathy to the bank

        Parameters
        ----------
        xyzfile: str
            xyz file.
        id: str
            Id of the bank. If ``None``, it is guessed from the file name:
           ``xyzfile="/path/toto.xyz"`` ->  ``id="toto"``
        """
        if id is None:
            id = os.path.basename(xyzfile)
            if id.endswith('.xyz'):
                id = id[:-4]
        if id in self:
            if force:
                del self.cfg[id]
            else:
                raise KeyError('Bathymetry %s already exists'%id)
        self[id] = xyzfile
        self._clients[id].update(**kwargs)
        self._update_order_()

    def __iadd__(self, d):
        self.add(d)
        return self

    def remove(self, id):
        """Remove a bathy from the bank

        *id*: A bathy id or :class:`XYZBathyBankClient` instance.
        """
        del self[id]

    def __isub__(self, d):
        self.remove(id)
        return self

    def rename(self, id1, id2):
        """Rename a bathymetry client"""
        if isinstance(id1, XYZBathyBankClient):
            id1 = id1.id
        client = self[id1]
        client._id = id2
        self.cfg[id2] = self.cfg[id1]
        self._clients[id2] = client
        del self[id1]
        self.validate_cfg()

    def copy(self, id1, id2, **kwargs):
        """Copy bathy infos to another"""
        if isinstance(id1, XYZBathyBankClient):
            id1 = id1.id
        self[id1]
        self.cfg[id2] = self.cfg[id1].copy()
        self._clients[id2] = XYZBathyBankClient(self, id2, **kwargs)

    def load(self, id, **kwargs):
        """Load a bathymetry as a :class:`XYZBathy` instance"""
        self.check_id(id)
        return self._clients[id].load(**kwargs)

    @property
    def clients(self):
        """Get all clients of the bank as a list"""
        return list(self._clients.values())

    def select(self, xmin=None, xmax=None, ymin=None, ymax=None, load=True,
               margin=None, ordered=None):
        """Return a list of bathymetries relevant to the specified region

        - *margin*: Margin for loading xyz bathies:

            - if None, the whole bathy is loaded from the bank,
            - else it should be a value relative to the approximative resolution (see :meth:`XYZBathy.resol`)
        """
        res = []
        for client in self.clients:
            if ((None not in (xmin, client.xmax) and xmin > client.xmax) or
                (None not in (xmax, client.xmin) and xmax < client.xmin) or
                (None not in (ymin, client.ymax) and ymin > client.ymax) or
                (None not in (ymax, client.ymin) and ymax > client.ymin)):
                    continue
            if load:
                client = client.load(zone=(xmin, ymin, xmax, ymax),
                                     margin=margin)
            res.append(client)
        return res

    def _sort_(self):
        clients = self.clients
        clients.sort(key=itemgetter('zorder'))
        self._clients = OrderedDict([(client.id, client)
                                     for client in clients])

#    def set_order(self, order):
#        """Set a list of ids to define the order"""
#        if order is None:
#            del self.cfg.defaults()['order']
#        else:
#            self.cfg.defaults()['order'] = str([id for id in order if id in self.ids()])
#    def get_order(self,):
#        """Get a list of ids that define the order"""
#        return eval(self.cfg.defaults().get('order', 'None'))
#
    def _update_order_(self):
        return
        if self.get_order() is None: return
        self.cfg.defaults()['order'] = str([id for id in self.get_order() if id in self.ids()])

    def save(self, cfgfile=None):
        if cfgfile is None:
            cfgfile = self.cfgfile
        else:
            self.cfgfile = cfgfile
        if cfgfile is None:
            raise VACUMMError('You must provide a valid file name')
        f = open(cfgfile, 'w')
        self.cfg.write(f)
        f.close()

    def merger(self, **kwargs):
        """Return a :class:`XYZBathyMerger` instance using the current bank

        Keywords are passed to :meth:`select`
        """
        merger = XYZBathyMerger()
        merger += self.select(**kwargs)
        return merger

    def plot(self, **kwargs):
        """Plot current bathies using a :class:`XYZBathyMerger` instance"""

        return self.merger(**kwfilter(kwargs, 'merger')).plot(**kwargs)


class XYZBathyMerger(XYZMerger):
    """Mix different bathymetries"""
    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None, long_name=None):

        self._bank = XYZBathyBank()
        XYZMerger.__init__(self, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            long_name=long_name, units='m')
        self._XYZ = XYZBathy

    def _load_(self, d):
        # Bathy from the bank
        if d in self._bank.ids:
            return self._bank[id].load()
        # Classic cases
        return XYZMerger._load_(self, d)

    def from_bank(self, margin=5):
        """Append all bathymetries from the bank that cover the grid

        - *margin*: Relative margin for clipping datasets from the bank (see :meth:`~vacumm.misc.io.XYZ.clip` and :meth:`~vacumm.misc.io.XYZ.resol`)

        All parameters are passed to :meth:`XYZBathyBank.select`
        """
        self += self._bank.select(xmin=self._xmin, xmax=self._xmax,
            ymin=self._ymin, ymax=self._ymax, load=True, margin=margin)

    def plot(self, color=None, marker=None, mode='cluster', title='XYZ merger',
        m=True, show=True, colorbar=True, savefig=None, savefigs=None, **kwargs):
        XYZMerger.plot(self, color=color, marker=marker, mode=mode, title=title,
            m=m, show=show, colorbar=colorbar, savefig=savefig, savefigs=savefigs, **kwargs)
    plot.__doc__ = XYZMerger.plot.__doc__

    def bathy(self, grid=None):
        """Get a gridded bathymetry"""
        return self.xyz().bathy(grid)


class _GriddedBathyMasked_:
    _shoreline_mask = None
    _shoreline_margin = 2
    _shoreline_arg = None
    def __init__(self, maxvalue=0., fillvalue=0., shoreline=None):
#        self._masked_bathy = None
        self.set_shoreline(shoreline)
        self.set_maxvalue(maxvalue)
        self.set_fillvalue(fillvalue)

    def set_maxvalue(self, maxvalue):
        """Set the value over which bathy is masked.

        :Params:

            - **maxvalue**: If not set to a scalar, bathy is not masked in this way.
        """
        self._maxvalue = maxvalue

    def set_fillvalue(self, fillvalue):
        """Set filling value when bathy is masked.

        :Params:

            - **fillvalue**: If not set to a scalar, bathy is not filled when masked.
        """
        self._fillvalue = fillvalue

    def set_shoreline(self, shoreline, margin=None):
        """Set the shoreline that to create a mask

        :Params:

            - **shoreline**: ``'auto'``, :class:`~vacumm.bathy.ShoreLine` instance,
              argument to :func:`~vacumm.bathy.shoreline.get_shoreline`.
              If ``None`` or ``False``, shoreline is not used.

        .. note:: The mask is applied only if :meth:`masked_bathy` is called.
        """
        if shoreline=='auto' or shoreline is True:
            shoreline = get_best_res((self.get_lon(), self.get_lat()))
        elif shoreline is False or shoreline is None:
            self._shoreline_mask = None
            self._masked_bathy = None
            shoreline = None

        # Get the shoreline
        self._shoreline_arg = shoreline
        if margin is not None: self._shoreline_margin = margin

    def get_shoreline_mask(self):
        """Get the shoreline mask"""
        # We already have it
        if self._shoreline_mask is not None:
            return self._shoreline_mask
        # We skip
        if self._shoreline_arg is None: return

        # We compute it
        lon = self.get_lon()
        lat = self.get_lat()
        xres, yres = self.get_res()
        clip = [lon.getValue().min()-self._shoreline_margin*xres,
            lat.getValue().min()-self._shoreline_margin*yres,
            lon.getValue().max()+self._shoreline_margin*xres,
            lat.getValue().max()+self._shoreline_margin*yres]
        if not isinstance(self._shoreline_arg, ShoreLine):
            shoreline = get_shoreline(self._shoreline_arg, clip=clip)
            if shoreline is None: return
        else:
            shoreline = shoreline.clip(clip)

        # Get the polygons
        polys = shoreline.get_shapes()

        # Create the mask
        self._shoreline_mask = polygon_mask(self.get_grid(), polys)
        del polys
        return self._shoreline_mask

    def reset_shoreline(self):
        """Remove the shoreline"""
        self.set_shoreline(False)

    def masked_bathy(self, **kwargs):
        """Get a masked version of the bathymetry"""
        # Mask with value
        bathy = self._bathy.clone()
        if self._maxvalue is not None and N.isreal(self._maxvalue):
            bathy[:] = MV2.masked_greater(bathy, self._maxvalue)

        # Not already available
        mask = self.get_shoreline_mask()
        if mask is not None:
            bathy[:] = MV2.masked_where(mask, self._bathy, copy=0)

        # Fill
        if self._fillvalue is not None and N.isreal(self._fillvalue):
            bathy.set_fill_value(self._fillvalue)
#            bathy[:] = bathy.filled(self._fillvalue)
        return format_bathy(bathy, **kwargs)


def format_bathy(bathy, **kwargs):
    """Format a bathymetry variable"""
    kwargs = dict([(k, v) for k, v in kwargs.items() if v is not None])
    return format_var(bathy, 'bathy', **kwargs)


fmt_bathy = format_bathy


class GriddedBathy(_GriddedBathyMasked_):
    """Interface to gridded bathymetries

    :Usage:

        >>> b = GriddedBathy(var, shoreline='f')
        >>> b.plot_bathy(title='Bathymetry')
        >>> b.save('bathy.nc')
        >>> newvar = b.bathy(mask=True)
    """
    def __init__(self, var, maxvalue=0., fillvalue=0., shoreline=None):
        if var.getLongitude() is None:
            var.getAxis(1).designateLongitude()
        if var.getLatitude() is None:
            var.getAxis(1).designateLatitude()
        if var.getGrid() is None:
            var.setGrid(create_grid(var.getLongitude(), var.getLatitude()))
        self._bathy = var
        self._xres, self._yres = resol(self.get_grid())
        _GriddedBathyMasked_.__init__(self, shoreline=shoreline, fillvalue=fillvalue, maxvalue=maxvalue)
        self._lon2d, self._lat2d = meshgrid(self.get_lon(), self.get_lat())
#        self._lon2dp,self._lat2dp = meshbounds(self.get_lon(), self.get_lat())
        self._xxs = self._yys = None

    def get_lon(self):
        return self._bathy.getAxis(1)
    def get_lat(self):
        return self._bathy.getAxis(0)
    def get_grid(self):
        if self._bathy.getGrid() is None:
            self._bathy.setGrid(create_grid(self.get_lon(), self.get_lat()))
        return self._bathy.getGrid()
    def get_res(self):
        return self._xres, self._yres

    def plot(self, mask=True, id=None, long_name=None, **kwargs):
        """Plot using :func:`plot_bathy`"""
        if not mask:
            kwargs.setdefault('fillcontinents', False)
        return plot_bathy(self.bathy(mask, id=id, long_name=long_name), **kwargs)

    def regrid(self, grid, method='auto', mask=True, id=None, long_name=None, **kwargs):
        """Regrid bathy to another grid"""
        return regrid2d(self.bathy(mask=mask, id=id, long_name=long_name),
                        grid, method=method, **kwargs)

    def bathy(self, mask=True, **kwargs):
        """Get the bathymetry variable

        :Params:

            - *mask*: If True and if a shoreline is set
              (with :meth:`set_shoreline`), apply a mask due to this shoreline.
        """
        if mask:
            return self.masked_bathy(**kwargs)
        return fmt_bathy(self._bathy.clone(), **kwargs)


    def save(self, ncfile, mask=True, **kwargs):
        """Save bathy to a netcdf file

        :Params:

            - **ncfile**: Output netcdf file name
            - Other keywords are passed to :meth:`bathy`
        """
        bathy = self.bathy(mask, **kwargs)
        f = cdms2.open(ncfile)
        f.write(bathy)
        f.close()

class GriddedBathyMerger(GriddedMerger, GriddedBathy):
    """Merger of gridded variables

    :Usage:

    >>> merger = GriddedMerger(mygrid, long_name='My bathy')
    >>> merger += etopo2(lon=(-10,0),lat=(42,50))
    >>> merger += my_bathy
    >>> merger.set_shoreline('i')
    >>> merged_bathy = merger.bathy()
    """
    def __init__(self, grid, shoreline=None, id=None, long_name=None, maxvalue=0., fillvalue=0.):

        GriddedMerger.__init__(self, grid, id=id, long_name=long_name, units='m')
        _GriddedBathyMasked_.__init__(self, shoreline=shoreline, maxvalue=maxvalue, fillvalue=fillvalue)
        self._bathy = self._masked_bathy = None

    def _load_(self, var, method, mask=False):
        if not cdms2.isVariable(var) and hasattr(var, 'bathy') and callable(var.bathy):
            var = var.bathy(mask)
        return GriddedMerger._load_(self, var, method)

    def merge(self, res_ratio=.5, mask=True, id=None, long_name=None):
        """Merge all variables onto the grid and apply final mask"""
        # Merge
        self._bathy = GriddedMerger.merge(self, res_ratio=res_ratio)

        # Mask
        masked_bathy = GriddedBathy.bathy(self, mask=mask, id=id, long_name=long_name)

        # Finalize
        del self._bathy, self._masked_bathy
        self._bathy = self._masked_bathy = None
        return masked_bathy

    def bathy(self, *args, **kwargs):
        """Shortcut to :meth:`merge`"""
        return self.merge(*args, **kwargs)

    def plot(self, **kwargs):
        """Plot merged bathy

        .. seealso:: :meth:`GriddedBathy.plot`
        """
        if self._bathy is None:
            self.merge()
        return GriddedBathy.plot(self, **kwargs)


class NcGriddedBathyError(Exception):
    pass
class NcGriddedBathy(GriddedBathy):
    """Get a gridded bathymetry from file

    :Params:

        - **lon/lat**, optional: Longitude and latitude selection.
        - **name**, optional: It is either

            - the name of bathymetry within the list given by :func:`print_bathy_list`
            - the name of netcdf file.
        - **i/jmargin**: Add or suppress a marginal grid points to the selected area (integer).

    :Example:

        >>> b = NcGriddedBathy(lon=(-6,-4), lat=(48,49), maxvalue=None)
        >>> b.save('bathy.nc')
        >>> b.regrid(mygrid)
        >>> b.plot(savefig='bathy.png')
        >>> b2 = NcGriddedBathy(name='bathy.nc')

    """


    def __init__(self, lon=None, lat=None, name=None, cfgfile=None, reverse=False,
        varname=None, lonname=None, latname=None, shoreline=None, maxvalue=0., fillvalue=0.,
        imargin=0, jmargin=0, **kwargs):
        # Source
        if name and os.path.exists(name): # direct from netcf file

            self.filename = name
            self.bathyname = None
            self.varname = varname
            self.lonname = lonname
            self.latname = latname

        else: # from config file

            bathyname = kwargs.get('bathyname', name)
            bathies = bathy_list(cfgfile=cfgfile)
            if bathyname is None:
                if 'etopo2' in bathies:
                    bathyname = 'etopo2'
                else:
                    kk = list(bathies.keys())
#                    if not len(kk):
#                        raise NcGriddedBathyError('No valid bathymetry found')
                    bathyname = kk[0]
            self.filename = _check_bathy_file_(bathyname)
            mybathy = bathies[bathyname]
            self.lonname = lonname or mybathy.get('lon', None)
            self.latname = latname or mybathy.get('lat', None)
            self.varname = varname or mybathy.get('var', None)
            self.bathyname = bathyname

        # Read
        try:
            f = cdms2.open(self.filename)
        except:
            raise NcGriddedBathyError('Bathymetry file not found: %s'%self.filename)
        # - var name
        if self.varname is None:
            for varname in f.listvariables():
                if varname.startswith('bounds_'): continue
                if len(f[varname].shape)==2:
                    self.varname = varname
                    break
            else:
                raise NcGriddedBathyError('Cannot find a valid 2D variable in bathy file: '+self.filename)
        elif self.varname not in f.listvariables():
            raise NcGriddedBathyError('Variable not in bathy file: '+self.varname)
        elif len(f[self.varname].shape)!=2:
            raise NcGriddedBathyError('Bathy variable not 2D: '+self.varname)
        # - dim names
        if self.lonname is None:
            self.lonname = f[self.varname].getAxis(1).id
        elif self.lonname in f.listdimension():
            if f[self.varname].getAxis(1).id != self.lonname:
                print('Bad longitude name "%s" for variable "%s". Skipping...'%(self.lonname, self.varname))
            else:
                cdms2.axis.longitude_aliases.append(self.lonname)
        if self.latname is None:
            self.latname = f[self.varname].getAxis(0).id
        elif self.latname in f.listdimension():
            if f[self.varname].getAxis(0).id != self.latname:
                print('Bad latitude name "%s" for variable "%s". Skipping...'%(self.latname, self.varname))
            else:
                cdms2.axis.latitude_aliases.append(self.latname)
        # - select
        f[self.varname].getAxis(0).designateLatitude()
        f[self.varname].getAxis(1).designateLongitude()
        select = {}
        if lon: select['lon'] = lon
        if lat: select['lat'] = lat
        # - read
        if 'lon' in select and isinstance(imargin, int) and imargin:
            i, j, k = f[self.varname].getLongitude().mapIntervalExt(select['lon'])
            select['lon'] = slice(max(i-k*imargin, 0), j+k*imargin, k)
        if 'lat' in select and isinstance(jmargin, int) and jmargin:
            i, j, k = f[self.varname].getLatitude().mapIntervalExt(select['lat'])
            select['lat'] = slice(max(i-k*jmargin, 0), j+k*jmargin, k)
        bathy = f(self.varname, **select)
        f.close()
        if self.bathyname is not None:
            bathy.long_name = '%s bathymetry'%self.bathyname.upper()

        # Reverse
        if reverse:
            bathy[:] *= -1

        # Finish initialisation
        GriddedBathy.__init__(self, bathy, shoreline=shoreline, maxvalue=maxvalue, fillvalue=fillvalue)




