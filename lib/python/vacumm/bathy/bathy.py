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
import cdms2,MV2
from genutil import minmax
import numpy as N
MV = MV2
cdms = cdms2

from vacumm.misc.axml.base import Base as XMLBase, def_setget
from vacumm.misc.io import XYZ, XYZMerger

from string import Template
from ConfigParser import ConfigParser, SafeConfigParser
import os.path, operator
from matplotlib.ticker import Formatter
import pylab as P
from warnings import warn
import pickle

from vacumm.config import check_data_file, get_config_sections, get_config_value

__all__ = ['print_bathy_list', 'Bathy', 'plot_bathy', 'fmt_bathy',  'bathy_list',
    'XYZBathy', 'XYZBathyBankClient', 'XYZBathyBank', 'XYZBathyMerger', 'GriddedBathy',
    'GriddedBathyMerger', 'NcGriddedBathy', 'NcGriddedBathyError',
#    'get_cfgfiles_gridded'
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
    sections, cfg = get_config_sections(cfg=cfgfile, parent_section=__name__,
        parent_option='cfgfile_gridded', getcfg=True)

#    if not cfgfile: cfgfile = get_cfgfiles_gridded() # File name
#    dvars = get_dir_dict(__names__) # Various directory names
#    bathys = {}
#    f = ConfigParser()
#    f.read(cfgfile)


    # Loop on section
    bathies = {}
    for section in sections:
        status = _check_bathy_file_(section, avail=True)
#        if not os.path.exists(f.get(section, 'file', vars=dvars)): continue
        if status<=minstatus: continue
        bathies[section] = {'status':status}
        for key in 'file','file_url','file_license','var','lon','lat':
            if cfg.has_option(section, key):
                bathies[section][key] = get_config_value(section, key, cfg=cfg)

    # Found any bathy?
    if name is not None:
        if not len(bathies): raise NcGriddedBathyError('No bathy are available')
        if not name in bathies:
            raise NcGriddedBathyError('Wrong bathymetry name: "%s". Please choose one of these: %s.'%
                (name, ', '.join(bathies.keys())))
        return bathies[name]
    return bathies

def _check_bathy_file_(bathyname, **kwargs):
    """Make sure to have a bathy file (see :func:`~vacumm.config.check_data_file`)"""
    return check_data_file(bathyname, 'file', parent_section=__name__,
        parent_option='cfgfile_gridded', **kwargs)


#def get_cfgfiles_gridded():
#    """Name of the configuration file containing the list of netcdf gridded bathy files"""
#    return get_config_value(__name__, 'cfgfile_gridded')

def print_bathy_list(cfgfile=None):
    """Print available bathymetries"""
#    if not cfgfile: cfgfile = get_cfgfiles_gridded()
    print 'List of bathymetries:'
    for bathyname,vv in bathy_list(cfgfile=cfgfile).items():
        print '\n[%s]' % bathyname
        print '  status: '+['Not available (need user input)',
            'Downloadable', 'On disk'][vv['status']]
        for key in 'file','file_url','file_license','var','lon','lat':
            if key in vv:
                print '  - %s: %s' % (key,vv[key])
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
        cdms.axis.longitude_aliases.append(self._lon_name)
        cdms.axis.latitude_aliases.append(self._lat_name)
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
        kwargs = {self._lon_name:self._lon_range,self._lat_name:self._lat_range}

        f = cdms.open(self.file_name)
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
        if not self.bathy.attributes.has_key('units'):
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
        self._lon2dp,self._lat2dp = axes2d(lon_axis.getValue(),lat_axis.getValue(),bounds=True,numeric=True)
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
        f = cdms.open(file,'w')
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
        for i in xrange(len(bathy1d)):
            f.write((fmt+'\n') % (self._lon2d.flat[i],self._lat2d.flat[i],bathy1d[i]))
        f.close()
        del bathy1d
        print 'Bathy wrote to ascii file',file


def plot_bathy(bathy, shadow=True, contour=True, shadow_stretch=1., shadow_shapiro=False,
    show=True, shadow_alpha=1., shadow_black=.3, white_deep=False, nmax=30,m=None, alpha=1.,
    zmin=None, zmax=None, **kwargs):
    """Plot a bathymetry

    - *lon*: Longitude range.
    - *lat*: Latitude range.
    - *show*:Display the figure [default: True]
    - *pcolor*: Use pcolor instead of contour [default: False]
    - *contour*: Add line contours [default: True]
    - *shadow*:Plot south-west shadows instead of filled contours.
    - *nmax*: Max number of levels for contours [default: 30]
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
        lon2d,lat2d = meshgrid(get_axis(bathy, -1).getValue(),get_axis(bathy, -2).getValue())

    # Masking
    if 'maxdep' in kwargs:
        zmin = -maxdep
    if 'maxalt' in kwargs:
        zmax = maxalt
    if zmin is not None:
        bathy[:] = MV2.masked_less(bathy, zmin)
    if zmax is not None:
        bathy[:] = MV2.masked_greater(bathy, zmax)

    # Default arguments for map
    if hasattr(bathy, 'long_name'):
        kwargs.setdefault('title',bathy.long_name)
    if not kwargs.has_key('cmap'):
        vmin, vmax = minmax(bathy)
#        print 'cmap topo', vmin, vmax
        kwargs['cmap'] = auto_cmap_topo((kwargs.get('vmin', vmin), kwargs.get('vmax', vmax)))
#       kwargs.setdefault('ticklabel_size','smaller')
    kwargs.setdefault('clabel_fontsize', 8)
    kwargs.setdefault('clabel_alpha',.7*alpha)
    kwargs.setdefault('clabel_glow_alpha', kwargs['clabel_alpha'])
    kwargs.setdefault('fill', 'contourf')
    kwargs['nmax'] = nmax
    kwargs['show'] = False
    kwargs['contour'] = contour
    if shadow: kwargs.setdefault('alpha',.5*alpha)
    kwargs.setdefault('projection', 'merc')
    kwargs.setdefault('fmt', BathyFormatter())
    kwargs.setdefault('colorbar_format', BathyFormatter())
    kwargs.setdefault('units', False)
    kwargs.setdefault('levels_mode','normal')
    kwargs.setdefault('bgcolor', '0.8')
    kwargs.setdefault('contour_linestyle', '-')
    savefig = kwargs.pop('savefig', None)
    kwsavefig = kwfilter(kwargs, 'savefig_')


    # White contour when dark
    if contour and white_deep:
        levels = auto_scale(bathy,nmax=nmax)
        colors = []
        nlevel = len(levels)
        for i in xrange(nlevel):
            if i < nlevel/2:
                colors.append('w')
            else:
                colors.append('k')
        kwargs.setdefault('contour_colors',tuple(colors))

    # Call to map
    m = map2(bathy, m=m, **kwargs)

    # Add shadow
    if shadow:

        # Filter
        data = MV.array(bathy,'f',fill_value=0.)
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
        pp = m.map.pcolormesh(xxs, yys, grdn,cmap=cmap)#P.get_cmap('gist_yarg'))
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
        if self.signed and val<0:
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


    def togrid(self, grid=None, mask=None, cgrid=False,  proj=None, **kwargs):


        # Shoreline for masking
        if mask is not None and mask is not False and mask is not MV2.nomask and grid is not None:

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
            if mask=='auto' and False:
                if grid is None:
                    grid = self.get_grid(**kwfilter(kwargs, 'grid'))
                mask = get_best(grid)

            # Direct handling
            if isinstance(mask, (str, Shapes)):
                # Get mask
                mask = get_shoreline(mask, clip=clip)

        # Normal case
        return XYZ.togrid(self, grid=grid, mask=mask, cgrid=cgrid,  proj=proj, **kwargs)
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
    _atts = (('xyzfile', ''),
        ('long_name', '', None),
        ('xmin', 'float', None),
        ('xmax', 'float', None),
        ('ymin', 'float', None),
        ('ymax', 'float', None),
        ('transp', 'boolean', True),
        )
    def __init__(self, bank, id, **kwargs):
        self.bank = bank
        self.cfg = bank.cfg
        assert self.cfg.has_section(id), 'No such bathy id "%s" in bathy bank'%id
        self.id = id
        assert self.cfg.has_option(id, 'xyzfile'), 'No "xyzfile" defined for bathy "%s"'%id
        self._xyz = None
        for att, val in kwargs.items():
            setattr(self, att, val)
        self._check_()
    def _atts_(self):
        return [att[0] for att in self._atts]

    def __str__(self):
        info = '[%s]\n'%self.id
        for att in self._atts_():
            info += '%s: %s\n'%(att, getattr(self, att))
        return info

    def _get_(self, att):
        idx = self._atts_().index(att)
        get = getattr(self.cfg,'get'+self._atts[idx][1])
        if not self.cfg.has_option(self.id, att):
            return self._atts[idx][2]
        return get(self.id, att)
    def _set_(self, att, val):
        self.cfg.set(self.id, att, str(val))

    def __setattr__(self, att, val):
        if att in self._atts_():
            self._set_(att, val)
            if att == 'xyzfile':
                del self._xyz
                self._xyz = None
                for att in 'xmin', 'xmax', 'ymin', 'ymax':
                    self.cfg.remove_option(self.id, att)
            self._check_(att=att)
        else:
            if att == 'id' and self.__dict__.has_key('id'):
                # Rename section
                assert not self.cfg.has_section(att), 'Bathy id "%s" already exists!'%att
                self.cfg.add_section(val)
                for opt in self.cfg.options(self.id):
                    self.cfg.set(val, opt, self.cfg.get(self.id, opt))
                self.cfg.remove_section(self.id)
                self.bank._clients[val] = self
                self.bank._clients.pop(self.id)
            elif id == 'cfg':
                raise AttributeError, "You can't change the cfg file"
            self.__dict__[att] = val
        self.save()
    def __getattr__(self, att):
        if att in self._atts_():
            return self._get_(att)
        else:
            return self.__dict__[att]

    def _check_(self, att=None):
        """Check that infos are consistent"""
        if not os.path.exists(self.xyzfile):
            warn('%s xyz bathy file does not exist'%self.xyzfile)
        else:
            # Fill what is missing
            if att is not None:
                atts = [att]
            else:
                atts = self._atts_()[1:]
            for att in atts:
                if self._get_(att) is not None: continue
                if self._xyz is not None:
                    xyz = self._xyz
                else:
                    xyz = self._xyz = XYZBathy(self.xyzfile)
                if att.startswith('x') or att.startswith('y'):
                    self._set_(att, getattr(xyz, att))
#               elif xyz.long_name is not None:
#                   self.long_name = xyz.long_name
#               elif xyz.units is not None:
#                   self.units = xyz.units
        self.save()

    def save(self):
        """Save the bank"""
        self.bank.save()

    def load(self, zone=None, margin=None):
        """Load the bathymetry as a :class:`XYZBathy` instance

        .. seealso::

            :meth:`~vacumm.misc.io.XYZ.clip`
        """
        if not os.path.exists(self.xyzfile):
            warn('%s xyz bathy file does not exist'%self.xyzfile)
        else:
            xyz = XYZBathy(self.xyzfile, long_name=self.long_name, transp=self.transp, id=self.id)
            if zone is not None and zone != (None, )*4:
                xyz = xyz.clip(zone=zone, margin=margin)
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
    def __init__(self, cfgfile=None):

        # Guess file name
        if cfgfile is None:
            f = ConfigParser()
            f.read(get_config_value(__name__, 'cfgfile_xyz'))
            if f.has_option('banks', 'main'):
                cfgfile = f.get('banks', 'main')
            #else:
            #    raise Exception, 'Cannot read config file to guess bank file name'
        # Read file
        self.cfg = SafeConfigParser()
        self.cfgfile = cfgfile
        if cfgfile and os.path.exists(cfgfile): self.cfg.read(cfgfile)
        self._clients = {}
        for id in self.cfg.sections():
            self._clients[id] = XYZBathyBankClient(self, id)
        self._order = None

    def ids(self):
        """Return the list of bathy ids"""
        return self.cfg.sections()

    def __str__(self):
        if not len(self):
            return 'No bathymetry'
        infos=[]
        for id in self.ids():
            infos.append(str(self._clients[id]))
        return '\n'.join(infos)
    def __len__(self):
        return len(self.cfg.sections())

    def add(self, xyzfile, id=None, force=False, **kwargs):
        """Add a bathy to the bank

        - **xyzfile**: xyz file.
        - *id*: Id of the bank. If ``None``, it is guessed from the file name:
           ``xyzfile="/path/toto.xyz"`` ->  ``id="toto"``
        """
        if id is None:
            id = os.path.basename(xyzfile)
            if id.endswith('.xyz'): id = id[:-4]
        if self.cfg.has_section(id):
            if force:
                self.cfg.remove_section(id)
            else:
                raise KeyError, 'Bathymetry %s already exists'%id
        self.cfg.add_section(id)
        self.cfg.set(id, 'xyzfile', xyzfile)
        self._clients[id] = XYZBathyBankClient(self, id, **kwargs)
#       self._save_()
        self._update_order_()
    def __iadd__(self, d):
        self.add(d)
        return self

    def remove(self, id):
        """Remove a bathy from the bank

        *id*: A bathy id or :class:`XYZBathyBankClient` instance.
        """
        if isinstance(id, XYZBathyBankClient):
            id = id.id
        if not self.cfg.has_section(id):
            warn('No such bathymetry %s'%id)
            return
        self.cfg.remove_section(id)
        del self._clients[id]
        self._update_order_()
#       self._save_()
    def __isub__(self, d):
        self.remove(id)
        return self

    def copy(self, id1, id2):
        """Copy bathy infos to another"""
        assert self.cfg.has_section(id1), 'No such bathymetry %s'%id1
        if not self.cfg.has_section(id2):
            self.cfg.add_section(id2)
        for opt in self.cfg.options(id1):
            self.cfg.set(id2, self.cfg.get(id1, opt))
        self._clients[id2] = XYZBathyBankClient(self, id2)
        self._update_order_()
##      self._save_()

    def __getitem__(self, id):
        assert self.cfg.has_section(id), 'No such bathymetry %s'%id
        return self._clients[id]
    def __setitem__(self, id, xyzfile):
        self.add(id, xyzfile)
    def __deltiem__(self, id):
        self.remove(id)

    def load(self, id):
        """Load a bathymetry as a :class:`XYZBathy` instance"""
        assert self.cfg.has_section(id), 'No such bathymetry '+id
        return self._clients[id].load()

    def list(self):
        """Get all clients of the bank as a list"""
        return self._clients.values()

    def select(self, xmin=None, xmax=None, ymin=None, ymax=None, load=True, margin=None, ordered=None):
        """Return a list of bathymetries relevant to the specified region

        - *margin*: Margin for loading xyz bathies:

            - if None, the whole bathy is loaded from the bank,
            - else it should be a value relative to the approximative resolution (see :meth:`XYZBathy.resol`)
        """
        res = []
        for id in self.ids():
            b = self._clients[id]
            if   (None not in (xmin, b.xmax) and xmin > b.xmax) or \
                (None not in (xmax, b.xmin) and xmax < b.xmin) or \
                (None not in (ymin, b.ymax) and ymin > b.ymax) or \
                (None not in (ymax, b.ymin) and ymax > b.ymin): continue
            if load:
                b = b.load(zone=(xmin, ymin, xmax, ymax), margin=margin)
            res.append(b)
        if self.get_order() is not None and ordered != False:
            ids = [b.id for b in res]
            order = [id for id in self.get_order() if id in ids]
            res.sort(cmp=lambda b1, b2: cmp(self.get_order().index(b1), self.get_order().index(b2)))
        return res

    def set_order(self, order):
        """Set a list of ids to define the order"""
        if order is None:
            del self.cfg.defaults()['order']
        else:
            self.cfg.defaults()['order'] = str([id for id in order if id in self.ids()])
    def get_order(self,):
        """Get a list of ids that define the order"""
        return eval(self.cfg.defaults().get('order', 'None'))

    def _update_order_(self):
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
        if d in self._bank.ids():
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
            shoreline = get_bestres((self.get_lon(), self.get_lat()))
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

    def masked_bathy(self, id=None, long_name=None):
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
            bathy[:] = bathy.filled(self._fillvalue)
        return fmt_bathy(bathy, id=id, long_name=long_name)

def fmt_bathy(bathy, id=None, long_name=None):
    """Format a bathymetry variable"""
    if id is not None: bathy.id = id
    if bathy.id.startswith('variable_'): bathy.id = 'bathy'
    if long_name is not None: bathy.long_name = long_name
    if not hasattr(bathy, 'long_name'): bathy.long_name = 'Bathymetry'
    bathy.units = 'm'
    return bathy



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
        return regrid2d(self.bathy(mask=mask, id=id, long_name=long_name), grid, method=method, **kwargs)

    def bathy(self, mask=True, id=None, long_name=None):
        """Get the bathymetry variable

        :Params:

            - *mask*: If True and if a shoreline is set (with :meth:`set_shoreline`), apply a mask due to this shoreline.
        """
        if mask:
            return self.masked_bathy(id=id, long_name=long_name)
        return fmt_bathy(self._bathy.clone(), id=id, long_name=long_name)


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

from vacumm.misc.grid.regridding import GriddedMerger
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
                    kk = bathies.keys()
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
                print 'Bad longitude name "%s" for variable "%s". Skipping...'%(self.lonname, self.varname)
            else:
                cdms.axis.longitude_aliases.append(self.lonname)
        if self.latname is None:
            self.latname = f[self.varname].getAxis(0).id
        elif self.latname in f.listdimension():
            if f[self.varname].getAxis(0).id != self.latname:
                print 'Bad latitude name "%s" for variable "%s". Skipping...'%(self.latname, self.varname)
            else:
                cdms.axis.latitude_aliases.append(self.latname)
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




######################################################################
from shorelines import get_best, get_shoreline, ShoreLine, get_bestres
from vacumm.misc.grid import get_xy, axes2d, resol, create_grid, get_grid, meshgrid, \
    meshbounds, get_axis
from vacumm.misc.grid.regridding import griddata, refine, regrid2d
from vacumm.misc.grid.masking import polygon_mask, GetLakes
from vacumm.misc.misc import kwfilter, auto_scale,geo_scale,deplab
from vacumm.misc.color import cmap_jete, cmap_bathy,bistre, cmap_custom, auto_cmap_topo
from vacumm.misc.plot import map2
from vacumm.misc.filters import norm_atan,deriv2d,shapiro2d
from vacumm.misc.io import Shapes


#b = NcGriddedBathy((9.2, 14.2), (-8.3, -3.4))
#b.plot(mask=True)
#
