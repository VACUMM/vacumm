# -*- coding: utf8 -*-
"""
Information about seaports and stations

StationInfo : A class to retrieve alot of information about seaports and stations.
MeanSeaLevel : A class dedicated to mean sea level along french coast (more complete than StationInfo).
"""
# Copyright or © or Copr. Actimar/IFREMER (2011-2015)
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

import os
#dirname = os.path.dirname(__file__)
#cfgfiles = dict(
#    default = os.path.join(dirname, 'config.cfg'),
#    site = os.path.join(dirname, 'site.cfg'),
#)
#def _get_option_(key,  **kwargs):
#    f = ConfigParser()
#    f.read([cfgfiles['default'], cfgfiles['site']])
#    if f.has_option(__name__, key):
#        return f.get(__name__, key)
from vacumm.config import get_config_value

from _geoslib import Point
import numpy as N

from vacumm.bathy.bathy import XYZBathy
from ConfigParser import ConfigParser
class StationInfoError(Exception):
    pass
class _StationPlot_(object):
    def plot(self,lon=None,lat=None,loc='upper left',
             fontsize=11.,alone=False,color='blue', show=True):
        """Show the current station on a map

        - *lon/lat*: Longitude and latitude plot ranges
        - *fontsize*: Font size of the station name [default:15.]
        - *color*: Color of name and symbol [default: 'blue']
        - *alone*: Do not plot other stations around [default: False]
        - Other keywords are passed to :func:`~vacumm.misc.plot.map`.
        """
        if self.nom is None:
            print "Aucune station n'est actuellement definie"
            return None

        import pylab as P
        from vacumm.misc.color import bistre
        from vacumm.misc import auto_scale
        from vacumm.plot import map
        A = N.array

        dlon = 2.8
        dlat = 2.
        if lon is None:
            lon = A([-dlon/2.,dlon/2.])+self.longitude
        if lat is None:
            lat = A([-dlat/2.,dlat/2.])+self.latitude
        mlon, mlat = m([self.longitude,],[self.latitude,])

        # Map
        kwargs.setdefault('proj', 'lcc')
        m = m(lon=lon, lat=lat, title='Station information', show=False, **kwargs)

        # Plot other stations
        if isinstance(self, StationInfo) and not alone:
            for station in self._stations:
                if station['nom'] != self.nom and \
                       station['longitude'] > min(lon) and \
                       station['longitude'] < max(lon) and \
                       station['latitude'] > min(lat) and \
                       station['latitude'] < max(lat):
                    long, lati = m([station['longitude'],],[station['latitude']])
                    m.plot(long,lati,'o',color='k',   label = '_nolegend_')
                    P.text(long,lati,station['nom']+'\n',va='bottom',ha='center',fontsize=fontsize*0.6)

        # Plot our station and add info
        long, lati = m([self.longitude,],[self.latitude,])
        ll = m.plot(mlon,mlat,'o',color=color,markersize=10.)
        P.text(mlon, mlat,self.nom+'\n',
             va='bottom',ha='center',fontsize=fontsize,color=color)
#       atts = ['longitude','latitude']
#       atts.extend(self._ids)
#       atts.extend(self._numerics)
        mylabel = 'Details of the station:\n'
        for att in self.attributes:
            if att not in ['zero_hydro']:
                val = getattr(self,att)
                if att in ['longitude','latitude']:
                    val = '%5.2f' % val
                elif att in self._ids:
                    att = 'ID '+att.upper()
                else:
                    att = att.upper()
                    val = '%g m' % val
                if val is not None:
                    mylabel += '   %s: %s\n' % (att,str(val))
        P.legend((ll,),(mylabel[:-1],),loc=loc,shadow=True,numpoints=1)

        if savefig:
            P.savefig(savefig)
        elif savefigs:
            savefigs(savefigs)
        if show:
            P.show()

        def show(self, *args, **kwargs):
            """Shortcut to :meth:`plot`"""
            return self.plot(*args, **kwargs)




class StationInfo(_StationPlot_):
    """Finding information about tidal stations

      This class helps finding information on
      ports and tide jauge measurement stations.
      It attempts to search for a station using several
      possible criteria, and returns attributes such as
      ids, positions, mean sea level, etc (listed using
      method :meth:`attributes`).
      See help on :meth:`search` method.

    :Params:

        - **nom**: name or regexp to search for.
        - **regexp**:
        - **file**: File in wich to search for info
        - **verbose**: Verbose mode
        - All other keywords are passed to :func:`search` method after initialization.
        """
    def __init__(self,nom=None, regexp=True, verbose=True,
                 file=None, *args, **kwargs):
        self.nom = None

        # File name
        if file is None:
#            file = _get_option_('ports_file')
            file = get_config_value(__name__, 'ports_file', ispath=True)
        if file is not None:
#            if not os.path.isabs(file):
#                file = os.path.abspath(os.path.join(dirname,  '../../../../data', file))
            if not os.path.exists(file):
                raise StationInfoError('File of ports not found: '+file)
        else:
            raise StationInfoError('No valid file of ports found')

        # Chargement du fichier
        self.set_file(file)

        # Send optional argument to search()
        if nom is not None or len(kwargs) or len(args):
            self.find(nom=nom,regexp=regexp,verbose=verbose,*args,**kwargs)


    def set_file(self,file,**kwargs):
        """ Set the current file and load it """
        import os
        self._loaded = False
        if not os.path.exists(file):
            print 'Fichier de station introuvable : '+file
        self._file = file
        self._load_file(**kwargs)


    def _load_file(self,**kwargs):

        # Open information file
        if kwargs.has_key('verbose'):
            if verbose:
                print 'Lecture de ',self._file
        f = open(self._file)

        # Parse general attributes
        self._headers = {}
        for line in f:
            sline = line[:-1].split(':')
            attribute = sline[0].strip().lower()
            definition = ':'.join(sline[1:]).strip()
            self._headers[attribute] = definition
            if attribute == 'zone':
                break

        # Parse ids and numerics
        self._ids = []
        self._numerics = {}
        for line in f:
            sline = line[:-1].split(':')
            attribute = sline[0].strip().lower()
            definition = ':'.join(sline[1:]).strip()
            if attribute == 'url_zone' :
                break
            else:
                setattr(self,attribute,None)
                if not definition.find('Identif'):
                    self._ids.append(attribute)
                else:
                    self._numerics[attribute] = definition

        # Parse zones
        self._zones = {}
        for line in f:
            sline = line[:-1].split(':')
            zone = sline[0].strip()[4:]
            description = ':'.join(sline[1:]).strip()
            if description != 'Gironde':
                self._zones[zone] = definition % description
            else:
                self._zones[zone] = description
                break

        # Header line
        line = f.next()
        attributes = line[27:-1].lower().split()

        # Loop on stations
        self._stations = []
        for line in f:

            # Name
            dic = {}
            dic['nom'] = dic['name'] = line[0:27].strip()

            # Create a list of values following each attribute
            sline = line[27:-1].split()
            ilon = len(self._ids)
            for i in 0,1:
                sline[ilon+i] = ' '.join(sline[ilon+i:ilon+i+2])
                del sline[ilon+i+1]
            if len(attributes) != len(sline):
                sline[len(attributes)-1] = ' '.join(sline[len(attributes)-1:])
                del sline[len(attributes):]

            # Ids
            for iatt in xrange(len(attributes)):
                value = sline[iatt]
                # Coordinates
                if attributes[iatt] in ['longitude','latitude']:
                    value = self._str2deg(value)
                # Undefined values
                elif value == 'None':
                    value = None
                # Url of zone
                elif attributes[iatt] == 'zone':
                    value = self._zones[value.lower()]
                # Zero hydro (!= None)
                elif attributes[iatt] == 'zero_hydro':
                    value = float(value.replace(',','.').split()[0])
                # Centimeters to meters
                elif len(value) < 5 and value.isdigit():
                    value = float(value) * 0.01
                # Store in a dictionary
                dic[attributes[iatt]] = value

            # Append dictionary to the station list
            self._stations.append(dic)
        f.close()
        self._loaded = True
        self._atts = attributes
    def attributes(self):
        """Return a list of attributes"""
        return self._atts


    def file(self):
        """ Return the file path where information are stored """
        return self._file

    def is_set(self):
        """ Return if a station is set """
        return self._loaded and self.nom is not None


    def _search(self,nom=None,regexp=True,nmax=5,**kwargs):
        """ Generic search within stations """

        stations = []
        if nom is not None:
            # Search within names using regular expression of strict comparisons
            nom = nom.strip()
            if regexp:
                # Regular expression
                import re
                this_re = re.compile(nom, re.I)
                for station in self._stations:
                    if this_re.search(station['nom']) is not None:
                        stations.append(station)
                        if len(stations) == nmax:
                            break
            else:
                # Strict equality
                for station in self._stations:
                    if station['nom'].lower() == nom.lower():
                        stations.append(station)
                        if len(stations) == nmax:
                            break

        else:
            # Search using other specific arguments
            for key,val in kwargs.items():

                if key in self._ids:
                    # Use ids
                    for station in self._stations:
                        if station in stations:
                            continue
                        if station[key] is not None and \
                               station[key].lower() == val.lower():
                            stations.append(station)
                            break

                elif key == 'position' and len(val) == 2:
                    # Find the closest station
                    import vacumm.misc as M,math
                    distances = []
                    x0 = deg2m(*val)
                    y0 = deg2m(val[1])
                    for station in self._stations:
                        if station in stations:
                            continue
                        x = deg2m(station['longitude'], station['latitude'])
                        y = deg2m(station['latitude'])
                        distances.append(math.sqrt((x-x0)**2+(y-y0)**2))
                    stations.append(self._stations[N.argmin(distances)])

                if len(stations) == nmax:
                    break

        if len(stations):
            if nmax == 1:
                return stations[0]
            else:
                return stations
        else:
            return None

    def search(self,nom=None,regexp=True,nmax=5,**kwargs):
        """Search for stations using several possible criteria and display results.

        - *nom*: Search within station names using this pattern
        - *regexp*: If True, use regular expression for search within names [default: True]
        - *position*: A two-element iterable ([lon,lat]) for searching for the closest station to this positions
        - *nmax*: Maximal number of stations displayed [default: 5]
        """
        stations = self._search(nom=nom,regexp=regexp,nmax=nmax,**kwargs)

        if stations is None:
            print 'Aucune station trouvee verifiant vos criteres'
        else:
            print 'Liste des stations trouvees (max %i) :' % nmax
            for station in stations:
                self._print_one(station)
                print 'Definition des termes accessible avec definitions()'


    def find(self,nom=None,regexp=True,verbose=False,*args,**kwargs):
        """Search for a station using several possible criteria and load it.

        - *nom*: Search within station names using this pattern
        - *regexp*: If True, use regular expression for search within names [default: True]
        - *position*: A two-element iterable ([lon,lat]) for searching for the closest station to this positions
        - *verbose*: Verbose mode [default: True]
        """
        station = self._search(nom=nom,regexp=regexp,nmax=1,*args,**kwargs)
        if station is None:
            if verbose:
                print 'Aucune station trouvee verifiant vos criteres'
        else:
            if verbose:
                print 'Chargement de la station suivante :'
                self._print_one(station)
                print 'Definition des termes accessible avec definitions()'
            self._set(station)

        return Station(station)


    def _str2deg(self,str):
        sstr = str.split()
        deg = float(sstr[0])
        deg += float(sstr[1][0:2])/60.
        deg += float(sstr[1][3:5])/3600.
        if sstr[1][-1].lower() in ['s','w']:
            deg = -deg
        return deg


    def _print_one(self,station):
        """ Print info for one station"""
        import types,vacumm.misc as M
        # Name and position
        self._print_one_arg('Nom',station['nom'])
        self._print_one_arg('Position','%s / %s' % \
          (M.lonlab(station['longitude']),M.latlab(station['latitude'])))
        self._print_one_arg('Zone',station['zone'])
        # Ids
        for id in self._ids:
            self._print_one_arg(id.upper(),station[id])
        # Properties
        keys = self._numerics.keys()
        keys.sort()
        for att in keys:
            if station[att] is not None:
                if type(station[att]) is types.StringType:
                    fmt = '%s'
                else:
                    fmt = '%g'
                self._print_one_arg(att.upper(),fmt % station[att])


    def _print_one_arg(self,att,val):
        if val is not None:
            print '  '+att.ljust(10)+' : '+val.encode('utf8')


    def _set(self,station):
        """ Set a station """
        # Internal data
        for att,val in station.items():
            setattr(self,att,val)


    def info(self):
        """ Show all available information about the current station """
        if not self.is_set():
            print "Aucune station n'est actuellement definie"
        else:
            print 'Station actuelle :'
            self._print_one(self.get_dict())


    def get_dict(self):
        """ Get the current station as a dictionnary """
        if self.nom is None:
            print "Aucune station n'est actuellement definie"
            return None
        station = {}
        for att in self._headers.keys():
            station[att] = getattr(self,att)
        for id in self._ids:
            station[id] = getattr(self,id)
        for att in self._numerics.keys():
            station[att] = getattr(self,att)
        return station


    def definitions(self):
        """ Print out the definition of all terms """
        print 'Definition des termes :'
        headers = ['nom','longitude','latitude','zone']
        for att in headers:
            self._print_one_arg(att,self._headers[att])
        for id in self._ids:
            self._print_one_arg(id.upper(),id.upper())
        for att in self._numerics.keys():
            self._print_one_arg(att.upper(),self._numerics[att])

class Station(dict, _StationPlot_):
    """A station with its info as a result from :class:`StationInfo`

    :Example:

    >>> from vacumm import StationInfo
    >>> station = StationInfo().find('Brest', nmax=1)
    >>> print station.longitude
    """

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self._atts = self.keys()
        for key, value in self.items():
            setattr(self, key, value)





#class _MeanSeaLevelData(object):
#    """Class for convenient sealevel data use"""
#    def __init__(self, data):
#        self._data = data
#    def __getattr__(self, att):
#        if att == 'data': return self._data
#        if att == 'names': return [d[0] for d in self._data]
#        if att in ['x', 'lon']:
#            return N.asarray([nxyz[1] for nxyz in self._data])
#        if att in ['y', 'lat']:
#            return N.asarray([nxyz[2] for nxyz in self._data])
#        if att in ['z', 'sealevel']:
#            return N.asarray([nxyz[3] for nxyz in self._data])
#        if att in ['xy', 'coords']:
#            return N.asarray([nxyz[1:3] for nxyz in self._data]).transpose()
#        if att in ['xyz', 'block']:
#            return N.asarray([nxyz[1:] for nxyz in self._data]).transpose()
#        return self.__dict__[att]
#    def __len__(self):
#        return len(self._data)
#
#
#class _MeanSeaLevel(_MeanSeaLevelData):
#    """Retrieve mean sea level at stations
#
#    Information can be retrieved thanks to the following properties:
#        .data: List of internal data, each element being [name, x, y, z]
#        .names: List of seaport and station names
#        .x or .lon: longitudes array (n,)
#        .y or .lat: latitudes array (n,)
#        .xy or coords: array of coordinates (2,n)
#        .z or .sealevel: mean sea level array (n,)
#        .xyz or .block: array combining coordinates and sea levels (3,n)
#
#    :Example:
#
#    >>> from vacumm.tide.station_info import MeanSeaLevel
#    >>> msl = MeanSeaLevel()
#    >>> for i in xrange(len(msl)):
#    >>>     print msl.names[i], msl.x[i], msl.y[i]
#    >>> print len(msl.clip((-4,44.,0,48)))
#    """
#    def __init__(self,  file=None):
#
#        data = []
#        f = open(meansealevel_file)
#        for line in f:
#            name = unicode(line[:29].strip(), 'utf8')
#            y, x, z = [float(v) for v in line[29:].split()]
#            data.append([name, x, y, z])
#        f.close()
#        _MeanSeaLevelData.__init__(self, data)
#    def clip(self, zone=None):
#        """Get data within specific bounds
#
#        - *zone*: A tuple of (xmin,ymin,xmax,ymax), a Polygon instance or an array to build a polygon (see vacumm.misc.grid.masking.polygons)
#        """
#        if zone is None: return self
#        zone = polygons([zone])[0]
#        mydata =  [d for d in self._data if Point((d[1], d[2])).within(zone)]
#        return _MeanSeaLevelData(mydata)
#
class MeanSeaLevelError(Exception):
    pass
class MeanSeaLevel(XYZBathy):
    """Retrieve mean sea level at stations

    Information can be retrieved thanks to the following properties:
        .data: List of internal data, each element being [name, x, y, z]
        .names: List of seaport and station names
        .x or .lon: longitudes array (n,)
        .y or .lat: latitudes array (n,)
        .xy or coords: array of coordinates (2,n)
        .z or .sealevel: mean sea level array (n,)
        .xyz or .block: array combining coordinates and sea levels (3,n)

    :Example:

    >>> from vacumm.tide.station_info import MeanSeaLevel
    >>> msl = MeanSeaLevel()
    >>> for i in xrange(len(msl)):
    >>>     print msl.names[i], msl.x[i], msl.y[i]
    >>> print len(msl.clip((-4,44.,0,48)))

    .. seealso::

        :class:`~vacumm.misc.io.XYZ`
    """
    def __init__(self, msl=None, names=None, **kwargs):
        if msl is None:
            msl = get_config_value(__name__, 'meansealevel_file')
#            msl = _get_option_('meansealevel_file')
            if msl is None:
                raise MeanSeaLevel("Can't find a default file of mean sea levels")
        if isinstance(msl, str): # file
            if not os.path.isabs(msl):
                msl = os.path.abspath(os.path.join(dirname,  '../../../../data', msl))
            if not os.path.exists(msl):
                raise MeanSeaLevel("File of mean sea levels not found: %s"%msl)
            data = []
            names = []
            xx = [] ; yy = [] ; zz = []
            f = open(msl)
            for line in f:
                name = unicode(line[:29].strip(), 'utf8')
                y, x, z = [float(v) for v in line[29:].split()]
                names.append(name)
                xx.append(x)
                yy.append(y)
                zz.append(z)
            f.close()
            msl = (xx, yy, zz)
        elif isinstance(msl, MeanSeaLevel): # instance
            names = list(msl.names())

        kwargs.setdefault('units', 'm')
        kwargs.setdefault('long_name', 'Mean sea level')
        XYZBathy.__init__(self, msl, **kwargs)
        self._names = names

    def get_names(self, mask=True):
        """Get :attr:`names`"""
        if self._names is None: return
        return self._filter_(self._names, mask)
    names = property(get_names, doc='Get name of valid stations')

    def clip(self, zone=None, inverse=False, **kwargs):
        """Geographical selection of part of the data

        - **zone**: (xmin,ymin,xmax,ymax) or a complex polygon (see :func:`~vacumm.misc.grid.masking.polygons`).
        - *inverse*: Inverse the selection.
        """
        if zone is None:
            return self
        mask = self._clip_mask_(zone, inverse)
        names=  self.names
        names = [names[i] for i, m in enumerate(mask) if m]
        return MeanSeaLevel(XYZ.clip(self, zone=mask, mask=True).xyz, names=names)



from vacumm.misc.grid.masking import polygons
from vacumm.misc.phys.units import deg2m
