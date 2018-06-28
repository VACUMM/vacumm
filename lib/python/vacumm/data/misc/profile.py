#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
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

__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__date__ = '2015-02-27'
__doc__ = '''
This module provides vertical profiles managers.

It is intented to be used with the two main goals:

- convert various profiles, eventually merged, to a standardized format
- operate on standardized profiles

The main classes are:

- :class:`Profile`: mainly used internally to store variables of one profile
- :class:`AbstractProfiles`: abstract base class used to handle a list of :class:`Profile`,
  it define what Profiles or other subclasses must implement to be used by other classes
  like :class:`ProfilesMerger`.
- :class:`Profiles`: concrete base class which implements the :class:`AbstractProfiles`:
  required methods. This class define a :meth:`Profiles.save` method which write a standardized
  NetCDF profiles file. It also define some utility methods like plotting utilities.
- :class:`ProfilesDuplicatesFilter`: a class used by default by :class:`ProfilesMerger`:
- :class:`ProfilesMerger`: split Profiles into unique / duplicate Profiles
- :class:`ProfilesDataset`:

Two classes may be used as a program entry point using their main class method:

- :func:`ProfilesMerger.main`

    - load profiles
    - save two profiles files: filtered (unique) and rejected (duplicates)

- :func:`ProfilesDataset.main`:

    - load a standardized profiles file
    - plot distribution, histogram and profile variables data
    - convert into a standardized profiles file format
'''

# ==============================================================================
# NOTES
# ==============================================================================
#
# * Some profile time value are invalid (e.g. ['1' '9' '9' '8' '0' '6' '0' '2' '9' '9' '9' '9' '0' '0']
#   in MANGA_SISMER_PR_CT.nc).
#   In this case, the time component will be set to 12H
#   Then when searching for duplicates, the first matching duplicate being in the
#   time range [T - N hours, T + N hours] will replace the time corrupted profile.
#   (with T the corrupted profile datetime and N a configurable number of hours)
# * Operations on masked values throws numpy warnings, this can be ignored with the following instruction
#     numpy.seterr(divide='ignore', over='ignore')
#
# ==============================================================================
# TODO
# ==============================================================================
#
# ~ Add quality flag filter (multiple checks: on depth, time, pos and variable data)
# x Add platform_code variable and check
# * Add threshold in lat,lon duplicate test ?
# x Replace time corrupted profile with valid duplicate if any.
# x Duplicates file must also contain the profiles of valid file
# * Depth selection
# * Check output data values (precision)
# * Check loaded depth/pressure units for seawater conversions (m/db)
# ~ Standardize outputs (names and attributes)
# * Do not convert depth at loading time, let level data abstract until usage
# * Defer in memory data loads when possible (keep FileVariable object
#   => imply keeping files opened : warn about limitation on number of opened files)
# * Track side effect of data manipulations for an api usage
# * Opt: add util methods:
#   - Profiles.get_profile(profile_slice, level_slice)
#   - Profile.get_profile(level_slice)
#   - Variables unit checks
#   - Export variable mapping functionnality to a dedicated module
# ==============================================================================

import copy, datetime, glob, optparse, os, pprint, sys, time, warnings

import cdms2, cdtime, MV2, numpy, pylab
from _geoslib import Point, Polygon


from vacumm.data.misc.dataset import Dataset, POST_PLOT_ARGS, OceanDataset
from vacumm.misc.atime import create_time, ch_units, comptime, datetime as adatetime
from vacumm.misc.atime import is_in_range as is_time_in_range, reltime, Intervals, round_date, add as add_date
from vacumm.misc.axes import create, create_time, create_dep, create_lat, create_lon
from vacumm.misc.bases import Object
from vacumm.misc.color import cmap_magic, StepsNorm
from vacumm.misc.config import ConfigManager
from vacumm.misc.exception import getDetailedExceptionInfo
from vacumm.misc.file import find, efind
from vacumm.misc.grid import meshweights
from vacumm.misc.grid.masking import masked_polygon#, polygons
from vacumm.misc.grid.regridding import interp1d
from vacumm.misc.io import ncget_var, netcdf3, netcdf4
from vacumm.misc.log import Logger
from vacumm.misc.misc import kwfilter, is_iterable
from vacumm.misc.plot import bar2, curve2, map2

try:
    from vacumm.diag.thermdyn import sw_depth, sw_pres, sw_dens
except ImportError:
    warnings.warn('Failed to import seawater fonctions, depth<=>pressure conversions not available')

# ==============================================================================


class Profile(Object):
    '''
    A profile variables container.
    Store variables in a dictionnary with variables identifiers as keys, cdms variables as values.
    Variables must share the same profile properties: platform code, datetime, position and depth.

    **NOTE: At present, this loads (copy) all data in memory ..!**

    :Required named argmuments:
            - **platform_code**: a string
            - **datetime**: datetime.datetime
            - **latitude**: float
            - **longitude**: float
            - **datetime_corrupted**: boolean indicating if the time component is corrupted (arbitrary set)
            - **depth**: cdms variable with shape (level)
            - **variables**: cdms variables dictionnary, variables with shape (level)

    '''
    __properties = ('platform_code', 'datetime', 'latitude', 'longitude', 'datetime_corrupted', 'depth', 'variables', 'parent', 'index', 'filename')

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        # TODO:
        # - Remove slot mecanism if it does not optimize memory usage
        # - Store depth as axis
        for a in self.__properties:
            setattr(self, a, kwargs.get(a, None))
        if not isinstance(self.platform_code, basestring):
            self.platform_code = ''.join(numpy.array(self.platform_code, dtype='c'))
        self.platform_code = self.platform_code.strip()
        if not hasattr(self.depth, 'shape'):
            self.depth = numpy.ma.array(self.depth).reshape((1,))
        if len(self.depth.shape) != 1:
            raise ValueError('Invalid depth shape, expecting 1D array, got %s'%(self.depth.shape,))
        self.level = cdms2.createAxis(range(len(self.depth)), id='level')
        self.depth = cdms2.createVariable(self.depth, id='depth', axes=(self.level,))
        for n,v in self.variables.items():
            if not hasattr(v, 'shape') and self.depth.shape[0] == 1:
                v = numpy.ma.array(v).reshape((1,))
            if v.shape != self.depth.shape:
                raise ValueError('Variable %s%s does not match expected shape: %s'%(n, v.shape, self.depth.shape))
            v = cdms2.createVariable(v, id=n, axes=(self.level,))
            #v = cdms2.createVariable(v, id=n)
            self.variables[n] = v
        self.__str__ = self.__newstr__

    # Replace __str__ once all __init__ actions are completed
    def __newstr__(self):
        return '<%s %s, depth shape: %s, variables: %s'%(
            self.__class__.__name__,
            ', '.join(('%s: %s'%(a, getattr(self, a, None)) for a in self.__properties[:5])), # update this on properties changes
            self.depth.shape, self.variables.keys())

    def get_depth_min(self):
        '''Return the minimum depth value.'''
        return float(MV2.min(self.depth))

    def get_depth_max(self):
        '''Return the maximum depth value.'''
        return float(MV2.max(self.depth))

    def get_depth_range(self):
        '''Return the depth range (min,max).'''
        return self.get_depth_min(), get_depth_max()

    def plot(self, variables=None, *args, **kwargs):
        '''
        Produce a 1D plot of the profile variable with depth as Y axis.
        :Params:
            - **variables**: optional list of strings used to restrict the variable to display
            - **args, kwargs**: passed to the underlying plot function
        '''
        if not variables: variables = self.variables.keys()
        if isinstance(variables, basestring): variables = (variables,)
        save = kwargs.pop('save', False)
        gobjs = []
        for i,varid in enumerate(variables):
            var = self.variables[varid]
            ckwargs = kwargs.copy()
            ckwargs.setdefault('title', '%s\n time: %s\nlatitude: %s, longitude: %s'%(var.id, self.datetime, self.latitude, self.longitude))
            ckwargs.update(dict(show=False, figure=i))
            gobjs.append(curve2(var, *args, **ckwargs))

            if save: pylab.savefig('profile_%s.png'%(varid))

        if kwargs.get('show', True) and gobjs:
            pylab.show()


# ==============================================================================


class AbstractProfiles(Object, list):
    '''
    A Profile objects collection.
    '''

    # TODO: Add NotImplementedError messages ?
    def __init__(self, *args, **kwargs):
        Object.__init__(self, **kwargs)
        list.__init__(self, *args)

    def __str__(self):
        return '<%s size: %s>'%(self.__class__.__name__, len(self))
    __repr__ = __str__

    @classmethod
    def get_type(cls, *a, **k):
        '''Must be defined by implementing subclasses.
        This method must return a string identifier of the profile type'''
        raise NotImplementedError('%s does not define a get_type method !'%(cls))
    @classmethod
    def get_variables_map(cls):
        '''Must be defined by implementing subclasses.
        This method must return a dict with variables identifiers strings as keys,
        list of variables candidates strings as values.
        '''
        raise NotImplementedError('')
    @classmethod
    def can_handle_dataset(cls, *a, **k):
        '''Must be defined by implementing subclasses.
        Return a boolean indicating if a file can be handled by the concrete subclass.'''
        raise NotImplementedError('')

    def get_default_variables(self, *a, **k):
        '''Must be defined by implementing subclasses.
        Return the list of strings indicating the default variables to be processed.'''
        raise NotImplementedError('')

    def load_variable(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_platform_code(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_filename(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_time(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_longitude(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_latitude(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_depth(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')
    def load_pressure(self, *a, **k):
        '''Must be defined by implementing subclasses'''
        raise NotImplementedError('')

    @classmethod
    def get_types(cls):
        '''Return a dict of known profiles types as keys, implementing classes as values.

        Available profiles types are based on this module AbstractProfiles implementation classes,
        with types determined by the get_type method.

        .. todo::
            - add a way to plug external profiles types (from other modules)

        '''
        types = {}
        for n,c in globals().items():
            try:
                if isinstance(c, type) and issubclass(c, cls):
                    types[c.get_type()] = c
            except NotImplementedError, e: self.warning('%s: %s', e.__class__.__name__, e)
            except Exception, e: self.warning(e)
        return types

    @classmethod
    def factory(cls, *args, **kwargs):
        '''Create a Profiles object from nothing, a profile type, a profile dataset and/or profiles specification.

        Examples:
            >>> factory(type='PROFILE_TYPE_ID')
            >>> factory(dataset='/path/to/dataset')
            >>> factory(dataset='/path/to/dataset', type='PROFILE_TYPE')
            >>> factory(dataset='/path/to/dataset', load=True)
            >>> factory(dataset='/path/to/dataset', type='PROFILE_TYPE', load=True)

        The dataset argument, if provided, is the dataset to use, either as a uri string, or CdmsFile.
        The type argument, if provided, is the profiles type identifier as returned by get_type.
        The load argument, if provided, indicate to load the provided file path.
        '''
        ptype = kwargs.pop('type', '').strip().upper()
        pdataset = kwargs.pop('dataset', '')
        safe = kwargs.pop('safe', False)
        load = kwargs.pop('load', False)
        types = cls.get_types()
        p = None
        if not ptype and pdataset:
            for t,c in types.items():
                if c.can_handle_dataset(pdataset):
                    ptype = t
                    break
        if ptype in types:
            p = types[ptype](*args, **kwargs)
        else:
            if safe: p = cls(*args, **kwargs)
            else: raise ValueError('Unsupported profiles type: "%s"'%(ptype))
        if pdataset and load: p.load(pdataset)
        return p


# ==============================================================================

class Profiles(AbstractProfiles):
    '''
    Implementation of AbstractProfiles.

    Define the load, get*, select and save features.

    This class contains a variable_map attribute which is a dict with
    variables identifiers strings as keys,  list of variables candidates
    strings as values. See load_variable for this mapping mecanism.
    Subclasses may modify this mapping attribute.

    :Params:
        - **profiles**: a list of Profile objects or Profiles instance or profiles file string.
        - **spec**: see below
        - **the spec attributes below can also be passed directly as keyword arguments**

    If the spec attribute is provided, it must be an object (usually another Profiles instance) with the following attributes:
        - **variables_map**: the variable mapping, default to :func:`get_variables_map`
        - **variables**: identifiers of variables to be loaded, default is from :func:`get_default_variables`
        - **qualities**: quality codes filter, default is from :func:`get_default_variables`
        - **time_units**: time unit for internal time axis, default is 'seconds since 1900-01-01 00:00:00'

    :Note:
        - When Profile or Profiles instances are given, their nested Profile objects are referenced, note copied.

    :Examples:
        >>> import glob
        >>> import vacumm.misc.axes as A
        >>> import vacumm.data.misc.profile as P
        >>> p = P.Profiles(
                profiles=glob.glob('argo-profiles-*.nc'),
                logger_level='debug',
                variables_map={
                    'time':('TIME',),
                    'depth':('DEPH',),
                    'latitude':('LATITUDE',),
                    'longitude':('LONGITUDE',),
                    'temperature':('TEMP',),
                    'salinity':('PSAL',),
                },
                variables=('temperature','salinity')
            )
        >>> time = p.get_time()
        >>> print 'time:', A.create_time(time, units=time.units).asComponentTime()
        >>> print 'depth:', p.describe(p.get_depth(), stats=True)
        >>> print 'latitude:', p.describe(p.get_latitude(), stats=True)
        >>> print 'longitude:', p.describe(p.get_longitude(), stats=True)
        >>> temp = p.get_variable('temperature')
        >>> print 'temperature:', p.describe(temp, stats=True)
        >>> psal = p.get_variable('salinity')
        >>> print 'salinity:', p.describe(psal, stats=True)

    .. todo::
        - fix load of single depth profile file
        - check depth quality test
    '''

    # See docstring above.
    default_variables_map = {
        'platform_code':['platform_code'],
        #'time':['time', 'station_date_time'],
        'time':['time', 'juld'],
        'time_qc':['time_qc', 'juld_qc'],
        'latitude':['latitude'],
        'longitude':['longitude'],
        'position_qc':['position_qc'],
        'depth':['deph', 'depth'],
        'depth_qc':['depth_qc', 'pres_qc'], # ok ?
        'pressure':['pressure', 'pres'],
        'temperature':['temperature', 'temp'],
        'salinity':['salinity', 'psal'],
    }

    default_time_units = 'seconds since 1900-01-01 00:00:00'

    def __init__(self, profiles=[], **kwargs):
        self.variables_map = self.get_variables_map()
        self.variables = self.get_default_variables()
        self.qualities = self.get_default_qualities()
        # Cloning
        # Note: A dedicated copy method should exists to copy the contained
        #       Profile objects (and their variables data) then add a kwarg
        #       to apply this here
        spec = kwargs.pop('spec', isinstance(profiles, Profiles) and profiles or None)
        if spec:
            for a in 'variables_map', 'variables', 'qualities', 'time_units':
                setattr(self, a, copy.copy(getattr(spec, a)))
        datasets = filter(lambda p: not isinstance(p, Profile), profiles)
        profiles = filter(lambda p: isinstance(p, Profile), profiles)
        # Apply args
        vmap = kwargs.pop('variables_map', None)
        if vmap:
            for k,v in vmap.iteritems():
                if not isinstance(v, (list,tuple)): v = (v,)
                if k in self.variables_map:
                    self.variables_map[k] = list(v) + list(self.variables_map[k])
                else:
                    self.variables_map[k] = v
        self.variables = kwargs.pop('variables', self.variables)
        self.qualities = kwargs.pop('qualities', self.qualities)
        self.timerange = kwargs.pop('timerange', [])
        self.lonrange = kwargs.pop('lonrange', [])
        self.latrange = kwargs.pop('latrange', [])
        self.time_units = kwargs.pop('time_units', self.default_time_units)
        AbstractProfiles.__init__(self, profiles, **kwargs)
        self.verbose('Initialize %s', self.__class__)
        self.verbose('Requested variables identifiers: %s', self.variables)
        self.verbose('Valid quality flags: %s', self.qualities)
        self.verbose('Time range filter: %s', self.timerange)
        self.verbose('Longitude range filter: %s', self.lonrange)
        self.verbose('Latitude range filter: %s', self.latrange)
        if datasets:
            self.load(datasets)

    def summary(self):
        '''
        Return a summary string of nested profiles.
        '''
        info = '%s\n  size: %s'%(self.__class__.__name__, len(self))
        if len(self):
            info += '\n  variables: %s\n  qualities: %s\n  time: %s\n  depth: %s\n  latitude: %s\n  longitude: %s'%(
            #numpy.unique(numpy.concatenate([p.variables.keys() for p in self]))[:],
            self.variables, self.qualities,
            self.get_time_range(), self.get_depth_range(), self.get_lat_range(), self.get_lon_range())
        return info

    def get_time_min(self):
        '''Get the minimum date of nested profiles as a datetime object'''
        if len(self): return min(p.datetime for p in self)

    def get_time_max(self):
        '''Get the maximum date of nested profiles as a datetime object'''
        if len(self): return max(p.datetime for p in self)

    def get_time_range(self):
        '''Get time min and max'''
        return self.get_time_min(), self.get_time_max()

    def get_depth_min(self):
        '''Get the minimum depth value'''
        if len(self): return min(p.get_depth_min() for p in self)

    def get_depth_max(self):
        '''Get the maxnimum depth value'''
        if len(self): return max(p.get_depth_max() for p in self)

    def get_depth_range(self):
        '''Get depth min and max'''
        return self.get_depth_min(), self.get_depth_max()

    def get_lat_min(self):
        '''Get the minimum latitude value'''
        if len(self): return min(p.latitude for p in self)

    def get_lat_max(self):
        '''Get the maximum latitude value'''
        if len(self): return max(p.latitude for p in self)

    def get_lat_range(self):
        '''Get lat min and max'''
        return self.get_lat_min(), self.get_lat_max()

    def get_lon_min(self):
        '''Get the minimum longitude value'''
        if len(self): return min(p.longitude for p in self)

    def get_lon_max(self):
        '''Get the maximum longitude value'''
        if len(self): return max(p.longitude for p in self)

    def get_lon_range(self):
        '''Get lon min and max'''
        return self.get_lon_min(), self.get_lon_max()

    @classmethod
    def get_type(cls):
        '''Return the profiles type identifier, used by :func:`get_types`'''
        return ''

    @classmethod
    def can_handle_dataset(cls, dataset):
        '''Return True if this class can handle the profiles of the given dataset.
        If dataset is a AbstractProfiles object, check that types matches.
        Else check that the type identifier is contained in the dataset file path, case insensitive (e.g. '/path/to/profiles_pr_ba_xxxx.nc' and type is 'pr_ba')
        '''
        t = cls.get_type()
        if isinstance(dataset, AbstractProfiles): return dataset.get_type() == t
        uri = isinstance(dataset, cdms2.dataset.CdmsFile) and dataset.id or dataset
        return t and t.lower() in uri.lower()

    @classmethod
    def get_variables_map(cls):
        '''Return a copy of variable_map
        '''
        return copy.deepcopy(cls.default_variables_map)

    @classmethod
    def get_default_variables(self):
        '''Return the variables identifiers this class should load by default'''
        return ()#'temperature', 'salinity')

    @classmethod
    def get_default_qualities(self):
        '''Return the quality code this class would filter by default'''
        return ()#'1', '2')

    def sort(self, cmp=None, *args, **kwargs):
        '''
        Sort profiles.
        :Params:
            - **cmp**: If a string, use a predifined sort behavior in:
                - 'time': Sort the nested profiles based on their datetime
        Other parameters passed to the original sort method
        '''
        import __builtin__
        if isinstance(cmp, basestring):
            cmp = cmp.strip().lower()
            if cmp == 'time':
                AbstractProfiles.sort(self, cmp=lambda p1, p2: __builtin__.cmp(p1.datetime, p2.datetime), *args, **kwargs)
            else:
                raise ValueError('Invalid sort method: "%s"'%cmp)
        else:
            AbstractProfiles.sort(self, cmp=cmp, *args, **kwargs)

    def load_variable(self, dataset, varid, *args, **kwargs):
        '''Load a variable

        :Params:
            - **dataset**: a CdmsFile object
            - **varid**: variable identifier
            - **as_axis**: if True, convert the loaded (1D) variable into a csdms axis
        '''
        # Variables candidates are case insensitive
        pfx, sfx = kwargs.get('prefix', ''), kwargs.get('suffix', '')
        varnames = self.variables_map.get(varid, varid)
        if isinstance(varnames, basestring): varnames = (varnames,)
        varnames = map(lambda n: '%s%s%s'%(pfx, n, sfx), varnames)
        var = None
        for varname in varnames:
            for n in dataset.variables.keys():
                if n.lower() == varname.lower():
                    var = dataset.getVariable(n)
                    break
            if not var:
                for n in dataset.axes.keys():
                    if n.lower() == varname.lower():
                        var = dataset.getAxis(n)
                        break
            else: break
        if var:
            self.info('Loaded %s%s%s from %s', pfx, varid, sfx, self.describe(var))
            if kwargs.get('as_axis', None):
                attrs = var.attributes.copy()
                var = cdms2.createAxis(var, id=var.id)
                for a,v in attrs.items(): setattr(var, a, v)
            #if kwargs.get('as_variable', None):
            else:
                var = cdms2.createVariable(var, id=var.id, attributes=var.attributes.copy())
            return var
        self.verbose('Could not load %s%s%s: No matching variable using varnames: %s', pfx, varid, sfx, varnames)
        self.debug('variable map: %s', self.variables_map)
        return None

    def load_platform_code(self, dataset, *args, **kwargs):
        '''Load and return the platform code variable (using :func:`load_variable`)'''
        v = self.load_variable(dataset, 'platform_code', *args, **kwargs)
        if v is not None: return v
        # Try getting platform code from global attributes
        for a in ('wmo_platform_code', 'platform_code'):
            if a in dataset.attributes:
                platcode = dataset.attributes[a]
                if not platcode.strip(): continue
                self.info('Loaded platform code from global attribute %s: %s', a, platcode)
                return platcode
        # Try getting platform code from filename last '_' delimited component
        if hasattr(dataset, 'id'):
            platcode = os.path.splitext(os.path.basename(dataset.id))[0].rsplit('_', 1)[-1]
            self.info('Loaded platform code from file name: %s', platcode)
            return platcode

    def load_filename(self,dataset):
        '''Load and return the filename'''
        filevar=os.path.basename(dataset.id)
        return filevar

    def load_time(self, dataset, *args, **kwargs):
        '''Load and return the time variable as an axis (using :func:`load_variable`)

        An extra keyword arguments getcflag indicates if the date corrupted flag must also be returned

        '''
        kwargs['as_axis'] = True
        timevar = self.load_variable(dataset, 'time', *args, **kwargs)
        #if timevar.dtype.kind in 'Sa':
        if timevar.typecode() in 'c':
            def parsetime(t):
                t = ''.join(t)
                d = datetime.datetime.strptime(t[:8], '%Y%m%d')
                hour, minute, second = int(t[8:10]), int(t[10:12]), int(t[12:14])
                # cf: time corruption flag
                cf = False
                if hour > 23:   cf = True; hour = 12
                if minute > 59: cf = True; minute = 0
                if second > 59: cf = True; second = 0
                ct = cdtime.comptime(d.year, d.month, d.day, hour, minute, second)
                return ct, cf
            times, corrflags = numpy.array(map(parsetime, timevar)).transpose()
            timevar = create_time(times, units=self.time_units)
        else: corrflags = None
        timevar = ch_units(timevar, self.time_units)
        if kwargs.get('getcflag', False): return timevar, corrflags
        else: return timevar

    def load_time_quality(self, dataset, *args, **kwargs):
        '''Load and return the time quality variable (using :func:`load_variable`)'''
        return self.load_variable(dataset, 'time_qc', *args, **kwargs)

    def load_longitude(self, dataset, *args, **kwargs):
        '''Load and return the longitude variable (using :func:`load_variable`)'''
        return self.load_variable(dataset, 'longitude', *args, **kwargs)

    def load_latitude(self, dataset, *args, **kwargs):
        '''Load and return the latitude variable (using :func:`load_variable`)'''
        return self.load_variable(dataset, 'latitude', *args, **kwargs)

    def load_position_quality(self, dataset, *args, **kwargs):
        '''Load and return the position quality variable (using :func:`load_variable`)'''
        return self.load_variable(dataset, 'position_qc', *args, **kwargs)

    def load_depth(self, dataset, *args, **kwargs):
        '''Load and return the depth variable (using :func:`load_variable`)

        :Params:
            - **tryconv**: if true, try to convert pressure to depth if depth cannot be loaded
            - **frompres**: force pressure to depth conversion

        '''
        frompres, tryconv = kwargs.pop('frompres', False), kwargs.pop('tryconv', True)
        if frompres:
            tryconv = True
            depth = None
        else:
            depth = self.load_variable(dataset, 'depth', *args, **kwargs)
        if depth is None and tryconv:
            if not frompres:
                self.warning('No depth found, now trying to convert pressure to depth')
            lat = self.load_latitude(dataset)
            if lat is None:
                self.error('Latitude is needed for this conversion')
                return None
            pres = self.load_pressure(dataset, tryconv=False)
            if pres is None:
                self.error('Pressure is needed for this conversion')
                return None
            self.debug('Latitude: %s', self.describe(lat))
            self.debug('Pressure: %s', self.describe(pres))
            # Get negative depth values
            # NOTE: pressure in decibar, depth in meters, latitude in decimal degrees [-90,+90]
            depth = -1. * sw_depth(pres[:].swapaxes(0,1), lat[:]).swapaxes(0,1)
            # TODO: add explicit named axes
            proax = cdms2.createAxis(xrange(depth.shape[0]), id="profile")
            depax = cdms2.createAxis(xrange(depth.shape[1]), id="depth")
            depth = cdms2.createVariable(depth, id='depth', axes=(proax, depax))
            self.verbose('Computed depth: %s', self.describe(depth))
        if depth is not None:
            # Check negative depth values
            nmod = 0
            for ip in xrange(len(depth) if len(depth.shape) == 1 else depth.shape[0]):
                d = depth if len(depth.shape) == 1 else depth[ip]
                mi, ma, av = MV2.min(d), MV2.max(d), MV2.average(d)
                if ma > 0:
                    self.verbose('Loaded depth[%s] contains positive values, expecting positive up (min: %s, max: %s, avg: %s)', ip, mi, ma, av)
                if av > 0:
                    self.verbose('Loaded depth[%s] average is positive, converting to negative values', ip)
                    if len(depth.shape) == 1: depth = -1 * d
                    else: depth[ip] = -1 * d
                    nmod += 1
            if nmod:
                self.warning('Depth of %s profiles were converted to negative values', nmod)
        return depth

    def load_depth_quality(self, dataset, *args, **kwargs):
        '''Load and return the depth quality variable (using variable_map)'''
        return self.load_variable(dataset, 'depth_qc', *args, **kwargs)

    def load_pressure(self, dataset, *args, **kwargs):
        '''Load and return the pressure variable (using :func:`load_variable`)

        :Params:
            - **tryconv**: if true, try to convert depth to pressure if depth cannot be loaded
            - **frompres**: force depth to pressure conversion

        '''
        fromdepth, tryconv = kwargs.pop('fromdepth', False), kwargs.pop('tryconv', True)
        if fromdepth:
            tryconv = True
            pres = None
        else:
            pres = self.load_variable(dataset, 'pressure', *args, **kwargs)
        if pres is None and tryconv:
            if not fromdepth:
                self.warning('No pressure found, now trying to convert depth to pressure')
            lat = self.load_latitude(dataset)
            if lat is None:
                self.error('Latitude is needed for this conversion')
                return None
            depth = self.load_depth(dataset, tryconv=False)
            if depth is None:
                self.error('Depth is needed for this conversion')
                return None
            self.debug('Latitude: %s', self.describe(lat))
            self.debug('Depth: %s', self.describe(depth))
            # NOTE: pressure in decibar, depth in meters, latitude in decimal degrees [-90,+90]
            pres = sw_pres(depth[:].swapaxes(0,1), lat[:]).swapaxes(0,1)
            # TODO: add explicit named axes
            proax = cdms2.createAxis(xrange(depth.shape[0]), id="profile")
            depax = cdms2.createAxis(xrange(depth.shape[1]), id="depth")
            pres = cdms2.createVariable(pres, id='pressure', axes=(proax, depax))
            self.verbose('Computed pressure: %s', self.describe(pres))
        return pres

    def load(self, dataset, **kwargs):
        '''Load a dataset.

        This loads a dataset or list of datasets and process quality filtering.

        :Params:
            - **dataset**: instance or list of instance of string (path) or CdmsFile
            - **variables**: variable identifiers to be loaded, defaults to self.variables
            - **qualities**: qualitie codes for filtering, defaults to self.qualities
            - **timerange**: if present, a time range for accepted profiles, defaults to self.timerange
            - **lonrange**: if present, a longitude range for accepted profiles, defaults to self.lonrange
            - **latrange**: if present, a latitude range for accepted profiles, defaults to self.latrange

        :Note: **this load all data into memory !**

        '''
        if isinstance(dataset, (list,tuple)): # replace this by is_iterable
            self.notice('%s loading %s profiles', self.__class__.__name__, len(dataset))
            for d in dataset:
                p = self.factory(dataset=d, load=False, spec=self, logger_config=self)# safe=True ?
                p.load(d, **kwargs)
                self.extend(p)
            return
        self.notice('%s loading: %s', self.__class__.__name__, dataset)
        varids = kwargs.get('variables', self.variables)
        qualities = kwargs.get('qualities', self.qualities)
        timerange = kwargs.get('timerange', self.timerange)
        lonrange = kwargs.get('lonrange', self.lonrange)
        latrange = kwargs.get('latrange', self.latrange)
        qtype = 'c' # int # lambda o:str(o).strip()
        #qualities = map(qtype, qualities)
        qualities = MV2z.array(qualities, dtype=qtype).tolist()
        # TODO: quality types in vars and in self ???
        self.verbose('Valid quality flags: %s', qualities)
        self.verbose('Requested variables identifiers: %s', varids)
        self.verbose('Time range filter: %s', timerange)
        self.verbose('Longitude range filter: %s', lonrange)
        self.verbose('Latitude range filter: %s', latrange)
        dsfile = None
        try:
            if isinstance(dataset, basestring):
                dsfile = cdms2.open(dataset)
            # Load dimensions or associated variables
            platvar = self.load_platform_code(dsfile)
            filevar = self.load_filename(dsfile)
            timevar, timecflag = self.load_time(dsfile, getcflag=True)
            timeqcvar = self.load_time_quality(dsfile)
            latvar = self.load_latitude(dsfile)
            lonvar = self.load_longitude(dsfile)
            posqcvar = self.load_position_quality(dsfile)
            depvar = self.load_depth(dsfile)
            depqcvar = self.load_depth_quality(dsfile)
            deprank = depvar.rank()
            # TODO:
            # - add other checks (var shapes, qc shapes)
            # - add variables standardisation in load_variable ?
            if platvar is None:
                self.error('No platform code variable found')
            if deprank > 2:
                raise ValueError('Unsupported depth rank: %s'%(depvar.rank()))
            npro, ndep = timevar.shape[0], depvar.shape[-1]
            varshape = (npro, ndep)
            variables = {}
            # Load variables
            for varid in varids:
                var = self.load_variable(dsfile, varid)
                if var is None:
                    #self.warning('Could not find %s ', varid)
                    continue
                # Check shape
                if var.shape != varshape:
                    self.error('Variable %s%s mismatch expected shape: %s', var.id, var.shape, varshape)
                    continue
                #self.verbose('Found %s: %s', varid, self.describe(var))
                variables[varid] = var
            #
            if len(qualities) and timeqcvar is None: self.warning('Could not perform time quality tests: quality flag not available')
            else: timeqcvar = MV2.array(timeqcvar, dtype=qtype)
            if len(qualities) and posqcvar is None: self.warning('Could not perform position quality tests: quality flag not available')
            else: posqcvar = MV2.array(posqcvar, dtype=qtype)
            if len(qualities) and depqcvar is None: self.warning('Could not perform depth quality tests: quality flag not available')
            else: depqcvar = MV2.array(depqcvar, dtype=qtype)
            timeqcerr = timeqcexc = posqcerr = posqcexc = depqcerr = depqcexc = 0
            # Load profiles timesteps
            npro = timevar.shape[0]
            for ipro in xrange(npro):
                self.debug('Processing profile[%s]', ipro)

                ptime = adatetime(cdtime.reltime(timevar[ipro], self.time_units))

                if timerange and not adatetime(timerange[0]) <= ptime <= adatetime(timerange[1]):
                    self.verbose('Excluding profile no %s: longitude %s not in range', ipro, lonvar[ipro])
                    continue
                if lonrange and not lonrange[0] <= lonvar[ipro] <= lonrange[1]:
                    self.verbose('Excluding profile no %s: longitude %s not in range', ipro, lonvar[ipro])
                    continue
                if latrange and not latrange[0] <= latvar[ipro] <= latrange[1]:
                    self.verbose('Excluding profile no %s: latitude %s not in range', ipro, latvar[ipro])
                    continue

                if not len(qualities) or timeqcvar is None: pass
                elif MV2.is_masked(timeqcvar[ipro]):
                    self.debug('Could not perform time quality test of profile no %s: quality flag not set', ipro)
                    timeqcerr += 1
                else:
                    self.debug('Checking time quality flag of profile no %s', ipro)
                    q = timeqcvar[ipro]
                    if q not in qualities:
                        self.verbose('Excluding profile no %s: time quality=%s', ipro, q)
                        timeqcexc += 1
                        continue

                if not len(qualities) or posqcvar is None: pass
                elif MV2.is_masked(posqcvar[ipro]):
                    self.debug('Could not perform position quality test of profil[%s]: quality flag not set', ipro)
                    posqcerr += 1
                else:
                    self.debug('Checking position quality flag of profile[%s]', ipro)
                    q = posqcvar[ipro]
                    if q not in qualities:
                        self.verbose('Excluding profile no %s: position quality=%s', ipro, q)
                        posqcexc += 1
                        continue

                if deprank == 1: pdepth = depvar
                else: pdepth = depvar[ipro]
                if not hasattr(pdepth, 'shape'):
                    pdepth = MV2.array(pdepth, fill_value=depvar.getMissing())
                pvars = {}
                for varid,var in variables.items():
                    pvars[varid] = var[ipro]

                if timecflag is not None: cftime = timecflag[ipro]
                else: cftime = False

                if isinstance(platvar, basestring): platcode = platvar
                elif platvar is not None: platcode = platvar[ipro]
                else: platcode = ''

                # Determine masked data based on the depth variable
                # This is done to optimize the memory/disk usage
                if pdepth.mask.size > 1:
                    orgshape = pdepth.shape
                    valid_indices = numpy.where(~pdepth.mask)[0]
                    if not len(valid_indices):
                        self.verbose('Excluding profile[%s]: depth is entirely masked', ipro)
                        continue
                    pdepth = pdepth.take(valid_indices)
                    for varid,var in pvars.items():
                        pvars[varid] = var.take(valid_indices)
                    self.debug('Shape of profile[%s] after depth optimization: %s => %s', ipro, orgshape, pdepth.shape)

##############################
                if False: # meth1
##############################
                    if not qualities or depqcvar is None: pass
                    else:
                        self.debug('Checking depth quality flag of profile no %s', ipro)
                        if depqcvar.rank() == 1: pdepqc = depqcvar
                        else: pdepqc = depqcvar[ipro]
                        valid_indexes = []
                        for idep,q in enumerate(pdepqc):
                            if MV2.is_masked(q):
                                self.debug('Could not perform depth no %s quality test of profile no %s: quality flag not set', idep, ipro)
                                depqcerr += 1
                            else:
                                if q in qualities: valid_indexes.append(idep)
                                else:
                                    self.verbose('Excluding depth no %s of profile no %s: depth quality=%s', idep, ipro, q)
                        pdepth = pdepth.take(valid_indexes)
                        for varid,var in variables.items():
                            pvars[varid] = pvars[varid].take(valid_indexes)
##############################
# OR meth2
##############################
#                    if not len(qualities) or depqcvar is None: pass
#                    else:
#                        if depqcvar.rank() == 1: pdepqc = depqcvar
#                        else: pdepqc = depqcvar[ipro]
#                        valid_indexes = []
#                        for qual in qualities:
#                            valid_indexes.extend(numpy.where(pdepqc==qual)[0])
#                        print pdepqc.shape
#                        pdepth = pdepth.take(valid_indexes)
#                        print pdepqc.shape
#                        for varid,var in variables.items():
#                            pvars[varid] = pvars[varid].take(valid_indexes)
##############################
# PRES_QC unknown format in /home11/caparmor/mars/VALID/DATA/MANGA/HYDRO_SISMER/MANGA_SISMER_PR_BO.nc
##############################

                profile = Profile(
                    parent=dsfile, index=ipro,
                    platform_code=platcode,
                    filename=filevar,
                    datetime=ptime, datetime_corrupted=cftime,
                    latitude=latvar[ipro],
                    longitude=lonvar[ipro],
                    depth=pdepth,
                    variables=pvars,
                    logger_config=self)

                self.append(profile)

                #self.debug('Adding profile: %s', profile)

            if qualities:
                if timeqcerr: self.warning('Ignored masked time quality flag of %s profile%s', timeqcerr, ('','s')[timeqcerr>1])
                if posqcerr: self.warning('Ignored masked position quality flag of %s profile%s', posqcerr, ('','s')[posqcerr>1])
                if depqcerr: self.warning('Ignored masked depth quality flag of %s depth%s', depqcerr, ('','s')[depqcerr>1])
                if timeqcexc: self.warning('Excluded time quality flag of %s profile%s', timeqcexc, ('','s')[timeqcexc>1])
                if posqcexc: self.warning('Excluded position quality flag of %s profile%s', posqcexc, ('','s')[posqcexc>1])
                if depqcexc: self.warning('Excluded depth quality flag of %s depth%s', depqcexc, ('','s')[depqcexc>1])
            n = len(self)
            self.notice('Loaded %s profile%s', n, ('','s')[n>1])

        except Exception, e:
            self.exception('Error loading %s', dataset)
        finally:
            if dsfile: dsfile.close()

    def get_profile_axis(self):
        '''Return the profile cdms axis of the nested profiles'''
        #indexes = numpy.array(range(len(self)))
        indexes = numpy.array([p.index for p in self])
        pro = cdms2.createAxis(indexes, id='profile')
        pro.description = 'Index of profile in original file'
        return pro

    def get_level_axis(self):
        '''Return the level cdms axis of the nested profiles'''
        npro, ndep = len(self), 0
        if npro:
            ndep = max([p.depth.shape[0] for i,p in enumerate(self)])
        depth = cdms2.createAxis(numpy.array(range(ndep)), id='level')
        depth.standard_name = 'level'
        depth.long_name = 'Level'
        depth.axis = 'Z'
        return depth

    def get_platform_code(self):
        '''Return the platform code cdms variable of the nested profiles, if suitable, None otherwise'''
        platmaxlen = max([len(p.platform_code) for i,p in enumerate(self)])
        if platmaxlen:
            platax = cdms2.createAxis(numpy.array(range(platmaxlen)), id='platform_code_string')
            platarr = numpy.array([numpy.array(p.platform_code.ljust(platmaxlen) or ' '*platmaxlen, dtype='c') for p in self], dtype='c')
            platvar = cdms2.createVariable(platarr, id='platform_code', typecode='c', axes=(self.get_profile_axis(), platax), attributes = {})
            return platvar

    def get_origin_filename(self):
        '''Return the filename cdms variable of the nested profiles '''
        filemaxlen = max([len(p.filename) for i,p in enumerate(self)])
        if filemaxlen:
            fileax = cdms2.createAxis(numpy.array(range(filemaxlen)), id='filename_string')
            filearr = numpy.array([numpy.array(p.filename.ljust(filemaxlen) or ' '*filemaxlen, dtype='c') for p in self], dtype ='c')
            filevar = cdms2.createVariable(filearr, id='filename', typecode='c', axes=(self.get_profile_axis(), fileax), attributes = {})
            filevar.description = 'File name of the original file'
            return filevar

    def get_time(self):
        '''Return the time cdms variable of the nested profiles'''
        time = cdms2.createVariable(numpy.array([reltime(p.datetime, self.time_units).value for p in self]),
            id='time', axes=(self.get_profile_axis(),), attributes = {
            'standard_name':'time', 'long_name':'Time', 'axis':'T',
            'units':self.time_units})
        return time

    def get_latitude(self):
        '''Return the latitude cdms variable of the nested profiles'''
        lat = cdms2.createVariable(numpy.array([p.latitude for p in self]),
            id='latitude', axes=(self.get_profile_axis(),), attributes = {
            'standard_name':'latitude', 'long_name':'Latitude', 'axis':'Y',
            'units':'degree_north', 'valid_min':-90.0, 'valid_max':90.0})
        return lat

    def get_longitude(self):
        '''Return the longitude cdms variable of the nested profiles'''
        lon = cdms2.createVariable(numpy.array([p.longitude for p in self]),
            id='longitude', axes=(self.get_profile_axis(),), attributes = {
            'standard_name':'longitude', 'long_name':'Longitude', 'axis':'X',
            'units':'degree_east', 'valid_min':-180.0, 'valid_max':180.0})
        return lon

    def get_depth(self):
        '''Return the depth cdms variable of the nested profiles'''
        npro, ndep = len(self), 0
        if npro:
            ndep = max([p.depth.shape[0] for i,p in enumerate(self)])
        depth = MV2.zeros((npro, ndep))
        depth = MV2.masked_where(depth == 0, depth)
        for i,p in enumerate(self):
            pdep = p.depth[:]
            depth[i][0:pdep.shape[0]] = pdep
            depth.mask[i][0:pdep.shape[0]] = pdep.mask
        depthvar = cdms2.createVariable(depth, id='depth', axes=(self.get_profile_axis(), self.get_level_axis()), attributes = {
            'standard_name':'depth', 'long_name':'Depth', 'axis':'Z',
            'units':'m', 'positive':'down', 'valid_min':MV2.min(depth) if len(depth) else None, 'valid_max':MV2.max(depth) if len(depth) else None})
        return depthvar

    def get_variable(self, varname):
        '''Return the varname (identifier) cdms variable of the nested profiles'''
        npro, ndep = len(self), 0
        if npro:
            ndep = max([p.depth.shape[0] for i,p in enumerate(self)])
        varattrs = {
            'temperature':{'long_name':'sea water temperature', 'standard_name':'sea_water_temperature', 'units':'degree_Celsius'},
            'salinity':{'long_name':'sea water salinity', 'standard_name':'sea_water_salinity', 'units':'PSU'},
        }
        emptyvar = MV2.zeros((ndep,))
        emptyvar = MV2.masked_where(emptyvar == 0, emptyvar)
        data = MV2.zeros((npro, ndep))
        data = MV2.masked_where(data == 0, data)
        for i,p in enumerate(self):
            a = p.variables.get(varname, emptyvar)
            data[i][0:a.shape[0]] = a
            data.mask[i][0:a.shape[0]] = a.mask
        var = cdms2.createVariable(data, id=varname,
            axes=(self.get_profile_axis(), self.get_level_axis()),
            attributes=varattrs.get(varname, {}))
        return var

    def save(self, filepath, **kwargs):
        '''
        Save profiles in a netcdf file.
        :Params:
            - **filepath**: output netcdf file path
            - **variables**: list of variable identifiers to be saved, default is self.variables
        '''
        variables = kwargs.get('variables', self.variables)
        self.info('Saving profiles in: %s', filepath)
        self.verbose('Requested variables identifiers: %s', variables)
        dataset = None
        try:
            npro, ndep = len(self), 0
            if npro:
                ndep = max([p.depth.shape[0] for i,p in enumerate(self)])
            if not npro:
                self.warning('There is no profile to save, aborting')
                return
            if not ndep:
                self.warning('There is no profile depth to save, aborting')
                return
            self.info('Number of profiles: %s', npro)
            self.info('Number of levels:   %s', ndep)
            dataset = cdms2.open(filepath, 'w')
            tim, lat, lon, plat, filename = self.get_time(), self.get_latitude(), self.get_longitude(), self.get_platform_code(), self.get_origin_filename()
            dataset.write(tim)
            if filename is None: self.warning('No filename variable')
            else:
                try: dataset.write(filename)
                except: self.exception('Cannot write filename variable')
            dataset.write(lat)
            dataset.write(lon)
            if plat is None: self.warning('No platform code variable')
            else:
                try: dataset.write(plat)
                except: self.exception('Cannot write platform variable')
            dataset.write(self.get_depth())
            for varname in variables:
                self.info('Writing variable %s', varname)
                dataset.write(self.get_variable(varname))
            creation_date = time.strftime('%Y-%m-%dT%H:%M:%SZ')
            dataset.software_version = '%s %s'%(os.path.basename(sys.argv[0]), __date__)
            dataset.creation_date = creation_date
            dataset.history = '%s : Creation'%(creation_date)
            dataset.title = 'Profiles'
            dataset.southernmost_latitude = MV2.min(lat)
            dataset.northernmost_latitude = MV2.max(lat)
            dataset.westernmost_longitude = MV2.min(lon)
            dataset.easternmost_longitude = MV2.max(lon)
            return True
        except Exception, e:
            self.exception('Failed to save file %s', filepath)
        finally:
            if dataset: dataset.close()
        return False

    def select(self, variables=None, time=None, polys=None, **kwargs):
        '''
        Return a new Profiles instance (profiles are still referenced) of selected profiles.

        :Params:
            -**variables**: variables identifiers filtering
            -**time**: time filtering, see :func:`is_time_in_range`
            -**polys**: list of polygons (_geoslib.Polygon) for position filtering
        '''
        profiles = self.__class__(self)
        del profiles[:]
        if not variables: variables = self.variables
        if isinstance(variables, basestring): variables = (variables,)
        if polys and not isinstance(polys, (list, tuple)): polys = (polys,)
        for ipro,pro in enumerate(self):
            # Check at least one variable matches
            if variables:
                for v in pro.variables.keys():
                    if v in variables:
                        break
                else:
                    self.debug('profile %s variables %s match %s', ipro, pro.variables.keys(), variables)
                    continue
            # Check time matches
            if time and not is_time_in_range(pro.datetime, time):
                self.debug('profile %s datetime %s not in range %s', ipro, pro.datetime, time)
                continue
            # Check at least on polygon matches
            if polys:
                for ipoly,poly in enumerate(polys):
                    if Point((pro.longitude, pro.latitude)).within(poly):
                        self.debug('profile %s position %s is in poly %s', ipro, (pro.longitude, pro.latitude), ipoly)
                        break
                else:
                    self.debug('profile %s position %s not in poly %s', ipro, (pro.longitude, pro.latitude), ipoly)
                    continue
            profiles.append(pro)
        self.info('Profiles after selection:\n%s', profiles.summary())
        return profiles



# ==============================================================================

class ProfilesWithDepthFromPres(Profiles):
    '''
    A special Profiles type which forces trying to load depth from pressure first.
    '''
    def load_depth(self, dataset, *args, **kwargs):
        # Force using pressure first if available
        kwargs['frompres'] = True
        d = Profiles.load_depth(self, dataset, *args, **kwargs)
        if d is not None: return d
        self.verbose('A depth from pressure conversion was expected but no pressure was found. Now trying to load depth')
        kwargs['frompres'] = False
        return Profiles.load_depth(self, dataset, *args, **kwargs)


class Profiles_PR_BA(ProfilesWithDepthFromPres):
    '''Profils BATHY provenant du GTS'''
    @classmethod
    def get_type(cls): return 'PR_BA'

    @classmethod
    def get_variables_map(cls):
        m = Profiles.get_variables_map()
        # We need to insert DEPH at the start of the depth names list
        # to ensure this will be the first name to be checked.
        # In these profiles, the DEPTH variable is not the one we are looking for !
        # NOTE: this is now the default behavior, see Profiles class
        m['depth'].insert(0, 'DEPH')
        return m


class Profiles_PR_BO(ProfilesWithDepthFromPres):
    '''Profils Bouteilles'''
    @classmethod
    def get_type(cls): return 'PR_BO'

    @classmethod
    def get_variables_map(cls):
        m = Profiles.get_variables_map()
        m['depth_qc'].insert(0, 'PROFILE_PRES_QC')
        return m


class Profiles_PR_CT(ProfilesWithDepthFromPres):
    '''Profils CTD'''
    @classmethod
    def get_type(cls): return 'PR_CT'


class Profiles_PR_PF(ProfilesWithDepthFromPres):
    '''Profils flotteurs Argo'''
    @classmethod
    def get_type(cls): return 'PR_PF'

    @classmethod
    def get_variables_map(cls):
        m = Profiles.get_variables_map()
        # We need to remove DEPTH because this kind of profiles do not really
        # provide it, this will cause to fall into pressure to depth conversion
        # NOTE: now this is done by subclassing ProfilesWithDepthFromPres
        m['depth'] = []
        return m


class Profiles_PR_RE(Profiles):
    '''Profils Recopesca'''
    @classmethod
    def get_type(cls): return 'PR_RE'


class Profiles_PR_TE(ProfilesWithDepthFromPres):
    '''Profils TESAC provenant du GTS'''
    @classmethod
    def get_type(cls): return 'PR_TE'

class Profiles_PR_XB(Profiles):
    '''Profils XBT ou XCTD'''
    @classmethod
    def get_type(cls): return 'PR_XB'


class Profiles_Selmed(ProfilesWithDepthFromPres):
    '''Profils utilises pour BOBYCLIM'''
    @classmethod
    def get_type(cls): return 'selmed'

    @classmethod
    def get_variables_map(cls):
        m = Profiles.get_variables_map()
        # We need to remove DEPTH because this kind of profiles do not really
        # provide it, this will cause to fall into pressure to depth conversion
        # NOTE: now this is done by subclassing ProfilesWithDepthFromPres
        m['depth'] = []
        m['time'] = ('STATION_DATE_TIME',)
        m['time_qc'] = []
        m['pressure'] = ('PRES',)
        m['latitude'] = ('LATITUDE',)
        m['longitude'] = ('LONGITUDE',)
        m['temperature'] = ('TEMP',)
        m['salinity'] = ('PSAL',)
        return m

class Profiles_ArgoGeo(ProfilesWithDepthFromPres):
    '''Profils ARGO - Geo'''
    @classmethod
    def get_type(cls): return '_prof.nc'

    @classmethod
    def get_variables_map(cls):
        m = Profiles.get_variables_map()
        # We need to remove DEPTH because this kind of profiles do not really
        # provide it, this will cause to fall into pressure to depth conversion
        # NOTE: now this is done by subclassing ProfilesWithDepthFromPres
        m['depth'] = []
        m['platform_code'] = ('PLATFORM_NUMBER',)
        m['time'] = ('JULD',)
        m['time_qc'] = ('JULD_QC',)
        m['pressure'] = ('PRES_ADJUSTED',)
        m['depth_qc'] = ('PRES_QC',)
        m['latitude'] = ('LATITUDE',)
        m['longitude'] = ('LONGITUDE',)
        m['position_qc'] = ('POSITION_QC',)
        m['temperature'] = ('TEMP_ADJUSTED',)
        m['salinity'] = ('PSAL_ADJUSTED',)
        return m


# ==============================================================================

class ProfilesFilter(Object):
    def filter(self, *args, **kwargs):
        '''
        Must be defined by implementing subclasses
        This method must return two Profiles instances (filtered and rejected)
        '''
        raise NotImplementedError('')

class ProfilesDuplicatesFilter(ProfilesFilter):
    '''
    Split a profiles list into two profiles list:
        - a list of unique profile (filtered)
        - a list of duplicated profiles (rejected)

    The duplication criterias are based on the equality of the profile's attributes:
        - platform code
        - latitude and longitude
        - datetime

    The profiles datetime must be flagged with a **datetime_corrupted** attribute
    indicating if the time component of the datetime atribute is reliable (False)
    or has been arbitrary set to 12H (True).
    In this case, two profiles are duplicates if the timedelta between their datetimes
    is less than **maxtimedelta**.
    '''
    def __init__(self, maxtimedelta=60*60*12, **kwargs):
        '''
        :Params:
            - **maxtimedelta**: number of seconds or timedelta representing the difference
              between two profiles dates from which they are not considered duplicates.
              This is only used when profiles have their datetime_corrupted flagged to True.
        '''
        ProfilesFilter.__init__(self, **kwargs)
        if not isinstance(maxtimedelta, datetime.timedelta):
            maxtimedelta = datetime.timedelta(seconds=maxtimedelta)
        self.maxtimedelta = maxtimedelta

    def is_duplicate(self, pro1, pro2):
        '''
        Return a boolean indicating if pro1 and pro2 are considered duplicates.
        '''
        # Check some attributes equality
        if pro1.latitude != pro2.latitude \
        or pro1.longitude != pro2.longitude \
        or pro1.platform_code.strip().upper() != pro2.platform_code.strip().upper():
            return False
        # If one of the datetimes is corrupted, equality is based on the dates timedelta
        if pro1.datetime_corrupted or pro2.datetime_corrupted:
            if abs(pro1.datetime - pro2.datetime) > self.maxtimedelta:
                return False
        # Else check strict equality
        else:
            if pro1.datetime != pro2.datetime: return False
        # Here we are, the two profiles are considered duplicates
        return True

    def filter(self, profiles):
        '''
        Filter profiles, split them into two list of filtered (unique) and rejected (duplicates).
        :Params:
            - **profiles**: A Profiles instance to be filtered
        :Return:
            - **filtered**: A Profiles instance
            - **rejected**: A Profiles instance
        :Note:
            The returned profiles are initialized with the same input profiles specification.
        '''
        # Build two list of profile indexes, unique and duplicates
        filtered_indexes, filtered, rejected_indexes, rejected = [], [], [], []
        # Loop on profiles
        for icur,curpro in enumerate(profiles):
            replaced = False
            # Loop on other profiles remove duplicates
            # and to find a replacement candidate if needed
            for i,pro in enumerate(profiles):
                # Check they don't refer to the same object
                if id(curpro) != id(pro):
                    # If curpro and pro are considered duplicates
                    if self.is_duplicate(curpro, pro):
                        # If curpro has a corrupted datetime and pro has a valid datetime and curpro has not been replaced
                        if curpro.datetime_corrupted and not pro.datetime_corrupted and not replaced:
                            # Add pro as a replacement of curpro, rejecting curpro
                            filtered_indexes.append(i)
                            rejected_indexes.append(icur)
                            # Also add curpro
                            rejected_indexes.append(i)
                            replaced = True
                        # Else curpro's datetime is valid or has been replaced
                        else:
                            rejected_indexes.append(i)
            # If curpro has not been replaced because it has a valid datetime or no replacement was found
            if not replaced:
                filtered_indexes.append(icur)
        filtered_indexes = numpy.unique(filtered_indexes)
        rejected_indexes = numpy.unique(rejected_indexes)
        if len(filtered_indexes): filtered = numpy.array(profiles).take(filtered_indexes)[:]
        if len(rejected_indexes): rejected = numpy.array(profiles).take(rejected_indexes)[:]
        return Profiles(filtered, spec=profiles), Profiles(rejected, spec=profiles)


# ==============================================================================


class ProfilesMerger(Object):
    '''
    Hold a Profiles object which must be populated with the load method and then merged
    with the merge method.
    '''
    def __init__(self, spec=None, **kwargs):
        Object.__init__(self, **kwargs)
        self.spec = spec if spec else {}
        self.profiles = Profiles(**self.spec)

    def load(self, dataset, **kwargs):
        '''
        Add profile(s) to the internal instance.
        :Params:
            - **dataset**: can be single or list of objects accepted by :meth:`AbstractProfiles.factory`
            - **kwargs**: additionnal arguments merged with self.spec and passed to :meth:`AbstractProfiles.factory`
        '''
        # This dataset collection detection is also available in Profiles, but here we
        # want to load them separatly to exclude unloadable profiles without a crash
        if isinstance(dataset, (list,tuple)):
            for d in dataset:
                self.load(d, **kwargs)
            return
        # self.psinfo()
        try:
            spec = self.spec.copy()
            spec.update(kwargs)
            profiles = Profiles.factory(dataset=dataset, logger_config=self, **spec)
            profiles.load(dataset)
            self.profiles.extend(profiles)
        except Exception:
            self.exception('Could not process %s', dataset)

    def merge(self, filter=None, filter_file=None, reject_file=None, filter_sort=None, reject_sort=None, **kwargs):
        '''
        Merge the profiles applying a filter originally designed to seperate duplicate profiles.
        See :class:`ProfilesDuplicatesFilter`.
        :Params:
            - **filter**: the filter class to be used
            - **filter_file**: optionnal, string, the filtered profiles output file path
            - **reject_file**: optionnal, string, the filtered profiles output file path
            - **filter_sort**: optional, string, see :meth:`Profiles.sort`
            - **reject_sort**: optional, string, see :meth:`Profiles.sort`
        '''
        if isinstance(filter, basestring):
            filter = globals().get(filter)
        if isinstance(filter, type):
            filter = filter()
        self.notice('Merging %s profiles with filter %s', len(self.profiles), filter.__class__.__name__)
        # self.psinfo()
        if filter:
            filtered, rejected = filter.filter(self.profiles)
        else:
            filtered, rejected = self.profiles , Profiles(self.spec)
        self.notice('Number of filtered profiles: %s', len(filtered))
        self.notice('Number of rejected profiles: %s', len(rejected))
        # self.psinfo()
        if filtered and filter_file:
            if filter_sort:
                self.notice('Sorting filtered profiles with method: %s', filter_sort)
                filtered.sort(filter_sort)
            self.notice('Saving filtered profiles in file: %s', filter_file)
            filtered.save(filter_file)
            # self.psinfo()
        if rejected and reject_file:
            if reject_sort:
                self.notice('Sorting rejected profiles with method: %s', reject_sort)
                rejected.sort(reject_sort)
            self.notice('Saving rejected profiles in file: %s', reject_file)
            rejected.save(reject_file)
            # self.psinfo()
        return filtered, rejected

    @classmethod
    def main(cls, args=None):
        '''
        Entry point of the bin/merge_profiles.py script. See merge_profiles.py --help.

        :Params:
            - **args**: passed to config.ConfigManager.opt_parse

        .. note::
            - You may pass the "args" argument to override command line arguments

        :Examples:

            >>> ProfilesMerger.main(args='--cfgfile myconf.cfg --verbose --ProfilesMerger-load-timerange "2010-01-01 2010-03-31"'.split())

        '''
        try:
            # Create the command line parser
            parser = optparse.OptionParser(
                description='Merge various vertical profiles into a netcdf file',
                usage='%prog [options] [profilesA.nc] [profilesB.nc] [...] ',
                version=__date__)
            parser.add_option('--netcdf3', action='store_true', dest='netcdf3', default=False,  help='Save files in NetCDF version 3')
            parser.add_option('--compress', action='store', dest='compress', default=3,  type=int, help='Set the compression level for NetCDF version 4 (default: %default, range [0,9], only used with --convert)')
            parser.add_option('--nonumpywarnings', action='store_true', dest='nonumpywarnings', default=False,  help='Hide all numpy warnings')
            # Add logging options
            Logger.add_optparser_options(parser)#, prefix='log')

            # Prepare the configuration
            # Specifications config file
            cfgini = '%s.ini'%os.path.splitext(os.path.realpath(__file__))[0]
            # Configuration section name
            cfgsec = cls.__name__
            # The configuration manager, using the specification file above
            cfgm = ConfigManager(cfgini)
            # Get default configuration
            cfgo = cfgm.defaults()
            # The manager parse the command line arguments and return a configuration
            cfgp = cfgm.opt_parse(parser, args=args)
            options, args = parser.values, parser.largs
            # If a configuration file is given, load it and override defaults
            if options.cfgfile:
                if not os.path.isfile(options.cfgfile):
                    cls.error('Configuration file does not exists: %s', options.cfgfile)
                else:
                    cls.notice('Loading configuration file: %s', options.cfgfile)
                    cfgf = cfgm.load(options.cfgfile)
                    cfgm.patch(cfgo, cfgf)
            # Last override from command line arguments
            cfgm.patch(cfgo, cfgp)
            # Apply logging command line arguments
            Logger.apply_class_optparser_options(options)
            # Show the final configuration
            cls.info('Loaded configuration:\n%s', pprint.pformat(cfgo.dict()))

            # Operations on masked values throws numpy warnings, this can be ignored with the following instruction
            if options.nonumpywarnings:
                numpy.seterr(all='ignore')
            # NetCDF format
            if options.netcdf3:
                cls.notice('Using on NetCDF version 3 writting mode')
                netcdf3()
            elif options.compress:
                cls.notice('Using on NetCDF version 4 writting mode with compression level %s', options.compress)
                netcdf4(options.compress)

            # Get the load section
            loadkw = cfgo[cfgsec].get('load', {})
            # Check that variable are requested (no sense to create a profile without data ?)
            if not loadkw['variables']:
                cls.warning('No variables were specified')

            # Create the merger with the load configuration (variables, qualities, time/lon/lat range, ...)
            merger = cls(spec=loadkw)

            # Load profiles specified as command line arguments
            merger.load(args)

            # Load profiles from explicit list of config file or command line option
            input_files = cfgo[cfgsec].get('input_files', [])
            if input_files:
                merger.load(input_files)

            # Load profiles with a search specified from config file or command line option
            findkw = cfgo[cfgsec].get('find', {})
            if findkw.get('path', None):
                regex, pattern = findkw.get('regex', None), findkw.get('pattern', None)
                if regex and pattern:
                    raise ValueError('Use regex or pattern option, not both')
                elif regex:
                    findkw.pop('pattern', None)
                    findfunc = efind
                else:
                    findkw.pop('regex', None)
                    findfunc = find
                cls.info('Searching files using method %s and options: %s', findfunc.__name__, findkw)
                files = findfunc(**findkw)
                cls.info('Found %s file%s', len(files), ('','s')[len(files)>1])
                merger.load(files)

            # Now merge the profiles with filters specified in configuration (e.g. duplicates filtering)
            merger.merge(**cfgo[cfgsec]['merge'])

            return 0
        except Exception, e:
            cls.exception('Fatal error: %s', e)
            return 1


# ==============================================================================


class ProfilesDataset(OceanDataset):
    '''
    Class for handling standardized profiles dataset as created by :meth:`Profiles.save`

    .. todo::
        - allow to load files with different depth size

    '''
    _platcodeid = ('platform_code',)
    _filenameid = ('filename',)
    _timeid = ('time',)
    _latid = ('latitude',)
    _lonid = ('longitude',)
    _depid = ('depth',)

    def get_axis(self, axname, select=None):
        # TODO: add dataset.get_axis by id behavior
        # TODO: make this method works with profiles dataset having different levels
        self.debug('Getting axis: %s, select: %s', axname, select)
        if not len(self.dataset):
            self.debug('No %s axis found, dataset is empty', axname)
            return None
        orgselect, (tmp,select) = select, self.get_selector(split=True, **select)
        seltime = select.pop('time', None)
        if seltime: seltime = map(comptime, seltime[:2])
        sellat = select.pop('latitude', select.pop('lat', None))
        if sellat: sellat = map(float, sellat[:2])
        sellon = select.pop('longitude', select.pop('lon', None))
        if sellon: sellon = map(float, sellon[:2])
        selpolys = select.pop('polygons', None)
        if selpolys and not is_iterable(selpolys, (list, tuple)): selpolys = (selpolys,)
        ax, allaxes = None, []
        for dataset in self.dataset:
            ax = dataset.getAxis(axname)
            if ax is None:
                self.warning('Axis %s not found in %s', axname, dataset.id)
                continue
            if axname == 'profile' and (seltime or sellat or sellon or selpolys):
                tmpax = []
                if seltime:
                    time = ncget_var(dataset, self._timeid)
                    time = create_time(time, units=time.units).asComponentTime()
                if sellat or selpolys:
                    lat = create_lat(ncget_var(dataset, self._latid))
                if sellon or selpolys:
                    lon = create_lon(ncget_var(dataset, self._lonid))
                for ipro in xrange(len(ax)):
                    # Check time matches
                    #if seltime and not (adatetime(seltime[0]) <= adatetime(time[ipro]) <= adatetime(seltime[1])): continue
                    if seltime and not (seltime[0] <= time[ipro] <= seltime[1]): continue
                    # Check position matches
                    if sellat and not sellat[0] <= lat[ipro] <= sellat[1]: continue
                    if sellon and not sellon[0] <= lon[ipro] <= sellon[1]: continue
                    # Check at least on polygon matches
                    if selpolys:
                        for ipoly,poly in enumerate(selpolys):
                            if Point((lon[ipro], lat[ipro])).within(poly):
                                self.debug('profile %s is in poly %s', ipro, ipoly)
                                break
                        else:
                            self.debug('profile %s not in poly %s', ipro, ipoly)
                            continue
                    # Selects pass, load the profile
                    tmpax.append(ax[ipro])
                if not tmpax: continue
                attrs = ax.attributes.copy()
                ax = cdms2.createAxis(tmpax, id=ax.id)
                for k,v in attrs.items(): setattr(ax, k, v)
            if ax is not None:
                allaxes.append(ax)
        if allaxes:
            refax = allaxes[0]
            attrs = refax.attributes.copy()
            ax = cdms2.createAxis(MV2.concatenate(allaxes), id=refax.id)
            for k,v in attrs.items(): setattr(ax, k, v)
            self.verbose('Loaded axis: %s', self.describe(ax))
            # self.psinfo()
            return ax
        self.warning('No data found for axis: %s, select: %s', axname, orgselect)

    def get_variable(self, varname, select=None, squeeze=False):
        # TODO: make this method works with profiles dataset having different levels
        self.debug('Getting variable: %s, select: %s', varname, select)
        if not len(self.dataset):
            self.debug('No %s variable found, dataset is empty', varname)
            return None
        if select is None: select = {}
        orgselect, (tmp,select) = select, self.get_selector(split=True, **select)
        seltime = select.pop('time', None)
        if seltime: seltime = map(comptime, seltime[:2])
        sellat = select.pop('latitude', select.pop('lat', None))
        if sellat: sellat = map(float, sellat[:2])
        sellon = select.pop('longitude', select.pop('lon', None))
        if sellon: sellon = map(float, sellon[:2])
        selpolys = select.pop('polygons', None)
        if selpolys and not is_iterable(selpolys, (list, tuple)): selpolys = (selpolys,)
        var, allvars = None, []
        for dataset in self.dataset:
            var = ncget_var(dataset, varname)
            if var is None:
                self.warning('Variable %s not found in %s', varname, dataset.id)
                continue
            # NOTE: data have to be loaded now (otherwise MV2.concatenate will not use the masks!)
            var = var(**select)
            if seltime or sellat or sellon or selpolys:
                axes = var.getAxisList()
                pro, tmppro, tmpvar = axes[0], [], []
                if seltime:
                    time = ncget_var(dataset, self._timeid)
                    time = create_time(time, units=time.units).asComponentTime()
                if sellat or selpolys:
                    lat = create_lat(ncget_var(dataset, self._latid))
                if sellon or selpolys:
                    lon = create_lon(ncget_var(dataset, self._lonid))
                for ipro in xrange(var.shape[0]):
                    # Check time matches
                    #if seltime and not (adatetime(seltime[0]) <= adatetime(time[ipro]) <= adatetime(seltime[1])): continue
                    if seltime and not (seltime[0] <= time[ipro] <= seltime[1]): continue
                    # Check position matches
                    if sellat and not sellat[0] <= lat[ipro] <= sellat[1]: continue
                    if sellon and not sellon[0] <= lon[ipro] <= sellon[1]: continue
                    # Check at least on polygon matches
                    if selpolys:
                        for ipoly,poly in enumerate(selpolys):
                            if Point((lon[ipro], lat[ipro])).within(poly):
                                self.debug('profile %s is in poly %s', ipro, ipoly)
                                break
                        else:
                            self.debug('profile %s not in poly %s', ipro, ipoly)
                            continue
                    # Selects pass, load the profile
                    tmppro.append(pro[ipro])
                    tmpvar.append(var[ipro])
                if not tmppro: continue
                newpro = cdms2.createAxis(tmppro, id=pro.id)
                for k,v in pro.attributes.items(): setattr(newpro, k, v)
                var = cdms2.createVariable(tmpvar, id=var.id, axes=[newpro]+axes[1:], attributes=var.attributes.copy())
            if var is not None:
                allvars.append(var)
        if allvars:
            refvar = allvars[0]
            attributes = refvar.attributes.copy()
            axis = cdms2.createAxis(MV2.concatenate([v.getAxis(0) for v in allvars]),
                id=refvar.getAxisIds()[0])
            for k,v in refvar.getAxis(0).attributes.items(): setattr(axis, k, v)
            axes = [axis] + [a.clone() for a in refvar.getAxisList()[1:]]
            var = cdms2.createVariable(MV2.concatenate(allvars), id=refvar.id, attributes=attributes, axes=axes)
            self.verbose('Loaded variable: %s', self.describe(var))
            if squeeze:
                var = self.squeeze_variable(var)
            # self.psinfo()
            return var
        self.warning('No data found for variable: %s, select: %s', varname, orgselect)

    def get_platform_code(self, **kwargs):
        return self.get_variable(self._platcodeid, **kwargs)

    def get_origin_filename(self, **kwargs):
        return self.get_variable(self._filenameid, **kwargs)

    def get_time(self, **kwargs):
        return self.get_variable(self._timeid, **kwargs)

    def get_latitude(self, **kwargs):
        return self.get_variable(self._latid, **kwargs)

    def get_longitude(self, **kwargs):
        return self.get_variable(self._lonid, **kwargs)

    def get_depth(self, **kwargs):
        return self.get_variable(self._depid, **kwargs)

    # TODO: check if profile/gridded data can be handled directly in Dataset
    def get_layer(self, varname, depth, select=None):
        '''
        Get layer data for a specified depth.

        :Params:
            - **variable**: variable name
            - **depth**: layer depth
            - **select**: selector
        '''
        var = self.get_variable(varname, select=select)#, squeeze=True)#, order='-z')
        if var is None: raise Exception('Variable %s not found'%varname)
        if var.getOrder() not in ('-z'):
            raise ValueError('Invalid var order: %s, -z expected'%(var.getOrder()))
        lon = self.get_longitude(select=select)
        lat = self.get_latitude(select=select)
        dep = self.get_depth(select=select)#, squeeze=True)#, order='-z')
        if dep is None: raise Exception('Depth not found')
        if dep.getOrder() not in ('-z'):
            raise ValueError('Invalid depth order: %s, -z expected'%(dep.getOrder()))
        npro = var.shape[0]

        # Exclude entirely masked depth/var coordinates
        goodfct = lambda a: numpy.ma.where(numpy.ma.count(a, axis=1)>0, True, False)
        depgood, vargood = goodfct(dep), goodfct(var)
        good = depgood & vargood
        depcnt, varcnt = len(depgood.nonzero()[0]), len(vargood.nonzero()[0])
        totcnt = len(good.nonzero()[0])
        self.debug('Number of profiles with valid data:'
            '\n  depth:            %s/%s'
            '\n  %-18s%s/%s'
            '\n  total:            %s/%s', depcnt, npro, '%s:'%var.id, varcnt, npro, totcnt, npro)
        if not totcnt:
            self.error('No valid profiles, aborting')
            return
        igood = numpy.where(good)[0]
        takefct = lambda a: MV2.take(a, igood, axis=0)
        dep, var = takefct(dep), takefct(var)

        # Depth interpolation axis
        odep = create_dep(is_iterable(depth) and depth or [depth])
        ovar = interp1d(var, axo=odep, xmap=[0], xmapper=dep) # Required order: ...z
        ovar = cdms2.createVariable(ovar, id=var.id, axes=[var.getAxis(0), odep], attributes=var.attributes.copy())
        self.describe(ovar, odep)
        return lon, lat, ovar

    def get_hist(self, tstep, **kwargs):
        '''
        Get histogram data.

        :Params:
            - **tstep**: bars time coverage (n,units) (ex: 1,'month')
            - **kwargs**: passed to :func:`get_time`
        '''
        time = self.get_time(**kwargs)
        dtime = adatetime(create_time(time, units=time.units))
        npro = time.shape[0]
        hist, hist_time = [], []
        for itv in Intervals((round_date(min(dtime), tstep[1]), round_date(add_date(max(dtime), tstep[0], tstep[1]), tstep[1])), tstep):
            hist_time.append(adatetime(itv[0]))
            n = len(filter(lambda d: adatetime(itv[0]) <= d < adatetime(itv[1]), dtime))
            hist.append(n)
            self.verbose('Interval %s: %s', itv[:2], n)
        hist_time = create_time(hist_time)
        hist = cdms2.createVariable(hist, axes=(hist_time,), id='hist', attributes=dict())
        return hist

    def plot_layer(self, varname, depth, *args, **kwargs):
        '''
        Plot a layer for a specified depth.

        See :func:`get_layer`
        '''
        mp = kwargs.pop('map', None)
        mapkw = kwfilter(kwargs, 'map', default=dict())
        for k in POST_PLOT_ARGS: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        lon, lat, var = self.get_layer(varname, depth, **kwargs)
        var = self.squeeze_variable(var)
        if len(var.shape) > 1:
            self.verbose('Averaging layer along depth')
            var = MV2.average(var, axis=1)
        self.describe(var, lat, lon)
        vmin, vmax = numpy.min(var), numpy.max(var)
        if mp is None:
            mapkw.update(dict(label=self.__class__.__name__, show=False,
                lon=(min(lon)-0.5,max(lon)+0.5), lat=(min(lat)-0.5,max(lat)+0.5),
                vmin=vmin, vmax=vmax, cmap='magic'))# cmap=cmap))
            mp = map2(**mapkw)
        # Plot profiles scatter points
        scakw = kwfilter(kwargs, 'sca', defaults=dict(s=20, vmin=vmin, vmax=vmax, cmap=mp.cmap))
        sc = mp.map.scatter(lon, lat, c=var, **scakw)
        # Post plotting
        plotkw.setdefault('title', var.id)
        mp.post_plot(**plotkw)

    def plot_hist(self, *args, **kwargs):
        '''
        Produce a histogram plot.

        :Keyword arguments:
            - **plot_<keyword>**: are passed to the plot function :func:`~vacumm.misc.plot.bar2`

        Other arguments, see :func:`get_hist`
        '''
        plotkw = kwfilter(kwargs, 'plot', defaults=dict(show=False))
        hist = self.get_hist(*args, **kwargs)
        if not len(hist):
            self.warning('No hist data')
            return
        plotkw.setdefault('title', 'Histogram')
#        vmin, vmax = MV2.min(hist), MV2.max(hist)
#        vamp = vmax - vmin
#        pmin, pmax = 0, vmax + (vamp / 10 or 1)
#        plotkw.setdefault('ylim', (pmin,pmax))
        plotkw.setdefault('xlabel', 'time')
        plotkw.setdefault('ylabel', 'number of profiles')
        self.info(self.describe(hist))
        bar2(hist, **plotkw)

    def plot_dist(self, polygons=None, drawpro='x', drawpol=True, drawcnt=True, **kwargs):
        '''
        Produce a map of profiles distribution.

        The polygons variable determine the plot mode:
            - if None, only profiles positions are drawn
            - otherwise, distribution is represented using colors for each specified polygon

        :Params:
            - **polygons**: list of distribution polygon (of type _geoslib.Polygon). Activate colored distribution if not None.
            - **drawpro**: symbol identifier of profiles (profiles not plotted if None)
            - **drawpol**: draw polygons contours if True
            - **drawcnt**: show textual profiles counts if True

        :Keyword arguments:
            - **plot_<keyword>**: are passed to the plot function :func:`~vacumm.misc.plot.map2`

        '''
        plotkw = kwfilter(kwargs, 'plot', defaults=dict(show=False))
        if not polygons: polygons = tuple()
        elif not isinstance(polygons, (list, tuple)): polygons = (polygons,)
        polygons = tuple(Polygon(z.boundary) for z in polygons)
        time = self.get_time(**kwargs)
        if time is None:
            self.error('No dist data')
            return
        npro = time.shape[0]
        time = create_time(time, units=time.units).asComponentTime()
        lat = self.get_latitude(**kwargs)
        lon = self.get_longitude(**kwargs)
        lon_min, lat_min = MV2.min(lon), MV2.min(lat)
        lon_max, lat_max = MV2.max(lon), MV2.max(lat)
        if npro or len(polygons):
            if len(polygons):
                lon_min = min(lon_min, *(min(*(c[0] for c in p.boundary)) for p in polygons))
                lat_min = min(lat_min, *(min(*(c[1] for c in p.boundary)) for p in polygons))
                lon_max = max(lon_max, *(max(*(c[0] for c in p.boundary)) for p in polygons))
                lat_max = max(lat_max, *(max(*(c[1] for c in p.boundary)) for p in polygons))
            lon_off, lat_off = abs(lon_max - lon_min) / 10., abs(lat_max - lat_min) / 10.
            plotkw.setdefault('lon_min', float(lon_min - lon_off))
            plotkw.setdefault('lat_min', float(lat_min - lat_off))
            plotkw.setdefault('lon_max', float(lon_max + lon_off))
            plotkw.setdefault('lat_max', float(lat_max + lat_off))
        plotkw.setdefault('title', '%s'%('Distribution'))
        show = plotkw.get('show', False)
        plotkw.update(dict(show=False, axes_rect=(0.075, 0.05, 0.75, 0.925)))
        m = map2(**plotkw)
        polyscounts = [0] * len(polygons)
        for ipro in xrange(npro):
            pt = Point((lon[ipro], lat[ipro]))
            for ipoly,poly in enumerate(polygons):
                if pt.within(poly):
                    polyscounts[ipoly] += 1
            if drawpro:
                x, y = m(lon[ipro], lat[ipro])
                #t = 'PRO:%s'%ipro
                t = drawpro%vars()
                pylab.gca().text(x, y, t, color='black', fontsize=7.75, horizontalalignment='center', verticalalignment='center')
        maxcount = len(polygons) and max(max(polyscounts), 10) or 10
        nsteps = 10
        step = maxcount/nsteps
        levels = range(0, maxcount+step, step)
        from matplotlib import mpl
        import matplotlib.colorbar
        colors = ['#0040FF', '#0080FF', '#00A2FF', '#00FFFB', '#00FF98', '#5DFF00', '#BFFF00', '#F1FF00', '#FFDD00', '#FFAC00', '#FF1800', '#CE1300']
        colors = colors[-nsteps:]
        cmap = mpl.colors.ListedColormap(colors)
        norm = mpl.colors.Normalize(vmin=0, vmax=maxcount)
        alpha = 0.95
        bbox_polygon = Polygon(numpy.array([
            [plotkw['lon_min'], plotkw['lat_min']], [plotkw['lon_max'], plotkw['lat_min']],
            [plotkw['lon_max'], plotkw['lat_max']], [plotkw['lon_min'], plotkw['lat_max']]]))
        for ipoly,poly in enumerate(polygons):
#            if not poly.within(bbox_polygon):
#                self.verbose('Ignoring polygon %s, not within bbox', ipoly)
#                continue
            self.verbose('%s profiles within polygon %s', polyscounts[ipoly], ipoly)
            coords = poly.get_coords()
            x, y = m(coords[:,0], coords[:,1])
            color = cmap(float(polyscounts[ipoly]) / maxcount)
            if drawpol:
                pylab.gca().plot(x, y, color='black', alpha=alpha)
            if drawcnt:
                pylab.gca().fill(x, y, color=color, alpha=alpha)
                pylab.gca().text(sum(x)/len(coords), sum(y)/len(coords),
                    #'POL:%d\nCNT:%s'%(ipoly, polyscounts[ipoly]),
                    '%s'%(polyscounts[ipoly]),
                    color='black', fontsize=7.75, horizontalalignment='center', verticalalignment='center')
        if polygons:
            fig = pylab.gcf()
            l,b,w,h = 0.85, 0.1, 0.05, 0.8
            ax = fig.add_axes((l,b,w,h))
            cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=levels, alpha=alpha)
            cb.set_label('Number of profiles')
        if show:
            pylab.show()

    def plot_pro(self, varname, *index, **kwargs):
        '''
        Produce a 1D plot of the profiles with depth as Y axis.

        :Params:
            - **varname**: variable name
            - **index**: profile(s) indexe(s) (after optionnal selection)

        :Keyword arguments:
            - **select**: selector
            - **plot_<keyword>**: are passed to the plot function :func:`~vacumm.misc.plot.curve2`

        '''
        kwargs.setdefault('select', {})
        plotkw = kwfilter(kwargs, 'plot')
        #show = plotkw.get('plot', False)
        curves = [None] * len(index)
        for icur,ipro in enumerate(index):
            kwargs['select']['profile'] = slice(ipro,ipro+1)
            var = self.get_variable(varname, **kwargs)
            if var is None or not len(var):
                self.error('Variable %s not found or index %s out of range', varname, ipro)
                continue
            var = var[0]
            platcode = self.get_platform_code(**kwargs)
            platcode = ''.join(platcode[0]).strip() if platcode is not None else ''
            time = self.get_time(**kwargs)
            time = create_time(time, units=time.units).asComponentTime()[0]
            dep = self.get_depth(**kwargs)[0]
            lat = self.get_latitude(**kwargs)[0]
            lon = self.get_longitude(**kwargs)[0]
            v, d = var, dep
            # Exclude masked data based on the depth variable
            if d.mask.size > 1:
                orgshape = d.shape
                valid_indices = numpy.where(~d.mask)
                if not len(valid_indices) or not len(valid_indices[0]):
                    self.warning('No valid data for variable %s profile at index %s', var.id, ipro)
                    continue
                valid_indices = valid_indices[0]
                v, d = v.take(valid_indices), d.take(valid_indices)
            # Reverse depth
            if d[0] > d[-1]:
                self.warning('Reversing data and level to go from bottom to surface')
                v, d = v[::-1], d[::-1]
            d = d.tolist() # tmp fix
            self.describe(d)
            #l = cdms2.createAxis(range(d.shape[0]), id='level')
            l = cdms2.createAxis(d, id='level')
            l.designateLevel()
            v = cdms2.createVariable(v, id=var.id, axes=[l], attributes=var.attributes.copy())
            d = cdms2.createVariable(d, id=dep.id, axes=[l], attributes=dep.attributes.copy())
            self.describe((v, d, l))
            curkw = plotkw.copy()
            curkw.setdefault('order', 'zd')
            curkw.setdefault('title', 'Profile %s'%(var.id))
            curkw.setdefault('label', 'i: %(ipro)s, c: %(platcode)s, t: %Y-%m-%d %Hh, x: %(lon).3f, y: %(lat).3f')
            if curkw['label']: curkw['label'] = adatetime(time).strftime(curkw['label'])%vars()
            curkw['show'] = False
            curves[icur] = curve2(v, **curkw)
        #if show: pylab.show()
        return curves

    @classmethod
    def main(cls, *args, **kwargs):
        '''
        Entry point of the bin/profile.py script. See profile.py - -help.

        :Params:
            - **args** and **kwargs**: passed to optparse.OptionParser.parse_args

        .. note::
            - You may pass the "args" named argument to override command line arguments

        .. todo:
            - Use config.ConfigManager instead of optparse.OptionParser

        :Examples:

            >>> ProfilesDataset.main(args=('myprofiles.nc', '--pro', '0,1,2'))
        '''
        try:
            from optparse import IndentedHelpFormatter
            class IndentedHelpFormatterWithNL(IndentedHelpFormatter):
                def format_description(self, description):
                    if not description: return ""
                    return description
            parser = optparse.OptionParser(
                description=
                    '''\
Show profiles informations: distribution, historgram, profile variable data

Input file[s] must be standardized profiles NetCDF file produced by merge_profiles.py or by this
program convertion utility (see --convert option).

Examples:
  This will produce a map with the profiles position for the given variable:
    %prog profiles.nc -v temperature --dist
  This will produce a map with the profiles distribution for the given variable and set of polygons:
    %prog profiles.nc -v temperature -t 2000,2005 --hist --dist -P data/polygon_manga/polygon_manga.txt
                    ''',
                usage='%prog [options] [file1] [file2] [...]',
                version=__date__,
                formatter=IndentedHelpFormatterWithNL())

            parser.add_option('--cfgfile', action='store', dest='cfgfile', default=None, help='Configuration file')
            parser.add_option('--hist', action='store', dest='hist', default=None, metavar='n,units', help='Plot profiles histogram with interval n,units (ex: "1,month" , "2,weeks" , ...)')
            parser.add_option('--dist', action='store_true', dest='dist', default=None, help='Plot profiles distribution (position and, if polygons are specified, the colored count of profiles per polygon)')
            parser.add_option('--drawpro', action='store_true', dest='drawpro', default=None, help='Plot profiles position in distribution plot when using polygons')
            parser.add_option('--pro', action='store', dest='pro', default=None, help='Plot a profile at specified comma separated indexes. Ranges may also be specified using start:stop[:step] notation, resulting in [start,...,stop[ indexes. Example: 0,1,2:5,5:11:1 would plot profiles 0 to 10')
            parser.add_option('-t', '--time', action='store', dest='time', default=None, help='Restrict timerange to be processed, format is "datemin,datemax" with date format: "YYYY-mm-ddTHH:MM:SS"')
            parser.add_option('-b', '--bbox', action='store', dest='bbox', default=None, help='Set the map bounding box, format is lon_min,lat_min,lon_max,lat_max')
            parser.add_option('-v', '--variables', action='append', dest='variables', default=[], help='Variables to be processed, Aliases may be specified using comma separated names. Repeat this option for each required variable. Ex: -v temperature,TEMP -v salinity,PSAL')
            parser.add_option('-p', '--poly', action='append', dest='poly', default=[], help='Add a polygon for profile restriction and accounting, polygon format is x1,y1,x2,y2,x3,y3,x4,y4... At least 4 couples must be set, otherwise if 4 values are set, it is considered as lonmin,latmin,lonmax,latmax')
            parser.add_option('--polyfile', action='append', dest='polyfile', default=[], help='Specify a polygon file containing couples coordinates (in lat, lon order), with values separated by a whitespace or coma')
            parser.add_option('-P', '--polysfile', action='append', dest='polysfile', default=[], help='Specify a polygons file (one per line)')
            parser.add_option('-o', '--output', action='store', dest='output', default='%(plot)s.png', metavar='pattern',  help='Output files pattern (default: "%default" with plot replaced by the %(plot)s type)')
            parser.add_option('--convert', action='store', dest='convert', default=None, metavar='converted_profiles.nc', help='Convert given profile file(s) into a standard netcdf output file. Variables must be specified using the -v option.')
            parser.add_option('--save', action='store_true', dest='save', default=None, help='Save figures, you may also use the -o option')
            parser.add_option('--show', action='store_true', dest='show', default=None, help='Show figures (enabled by default if --save is not used)')
            parser.add_option('--netcdf3', action='store_true', dest='netcdf3', default=False,  help='Save files in NetCDF version 3 (only used with --convert)')
            parser.add_option('--compress', action='store', dest='compress', default=3,  type=int, help='Set the compression level for NetCDF version 4 (default: %default, range [0,9], only used with --convert)')
            parser.add_option('--nonumpywarnings', action='store_true', dest='nonumpywarnings', default=False,  help='Hide all numpy warnings')
            parser.add_option('--debug', action='store_true', dest='debug', default=False,  help='Debug mode')

            #Logger.add_optparser_options(parser)
            options, args = parser.parse_args(*args, **kwargs)
            #Logger.apply_class_optparser_options(options)

            # Operations on masked values throws numpy warnings, this can be ignored with the following instruction
            if options.nonumpywarnings:
                numpy.seterr(all='ignore')
            if options.cfgfile:
                Object.load_default_config(options.cfgfile, nested=True)
                cls.load_default_config(Object.get_default_config())
            if options.debug:
                for c in [Object,ProfilesDataset]+Profiles.get_types().values():
                    c.get_logger().set_level('debug')
                    c._log_obj_stats = True
            if not options.save:
                options.show = True
            select = {}
            # format variables to a list of variables candidates
            # ex: [['temperature','TEMP'], ['salinity','PSAL']]
            options.variables = map(lambda v: v.split(','), options.variables)
            if options.hist:
                options.hist = options.hist.split(',')
                options.hist[0] = int(options.hist[0])
            if options.time:
                select.update(time=options.time.split(','))
            if options.bbox:
                b = map(float, options.bbox.split(','))
                select.update(lat=(b[1],b[3]), lon=(b[0],b[2]))
            polys = []
            polys.extend(options.poly)
            for fp in options.polyfile:
                f = file(fp)
                polys.append(' '.join(l.strip() for l in f if l.strip() and not l.strip().startswith('#')))
                f.close()
            for fp in options.polysfile:
                f = file(fp)
                polys.extend(l.strip() for l in f if l.strip() and not l.strip().startswith('#'))
                f.close()
            if polys:
                select['polygons'] = []
                for ipoly,poly in enumerate(polys):
                    c = []
                    # two possible delimiters: whitespaces and ','
                    for s in poly.split(): c.extend(s.split(','))
                    c = [v.strip() for v in c if v.strip()]
                    c = map(float, c)
                    if len(c) == 4:
                        x1, y1, x2, y2 = c
                        p = Polygon(numpy.array([[x1, y1], [x2, y1], [x2, y2], [x1, y2]]))
                    elif len(c) >= 8 and len(c)%2 == 0:
                        p = Polygon(numpy.array(zip(c[::2],c[1::2])))
                    else:
                        cls.error('Invalid format in polygon %s', ipoly)
                        continue
                    select['polygons'].append(p)
                    cls.debug('Loaded polygon %s with %s coordinates'%(ipoly, len(select['polygons'][-1].boundary)))

            plot_kwargs = dict(select=select)

            if options.bbox:
                b = map(float, options.bbox.split(','))
                plot_kwargs.update(dict(plot_lon_min=b[0], plot_lat_min=b[1], plot_lon_max=b[2], plot_lat_max=b[3]))

            if options.convert:
                if options.netcdf3:
                    cls.notice('Using on NetCDF version 3 writting mode')
                    netcdf3()
                elif options.compress:
                    cls.notice('Using on NetCDF version 4 writting mode with compression level %s', options.compress)
                    netcdf4(options.compress)
                vmap = dict()
                for v in options.variables:
                    vmap[v[0]] = len(v) > 1 and v[1:] or v[0]
                profiles = Profiles(variables_map=vmap, variables=vmap.keys(), qualities=())
                for f in args: profiles.load(f)
                profiles.save(options.convert)
                if not options.variables:
                    cls.warning('No variable were specifed')
                return 0

            profiles = cls()
            if options.cfgfile:
                profiles.load_default_config(options.cfgfile, nested=True)
            profiles.load_dataset(args, append=True)

            if options.hist:
                hist_kwargs = plot_kwargs.copy()
                pylab.clf()
                cls.notice('Plotting histogram',)
                profiles.plot_hist(options.hist, **hist_kwargs)
                if options.save:
                    output = options.output%dict(plot='hist')
                    cls.notice('Saving %s', output)
                    pylab.savefig(output)
                if options.show: pylab.show()

            if options.dist:
                dist_kwargs = plot_kwargs.copy()
                dist_kwargs['drawpro'] = 'x' if not select.get('polygons', None) or options.drawpro else None
                dist_kwargs['drawcnt'] = True
                dist_kwargs['drawpol'] = True
                pylab.clf()
                cls.notice('Plotting distribution')
                profiles.plot_dist(polygons=select.get('polygons', None), **dist_kwargs)
                if options.save:
                    output = options.output%dict(plot='dist')
                    cls.notice('Saving %s', output)
                    pylab.savefig(output)
                if options.show: pylab.show()

            if options.pro:
                pro_kwargs = plot_kwargs.copy()
                pro_kwargs['plot_axes_rect'] = [0.1,0.4,0.8,0.5]
                indexes = []
                try:
                    for idx in options.pro.split(','):
                        if ':' in idx: indexes.extend(range(*map(lambda i: int(i) if len(i) else None, idx.split(':'))))
                        else: indexes.append(int(idx))
                except Exception,e: cls.error('Invalid profiles indexes format: %s: %s'%(options.pro, e))
                if not options.variables: cls.error('At least one variable option is required for profile plot')
                for varname in options.variables:
                    pylab.clf()
                    cls.notice('Plotting profiles of variable %s', varname)
                    curves = profiles.plot_pro(varname, *indexes, **pro_kwargs)
                    #pylab.legend(loc='best')
                    pylab.legend(bbox_to_anchor=(1.0,-0.15))
                    if options.save:
                        output = options.output%dict(plot='pro_'+varname)
                        cls.notice('Saving %s', output)
                        pylab.savefig(output)
                    if options.show: pylab.show()

            return 0
        except Exception, e:
            cls.exception('Fatal error: %s', e)
            return 1

    def get_stratification_data(self, select):
        '''Get stratification data

        :Params:
            - **select**: selector

        :Return:
            - time, lat, lon, depth, depthimax, deltadepth, temp, sal, dens, densmin, densmax, densmean with shape (profile,)

        '''
        prof = self.get_axis('profile', select=select)
        lev = self.get_axis('level', select=select)
        time = self.get_variable('time', select=select)
        lon = self.get_variable('longitude', select=select)
        lat = self.get_variable('latitude', select=select)
        depth = self.get_variable('depth', select=select)
        temp = self.get_variable('temperature', select=select)
        sal = self.get_variable('salinity', select=select)
        self.verbose('Initial stratification data:\n  %s', '\n  '.join(self.describe(o) for o in (prof,lev,time,lon,lat,depth,temp,sal)))

        # Set intial order for masked profiles exclusion
        temp, sal, depth = temp.reorder('-z'), sal.reorder('-z'), depth.reorder('-z')

        # Exclude entirely masked depth/temp/sal coordinates
        goodfct = lambda a: numpy.ma.where(numpy.ma.count(a, axis=1)>0, True, False)
        depgood, tempgood, salgood = goodfct(depth), goodfct(temp), goodfct(sal)
        good = depgood & tempgood & salgood
        depcnt, tempcnt, salcnt = len(depgood.nonzero()[0]), len(tempgood.nonzero()[0]), len(salgood.nonzero()[0])
        totcnt = len(good.nonzero()[0])
        self.debug('Number of profiles with valid data:'
            '\n  depth:        %s/%s'
            '\n  temperature:  %s/%s'
            '\n  salinity:     %s/%s'
            '\n  total:        %s/%s', depcnt, len(depth), tempcnt, len(depth), salcnt, len(depth), totcnt, len(depth))
        if not totcnt:
            self.error('No valid profiles, aborting')
            return
        igood = numpy.where(good)[0]
        takefct = lambda a: MV2.take(a, igood, axis=0)
        prof = takefct(prof)
        time, lat, lon = takefct(time), takefct(lat), takefct(lon)
        depth, temp, sal = takefct(depth), takefct(temp), takefct(sal)

        # Restore axes as tz as loaded above
        prof, lev = cdms2.createAxis(prof, id='profile'), cdms2.createAxis(lev, id='level')
        for v in time, lat, lon: v.setAxisList([prof])
        for v in depth, temp, sal: v.setAxisList([prof, lev])
        # Set correct order for calculation
        depth, temp, sal = depth.reorder('z-'), temp.reorder('z-'), sal.reorder('z-')

        self.debug('Selected data:\n  %s', '\n  '.join(self.describe(o) for o in (prof,lev,time,lon,lat,depth,temp,sal)))

        # Maximum depth indice
        depthimax = MV2.array((~temp.mask).sum(axis=0)-1, id='depthimax')
        # Compute thick for each level
        deltadepth = MV2.array(meshweights(depth, axis=0), id='deltadepth')
        # Adjust thick of last level
        #####
        try: ddiff = numpy.diff(depth, axis=0)
        except ValueError: # Error when depth shape is (.., 1) (one profile) ...???
            #self.exception('')
            ddiff = numpy.diff(numpy.array(depth), axis=0)
            #ddiff = numpy.swapaxes(ddiff, 0, 1)
        #####
        self.debug('ddiff: %s', self.describe(ddiff))
        for ip, k in numpy.ndenumerate(depthimax):
            deltadepth[k, ip[0]] = ddiff[k-1, ip[0]]

        # Compute pressure
        pres = MV2.array(sw_pres(depth, numpy.resize(lat, depth.shape)), id='pres')
        # Compute pressure density
        dens = MV2.array(sw_dens(sal, temp, pres), id='dens')
        # Compute pressure density min, max and mean
        densmin = MV2.min(dens, axis=0)
        densmin.id = 'densmin'
        densmax = MV2.max(dens, axis=0)
        densmax.id = 'densmax'
        densmean = MV2.average(dens, axis=0, weights=deltadepth)
        densmean.id = 'densmean'

        ret = time, lat, lon, depth, depthimax, deltadepth, temp, sal, dens, densmin, densmax, densmean
        self.verbose('Final stratification data:\n  %s', '\n  '.join(self.describe(o) for o in ret))
        return ret


    def get_mld(self, select):
        '''Get mixed layer depth

        MLD is computed for each selected profiles

        :Params:
            - **select**: selector

        :Return:
            - mld with shape (profile,)
            - lat with shape (profile,)
            - lon with shape (profile,)

        '''
        time, lat, lon, depth, depthimax, deltadepth, temp, sal, dens, densmin, densmax, densmean = self.get_stratification_data(select)
        # Maximum depth (max + thick * 0.5)
        H = numpy.array([depth[depthimax[i]][i] + deltadepth[depthimax[i]][i] * 0.5 for i in xrange(depth.shape[1])])
        self.debug('H:   %s', self.describe(H))
        # Compute Mixed Layer Depth
        mld = H*(densmax - densmean) / (densmax - densmin)
        mld.id = 'MLD'
        mld.units = 'm'
        mld.long_name = u'Profondeur de la couche de melange'
        self.verbose('MLD:   %s', self.describe(mld))
        return mld, lat, lon


    def get_ped(self, select):
        '''Get potential energy deficit

        PED is computed for each selected profiles

        :Params:
            - **select**: selector

        :Return:
            - ped with shape (profile,)
            - lat with shape (profile,)
            - lon with shape (profile,)

        '''
        time, lat, lon, depth, depthimax, deltadepth, temp, sal, dens, densmin, densmax, densmean = self.get_stratification_data(select)
        from vacumm.misc.phys.constants import g
        # Anomalie de densité
        danom = dens - densmean
        # Énergie potentielle disponible
        ape = danom * g
        ape *= deltadepth
        # Deficit
        ped = MV2.average(ape, axis=0, weights=deltadepth)
        ped.id = 'PED'
        ped.units = 'J.m^{-2}'
        ped.long_name = u"Potential energy deficit"
        self.verbose('PED:   %s', self.describe(ped))
        return ped, lat, lon



    def plot_mld(self, select, **kwargs):
        '''
        Produce mixed layer depth map.

        :Params: see :func:`get_mld`

        :Plot params:
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations

        '''
        mld, lat, lon = self.get_mld(select)
        # Plot
        #self.psinfo()
        vmin, vmax = numpy.min(mld), numpy.max(mld)
        m = kwargs.pop('map', None)
        c = kwargs.pop('cmap', 'magic')
        # Plot the map if not provided
        if m is None:
#            from matplotlib import mpl
#            import matplotlib.colorbar
#            nsteps = 10
#            step = (vmax-vmin)/nsteps
#            levels = range(0, vmax+step, step)
#            colors = ['#0040FF', '#0080FF', '#00A2FF', '#00FFFB', '#00FF98', '#5DFF00', '#BFFF00', '#F1FF00', '#FFDD00', '#FFAC00', '#FF1800', '#CE1300']
#            colors = colors[-nsteps:]
#            cmap = mpl.colors.ListedColormap(colors)
#            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#            alpha = 0.8
#            fig = pylab.gcf()
#            l,b,w,h = 0.85, 0.1, 0.05, 0.8
#            ax = fig.add_axes((l,b,w,h))
#            cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=levels, alpha=alpha)
#            cb.set_label('')

            mapkw = kwfilter(kwargs, 'map', defaults=dict(
            lon=(min(lon)-0.5,max(lon)+0.5), lat=(min(lat)-0.5,max(lat)+0.5),
            vmin=vmin, vmax=vmax, cmap=c))
            mapkw.update(show=False)
            m = map2(**mapkw)

        # Plot profiles scatter points
        scakw = kwfilter(kwargs, 'sca', defaults=dict(s=20, vmin=vmin, vmax=vmax, cmap=m.cmap))
        sc = m.map.scatter(lon, lat, c=mld, **scakw)
        # Post plotting
        plotkw = kwfilter(kwargs, 'plot')
        #plotkw.setdefault('title', '')
        m.post_plot(**plotkw)


    def plot_ped(self, select, **kwargs):
        '''
        Produce potential energy deficit map.

        :Params: see :func:`get_ped`

        :Plot params:
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations

        '''
        ped, lat, lon = self.get_ped(select)
        # Plot
        #self.psinfo()
        vmin, vmax = numpy.min(ped), numpy.max(ped)
        m = kwargs.pop('map', None)
        c = kwargs.pop('cmap', 'magic')
        # Plot the map if not provided
        if m is None:
#            from matplotlib import mpl
#            import matplotlib.colorbar
#            nsteps = 10
#            step = (vmax-vmin)/nsteps
#            levels = range(0, vmax+step, step)
#            colors = ['#0040FF', '#0080FF', '#00A2FF', '#00FFFB', '#00FF98', '#5DFF00', '#BFFF00', '#F1FF00', '#FFDD00', '#FFAC00', '#FF1800', '#CE1300']
#            colors = colors[-nsteps:]
#            cmap = mpl.colors.ListedColormap(colors)
#            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#            alpha = 0.8
#            fig = pylab.gcf()
#            l,b,w,h = 0.85, 0.1, 0.05, 0.8
#            ax = fig.add_axes((l,b,w,h))
#            cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=levels, alpha=alpha)
#            cb.set_label('')

            mapkw = kwfilter(kwargs, 'map', defaults=dict(
            lon=(min(lon)-0.5,max(lon)+0.5), lat=(min(lat)-0.5,max(lat)+0.5),
            vmin=vmin, vmax=vmax, cmap=c))
            mapkw.update(show=False)
            m = map2(**mapkw)

        # Plot profiles scatter points
        scakw = kwfilter(kwargs, 'sca', defaults=dict(s=20, vmin=vmin, vmax=vmax, cmap=m.cmap))
        sc = m.map.scatter(lon, lat, c=ped, **scakw)
        # Post plotting
        plotkw = kwfilter(kwargs, 'plot')
        #plotkw.setdefault('title', '')
        m.post_plot(**plotkw)




if __name__ == '__main__':
    Profiles.main()




