# -*- coding: utf8 -*-
"""
Provides a class to perform basic operations on marigraphic sea level data
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
import numpy as N,numpy.ma as MA,cdms2 as cdms,MV2 as MV,cdtime
from genutil import minmax
from genutil.statistics import std as STD

import pylab as P,types
import scipy.interpolate as S
import types
from operator import isNumberType
from warnings import warn

__all__ = ['Marigraph']

_docs = dict(
    time_range= "time range for data selection (like ('1999','2000-01-01','co')). If None, time range is reset [default: None]",
    tide_filter=    "Tide ilter name within [False,'demerliac','godin']. Default: False"
)


for key, val in _docs.items():
    ss = '*'*(1+val.startswith('+'))
    _docs[key] = '- %(ss)s%(key)s%(ss)s: %(val)s' % dict(ss=ss, key=key, val=val)
_docs['base'] = '\n\t\t\t'.join([_docs[key] for key in 'time_range', 'tide_filter'])

def _fmtdoc_(func):
#   try:
#       func.__doc__ = func.__doc__ % _docs
#   except:
#       pass
    return func



class Marigraph(object):
    """A marigraph class that support statistics and plots.

    You give as input a 1D cdms sea level variable with a suitable time,
    or a tide.StationInfo() object, or simply a 5-letter shom ID.
    The class will design an object with high and low tides,
    and pure tidal and without tidal signal.
    You can optionnally specify a time axis on wich the input signal
    is to be interpolated before any processing.
    You can also restrict the time range and perfom a running average.

    :Parameters:

        All parameters are passed to :meth:`set`

            - **data** : It can be either

                - A ``cdms2`` 1D array with a time axis.
                - A :class:`~vacumm.tide.StationInfo` instance.
                - The name of a station that can be search for using :class:`~vacumm.tide.StationInfo`.

            - *select*: time range for data selection (like ('1999','2000-01-01','co')). If None, time range is reset [default: None],
            - *global_anomaly*: First remove the anomaly of the whole signal before any processing [default: False]
            - *anomaly*: Remove the anomaly of the restricted sample before any processing [default: False]
            - *shom_anomaly*: Same but removed mean is get from observation at correspondind shom station [default: False]
%s

    :Example:

    >>> # Get data
    >>> from vacumm.tide import Marigraph
    >>> brest = Marigraph('BREST',range=('2006-03','2006-08'),hourly=True)
    >>> # Plot
    >>> brest.plot_sea_level()
    >>> brest.plot_low_high()
    >>> brest.plot_cotes()
    >>> # Change time range
    >>> brest.set_time_range(('2006-07','2006-08','co'))
    >>> # Change tide filter
    >>> brest.set_tide_filter('godin')
    >>> # Save
    >>> f = cdms.open('mytide.nc')
    >>> f.write(brest.cotes())
    >>> f.write(brest.high())
    >>> f.write(brest.low())
    >>> f.write(brest.zeros())
    >>> f.close()
    """

    # Initialization
    # --------------
    def __init__(self,data, tide_filter='demerliac',
        outside_std=3.5,bad_values=None,verbose=False, mean=None, **kwargs):

        # Data synchronisation
        self._data_funcs = ['_check_tide_filter_','_check_low_high_',]

        """Assign a data object to the marigraph

        :Parameters:

            - **data** : It can be either

                - A ``cdms2`` 1D array with a time axis.
                - A :class:`~vacumm.tide.StationInfo` instance.
                - The name of a station that can be search for using :class:`~vacumm.tide.StationInfo`.

            - *global_anomaly*: First remove the anomaly of the whole signal before any processing [default: False]
            - *anomaly*: Remove the anomaly of the restricted sample before any processing [default: False]
            - *shom_anomaly*: Same but removed mean is get from observation at correspondind shom station [default: False]
%s
        """
        self.clean()
        self._verbose = verbose
        if isinstance(data,StationInfo) or \
               (isinstance(data, str) and len(data) == 5):
            # Data is not a variable, but a StationInfo object or a string = shom id
            if isinstance(data, str):
                # We start from a shom id
                if verbose:
                    print 'Searching for SHOM station "'+data+'"...'
                    print '-'*80
                station = StationInfo(shom=data,**self._clean_kwargs(kwargs))
                if station is None:
                    print 'Not found'
                    raise StandardError
                if verbose: print '-'*80
                shom = data
                self.name = station.name
                self.longitude = station.longitude
                self.latitude = station.latitude
            else:
                # We already have a station info
                shom = data.shom
                if shom is None:
                    raise StandardError, 'Not a valid SHOM station'
                if data.nom is not None:
                    self.name = data.nom
                    self.longitude = data.longitude
                    self.latitude = data.latitude

            from vacumm.tide.sonel_mareg import get_slv_shom
            data = get_slv_shom(shom,**self._clean_kwargs(kwargs,bad=['tide_filter']))
            if mean == 'shom':
                mean = station.nm

        else:
            # Here we have cdms variable, so we just check it
            assert cdms.isVariable(data), \
                   'The marigraph data object is not a valid cdms variable'
            data = data.clone()

            # Check time axis
            taxis = data.getTime()
            if taxis is None:
                taxis = data.getAxis(0).designateTime()
            assert hasattr(taxis, 'units'), 'Time axis has no units'

            # Must be a vector
            if data.ndim > 1:
                # Reorder
                if not data.getOrder().endswith('t'):
                    data = data.reorder('...t')
                # Reshape to average
                tmp = MV2.reshape(data, (N.multiply.reduce(data.shape[:-1]), data.shape[-1]))
                data  = MV2.average(tmp, axis=0)
                del temp

            # Get some attributes
            for att in ['name','longitude','latitude']:
                latt = 'station_'+att
                if data.attributes.has_key(latt):
                    setattr(self,att,getattr(data,latt))

        # Remove extrem values
        if outside_std not in [False,None,0.]:
            std = data.std()
            if mean is not None:
                tmpmean = mean
            else:
                tmpmean = float(data.mean())
            data[:] = MV.masked_outside(data,tmpmean-outside_std*std,   tmpmean+outside_std*std)

        # Attributes
        self.data['base'] = data
        if not self.data['base'].attributes.has_key('units'):
            self.data['base'].units = 'm'
        if not self.data['base'].attributes.has_key('long_name'):
            self.data['base'].units = 'Sea level'

        # Real mean
        if mean is None:
            mean = data.mean()
        self._mean = float(mean)
        self.data['anom'] = data
        self.data['anom'][:] -= mean
        self.data['anom'].id += '_anom'
        self.data['anom'].long_name = 'Anomaly of '+self.data['anom'].long_name

        self.tide_filter = tide_filter




    # Apply tide filters
    # ------------------
    def mean(self):
        """Mean sea level"""
        return self._mean

    def anomaly(self):
        """Sea level anomaly"""
        return self._data['anom']

    @_fmtdoc_
    def set_tide_filter(self,tide_filter):
        """Set the tide filter

        :Parameters:
            %(tide_filter)s
        """
        self._check_tide_filter_(tide_filter)

    def _check_tide_filter_(self, tide_filter=None, **kwargs):

        # Do we have something?
        need = self.data['tide'] is None

        # A filter is specified
        if tide_filter is not None:
            if tide_filter is True:
                # Default filter
                tide_filter = 'demerliac'
            else:
                # Valid filter?
                valid_filters = ['demerliac','godin',False]
                assert tide_filter in valid_filters,  "Wrong filter name. Must be within %s. Please use set_tide_filter()) method."%str(valid_filters)
            # Change?
            need = need or (self.tide_filter != tide_filter)

        # Something changed?
        if not need: return
        if self.tide_filter is False:
            # Reset
            self.data['tide'] = self.data['base']
            self.data['cotes'] = self.data['base'].clone()
            self.data['cotes'][:] = 0.
        else:
            # Apply filter
            filter_func = getattr(filters, tide_filter)
            self.data['cotes'],self.data['tide'] = filter_func(self.data['anom'],get_tide=True)
            if self._verbose:
                print 'Computed tide signal using %s filter' % tide_filter
        self.tide_filter = tide_filter


    # Low and high tides
    # ------------------
    def _check_low_high_(self, ref, **kwargs):
        """Computes high and low tides from a sea level anomaly. It defines following variables: high,low
        """
        # Check changes
        if self.data['lows'] is not None: return

        # Reference
        self._lowhigh_ref = ref
        if ref in ['demerliac', 'godin']:
            if self.tide_filter != ref:
                ref, tide = getattr(filters, ref)(self.data['base'])
            else:
                ref = self.data['cotes']

        # Call to extema
        self.data['lows'], self.data['highs'] = filters.extrema(self.data['anom'], ref=ref, **kwargs)
        self.data['zeros'] = filters.zeros(self.data['anom'], ref=ref, **kwargs)
        for ex in 'lows', 'highs', 'zeros':
            self.data[ex][:] += self._mean

        if self._verbose:
            print 'Computed low and high tides'



    # Check changes
    # -------------


    #### Get data
    def sea_level(self):
        """Sea level"""
        return self.data['base']

    @_fmtdoc_
    def tide(self, anom=True, **kwargs):
        """Return the tide signal: this signal is obtained by filtering the original sea level signal
          using a demerliac or a godin filer.

        %(ref)s

        .. seealso:

            :meth:`set_tide_filter`, :meth:`cotes`, :mod:`vacumm.tide.filters`
        """

        self._check_tide_filter_(**kwargs)
        data = self.data['tide']
        if not anom:
            data[:] += self._mean
        return data

    @_fmtdoc_
    def cotes(self, anom=True, **kwargs):
        """Return the sea level without the tide signal (surcotes and decotes). You must specify the filter before using this method.

        %(ref)s

        .. seealso:

            :meth:`set_tide_filter`, :meth:`tide`, :mod:`vacumm.tide.filters`
        """
        self._check_tide_filter_(**kwargs)
        data = self.data['cotes']
        if not anom:
            data[:] += self._mean
        return data

    @_fmtdoc_
    def lows(self, ref='mean', **kwargs):
        """Low tides

        %(ref)s

        .. seealso:
            :meth:`high`, :meth:`zeros`, :mod:`vacumm.tide.filters`
        """
        self._check_low_high_(ref, **kwargs)
        return self.data['lows']

    @_fmtdoc_
    def highs(self, ref='mean',**kwargs):
        """High tides

        %(ref)s

        .. seealso:
            :meth:`low`, :meth:`zeros`, :mod:`vacumm.tide.filters`
        """
        self._check_low_high_(ref, **kwargs)
        return self.data['highs']


    @_fmtdoc_
    def zeros(self, ref='mean',**kwargs):
        """Zeros of sea level anomaly

        %(ref)s

        .. seealso:
            :meth:`low`, :meth:`high`, :mod:`vacumm.tide.filters`
        """
        self._check_low_high_(ref, **kwargs)
        return self.data['zeros']


    #### Misc methods

    def clean(self):
        """ Clean out the marigraph object """
        # Data
        self.data = {}
        for d in 'base', 'tide', 'cotes', 'lows', 'highs', 'zeros':
            self.data[d] = None

        # Filters
        self.tide_filter = False

        # Misc
        self.name = None
        self.longitude = None
        self.latitude = None
        self.mean_sea_level = -999.
        self.global_mean_sea_level = -999.

    def _clean_kwargs(self,kwargs,bad = ['tide_filter']):
        kwargs = kwargs.copy()
        for kw in bad:
            if kwargs.has_key(kw):
                del kwargs[kw]
        return kwargs

    def name(self,name=None):
        """ Set or return the station name """
        if name is None:
            return self.name
        else:
            self.name = name

    def longitude(self,val=None):
        """ Set or return the station longitude """
        if val is None:
            return self.longitude
        else:
            self.longitude = val

    def latitude(self,val=None):
        """ Set or return the station latitude """
        if val is None:
            return self.latitude
        else:
            self.latitude = val



    #### Plots

    def plot(self, var=None, orig=True, tide=True, cotes=True, highs=True, lows=True, zeros=True, legend=True, title=None, savefig=None, savefigs=None, show=True, marker='o', **kwargs):

        # Which variable?
        vtypes = 'orig', 'tide', 'highs', 'lows', 'zeros', 'cotes'
        if var is not None:
            assert var in vtypes
            for vt in vtypes:
                exec "%s = %s"%(vt, vt==var)

        # Plot params
        # - specific
        kwplot = {}
        for vt in vtypes:
            kwplot[vt] = kwfilter(kwargs, vt)
            kwplot[vt].update(show=False)
            kwplot[vt].update(title=False)
        # - global
        nt = self.data['base'].shape[0]
        times = T.mpl(self.data['base'][0:nt:nt-1].getTime())
        kwargs.setdefault('xmin', times[0])
        kwargs.setdefault('xmax', times[1])
        # - complete
        for vt in vtypes:
            tmp = kwargs.copy()
            tmp.update(kwplot[vt])
            kwplot[vt] = tmp
        anom = not (orig or highs or lows or zeros)

        # Original signal
        if orig:
            kwplot['orig'].setdefault('color', '#888888')
            kwplot['orig'].setdefault('zorder', '10')
            curve(self.sea_level(), **kwplot['orig'])

        # High tides
        if highs:
            kwplot['highs'].setdefault('color', 'r')
            kwplot['highs'].setdefault('linewidth', 0)
            kwplot['highs'].setdefault('markersize', 5)
            kwplot['highs'].setdefault('zorder', '12')
            curve(self.highs(), marker, **kwplot['highs'])

        # Low tides
        if lows:
            kwplot['lows'].setdefault('color', 'b')
            kwplot['lows'].setdefault('linewidth', 0)
            kwplot['lows'].setdefault('markersize', 5)
            kwplot['lows'].setdefault('zorder', '12')
            curve(self.lows(), marker, **kwplot['lows'])

        # Zeros
        if zeros:
            kwplot['zeros'].setdefault('color', 'k')
            kwplot['zeros'].setdefault('linewidth', 0)
            kwplot['zeros'].setdefault('markersize', 5)
            kwplot['zeros'].setdefault('zorder', '12')
            curve(self.zeros(), marker, **kwplot['zeros'])

        # Tidal signal
        if tide:
            kwplot['orig'].setdefault('color', 'k')
            kwplot['orig'].setdefault('zorder', '11')
            curve(self.tide(anom=anom), **kwplot['tide'])

        # Surcote/decotes
        if cotes:
            kwplot['cotes'].setdefault('color', 'g')
            kwplot['cotes'].setdefault('zorder', '12')
            curve(self.cotes(anom=anom), **kwplot['cotes'])

        # Misc
        if title is None:
            if var is not None:
                title = var.title()
            else:
                title = 'Marigraph'
        if title:
            P.title(title)
        if savefig:
            P.savefig(savefig,  **kwfilter(kwargs, 'savefig'))
        if savefigs:
            Savefigs(savefigs,  **kwfilter(kwargs, 'savefigs'))
        if show:
            P.show()







from station_info import StationInfo
from vacumm.misc.grid import bounds1d
import vacumm.misc.atime as T, filters
from vacumm.misc import kwfilter
from vacumm.misc.plot import curve2 as curve, savefigs as Savefigs

