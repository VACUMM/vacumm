# -*- coding: utf8 -*-
"""Utilities to handle shorelines"""
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

import numpy as N

from vacumm import VACUMMError, VACUMM_CFG
from vacumm.config import get_data_file
from vacumm.misc.grid import get_xy
from vacumm.misc.basemap import gshhs_reslist, gshhs_autores
from vacumm.misc.core_plot import register_plot_coastline
from vacumm.misc.sdata import Shapes, GSHHSBM

_shoreline_list = ['Histolitt', 'GSHHS_SF', 'GSHHS','GSHHSBM']
__all__ = _shoreline_list+['ShoreLine', 'get_best', 'get_shoreline',
                           'list_shorelines', 'get_best_res',
                           'Histolitt', 'GSHHS', 'GSHHS_SF', 'GSHHSSF',
                           'ShoreLineError', 'get_best_shoreline',
                           'get_bestres',]



class ShoreLineError(VACUMMError):
    pass

def _shorelines_list_(name=None):
    # Get list a dict
    shorelines = VACUMM_CFG[__name__]['xyz']

    # Check requested name
    if name is not None:
        if name not in shorelines:
            raise ShoreLineError('Shoreline not available: %s. Please use list_shorelines() to print the list of available shorelines.'%name)
        return shorelines[name]

    return shorelines

def list_shorelines():
    """List available shorelines"""
    print('Available shorelines:')
    shorelines = _shorelines_list_()
    if not shorelines:
        print('None available')
    else:
        for name, specs in shorelines.items():
            print('- {} ({})'.format(name, specs['file']))

class _CoastalBathy_:
    _factor = 2.
    def bathy(self, factor=None, ext=True, nl=True, xyz=True, **kwargs):
        """Get bathymetry at coast from mean sea level at stations

        :Params:

            - *factor*: Multiplicative correction factor
            - *xyz*: If True, return a :class:`XYZ` instance instead of a :mod:`numpy` array
            - *ext*: Allow extrapolation
            - *nl*: Nonlinear interpolation

        :Return: A :class:`~vacumm.misc.io.XYZ` instance of sea level
        """
        # Mean sea levels at stations
#       msl = MeanSeaLevel(sel=self.zone())
        try:
            from vacumm.tide.station_info import MeanSeaLevel
        except:
            raise ImportError('Cannot guess bathymetry along shoreline because module "vacumm.tide.station_info" not available')
        msl = MeanSeaLevel()
        # Interpolate to coastal positions
        cb = msl.interp(self, ext=ext, nl=nl, xyz=xyz, **kwargs)
        # Rectification
        if factor is None: factor = self.factor
        cb *= factor

        return cb

    def xyz(self, *args, **kwargs):
        """Shortcut to :meth:`bathy`"""
        return self.bathy(*args, **kwargs)

    def get_factor(self):
        """Get factor to apply to the sea level on the shoreline to convert to bathymetic 'depth'"""
        return self._factor

    def set_factor(self, factor):
        """Set factor to apply to the sea level on the shoreline to convert to bathymetic 'depth'"""
        self._factor = factor
    factor = property(get_factor, set_factor, doc='Factor to apply to the sea level on the shoreline to convert to bathymetic "depth"')

class _PolyShapes_:

    def greatest_polygon(self):
        """Get the greatest polygon (only works with polygons!)"""
        poly = None
        for p in self._shapes:
            if poly is None or p.area() > poly.area():
                poly = p
        return poly

    def plot(self, select=None, ax=None, m='auto', fill=None, color=None,
             edgecolor=None, linewidth=None, show=True, **kwargs):
        return Shapes.plot(self, select=select, ax=ax, fill=fill, m=m,
                           color=color, edgecolor=edgecolor,
                           linewidth=linewidth,
                           show=show, **kwargs)
    plot.__doc__ = Shapes.plot.__doc__



class ShoreLine(_CoastalBathy_, _PolyShapes_,  Shapes):
    """Version of :class:`vacumm.misc.io.Shapes` dedicated to shorelines"""
    _extent = None
    _name = None

    def __init__(self, input=None, *args, **kwargs):
        input = self._check_input_(input)
        Shapes.__init__(self, input, *args, **kwargs)

    def _check_input_(self, input):
        if input is None and self._name is not None:
            input = get_data_file(VACUMM_CFG[__name__]['shapefiles'][self._name],
                                  raiseerr=False)
            if input is None:
                raise VACUMMError("Can't find shoreline file")
            if input.endswith('.shp'):
                input = input[:-4]
        return input

    @classmethod
    def contains(cls, lon, lat):
        """Is this shoreline contains this point?"""
        if cls._extent is None: return True
        xmin, ymin, xmax, ymax = cls._extent
        lon = N.asarray(lon)
        lat = N.asarray(lat)
        return lon.max()>xmin and lon.min()<xmax and lat.max()>ymin and lat.min()<ymax

    embed = contains

    @classmethod
    def avail(cls):
        """Is this shoreline available locally or for download?"""
        if cls._name is None: return True
        return get_data_file(VACUMM_CFG[__file__]['shapefiles'][cls._name], avail=True)
#        return 'shapefile_%s'%cls._name in _shorelines_list_()


class Histolitt(ShoreLine):
    """Shoreline of France from SHOM/IGN at 1/25000from shapefile of Polygons covering metropolitan France"""
    _factor = 1.
    _mres = 18.5532411
    _extent = [-5.1, 41.3, 9.7, 51.12]
    _name = 'histolitt'


class Histolitt2(ShoreLine):
    """Shoreline of France from SHOM/IGN V2 at 1/25000from shapefile of Polygons covering metropolitan France"""
    _factor = 1.
    _mres = 18.5532411
    _extent = [-5.1, 41.3, 9.7, 51.12]
    _name = 'histolitt2'

    def __init__(self, input=None, clip=True, *args, **kwargs):
        ShoreLine.__init__(self, input, clip=False, *args, **kwargs)
        if clip is not None and clip is not False:
            self.clip(clip, copy=False)


class GSHHSSF(ShoreLine):
    """Fine world shoreline from USGS shapefile
    .. warning::

        HUGE! Please use :class:`GSHHS` instead
    """
    _factor = 1.
    _name = 'gshhs'


GSHHS_SF = GSHHSSF  # Compat


class GSHHS(_CoastalBathy_, _PolyShapes_, GSHHSBM):
    _factor = 1.


def get_best_res(gg):
    """Get the best shoreline resolution as letter, for my grid or (lon,lat)"""
    # TODO: get_best_res must be based on meta data!
    lon, lat = get_xy(gg, num=True)
    testresol = ((lon.max()-lon.min())+(lat.max()-lat.min()))/2.0
    if testresol < 1. and Histolitt.avail() and Histolitt.contains(lon, lat):
        return 's'
    return gshhs_autores(lon.min(), lon.max(), lat.min(), lat.max())


get_bestres = get_best_res  # compat


def get_best_shoreline(gg, **kwargs):
    """Get the best shoreline (:class:`GSHHS` or :class:`Histolitt`)
    for my grid or (lon,lat)"""
    return get_shoreline(get_bestres(gg), clip=gg, **kwargs)


get_best = get_best_shoreline  # compat


def get_shoreline(arg, *args, **kwargs):
    """Get a valid ShoreLine instance according to arg.

    :Params:

    - **arg**: it can be either:

        - A :class`~vacumm.misc.io.Shapes` instance (other args are ignored)
        - A string within the folowing lists:

            - for GSHHS resolutions: %s
            - for other shorelines: %s
    """
    if isinstance(arg, Shapes):
        return arg
    if isinstance(arg, str):
        if arg in gshhs_reslist:
            return GSHHS(arg, *args, **kwargs)
        elif arg == 's':  #FIXME: must work with registered shorelines
            if not Histolitt.avail():
                raise ShoreLineError('Shoreline not available: histolitt')
            return Histolitt(*args, **kwargs)
        elif arg in _shoreline_list:
            return eval(arg)(*args, **kwargs)
    return
#   raise TypeError, 'cannot guess a ShoreLine from arg:', arg
get_shoreline.__doc__ =  get_shoreline.__doc__% (gshhs_reslist, _shoreline_list)

register_plot_coastline(Histolitt, 's')
