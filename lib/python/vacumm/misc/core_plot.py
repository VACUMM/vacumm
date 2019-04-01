# -*- coding: utf8 -*-
"""Classes for all plots"""
# Copyright or © or Copr. Actimar/IFREMER (2012-2018)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is a computer program whose purpose is to [describe
# functionalities and technical features of your software].
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
import re
from operator import isNumberType
from warnings import warn

import matplotlib.dates
import matplotlib.pyplot as P
import matplotlib.transforms as mtransforms
import numpy as N, MV2, cdms2
from matplotlib.artist import Artist
from matplotlib.axes import Subplot, Axes
from matplotlib.axis import YAxis
from matplotlib.colors import Colormap, Normalize
from matplotlib.dates import (DateFormatter, MonthLocator, WeekdayLocator,
    YearLocator,DayLocator, HourLocator, MinuteLocator, SecondLocator,
    MONDAY, WEEKLY, YEARLY, MONTHLY,
    AutoDateLocator, AutoDateFormatter, MO, DAILY, HOURLY, num2date,
    MINUTELY, SECONDLY)
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.path import Path
from matplotlib.patheffects import Normal
from matplotlib.text import Text
from matplotlib.ticker import (FormatStrFormatter, Formatter, Locator,
    NullLocator, AutoMinorLocator, AutoLocator)
from matplotlib.transforms import offset_copy
from mpl_toolkits.basemap import Basemap, _setlatlab, _setlonlab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection, LineCollection

from ._ext_plot import (DropShadowFilter, FilteredArtistList, GrowFilter,
    LightFilter)
from .misc import (kwfilter, dict_aliases, geo_scale, lonlab, latlab, deplab, cp_atts,
    auto_scale, dict_check_defaults, basic_auto_scale, dict_copy_items,
    dict_merge, phaselab, set_atts)
from .atime import mpl, strftime, is_numtime, numtime
from .axes import (check_axes, axis_type, set_order, merge_orders,
    check_order, order_match, isaxis, get_axis_type)
from .color import (get_cmap, RGB, land,
    RGBA, change_luminosity, change_saturation, pastelise,
    CMAP_POSITIVE, CMAP_NEGATIVE, CMAP_SYMETRIC)
from .docstrings import docfiller
from .filters import generic2d
from .grid import get_axis, meshbounds, meshgrid
from .grid.masking import resol_mask
from .grid.regridding import shift1d
from .phys.units import deg2m, tometric, m2deg
from .remote import OutputWorkFile
from ..config import get_config_value
import color


MA = N.ma



__all__ = ['PlotError','Plot', 'Plot1D', 'Curve', 'Bar', 'Stick',
    'Plot2D', 'Map', 'Hov', 'QuiverKey',
    'ScalarMappable','AutoDateFormatter2', 'AutoDateLocator2',
    'AutoDateMinorLocator', 'AutoDualDateFormatter', 'DualDateFormatter',
    'MinuteLabel', 'Section', 'twinxy', 'DepthFormatter',
    'AutoDegreesMinutesFormatter', 'AutoDegreesMinutesLocator']


#: Aliases for argisimage services hosted by the arcgis server
ARCGISIMAGE_ALIASES = dict(
    esriimagery="ESRI_Imagery_World_2D",
    esristreet="ESRI_StreetMap_World_2D",
    esristreetmap="ESRI_StreetMap_World_2D",
    natgeo="NatGeo_World_Map",
    ngstopo="NGS_Topo_US_2D",
    usatopo="USA_Topo_Maps",
    ocean="Ocean_Basemap",
    topo="World_Topo_Map",
    shaded="World_Shaded_Relief",
    physical="World_Physical_Map",
    imagery="World_Imagery",
    street="World_Street_Map",
    streetmap="World_Street_Map",
    terrain="World_Terrain_Base",
)


class PlotError(Exception):
    pass

class Plot(object):
    """Base class for all plots

    :Generic params:

        - **load_data**, optional: Load data.
        - **pre_plot**, optional: Initialize the plot (preprocessing).
        - **plot**, optional: Plot data.
        - **post_plot**, optional: Finalize plot.
        - **Data loading**: See :meth:`load_data`, :meth:`_set_axes_`, :meth:`_check_order_` .
        - **Plot initialisation**: see  :meth:`pre_plot`.
        - **Plot**: see :meth:`plot`
        - **Plot finalization**: see :meth:`post_plot`


    :Attribute params:


        - **long_name**, optional: Force the :attr:`~vacumm.misc.core_plot.Plot.long_name`
          attribute used in :attr:`~vacumm.misc.core_plot.Plot.title`.
          See also param ``x/ylong_name``.
        - **x/ylong_name**, optional: Same as :attr:`~vacumm.misc.core_plot.Plot.long_name`
          but for X and Y axes (:attr:`~vacumm.misc.core_plot.Plot.xlong_name`
          or :attr:`~vacumm.misc.core_plot.Plot.ylong_name`). It refers only to
          axes for 2D plots, and potentially to data for 1D plots.
        - **units**, optional: Force the :attr:`~vacumm.misc.core_plot.Plot.units`
          attribute used in labels (:attr:`~vacumm.misc.core_plot.Plot.label`,
          :attr:`~vacumm.misc.core_plot.Plot.xlabel` or
          :attr:`~vacumm.misc.core_plot.Plot.ylabel`).
        - **latex_units**, optional: Interpret units with latex after defining
          the  :attr:`~vacumm.misc.core_plot.Plot.latex_units` attribute.
          Alternatively, you can simply specify units enclosed with ``"$"``.
        - **x/yunits**, optional: Same as :attr:`~vacumm.misc.core_plot.Plot.units`
          but for X and Y axes (:attr:`~vacumm.misc.core_plot.Plot.xunits`
          or :attr:`~vacumm.misc.core_plot.Plot.yunits`). It refers only to
          axes for 2D plots, and potentially to data for 1D plots.
          See also param ``x/yunits``.
        - **x/ymin/max**, optional: Force min and max along X and Y by setting
          attributes :attr:`~vacumm.misc.core_plot.Plot.xmin`,
          :attr:`~vacumm.misc.core_plot.Plot.xmax`,
          :attr:`~vacumm.misc.core_plot.Plot.ymin` or
          :attr:`~vacumm.misc.core_plot.Plot.ymax`.
        - **vmin/max**, optional: Force min and max value for data by setting
          attributes :attr:`~vacumm.misc.core_plot.Plot.vmin`,
          :attr:`~vacumm.misc.core_plot.Plot.vmax`.
          This may be equivalent to
          to set X or Y extrema for 1D plots.
        - **x/ylabel**, optional: Force label used for X and Y axes
          by setting attributes :attr:`~vacumm.misc.core_plot.Plot.xlabel` or
          :attr:`~vacumm.misc.core_plot.Plot.ylabel`.
        - **title**, optional: Force the title of the plot by setting
          :attr:`~vacumm.misc.core_plot.Plot.title`
        - **x/ymasked**, optional: Force the plot to fit all data
          positions along X and/or Y,
          even if data is missing, by setting attributes
          :attr:`~vacumm.misc.core_plot.Plot.xmasked`,
          :attr:`~vacumm.misc.core_plot.Plot.ymasked`
          or :attr:`~vacumm.misc.core_plot.Plot.xymasked`.
        - **anim**, optional: Create an animation for the current figure.

        Example:

        >>> myplot = Plot(data, xmin=3.5, title='My plot')
        >>> print myplot.xmin, myplot.title
        3.5 Myplot



        The followwing rules apply:

        - :attr:`title` defaults to :attr:`long_name`,
          which defaults to the ``long_name`` attribute of input data.
          You can use templates with all other attributes as replacement
          keys (example ``"%(long_name)s (max=%(vmax)g)"``) but
          :attr:`xlabel` and :attr`ylabel`.
        - :attr:`units` defaults to the ``units`` attribute of input data.
        - :attr:`xlong_name` defaults to :attr:`long_name` if X
          axis refers to data, otherwize to the ``long_name`` attribute
          of input second axis.
        - :attr:`ylong_name` defaults to :attr:`long_name` if Y
          axis refers to data, otherwize to the ``long_name`` attribute
          of input first axis.
        - :attr:`xunits` defaults to :attr:`units` if X
          axis refers to data, otherwize to the ``units`` attribute
          of input second axis.
        - :attr:`yunits` defaults to :attr:`units` if Y
          axis refers to data, otherwize to the ``units`` attribute
          of input first axis.
        - :attr:`xlabel` is empty if X axis refers to a spatial or
          temporal axis, else it defaults ``"%(xlong_name)s [%(xunits)s]"``.
          You can use templates with all other attributes as replacement
          keys (example ``"%(long_name)s (max=%(vmax)g)"``) but
          :attr:`title` and :attr`ylabel`.
        - :attr:`ylabel` is empty if Y axis refers to a spatial or
          temporal axis, else it defaults ``"%(ylong_name)s [%(yunits)s]"``.
          You can use templates with all other attributes as replacement
          keys (example ``"%(long_name)s (max=%(vmax)g)"``) but
          :attr:`title` and :attr`xlabel`.
        - :attr:`vmin`/:attr:`vmax` defaults to the min/max of all plotted data.
        - :attr:`xmin`/:attr:`xmax` defaults to :attr:`vmin`/:attr:`vmin`
          if X axis refers to data, else to the min/max of
          the second axis of input data.
        - :attr:`ymin`/:attr:`ymax` defaults to :attr:`vmin`/:attr:`vmin`
          if Y axis refers to data, else to the min/max of
          the first axis of input data.
        - :attr:`uvlat` and :attr:`uvscaler` are used to convert vectors
          (like speed) used by quiver plots.

        .. warning::

            These attributes may be used for the plotting process but may not be
            equal to their graphical counterpart. For instance, :attr:`xmin`
            may not be equal to ``matplotlib.pyplot.xlim()[0]``.


    :Generic tasks:

        #. Call to :meth:`load_attributes`.
        #. Call to :meth:`load_data` if ``load_data is True``.
        #. Call to :meth:`pre_plot` if ``pre_plot is True``.
        #. Call to :meth:`register`.
        #. Call to :meth:`load_attributes` for remaining attributes.
        #. Call to :meth:`plot` if ``plot is True``.
        #. Call to :meth:`post_plot` if ``plot is True and post_plot is True``.


    """
    rank = None
    _order = None
    _primary_attributes = ['long_name', 'xlong_name', 'ylong_name',
        'units', 'xunits', 'yunits', 'vmin', 'xmin', 'ymin', 'vmax',
        'xmax', 'ymax']
    _secondary_attributes = ['xlabel', 'ylabel', 'title']
    _special_attributes = ['uvlat', 'uvscaler', 'anim', 'xymasked',
        'xmasked', 'ymasked']
    _plot_status_none = 0
    _plot_status_plot = 1
    _plot_status_post = 2
    masked = True


    def __init__(self, data, load_data=True, pre_plot=True, plot=False, post_plot=False,
        **kwargs):
        # First inits
        self._gobjs = {}
        self.data = self.axes = None
#        self._kwargs = kwargs
        self._post_plotted = False
        self._finalize = []

        # Load data
        if load_data :
            self.load_data(data, **kwargs)

        # Load attributes: for current instance
        self.load_attributes(kwargs)

        # Initialize plot
        self.pre_plot(**kwargs)

        # Register
        self.register()

        # Load attributes: stored in axes instance
        self.load_attributes(kwargs)

        # Plot
        if plot:
            self.plot(**kwargs)

            # Post plot
            if post_plot:
                self.post_plot(**kwargs)

    def plot_status(self):
        """Guess if current instance has a plot"""
        if self.get_obj('plots') is None:
            return self._plot_status_none
        if self._post_plotted:
            return self._plot_status_post
        return self._plot_status_plot

#    def update(self, status=None, glob=True):
#        """Update what is needed of the plot"""
#        # Guess the status
#        if status is None: status = self.plot_status()
#        if not status: return
#
#        # Update
#        if glob: # Global update (all plotters)
#
#            for plotter in self.get_brothers():
#                plotter.update(status=status, glob=False)
#
#        else: # Update this plotter only
#
#            if status>=self._plot_status_plot: # Re-plot
#                self.plot(*self._pargs, **self._kwargs)
#
#                if status==self._plot_status_post: # Re-post_plot
#                    self.post_plot(**self._kwargs)

    def is_plotted(self):
        return self.plot_status()>=self._plot_status_plot
    def is_post_plotted(self):
        return self.plot_status()>=self._plot_status_post

    def load_attributes(self, items, select=None):
        """Load (selected) attributes from ``items``

        Attributes not found are set to ``None``.

        :Example:

            >>> items = dict(xmin=4., title='%(long_name)s')
            >>> myplot.load_attributes('xmin', 'units', units='degC', **items)
        """
        # User attributes
        for att in self._primary_attributes+self._secondary_attributes+self._special_attributes:
            if select and att not in select: continue
            value = items.get(att, None)
            if value is not None:
                setattr(self, att, value)
                if self.isset(att): del items[att] # clean out

    def register(self):
        """Register current instance into self.fig and self.axes"""
        if self.axes is None:
            raise PlotError('Cannot register plot instance since no MPL axes available')
        if self.get_axobj('plotters') is None or self not in self.get_axobj('plotters'):
            self.add_axobj('plotters', self)
            self.fig._vacumm = self.get_axobj()


    def load_data(self, data, **kwargs):
        """Load data and format data, and check rank.
        It finally calls :meth:`_check_order_`.

        :Params:

            - **data**: A single :mod:`cdms2` variable or a tuple of them,
              in the forms ``(m,)``, or ``(u,v)``, ``(m,u,v)``, where:

                - ``u``: X component of a vector.
                - ``v``: Y component of a vector.
                - ``m``: A scalar variable. If ``u`` and ``v`` are set,
                  it defaults to their modulus.

            - **order**, optional: See :meth:`_check_order_`
            - **transpose**, optional: See :meth:`_check_order_`
            - Keywords are passed to :meth:`_set_axes_` for
              axis subtitutions.

        :Attributes:

            The following attributes are defined:

            .. attribute:: data

                A 1- to 3-element tuple of :class:`MV2.array`
                in a form of similar to **data** above.

            .. attribute:: x

                The last axis of the first element of data,
                or ``None`` if X axis refer to data.

            .. attribute:: y

                The first axis of the first element of data,
                or ``None`` if Y axis refer to data.

            .. attribute:: mask

                The data mask.

        :See also:

            :meth:`_check_order_` :meth`_set_axes_`
        """
        self.raw_data = data
        self.data =  None
#        if data is None:
#            self._check_order_()
#            return

        # Arrows
        if data is not None:
            N.seterr(over='ignore')
            if isinstance(data, list):
                data = tuple(data)
            if not isinstance(data, tuple):
                self.data = [MV2.array(data, copy=1),]
            elif len(data) in [1,3]: # (mod,) or (mod,u,v)
                self.data = [MV2.array(v, copy=1) for v in data]
            elif len(data)==2:
                # Get vectors and modulus (u,v) -> (mod,u,v)
                vx,vy = [MV2.array(v, copy=1) for v in data]
                self.data = [MV2.sqrt(vx**2+vy**2),vx,vy]
                ln = []
                for vv in vx,vy:
                    if hasattr(vv,'units'):
                        self.data[0].units = vv.units
                        break
                    if hasattr(vv,'long_name'):
                        ln.append(vv.long_name)
                ln  = ' and '.join(ln)
                if len(ln) == 0:
                    ln = 'Modulus'
                else:
                    ln = 'Modulus of '+ln
                self.data[0].long_name = ln
                self.data[0].setAxisList(self.data[2].getAxisList())
                self.data[0].setGrid(self.data[2].getGrid())
                self.data[0].id = 'modulus'
            else:
                raise PlotError('You must provide no more than 3 arrays as input for '+self.__class__.__name__)

            # Check all variables and axes
            units = [None]*self.rank
            for ivar,var in enumerate(self.data):

                # Check rank
                if var.rank() != self.rank:
                    raise PlotError('Your variable must have a rank = %i (current rank is %i)' % (self.rank, var.rank()))

                # Guess axis types
                if '-' in var.getOrder():
                    check_axes(var)

                # Make time axes comptatible with matplotlib
                mpl(var, copy=True)

        # Store axes in x and y
        self._set_axes_(**kwargs)

        # Homogeneous axis units and long_names
        self._check_axis_atts_()

        # Check order
        self._check_order_(**kwargs)

        # Store  mask
        #self.x = get_axis(self.data[0], -1) if self.order[1]!='d' else None
        #self.y = get_axis(self.data[0], 0) if self.order[0]!='d' else None
        #self.x = self.data[0].getAxis(-1) if self.order[1]!='d' else None
        #self.y = self.data[0].getAxis(0) if self.order[0]!='d' else None
        self.mask = (N.ma.getmaskarray(self.data[0]) if self.has_data() else
                     N.ma.nomask)
        self.masked = self.mask.all()
        if self.masked and self.has_data():
            warn('All your data are masked')

        return self.data

    def _set_axes_(self, **kwargs):
        pass

    def _check_axis_atts_(self):
        """Check that X/Y units and long_name can be used if possible

        If current axes (self.x and self.y) have no units or long_name,
        it tries to get them from variables (self.data).

        :attr`x` and :attr`y` must have been set _set_axes_ method.
        """
        if not self.has_data():
            return
        for i, xy in enumerate(['y','x']):
            axis = getattr(self, xy)
            if axis is None: continue
            for att in 'units', 'long_name':
                if hasattr(axis, att): continue
                for dat in self.data:
                    if self.rank==1:
                        ax = dat.getAxis(0)
                    else:
                        ax = get_axis(dat, i, strict=False)#, geo=False)
                    if hasattr(ax, att):
                        setattr(axis, att, getattr(ax, att))
                        break


    def _check_order_(self, order=None, transpose=False, vertical=None, **kwargs):
        """Check the order of axes

        :Params:

            - **order**, optional: A string of length 2 that specify the physical type
              of plot axes. The first char refers to the Y axis, and the second
              to the X axis. It must contain one of the following characters:

                - ``x``: longitude
                - ``y``: latitude
                - ``z``: level
                - ``t``: time
                - ``-``: any of the latters
                - ``d``: data

              If char is upper-cased, its type is mandatory and must be found
              in input :mod:`cdms2` variables.
              It defines the :attr:`~vacumm.misc.core_plot.Plot.order` attribute.

            - **vertical**: Force data to be plotted along the vertical axis.
            - **transpose**: Transpose the plot axes
              (:attr:`~vacumm.misc.core_plot.Plot.order`).

        :Attributes:

            This methods defines the following attributes:

            .. attribute:: order

                Two-character string define the order of axis types.
                It is the sum of the :attr:`xtype` and :attr:`ytype`.
                Examples: ``"xt"``, ``z-``.

            .. attribute:: xtype

                Type of the X axis (see above).

            .. attribute:: ytype

                Type of the Y axis (see above).
        """
        if self._order is None:
            raise NotImplementedError
        if not self.has_data():
            if not isinstance(self._order, list):
                self.order = self._order.lower()
            else:
                self.order = ('d' if self.y is None else
                              get_axis_type(self.y, checkaxis=False))
                self.order = self.order + ('d' if self.x is None else
                              get_axis_type(self.x, checkaxis=False))
        else:

            # Allowed orders
            if not isinstance(self._order, list): self._order = [self._order]
            withd = 'd' in self._order[0]
            if order is not None:
                if not isinstance(order, list): order = [order]
                for oo in order:
                    if len(oo) != 2 or (withd and not 'd' in oo):
                        raise PlotError('You order must be a 2-element string with "d" inside')
                self._order = order

            # Loop on variables
            data_order = None
            for i, var in enumerate(self.data):

                self.data[i], this_order, _reordered = check_order(var,
                    self._order, reorder=True,
                    vertical=vertical, extended=True, getorder=True)
                if data_order is not None:
                    assert order_match(data_order, this_order), \
                        "Axis order incompatible: %s and %s"%(last_order, this_order)
                    data_order, _ = merge_orders(data_order, this_order)
                else:
                    data_order = this_order

                if i==0:
                    reordered = _reordered

            self.order = data_order

            # Transpose axes when reordered
            if reordered or (self.rank==1 and data_order[1]=='d'):
                self._transpose_axes_()

        # Transpose?
        if transpose: self._transpose_()

        self.ytype, self.xtype = self.order
        return self.order

    def _transpose_(self):
        """Transpose data and axes"""

        # Variables
        if self.rank==2:
            for ivar, var in enumerate(self.data):
                self.data[ivar] = MV2.transpose(var)
                del var

        # Axes
        self._transpose_axes_()

        # Order
        self.order = self.order[::-1]

    def _transpose_axes_(self):
        for xy in 'x','y':
            axis = getattr(self, xy)
            if axis is None: continue
            if len(axis.shape)==2:
                axis = MV2.transpose(axis)
            setattr(self, xy, axis)
        self.x, self.y = self.y, self.x

    def get(self, key, **kwargs):
        return getattr(self, 'get_'+key)(**kwargs)

    def get_xmask(self):
        """Get the data mask projected along X"""
        if self.mask is None: return
        if self.xtype == 'd': return self.mask
        xdata = self.get_xdata(masked=False)
        if xdata.ndim==self.mask.ndim==1:
            return self.mask
        return self.mask.min(axis=0)

    def get_ymask(self):
        """Get the data mask projected along Y"""
        if self.mask is None: return
        if self.ytype == 'd': return self.mask
        ydata = self.get_ydata(masked=False)
        if ydata.ndim==self.mask.ndim==1:
            return self.mask
        return self.mask.min(axis=-1)


    def get_data(self, scalar=False):
        """Get data as a tuple of :class:`~numpy.ma.core.MaskedArray`

        :See also: :meth:`get_xdata` :meth:`get_ydata` :attr:`uvscaler`
        """
        data = [var.asma() for var in self.data]
        if len(data)==3:
            uvscaler = self.uvscaler
            if uvscaler is not None:
                data = [data[0]]+list(uvscaler(data[1], data[2]))
            if scalar=='uv':
                data = data[1:]
        if scalar is not False:
            if scalar is True:
                scalar = 0
            if isinstance(scalar, int):
                return data[scalar]
        return tuple(data)

    def get_xdata(self, scalar=True, masked=False, bounds=False):
        """Get the numerical data associated with the X axis

        .. note::

            It can come from a physical axis or data
            depending on the axis type :attr:`xtype`.

        :Params:

            - **scalar**, optional: Set it to ``True`` to get data
              as a scalar array in case X axis refers to a tuple of data.
              If set to an int, it takes the element #scalar of this tuple.
            - **masked**, optional: If it is an axis (not data), values are
              masked with data mask.
            - **bounds**, optional: The data bounds (valid only of X
              is an axis).

        :See also: :meth:`get_ydata` :meth:`get_data`
        """
#        # Nothing
#        if not self.has_data(): return
        # Axis
        if self.x is not None:
            x = self.x[:] #.getValue()
            if cdms2.isVariable(x): x = x.asma()
            if N.ma.isMA(x): x = x.filled(x.max())
            if masked is None:
                masked = self.xmasked
            if masked:
                mask = self.get_xmask()
                if mask.shape!=x.shape:
                    mask = N.resize(mask, x.shape)
                x = N.ma.masked_array(x, mask=mask)
            if bounds:
                if not hasattr(self, '_xb'):
                    self._xb = meshcells(x)
                return self._xb
            return x
        # Data to plot
        if scalar is True:
            scalar = max(0, len(self.data)-2)
        return self.get_data(scalar)

    def get_ydata(self, scalar=True, masked=False, bounds=False):
        """Get the numerical data associated with the Y axis

        .. note::

            It can come from a physical axis or data
            depending on the axis type :attr:`ytype`.

        :Params:

            - **scalar**, optional: Set it to ``True`` to get data
              as a scalar array in case Y axis refers to a tuple of data.
              If set to an int, it takes the element #scalar of this tuple.
            - **masked**, optional: If it is an axis (not data), values are
              masked with data mask.
            - **bounds**, optional: The data bounds (valid only of Y
              is an axis).

        :See also: :meth:`get_xdata` :meth:`get_data`
        """
#        # Nothing
#        if not self.has_data(): return
        # Axis
        if self.y is not None:
            y = self.y[:] #.getValue()
            if cdms2.isVariable(y): y = y.asma()
            if N.ma.isMA(y): y = y.filled(y.max())
            if masked is None:
                masked = self.ymasked
            if masked:
                mask = self.get_ymask()
                if mask.shape!=y.shape:
                    mask = N.resize(mask, y.shape[::-1]).T
                y = N.ma.masked_array(y, mask=mask)
            if bounds:
                if not hasattr(self, '_yb'):
                    self._yb = meshcells(y)
                return self._yb
            return y
        # Data to plot
        if scalar is True:
            scalar = max(0, len(self.data)-1)
        return self.get_data(scalar)




    def pre_plot(self, axes=None, figure=None, figsize=None, subplot=None, twin=None,
            subplots_adjust=None, bgcolor=None, noframe=False, fullscreen=False,
            verbose=False, axes_host=False, axes_xoffset=0, elev=None, azim=None,
            **kwargs):
        """Initialize the plot

        :Tasks:

            #. Filter keyword parameters.
            #. Create the :class:`~matplotlib.fig.Figure` instance and
               store it into :attr:`fig`.
            #. Create the :class:`~matplotlib.axes.Axes` instance and
               store it into :attr:`axes`

        :Params:

            - **fig**, optional: Figure number.
            - **figsize**, optional: Initialize the figure with this size.
            - **axes**, optional: Use this axes object.
            - **subplot**, optional: Call to :func:`~matplotlib.pyplot.subplot` to create axes.
            - **subplots_adjust**, optional: Dictionary sent to :func:`~matplotlib.pyplot.subplots_adjust`.
              You can also use keyparams 'left', 'right', 'top', 'bottom', 'wspace', 'hspace'!
            - **top/bottom/left/right/wspace/hspace**, optional: Override ``subplots_adjust``.
            - **sa**, optional: Alias for subplots_adjust.
            - **twin**, optional: Use ``"x"`` or ``"y"`` or ``"xy"``  to make a copy of current
              X or Y axes (see :func:`matplotlib.pyplot.twinx`).
              You can also provide a dictionary : ``twin=dict(x=axes1, y=axes2)``.
            - **bgcolor**, optional: Background axis color.
            - **axes_rect**, optional: [left, bottom, width, height]
              in normalized (0,1) units to create axes using :func:`~matplotlib.pyplot.axes`.
            - **axes_<param>**, optional: <param> is passed to :func:`~matplotlib.pyplot.axes`.
            - **noframe**, optional: Suppress plot frames.
            - **fullscreen**, optional: Plot in full screen mode (thus, ``noframe==True``).
            - **verbose**, optional: Informs about errors with axes.


        :Attributes:

            .. attribute:: fig

                :class:`~matplotlib.fig.Figure` on which plots are drawn.

            .. attribute:: axes

                :class:`~matplotlib.axes.Axes` instance of the current plot.

        """
        # Keywords management
        figure = kwargs.pop('fig', figure)
        kwfig = kwfilter(kwargs,'figure')
        kwaxes = kwfilter(kwargs,'axes')
        subplots_adjust = kwargs.pop('sa', subplots_adjust)
        for adj in 'left', 'right', 'top', 'bottom', 'wspace', 'hspace':
            if kwargs.has_key(adj):
                if subplots_adjust is None: subplots_adjust = {}
                subplots_adjust[adj] = kwargs.pop(adj)
        if fullscreen:
            noframe = True
            subplot = None
            subplots_adjust = None
            default_axes_rect = [0,0,1,1]
        else:
            default_axes_rect = None
        axes_rect = kwaxes.pop('rect', default_axes_rect)
        if axes_rect is not None:
            axes_rect = [axes_rect]
        else:
            axes_rect = []
        if noframe:
            kwargs['figure_frameon'] = False
            kwargs['axes_frameon'] = False

        # Figure
        if figure == 'old':
            figure = None
        if isinstance(figure, Figure):
            self.fig = figure
        elif figure is True or figure == 'new':
            self.fig = P.figure(**kwfig)
        elif figure is not None:
            self.fig = P.figure(figure, **kwfig)
        else:
            self.fig = P.gcf()
        if figsize is not None:
            if not isinstance(figsize, tuple):
                figsize = (figsize, figsize)
            else:
                figsize = tuple(figsize)
            self.fig.set_size_inches(figsize, forward=True)
        if noframe:
            self.fig.set_frameon(False)
        if subplots_adjust is not None:
            self.fig.subplots_adjust(**subplots_adjust)

        # Axes
        if (axes_host and axes is None and subplot is None and not axes_rect
                and twin is None):
            subplot = 111
        if axes=="3d":
            kwaxes['projection'] = "3d"
            if subplot is None and not axes_rect:
                subplot = 111
            axes = None
        if axes is not None:
            self.axes = axes
            if self.axes.get_figure() != self.fig:
                if verbose: print 'Axes does not match figure'
                self.fig = self.axes.get_figure()
        elif subplot is not None:
            if axes_host:
#                from mpl_toolkits.axes_grid1 import host_subplot
                from mpl_toolkits.axes_grid1.parasite_axes import host_subplot_class_factory
                from mpl_toolkits.axisartist import Axes as AAxes
                host_subplot_class = host_subplot_class_factory(AAxes)
                if isinstance(subplot,(list,tuple)):
                    self.axes = host_subplot_class(self.fig, *subplot, **kwaxes)
                else:
                    self.axes = host_subplot_class(self.fig, subplot, **kwaxes)
                self.fig.add_subplot(self.axes)
            else:
                if isinstance(subplot,(list,tuple)):
                    self.axes = self.fig.add_subplot(*subplot,**kwaxes)
                else:
                    self.axes = self.fig.add_subplot(subplot,**kwaxes)
        elif axes_rect:
            self.axes = self.fig.add_axes(*axes_rect,**kwaxes)
        else:
            if isinstance(twin, str):
                self.axes = twinxy(twin, ax=self.axes)
            else:
                self.axes = P.gca()
        try: # Fails with host
            self.fig.sca(self.axes)
        except:
            pass
        if axes_xoffset and hasattr(self.axes, 'get_grid_helper'):
            new_fixed_axis = self.axes.get_grid_helper().new_fixed_axis
            if axes_xoffset>0:
                loc = 'right'
            else:
                loc = 'left'
            offset = (axes_xoffset, 0)
            self.axes.axis[loc] = new_fixed_axis(loc=loc, axes=self.axes, offset=offset)
            self.axes.axis[loc].toggle(all=True)
        self.is3d = isinstance(self.axes, Axes3D)
        if self.is3d:
            self.axes.view_init(elev=elev, azim=azim)
        if noframe:
            self.axes.set_frame_on(False)
        if bgcolor is not None:
            try:
                self.axes.set_facecolor(bgcolor)
            except AttributeError:
                self.axes.set_axis_bgcolor(bgcolor)

    def plot(self, **kwargs):
        """The main plot"""
        if self.plot_status()==2:
            self.axes.cla()

    def clear(self):
        """Clear axes from plotted objects"""
        objs = self.get_obj('plotted')
        if objs is None: return
        for obj in objs:
            self.remove(obj)
        if self.axes is not None:
            for container in self.axes.xaxis, self.axes.yaxis, self.axes:
                if hasattr(container, '_vacumm'):
                    del container._vacumm
        if hasattr(self, '_gobjs'):
            del self._gobjs

    def remove(self, objs):
        """Remove an graphical object from axes"""
        if isinstance(objs, tuple):
            objs = list(objs)
        if not isinstance(objs, list):
            objs = [objs]
        for obj in objs:
            for axobjs in self.axes.lines, self.axes.patches:
                if obj in axobjs:
                    axobjs.remove(obj)
            if self.axes is not None:
                for container in self.axes.xaxis, self.axes.yaxis, self.axes:
                    if hasattr(container, '_vacumm'):
                        for key, val in container._vacumm.iteritems():
                            if val is obj:
                                del container._vacumm[key]
            pp = self.get_obj('plotted')
            if pp is not None and obj in pp:
                pp.remove(obj)

    def cla(self):
        """Clear axes of everything

        See :func:`clear` and :func:`matplotlib.pyplot.cla`
        """
        self.clear()
        self.axes.cla()

    def clf(self):
        """Clear figure of everything

        See :func:`clear`, :func:`cla` and
        :func:`matplotlib.pyplot.clf`
        """
        self.cla()
        self.fig.clf()

    def close(self):
        """Close the current
        See :func:`clear`, :func:`cla` and :func:`clf`
        """
        self.clf()
        try:
            P.close(self.fig)
        except:
            pass


    def get_brothers(self, notme=False, mefirst=True, filter=False):
        """Return all :class:`Plot` instances that belongs to current axes

        :Params:

            - **notme**, optional: Do not include current object in the list.
            - **mefirst**, optional: Place me at the beginning of the list.
            - **filter**, optional: If callable, use it to filter out brothers.
        """
        brothers = self.get_axobj('plotters')
        brothers = [] if brothers is None else list(brothers)
        if self in brothers:
            if notme:
                brothers.remove(self)
            elif mefirst:
                brothers = [self]+brothers.pop(self)
        if callable(filter):
            brothers = [b for b in brothers if filter(b)]
        return brothers

    @classmethod
    def get_current(cls, axes=None):
        """Retreive an instance of this class if found to be plotted in currents axes

        :Params:

            - **axes**, optional: Check this axes instance instead of the current one.

        :Return: Last plotted instance, else ``None``

        :Example:

            >>> m = Map.get_current()
        """
        if axes is None:
            if P.get_fignums() and P.gcf().axes:
                axes = P.gca()
            else:
                return
        if not hasattr(axes, '_vacumm'): return
        for p in axes._vacumm.get('plotters', [])[::-1]:
            if isinstance(p, cls): return p



    def add_obj(self, gtype, obj, single=False):
        """Add a graphic object to the bank of current instance

        :Params:

            - **gtype**: A list (or a single element) of string keys
              to name the object.
            - **obj**: The object it self (may be a list).
            - **single**, optional: If ``True``, ``obj`` if store as is
              (i.e is not appended to existing store objects having the same name).

        :Return:

            The object added.

        :Example:

            >>> text_object = myplot.add_obj(['plotted', 'text', myplot.axes.text(10, 20, 'text'))
            >>> text_object = myplot.add_obj('colorbar', myplot.colorbar(), single=True)

        :See also:

            :meth:`set_obj` :meth:`get_obj`
        """
        # Store the object
        if not hasattr(self, '_gobjs'):
            self._gobjs = {}
        if not isinstance(gtype, list): gtype = [gtype]
        for gt in gtype:
            if single:
                 self._gobjs[gt] = obj
            else:
                if not self._gobjs.has_key(gt):
                    self._gobjs[gt] = []
                if isinstance(obj, list):
                    self._gobjs[gt].extend(obj)
                else:
                    self._gobjs[gt].append(obj)

        return obj

    def set_obj(self, gtype, obj):
        """Shortcut to :meth:`add_obj` called with ``single=True``"""
        return self.add_obj(gtype, obj, single=True)

    def get_obj(self, gtype):
        """Get a graphic object stored in the bank

        :Example:

            >>> myplot.get_obj('pcolor')[0].set_zorder(15)
            >>> myplot.get_obj('key').set_color('red')

        :Return:

            The object or ``None`` if not found.

        :See also:

            :meth:`add_obj`
        """
        if not hasattr(self, '_gobjs'):
            self._gobjs = {}

        # Dict key
        if isinstance(gtype, str):
            return self._gobjs.get(gtype, None)

        # Check type
        for oo in self._gobjs.values():
            if isinstance(oo, list):
                for o in oo:
                    if isinstance(oo, gtype):
                        return o
            elif isinstance(oo, gtype):
                return oo
    def del_obj(self, gtype):
        self._gobjs.pop(gtype, None)

    def __getitem__(self, gtype):
        return self.get_obj(gtype)

    def add_axobj(self, gtype, obj, single=False, axis=None):
        """Add a object to the bank of current :class:`matplotlib.axes.Axes` instance

        :Return:

            The object added.

        :Example:

            >>> text_object = myplot.add_axobj('vmin', 24.5)

        :See also:

            :meth:`get_axobj`
        """
        if self.axes is None: return
        container = self.axes if axis is None else getattr(self.axes, axis+'axis')
        if not hasattr(container, '_vacumm'):
            container._vacumm = {}
        if single:
             container._vacumm[gtype] = obj
        else:
            if not container._vacumm.has_key(gtype):
                container._vacumm[gtype] = []
            container._vacumm[gtype].append(obj)
        return obj

    def set_axobj(self, gtype, obj, axis=None):
        """Shortcut to :meth:`add_axobj` called with ``single=True``"""
        return self.add_axobj(gtype, obj, single=True, axis=axis)

    def get_axobj(self, gtype=None, axis=None, axes=None):
        """Get an object stored in the bank of current
        :class:`matplotlib.axes.Axes` instance

        :Params:

            - **gtype**, optional: Object type (name).
              If not set, all objects are returned.
            - **axis**, optional: If one of ``"x"`` or ``"y"``,
              get objects stored in current
              xaxis or yaxis instead if current axes instance.
            - **axes**, optional: Target axes, which defaults to

                #. attribute :attr:`axes`,
                #. result from :func:`matplotlib.pyplot.gca`.

        :Example:

            >>> myplot.get_axobj()
            >>> myplot.get_axobj('vmin')
            >>> myplot.get_axobj('hlitvs', axis='x')

            >>> Plot.get_axobj()

        :Return:

            The object or ``None`` if not found.

        :See also:

            :meth:`add_axobj`
        """
        if axes is None and hasattr(self, 'axes'): axes = self.axes
        if axes is None: return
            #if not P.get_fignums() or not P.gcf().axes: return
            #axes = P.gca()
        container = axes if axis is None else getattr(axes, axis+'axis')
        if container is None: return
        if not hasattr(container, '_vacumm'):
            container._vacumm = {}
        if gtype is None: return container._vacumm
        return container._vacumm.get(gtype, None)

    def del_axobj(self, gtype, axis=None):
        container = self.axes if axis is None else getattr(self.axes, axis+'axis')
        if container is None: return
        if container._vacumm.has_key(gtype):
            del container._vacumm[gtype]

    def isset(self, key):
        """Check if an attribute has been manually set different from ``None``


        :Example:

            >>> return myplot.iset('xmin')

        :See also: :meth:`get_obj` :meth:`get_axobj`
        """
        if key in ['units', 'long_name']:
            return getattr(self, key) is not None
        return self.get_axobj(key) is not None or self.get_obj(key) is not None


    def post_plot(self, grid=True, figtext=None, show=True,
        close=False, savefig=None, savefigs=None, title=None,
        fullscreen=False, anchor=None, autoresize=2, finalize=None,
        key=False, hlitvs=False, legend=False, tight_layout=False,
        param_label=None, **kwargs):
        """Finish plotting stuff (plot size, grid, texts, saves, etc)

        :Params:

            - **title**: Title of the figure [defaults to var.long_name or '']
            - **grid**: Plot the grid [default: True]
            - **grid_<param>**: <param> is passed to :func:`~matplotlib.pyplot.grid`
            - **hlitvs**: Add highlithing if time axis [default: False]
            - **figtext**: figtext Add text at a specified position on the
              figure. Example: figtext=[0,0,'text'] add a 'text' at the
              lower left corner, or simply figtext='text'.
            - **figtext_<param>**: <param> is passed to
              :func:`~matplotlib.pyplot.figtext`
            - **anchor**: Anchor of the axes (useful when resizing) in
              ['C', 'SW', 'S', 'SE', 'E', 'NE', 'N', 'NW', 'W'].
            - **legend**, optional: Draw the legend using :func:`~matplotlib.pyplot.legend`.
            - **legend_<param>**: <param> is passed to :func:`~matplotlib.pyplot.legend`
            - **show**: Display the figure [default: True]
            - **savefig**: Save the figure to this file.
            - **savefig_<param>**: <param> is passed to method :meth:`savefig`
              and finally to the matplotlib function :func:`~matplotlib.pyplot.savefig`.
            - **savefigs**: Save the figure into multiple formats using
              :func:`savefigs` and 'savefigs' as the prefix to the files.
            - **savefigs_<param>**: <param> is passed to :func:`savefigs`
            - **autoresize**: Auto resize the figure according axes (1 or True),
              axes+margins (2). If 0 or False, not resized [default: False=2].
            - **key**: Add a key (like 'a)') to the axes using add_key
              if different from None [default: None]
            - **key_<param>**: <param> is passed to :func:`add_key`
            - **param_label**: Add a param label to the figure using
              :meth:`add_param_label` if different from None [default: None]
            - **param_label_<param>**: <param> is passed to :meth:`add_param_label`
            - **close**: Close the figure at the end [default: False]
            - **title_<param>**: <param> is passed to :func:`~matplotlib.pyplot.title`
            - **logo_<param>**: <param> is passed to :func:`add_logo`
            - **tight_layout**: To make a tight layout one everything is plotted.
        """
        self._post_plotted = True

        # Format X and Y axis
        self.format_axes(**kwargs)

        # Filter kewords
        kw = {}
        for kwtype in ['grid', 'title', 'hlitvs', 'hldays', 'dayhl', 'finalize',
            'figtext', 'key', 'savefig', 'savefigs', 'show', 'legend',
            'tight_layout', 'param_label', 'autoresize']:
            kw[kwtype] = kwfilter(kwargs, kwtype+'_')
            if (kwtype in self._primary_attributes+self._secondary_attributes+
                    self._special_attributes and kw[kwtype].has_key(kwtype)):
                del kw[kwtype][kwtype]
        kwanim = kwfilter(kwargs, 'anim_', keep=True)
        kw['show'].update(**kwanim)
        kw['hlitvs'].update(kw['dayhl']) # compat
        kw['hlitvs'].update(kw['hldays']) # compat

        # Resize plot
        autoresized = self.autoresize(autoresize, **kw['autoresize'])

        # Anchor
        if anchor is None and autoresized:
            anchor='C'
        if anchor is not None:
            self.axes.set_anchor(anchor)

        # Grid
        self.grid(grid, **kw['grid'])

        # Highlight intervals
        if kwargs.pop('hldays', hlitvs):
            if kwargs.pop('hldays', False):
                kw['hlitvs'].setdefault['units'] = 'day'
            self.hlitvs(**kw['hlitvs'])

        # Title
        self.ptitle(title, **kw['title'])

        # Fig text
        self.figtext(figtext, **kw['figtext'])

        # Key of axes
        self.add_key(key, **kw['key'])

        # Params
        if param_label:
            self.add_param_label(param_label, **kw['param_label'])

        # Legend
        if legend:
            self.legend(**kw['legend'])

        # Tight layout
        if tight_layout:
            self.fig.tight_layout(**kw['tight_layout'])

        # User finalization
        if callable(finalize):
            dict_check_defaults(kw['finalize'], fig=self.fig, ax=self.axes)
            finalize(**kw['finalize'])

        # Save it
        self.savefig(savefig, **kw['savefig'])
        self.savefigs(savefigs, **kw['savefigs'])

        # Show or close
        if show:
            self.show(**kw['show'])

        if close:
            self.close()



    def format_axes(self, tz=None, nodate=False, date_rotation=None, date_fmt=None,
            date_locator=None, date_minor_locator=None, date_nominor=False,
            log=False, **kwargs):
        """Scale and format X and Y axes

        :Params:

            - **x/y/vskip**, optional: Skip axis formating.
            - **nodate**, optional: do not format as date.
            - **date_rotation**, optional: Rotate date labels.
            - **date_fmt**, optional: Date format (like ``"%s/%m/%Y"``).
            - **date_locator**, optional: Major locator (see :func:`setup_time_axis`).
            - **date_minor_locator**, optional: Minor locator (see :func:`setup_time_axis`).
            - **date_nominor**, optional: Do not plot minor localor.
            - **x/y/vmin/max**, optional: Force min/max of X or Y axis (defaults to :attr:`xmin`, etc).
            - **x/y/vlim**, optional: Force min/max of X or Y axis with `(min,max)`` like argument.
            - **x/y/vminmax**, optional: Minimal max value
              -> use this value if max is too low.
            - **x/y/vmaxmin**, optional: Maximal min value.
            - **x/yticks**, optional: Position of ticks.
            - **x/yfmt** (or **...format**, **...tickfmt**, **...tickformat**, optional: Format of ticks.
            - **x/yticklabels**, optional: Label of ticks.
            - **x/yhide**, optional: Hide labels.
            - **x/ynmax** (or **...nmax_ticks***), optional: Max number of ticks for some locators.
            - **x/y/vtitle** (or **..label**), optional: Title of the axis (defaults to :attr:`xlabel`, etc).
        """
        vkwargs = kwfilter(kwargs, 'v', copy=True, short=True)
        for iaxis, xy in enumerate(('x', 'y')):

            if kwargs.get(xy+'skip'):
                continue

            # Base
            axis = getattr(self.axes, xy+'axis')
            data = self.get(xy+'data')
            iorder = 1-iaxis

            # Keywords
            defaults = {}
            props = {}
            props['type'] = self.order[1-iaxis]
            vatts = ['min', 'max', 'minmax', 'maxmin', ('label', 'title'), 'units', ]
            aatts = ['strict', ('fmt', 'format', 'tickfmt', 'tickformat'), ('rotation', 'rotate'),
                'lim', 'hide', 'ticks', 'ticklabels','locator', 'minor_locator', 'formatter', 'minor_formatter',
                ('nmax', 'nmax_ticks'), 'scale', 'minutes']
            defaults = dict(minutes=True)#strict=True)
            if props['type']!='d':
                defaults['strict'] = True
            xykwargs = kwfilter(kwargs, xy, copy=1, short=True)
            for raw_att in aatts+vatts:

                # Merge aliases
                if isinstance(raw_att, tuple):
                    dict_aliases(xykwargs, raw_att)
                    att = raw_att[0]
                else:
                    att = raw_att

                # Default values
                # - base
                default = None
                # - get it from cdms variable
                if props['type']=='d' and raw_att in vatts:
                    kwvar = dict_aliases(vkwargs, raw_att)
                    if not xykwargs.has_key(att) and kwvar.has_key(att):
                        default = kwvar[att]
                # - special
                if defaults.has_key(att):
                    default = defaults[att]

                # Get attribute
                props[att] = xykwargs.get(att, default)

            props['label_kwargs'] = kwfilter(xykwargs, 'label_', copy=1)
            props['label_kwargs'].update(kwfilter(xykwargs, 'title_', copy=1))
            props['fmt_kwargs'] = kwfilter(xykwargs, 'fmt_', copy=1)
            props['locator_kwargs'] = kwfilter(xykwargs, 'locator_', copy=1)
            props['ticklabels_kwargs'] = kwfilter(xykwargs, 'ticklabels_', copy=1)
            if props['type']=='d' and log is True:
                props['locator'] = False
                props['minor_locator'] = False

            # Axis scale
            if props['scale'] is not None:
                getattr(self.axes, 'set_%sscale'%xy)(props['scale'])

            # Limits
            axmin, axmax = None, None
            # - wee clearly need strict limits
            if props['strict'] is None:
                props['strict'] = props['type'] == 't'
            if props['strict']:
                axmin = getattr(self, 'get_%smin'%xy)(glob=True)
                axmax = getattr(self, 'get_%smax'%xy)(glob=True)
            # - use of the form xlim/ylim
            if props['lim']:
                if isinstance(props['lim'],dict):
                    if props['lim'].has_key(xy+'min'):
                        axmin = props['lim'][xy+'min']
                    elif props['lim'].has_key('min'):
                        axmin = props['lim']['min']
                    if props['lim'].has_key(xy+'max'):
                        axmax = props['lim'][xy+'max']
                    elif props['lim'].has_key('max'):
                        axmax = props['lim']['max']
                else:
                    axmin, axmax = props['lim']
            # - we directly specify xmin/xmax...
            if props['min'] is not None:
                axmin = props['min']
            elif self.isset(xy+'min') or (self.order[iorder]=='d' and self.isset('vmin')):
                axmin = getattr(self, 'get_%smin'%xy)()
            if props['max'] is not None:
                axmax = props['max']
            elif self.isset(xy+'max') or (self.order[iorder]=='d' and self.isset('vmax')):
                axmax = getattr(self, 'get_%smax'%xy)()
            # - we put min and max in boundary
            if props['minmax'] is not None:
                if axmin is None:
                    axmin = getattr(self, 'get_%smin'%xy)(glob=True)
                axmin = min(axmin, props['minmax'])
            if props['maxmin'] is not None:
                if axmax is None:
                    axmax = getattr(self, 'get_%smax'%xy)(glob=True)
                axmax = max(axmax, props['maxmin'])
            # - ok, lets set it
            xylim_func = getattr(self.axes, 'set_%slim'%xy)
            if axmin is not None and axmin is not False:
                if props['type'] == 't' and (isinstance(axmin, basestring) or not N.isscalar(axmin)):
                    try:
                        axmin = mpl(axmin)
                    except:
                        raise PlotError("Can't convert axmin to date: %s"%axmin)
                xylim_func(**{xy+'min':_asnum_(axmin)})
            if axmax is not None and axmax is not False:
                if props['type'] == 't' and (isinstance(axmax, basestring) or not N.isscalar(axmax)):
                    try:
                        axmax = mpl(axmax)
                    except:
                        raise PlotError("Can't convert axmax to date: %s"%axmin)
                xylim_func(**{xy+'max':_asnum_(axmax)})

            # Ticks values and formats
            if props['type'] in ['x','y','z']: # Lon, lat or dep axis
                kwscale = dict(vmin=axmin, vmax=axmax, geo=props['minutes'])
                kwlf = props['fmt_kwargs']
                lab_val = props['ticks']
                if hasattr(lab_val, '__len__') and len(lab_val)==0:
                    lab_val = None
                if props['type'] in ['x','y']:
                    if lab_val is None:
                        lab_val = geo_scale(data, nmax=props['nmax'], **kwscale)
                    if props['type'] == 'x':
                        lab_func = lonlab
                    else:
                        lab_func = latlab
                    kwlf.setdefault('decimal', not props['minutes'])
                    kwlf.setdefault('auto_minutes', int(min(lab_val))!=int(max(lab_val)))
                else:
                    if lab_val is None:
                        lab_val = auto_scale(data, **kwscale)
                    lab_func = deplab
                if callable(props['fmt']):
                    lab_func = props['fmt']
                    props['fmt'] = None
                props['ticks'] = None #FIXME: bad choice (due to labels?)
                axis.set_ticks(lab_val)
                if not self.has_data(): # set min/max from ticks when no data available
                    if axmin is None :
                        axmin = lab_val[0]
                    if axmax is None:
                        axmax = lab_val[-1]
#                kwtf = {}
                if props['fmt'] is not None:
                    axis.set_ticklabels([props['fmt']%l for l in lab_val], **props['ticklabels_kwargs'])
                else:
#                    kwtf['fmt'] = props['fmt']
                    axis.set_ticklabels(lab_func(lab_val, **kwlf), **props['ticklabels_kwargs'])
                props['locator'] = None

            elif props['type'] == 't' and not nodate: # Time axis

                # Rotate dates
#                if date_rotation is None:
#                    if xy=='y':
#                        date_rotation = 0.
#                    else:
#                        date_rotation = 30.

                # Ok let's format
                trange = data.ptp()
                kwdate = kwfilter(kwargs,'date', copy=1)
                kwdate.setdefault('nmax_ticks', props['nmax'])
                kwdate.setdefault('auto', axmin is None or axmax is None)
                setup_time_axis(axis, rotation=date_rotation,fmt=date_fmt,locator=date_locator,
                    minor_locator=date_minor_locator,nominor=date_nominor,trange=trange,
                    nodual=date_rotation is not None, **kwdate)
                props['locator'] = None

            else: # Other axes
                if props['fmt'] is not None:
                # Numeric format
                    axis.set_major_formatter(FormatStrFormatter(props['fmt']))
            if props['locator']:
                locator = axis.set_major_locator(props['locator'])
                if props['locator_kwargs']:
                   locator.set_params(**props['locator_kwargs'])
            if props['minor_locator'] is not False and props['type'] not in 'xyzt':
                props['minor_locator'] = AutoMinorLocator()
            if props['minor_locator']:
                axis.set_minor_locator(props['minor_locator'])
            if props['formatter']:
                axis.set_major_formatter(props['formatter'])
            if props['minor_formatter']:
                axis.set_minor_formatter(props['minor_formatter'])
            if props['type'] != 't' or nodate:
                if props['rotation'] is not None:
                    P.setp(axis.get_ticklabels(), 'rotation', props['rotation'])
            if nodate: props['type'] = '-'
            if props['ticks'] is not None:
                axis.set_ticks(props['ticks'])
                if not self.has_data(): # set min/max from ticks when no data available
                    if axmin is None :
                        axmin = props['ticks'][0]
                    if axmax is None:
                        axmax = props['ticks'][-1]
            if props['ticklabels'] is not None:
                if props['ticklabels'] is False:
                    props['ticklabels'] = []
                axis.set_ticklabels(props['ticklabels'], **props['ticklabels_kwargs'])

            # Hiding
            if self._xyhide_(xy, props['hide']):
                P.setp(axis.get_ticklabels(), 'visible', not props['hide'])

            # Label
            elif props['label'] or getattr(self, xy+'label') and  axis.get_label().get_text()=='':
                if props['label'] is None:
                    props['label'] = getattr(self, xy+'label')
                if props['label'] is not False:
                    func = getattr(self.axes, 'set_%slabel'%xy)
                    func(props['label'], **props['label_kwargs'])

            # Now set again min/max
            if axmin is not None and axmax is not False:
                xylim_func(**{xy+'min':_asnum_(axmin)})
            if axmax is not None and axmax is not False:
                xylim_func(**{xy+'max':_asnum_(axmax)})

            # Properties
            if props['ticklabels_kwargs']:
                P.setp(axis.get_ticklabels(), **props['ticklabels_kwargs'])
            if props['label_kwargs']:
                P.setp(axis.get_label(), **props['label_kwargs'])

    def _xyhide_(self, xy, hide):
        if not isinstance(hide, basestring):
            return hide
        hide = hide.lower()
        if hide=='auto' or hide=='subplot':
            if not isinstance(self.axes, Subplot): return False
            if xy=='x':
                return not self.axes.is_last_row()
            return not self.axes.is_first_col()
        return False

    def autoresize(self, autoresize=True, minaspect=None):
        """Resize figure or axes to fit to data axes"""
        if (autoresize and self.axes.get_aspect() != 'auto' and
                isinstance(self.axes, Subplot) and
                self.axes.is_first_col() and self.axes.is_first_row() and
                self.axes.is_last_col() and self.axes.is_last_row()):

            r = self.axes.get_aspect() # dy/dx
            if r=='equal':
                r = 1.
            r *= self.axes.get_data_ratio()
    #       rect = self.axes.get_position(True)
    #       if autoresize == 2: r *= rect.width/rect.height

            sp = self.fig.subplotpars
            x0 = sp.left
            x1 = sp.right
            y0 = sp.bottom
            y1 = sp.top
            Dx = x1 - x0
            Dy = y1 - y0
            R = r * Dx / Dy

            if minaspect is not None:
                R = N.clip(R, minaspect, 1/minaspect)


            w, h = self.fig.get_size_inches()

            if autoresize=='x': # Resize x only, change the surface
                H = h
                W = H / R
            elif autoresize=='y': # Resize y only, change the surface
                W = w
                H = R * W
            else: # Resize both x and y without changing the surface
                a = 1.*w*h
                W = N.sqrt(a / R)
                H = r * W

            self.fig.set_size_inches(W, H ,forward=True)
            return True
        return False

    def _transform_(self, transform, default=None):
        return _transform_(transform, default, ax=self.axes, fig=self.fig)

    def get_xy(self, x, y, transform=None, xyscaler=None,
               default_transform=None, atleast_1d=False):
        """Convert (x,y) in data coordinates

        :Params:

            - **x/y**: Coordinates referenced to data, axes or figure.
            - **transform**, optional: Transform applied to coordinates.
              This either a :class:`matplotlib.transforms.Transform` or a string:
              ``"data"``, ``"axes"``, ``"figure"``.
            - **xyscaler**, optional: Converter of coordinates used when input
              coordinates are in data coordinates. It must be a callable,
              and it defaults to attribute :attr:`xyscaler` if existing.
              It converts for instance from degrees to meters for :class:`Map`
              instances. If equal to False, no conversion is performed.
        """
        transform = self._transform_(transform, default_transform)

        # From figure or axes coordinates
        if transform is self.fig.transFigure or transform is self.axes.transAxes:
            return tuple((transform-self.axes.transData).transform_point((x, y)))

        # From data coordinates (check xyscaler only)
        if (transform is self.axes.transData or transform in self.axes.transData._parents.values() \
            or str(self.axes.transData) in str(transform)) and xyscaler is not False: # Add transform from ie degrees to m
            if xyscaler is None and hasattr(self, 'xyscaler'): xyscaler = self.xyscaler
            x = _asnum_(x)
            y = _asnum_(y)
            if xyscaler:
                x, y = xyscaler(x, y)

        return x, y


    def get_transoffset(self, x, y, units='points', transform='data'):
        """Return a translation :class:`~maplotlib.transforms.Transform`

        It can be used for instance to plot an object with an
        offset with respect to its specified position.

        :Params:

            - **x/y**: Relative position.
            - **units**, optional: Units ("points", "inches", "pixels", ...)
            - **transform**, optional: Base transform for reference position.
              Choose for instance "data" or "axes".


        :Example:

            >>> o = Plot2D(data)
            >>> o.add_point(-4, 43)
            >>> t = o.get_transoffset(0, 10)
            >>> o.text(-4, 43, transform=t)

        :See also: :func:`~maplotlib.transforms.offset_copy`
        """
        transform = self._transform_(transform, 'data')
        return offset_copy(transform, fig=self.fig, x=x, y=y, units=units)

    def grid(self, b=True, **kwargs):
        """Add a grid to axes using :func:`~matplotlib.pyplot.grid`


        :Example:

            >>> myplot.grid(color='r')
            >>> myplot.grid(False)
        """
#        grid = self.get_axobj('grid')
#        print 'grid',grid
#        if grid is None:
        self.axes.grid(b=b, **kwargs)
#           print grid
#           self.set_axobj('grid', grid)
#           print self.get_axobj('grid')
#        print 'done'
#        return grid

    def add_figtext(self, *args, **kwargs):
        """Add text to the current figure using :func:`~matplotlib.pyplot.figtext`

        :Defaults:

            - Position: defaults to the top center.
            - Alignement: ``ha="center", va="top"``

        :Example:

            >>> myplot.figtext('Group of plots')
            >>> myplot.figtext(0.2, 0.92, 'My plots', color='b', ha='left', va='center')
        """
        # Arguments
        if len(args)==0: return
        if isinstance(args[0], (list, tuple)):
            args = args[0]
        if len(args) == 1:
            text = args[0]
            if text is None: return
            if not isinstance(text, dict):
                x = .5
                y = .98
            else:
                x = 0
                y = 0
        else:
            x, y, text = args

        # Keywords
        if isinstance(text, dict):
            text = ' '.join('%s=%s'%(key, text[key]) for key in sorted(text.keys()))
            defaults = dict(ha='left', va='bottom', size=9, color='.4', family='monospace')
        else:
            defaults = dict(ha='center', va='center', size='12')
        for key, val in defaults.items():
            kwargs.setdefault(key, val)

        # Plot
        return self.add_obj('figtext', self.fig.text(x, y, text, **kwargs))
    figtext = add_figtext # backward compat

    def add_text(self, x, y, text, transform='axes', shadow=False, glow=False,
        xyscaler=None, strip=True, **kwargs):
        """Add text to the plot axes

        :Params:

            - **x,y**: Coordinates of the text.
            - **text**: Text to plot.
            - **transform**, optional: Type of coordinates
              (like ``"axes"`` or ``"data"``).
            - **shadow**, optional: Add a droped shadow below the text
              (see :func:`add_shadow`).
            - **shadow_<param>**, optional: ``<param>`` is passed to :func:`add_shadow`.
            - **glow**, optional: Add a glow effect the text
              (see :func:`add_glow`).
            - **glow_<param>**, optional: ``<param>`` is passed to :func:`add_glow`.
            - Other keywords are passed to :func:`matplotlib.pyplot.text`.

        """
        # Keywords
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')

        # Coordinates transform
        transform = self._transform_(transform, 'axes')
        if transform not in [self.axes.transAxes, self.fig.transFigure]:
            x, y = self.get_xy(x, y, transform, xyscaler=xyscaler)

        # Plot
        if strip: text = text.strip()
        obj = self.add_axobj('text', self.axes.text(x, y, text, transform=transform, **kwargs))
        self.register_obj(obj, **kwargs)

        # Path effects
        if shadow: self.register_obj(self.add_shadow(obj, **kwsh), **kwargs)
        if glow: self.register_obj(self.add_glow(obj, **kwgl), **kwargs)

        return obj
    text = add_text # backward compat

    def add_water_mark(self, text, x=.5, y=.5, ha='center', va='center', size=20,
        color='k', alpha=.7, zorder=0, transform='axes', **kwargs):
        """Add a background text to the plot

        All arguments are passed to :meth:`add_text`.
        """
        return self.add_text(x, y, text, ha=ha, va=va, size=size, color=color,
            zorder=zorder, transform=transform, **kwargs)

    def add_time_label(self, x, y, mytime, fmt="%Y/%m/%d %H:%M", **kwargs):
        """Add a time label

        See :meth:`add_text` for other keywords
        """
        text = strftime(fmt, mytime)
        kwargs.setdefault('family', 'monospace')
        return self.add_text(x, y, text, **kwargs)

    def add_lon_label(self, x, y, mylon, fmt='%5g', **kwargs):
        """Add a longitude label

        See :meth:`add_text` for other keywords
        """
        text = lonlab(mylon, fmt=fmt)
        return self.add_text(x, y, text, **kwargs)

    def add_lat_label(self, x, y, mylat, fmt='%5g', **kwargs):
        """Add a longitude label

        See :meth:`add_text` for other keywords
        """
        text = latlab(mylat, fmt=fmt)
        return self.add_text(x, y, text, **kwargs)

    def add_left_label(self, text, **kwargs):
        """Add a text label to the left of the plot

        :See also: :func:`~vacumm.misc.plot.add_left_label`
        """
        kwargs['ax'] = self.axes
        return add_left_label(text, **kwargs)

    def add_right_label(self, text, **kwargs):
        """Add a text label to the right of the plot

        :See also: :func:`~vacumm.misc.plot.add_right_label`
        """
        kwargs['ax'] = self.axes
        return add_right_label(text, **kwargs)

    def add_top_label(self, text, **kwargs):
        """Add a text label to the top of the plot

        :See also: :func:`~vacumm.misc.plot.add_top_label`
        """
        kwargs['ax'] = self.axes
        return add_top_label(text, **kwargs)

    def add_bottom_label(self, text, **kwargs):
        """Add a text label to the bottom of the plot

        :See also: :func:`~vacumm.misc.plot.add_bottom_label`
        """
        kwargs['ax'] = self.axes
        return add_bottom_label(text, **kwargs)

    def add_key(self, key, **kwargs):
        """Add a key to specify the plot number

        See :func:`~vacumm.misc.plot.add_key` for more information.
        """
        if key is False: return
        from vacumm.misc.plot import add_key
        kwargs.update(fig=self.fig, axes=self.axes)
        return self.set_axobj('key', add_key(key, **kwargs))

    def add_param_label(self, text, **kwargs):
        """Add parameters description to the bottom/left of the figure using
        :func:`~vacumm.misc.plot.add_param_label`

        :Example:

            >>> c = curve2(sst, show=False)
            >>> c.add_param_label(dict(max=.23, kz=0.25), color='r')

        :Params:

            - **text**: Either a string or a dictionary.
            - See :func:`~vacumm.misc.plot.add_param_label` for other parameters
        """
        kwargs['fig'] = self.fig
        return add_param_label(text, **kwargs)

    def add_annotation(self, x, y, xtext, ytext, text='', xycoords='data',
            textcoords='offset points', arrowprops='->',
            shadow=False, glow=False,
            xyscaler=None, strip=True, **kwargs):
        """Add an annotation to the plot axes using :func:`matplotlib.pyplot.annotate`

        :Params:

            - **x,y**: Coordinates of the text.
            - **text**: Text to plot.
            - **xycoords/transform**, optional: Type of coordinates of point
              (like ``"axes"`` or ``"data"``).
            - **textcoords**, optional: Type of coordinates of text
              (like ``"axes"`` or ``"data"``).
            - **arrowprops**, optional: Dictionary of arrow properties or
              string thet defines the arrow style.
            - **shadow**, optional: Add a droped shadow below the text
              (see :func:`add_shadow`).
            - **shadow_<param>**, optional: ``<param>`` is passed to :func:`add_shadow`.
            - **glow**, optional: Add a glow effect the text
              (see :func:`add_glow`).
            - **glow_<param>**, optional: ``<param>`` is passed to :func:`add_glow`.
            - Other keywords are passed to :func:`matplotlib.pyplot.annotate`.

        """
        # Keywords
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')

        # Coordinates transform
        xycoords = kwargs.pop('transform', xycoords)
        xycoords = self._transform_(xycoords, 'data')
        if not isinstance(xycoords, basestring) and \
                xycoords not in [self.axes.transAxes, self.fig.transFigure]:
            x, y = self.get_xy(x, y, xycoords, xyscaler=xyscaler)
        if textcoords is None:
            textcoords = xycoords
        else:
            textcoords = self._transform_(textcoords, 'offset points')
        if not isinstance(textcoords, basestring) and \
                textcoords not in [self.axes.transAxes, self.fig.transFigure]:
            xtext, ytext = self.get_xy(xtext, ytext, textcoords, xyscaler=xyscaler)

        # Arrow properties
        if isinstance(arrowprops, basestring):
            arrowprops = dict(arrowstyle=arrowprops)
        elif not isinstance(arrowprops, dict):
            arrowprops = {}
        arrowprops.setdefault('arrowstyle', '->')

        # Plot
        if strip: text = text.strip()
        obj = self.axes.annotate(xy=(x, y), xytext=(xtext, ytext), s=text,
            xycoords=xycoords, textcoords=textcoords, arrowprops=arrowprops, **kwargs)
        obj = self.add_axobj('text', obj)
        self.register_obj(obj, **kwargs)

        # Path effects
        if shadow: self.register_obj(self.add_shadow(obj, **kwsh), **kwargs)
        if glow: self.register_obj(self.add_glow(obj, **kwgl), **kwargs)

        return obj


    def hlitvs(self, **kwargs):
        """Highlight intervals with grey/white background alternance

        See :func:`~vacumm.misc.plot.hlitvs` for more information.
        """
        ret = []
        for i, xy in enumerate('yx'):
            if self.order[i] == 't':
                cached = self.get_axobj('hlitvs', axis=xy)
                if cached is not None:
                    self.remove(cached)
                kwargs['axis'] = getattr(self.axes, xy+'axis')
                obj = self.set_axobj('hlitvs', hlitvs(**kwargs), axis=xy)
                self.register_obj(obj, **kwargs)
                return obj

    def hldays(self, **kwargs):
        """Alias for :

            >>> myplot.hlitv(units='day')
        """
        kwargs.setdefault('units', 'day')
        return self.hlitvs(**kwargs)

    def ptitle(self, title=None, force=None, **kwargs):
        """Add a title to the plot

        .. note::

            No title is added to the plot if a title already exists
            and the specified title is guessed (not hard set).

        :Params:

            - **title**: Title to add to plot.

                - A string: directly used.
                - ``True`` or ``None``: the :attr:`title` attribute is used.
                - ``False``: not title is plotted.

            - **force**, optional: If the title is already plotted,
              it is not overwritten, except if ``force is True``.
            - Other keywords are passed to :func:`matplotlib.pyplot.title`.

        """
        if title is None:
            title = self.title
            isset = self.isset('title')
        else:
            isset = True
        if force is None: force = isset
        if title is False: return
        if title is not None and (force or self.axes.get_title()==''):
            return self.set_axobj(title, self.axes.set_title(title, **kwargs))

    def legend(self, *args, **kwargs):
        """A simple call to the :meth:`matplotlib.axes.Axes.legend` method

        Arguments and keywords are passed to :meth:`~matplotlib.axes.Axes.legend`.

        Defaults values :

            - **loc**: ``"best"``
            - **shadow**: ``False``
            - **fancybox**: ``True``
            - **alpha** (applied to legend patch): ``0.5``
        """
        kwargs.setdefault('loc', 'best')
#        kwargs.setdefault('shadow', False)
#        kwargs.setdefault('fancybox', True)
        alpha = kwargs.pop('alpha', .5)
        zorder = kwargs.pop('zorder', None)
#        zorder = kwargs.pop('framealpha', alpha)
        leg = self.axes.legend(*args, **kwargs)
        if leg is not None:
            leg.legendPatch.set_alpha(alpha)
            self.set_axobj('legend', leg)
            if zorder:
                leg.set_zorder(zorder)
        return leg

    def _get_xykeys_(self, xy):
        """Get possible keys for interval selections along X or Y"""
        keys = [xy]
        ixy = 'yx'.index(xy)
        if self.order[ixy]=='x':
            keys.extend(['lon', 'longitude'])
        elif self.order[ixy]=='y':
            keys.extend(['lat', 'latitude'])
        elif self.order[ixy]=='z':
            keys.extend(['level', 'dep', 'depth'])
        elif self.order[ixy]=='t':
            keys.extend(['time'])
        elif self.order[ixy]=='d':
            keys.extend(['data', 'value'])
        return keys

    def _get_boxminmax_(self, box):
        """Get ``xmin,ymin,xmax,ymax`` from box specs

        Two cases:

        - ``box=xmin,ymin,xmax,ymax``
        - ``box=dict(x=(xmin,xmax),y=(xmin,xmax))`` or with for instance
          ``lon``, ``lat``, ``time`` depending on axis type.
        """
        if isinstance(box, dict):
            # Along X
            for key in self._get_xykeys_('x'):
                if key in box:
                    xmin, xmax = box[key][:2]
                    break
            else:
                xmin, xmax = self.xmin, self.xmax
            # Along Y
            for key in self._get_xykeys_('y'):
                if key in box:
                    ymin, ymax = box[key][:2]
                    break
            else:
                ymin, ymax = self.ymin, self.ymax
            box = xmin, ymin, xmax, ymax
        return box


    def add_box(self, box, zorder=150, shadow=False, glow=False, color='r',
            npts=10, xyscaler=None, fill=False, **kwargs):
        """Add a box to the plot using :meth:`matplotlib.pyplots.plot`

        :Params:

            - **box**: Box limits in the forms ``[xmin,ymin,xmax,ymax]``
              ``dict(x=(xmin,xmax),y=(xmin,xmax)``.
            - **color**, optional: Line color of the box.
            - **npts**, optional: Number of points per side
              (useful with special map projections).
            - **shadow**, optional: Add a droped shadow below the box
              (see :func:`add_shadow`).
            - **shadow_<param>**, optional: ``<param>`` is passed to :func:`add_shadow`.
            - **glow**, optional: Add a glow effect the box
              (see :func:`add_glow`).
            - **glow_<param>**, optional: ``<param>`` is passed to :func:`add_glow`.
            - Other keywords are passed to :func:`matplotlib.pyplot.plot`.

        :See also: :func:`matplotlib.pyplot.plot`
        """
        # Limits
        try:
            xmin, ymin, xmax, ymax = self._get_boxminmax_(box)
        except:
            raise PlotError('Box limits should be a list in the form [xmin,ymin,xmax,ymax]'
                ' or a dictionary in the form dict(x=(xmin,xmax),y=xmin,xmax): %s'%box)
        xmin = _asnum_(xmin)
        xmax = _asnum_(xmax)
        ymin = _asnum_(ymin)
        ymax = _asnum_(ymax)
        if xyscaler is None and hasattr(self, 'xyscaler'): xyscaler = self.xyscaler
        if xyscaler:
            xmin, ymin = xyscaler(xmin, ymin)
            xmax, ymax = xyscaler(xmax, ymax)

        # Params
        kwargs.update(ec=color, zorder=zorder)
        if not fill:
            kwargs.update(fc='none')
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        kwsh.setdefault('zorder', zorder)
        kwgl.setdefault('zorder', zorder)
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')

        # Sides
        X = N.linspace(xmin, xmax, npts),
        X += N.linspace(xmax, xmax, npts)[1:],
        X += N.linspace(xmax, xmin, npts)[1:],
        X += N.linspace(xmin, xmin, npts),
        Y = N.linspace(ymin, ymin, npts),
        Y += N.linspace(ymin, ymax, npts)[1:],
        Y += N.linspace(ymax, ymax, npts)[1:],
        Y += N.linspace(ymax, ymin, npts),

        # Plot
        o = self.axes.fill(N.concatenate(X), N.concatenate(Y), **kwargs)
        self.register_obj(o, **kwargs)

        # Effects
        if shadow: self.add_shadow(o,**kwsh)
        if glow: self.add_glow(o, **kwgl)


        return o

    def add_arrow(self, x, y, udata, vdata, zorder=150, polar=False, degrees=True,
            shadow=False, glow=False, quiverkey=False, xyscaler=None,
            color=False, **kwargs):
        """Add an arrow to the map using :func:`matplotlib.pyplot.quiver`

        :Params:

            - **x,y**: Coordinates of the position of the tail
            - **udata**: X or radial component of arrows.
            - **vdata**: Y or directional component of arrows.
            - **polar**, optional: Consider polar coordinates: ``(u, v) -> (rho, theta)``
            - **degrees**, optional: If True (default), trat ``theta`` as degrees, else radians.
            - **quiver_<param>**, optional: ``<param>`` is passed to :func:`matplotlib.pyplot.quiver`.
            - Other keywords are passed to :func:`matplotlib.pyplot.plot`.

        :See also: :func:`matplotlib.pyplot.scatter`
        """
        # Data for arrows
        udata = MV2.asarray(udata)
        vdata = MV2.asarray(vdata)
        if polar:
            u, v = udata,vdata
            m = u
            angle = (v*N.pi/180.) if degrees else v
            u = m * MV2.cos(angle)
            v = m * MV2.sin(angle)
            del angle
            if hasattr(m, 'units'):
                u.units = v.units = m.units
            if hasattr(m, 'long_name'):
                u.long_name = 'X component of '+m.long_name
                v.long_name = 'Y component of '+m.long_name
            # Data
            udata, vdata = u,v

        # Coordinates
        x = _asnum_(x, atleast_1d=True)
        y = _asnum_(y, atleast_1d=True)
        if xyscaler is None and hasattr(self, 'xyscaler'): xyscaler = self.xyscaler
        if callable(xyscaler):
            x, y = xyscaler(x, y)

        # Params
        kwargs.update(zorder=zorder)
        kwqv = kwfilter(kwargs,'quiver')
        kwqvk = kwfilter(kwargs,'quiverkey')
        #dict_copy_items(kwargs, [kwqv],'anim')
        kwargs = dict_merge(kwargs, kwqv)
        args = []

        # Color
        if color is True:
            color = N.ma.sqrt(u**2+v**2)
        if isinstance(color, (N.ndarray, list)):
            args.append(color)
        elif color is not False and color is not None:
            kwargs['color'] = color

        # Plot
        o = self.axes.quiver(x, y, udata, vdata, *args, **kwargs)
        self.register_obj(o, **kwargs)
        if quiverkey:
          self.quiverkey(o, **kwqvk)
        return o


    def quiverkey(self, qv, value, pos=(0.,1.02), text='%(value)g %(units)s',
            units=None, latex_units=None, value_mode=80, **kwargs):
        """Add a quiver key to the plot

        :Params:

            - **qv**: Results of :func:`~matplotlib.pyplot.quiver`.
            - **value**: Numeric value for key (used by text).
            - **pos**, optional: Position of key for arrow .
            - **text**, optional: Text or format with variables 'value' and 'units'.
            - **units**, optional: Units for key (used by text).
            - **latex_units**, optional: Interpret units using latex.
            - Extra keywords are passed to :func:`~matplotlib.pyplot.quiverkey`.
        """

        # Value
        value = get_quiverkey_value(value, mode=value_mode)

        # Text
        if units is None:
            units = self.quiverkey_units
        elif cdms2.isVariable(units) and  hasattr(units,'units'):
            units = units.units
        elif not isinstance(units, basestring):
            units = ''
        latex_units = kwargs.pop('tex', None)
        if latex_units is None: latex_units = self.latex_units
        if units is None:
            units = ''
        if latex_units and not self.is_latex(units):
            units = '$%s$'%units
        try:
            text = text % vars()
        except:
            text = '%(value)g' % value

        # Plot
        pos = kwargs.pop('loc', pos)
        qvk = self.axes.quiverkey(qv, pos[0], pos[1], value, text, **kwargs)
        return self.register_obj(qvk, 'quiverkey', **kwargs)

    def get_quiverkey_units(self):
        """Get :attr:`quiverkey_units`"""
        units = self.get_obj('quiverkey_units')
        if units is None and isinstance(self, Plot1D) and self.isset('units'):
            units = self.units
        if units is None:
            units = self.get_units(idata=[-2, -1])
        return units
    def set_quiverkey_units(self, value):
        """Set :attr:`quiverkey_units`"""
        self.set_axobj('quiverkey_units', value)
    def del_quiverkey_units(self):
        """Del :attr:`quiverkey_units`"""
        self.del_axobj('quiverkey_units')
    quiverkey_units = property(get_quiverkey_units, set_quiverkey_units,
        del_quiverkey_units, doc="Units used for quiverkey")



    def add_line(self, extents, zorder=150, shadow=False, glow=False, color='r',
        npts=10, xyscaler=None, **kwargs):
        """Add a line to the plot using :meth:`matplotlib.pyplots.plot`

        :Params:

            - **extents**: Extents in the forms ``[xmin,ymin,xmax,ymax]``
              ``dict(x=(xmin,xmax),y=xmin,xmax)``.
            - **color**, optional: Line color of the line.
            - **npts**, optional: Number of points per side
              (useful with special map projections).
            - **shadow**, optional: Add a droped shadow below the box
              (see :func:`add_shadow`).
            - **shadow_<param>**, optional: ``<param>`` is passed to :func:`add_shadow`.
            - **glow**, optional: Add a glow effect the box
              (see :func:`add_glow`).
            - **glow_<param>**, optional: ``<param>`` is passed to :func:`add_glow`.
            - Other keywords are passed to :func:`matplotlib.pyplot.plot`.

        :See also: :func:`matplotlib.pyplot.plot`
        """
        # Positions
        try:
            xmin, ymin, xmax, ymax = self._get_boxminmax_(extents)
        except:
            raise PlotError('Extents should be a list in the form [xmin,ymin,xmax,ymax]'
                ' or a dictionary in the form dict(x=(xmin,xmax),y=xmin,xmax): %s'%extents)
        xmin = _asnum_(xmin)
        xmax = _asnum_(xmax)
        ymin = _asnum_(ymin)
        ymax = _asnum_(ymax)
        if xyscaler is None and hasattr(self, 'xyscaler'): xyscaler = self.xyscaler
        if xyscaler:
            xmin, ymin = xyscaler(xmin, ymin)
            xmax, ymax = xyscaler(xmax, ymax)

        # Params
        kwargs.update(color=color, zorder=zorder)
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        kwsh.setdefault('zorder', zorder)
        kwgl.setdefault('zorder', zorder)
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')

        # Sides
        X = N.linspace(xmin, xmax, npts)
        Y = N.linspace(ymin, ymax, npts)

        # Plot
        o = self.axes.plot(X, Y, **kwargs)
        self.register_obj(o, **kwargs)

        # Effects
        if shadow: self.add_shadow(o, **kwsh)
        if glow: self.add_glow(o, **kwgl)

        return o

    def add_point(self, x, y, zorder=150, shadow=False, glow=False,
            color='r', size=20, xyscaler=None, **kwargs):
        """Add a point to the map using :meth:`matplotlib.pyplots.plot`

        :Params:

            - **x,y**: Coordinates.
            - **color**, optional: Line color of the point.
            - **shadow**, optional: Add a droped shadow below the box
              (see :func:`add_shadow`).
            - **shadow_<param>**, optional: ``<param>`` is passed
              to :func:`add_shadow`.
            - **glow**, optional: Add a glow effect the box
              (see :func:`add_glow`).
            - **glow_<param>**, optional: ``<param>`` is passed to :func:`add_glow`.
            - Other keywords are passed to :func:`matplotlib.pyplot.plot`.

        :See also: :func:`matplotlib.pyplot.scatter`
        """
        # Coordinates
        x = _asnum_(x, atleast_1d=True)
        y = _asnum_(y, atleast_1d=True)
        if xyscaler is None and hasattr(self, 'xyscaler'): xyscaler = self.xyscaler
        if callable(xyscaler):
            x, y = xyscaler(x, y)

        # Params
        kwargs.update(zorder=zorder)
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        size = kwargs.pop('s', size)
        if size is not None:
            kwargs['s'] = size
        color = kwargs.pop('c', color)
        if color is not None:
            kwargs['c'] = color
        #kwsh.setdefault('zorder', zorder-0.01)
        #kwgl.setdefault('zorder', zorder-0.01)

        # Plot
        o = self.axes.scatter(x, y, **kwargs)
        self.register_obj(o, **kwargs)

        # Effects
        if shadow: self.add_shadow(o, **kwsh)
        if glow: self.add_glow(o, **kwgl)

        return o

    def add_place(self, x, y, text, zorder=150, color='k',
            shadow=False, glow=False,
            text_offset=(0, 10), ha='center', va='center', **kwargs):
        """Place a point using :meth:`add_point` and a label using :meth:`add_text`

        :Examples:

            >>> m = map2(sst, show=False)
            >>> m.add_place(-5, 44, 'Buoy 654', text_offset=(20,0), text_ha='left',
                text_color='b', point_size=100, shadow=True)

        :Params:

            - **x/y**: Coordinates of the place in data units.
            - **text**: Name of the place.
            - **text_offset**, optional: Offset of the text in points with relative to
              coordinates.
            - **point_<param>**, optional: ``<param>`` is passed to :meth:`add_point`.
            - **text_<param>**, optional: ``<param>`` is passed to :meth:`add_text`.

        """
        kwpoint = kwfilter(kwargs, 'point_')
        kwtext = kwfilter(kwargs, 'text_')
        kwcom = dict(zorder=zorder, shadow=shadow, glow=glow, color=color)
        dict_check_defaults(kwpoint, **kwcom)
        dict_check_defaults(kwtext, ha=ha, va=va, weight='bold', **kwcom)
        tha = kwtext['ha']
        tva = kwtext['va']
        if tha=='auto':
            if text_offset[0]==0:
                tha = 'center'
            elif text_offset[0]>0:
                tha = 'left'
            else:
                tha = 'right'
        if tva=='auto':
            if text_offset[1]==0:
                tva = 'center'
            elif text_offset[1]>0:
                tva = 'bottom'
            else:
                tva = 'top'
        kwtext['ha'] = tha
        kwtext['va'] = tva
        p = self.add_point(x, y, **kwpoint)
        kwtext.setdefault('transform', self.get_transoffset(*text_offset))
        t = self.add_text(x, y, text, **kwtext)
        return p, t


    def add_lines(self, xx, yy, zorder=150, shadow=False, glow=False, color='r',
        xyscaler=None, closed=False, **kwargs):
        """Add lines to the plot using :meth:`matplotlib.axes.Axes.plot`

        :Params:


            - **xx/yy**: Coordinates (in degrees).
            - **color**, optional: Line color of the line.
            - **closed**, optional: Close the lines to form a polygon.
            - **shadow**, optional: Add a droped shadow below the box
              (see :func:`add_shadow`).
            - **shadow_<param>**, optional: ``<param>`` is passed to :func:`add_shadow`.
            - **glow**, optional: Add a glow effect the box
              (see :func:`add_glow`).
            - **glow_<param>**, optional: ``<param>`` is passed to :func:`add_glow`.
            - Other keywords are passed to :func:`matplotlib.pyplot.plot`.

        :See also: :func:`matplotlib.pyplot.plot`
        """
        # Positions
        xx = N.ma.ravel(_asnum_(xx, atleast_1d=True))
        yy = N.ma.ravel(_asnum_(yy, atleast_1d=True))
        if xx.size!=yy.size:
            raise PlotError('xx and yy must have the same number of elements (%i!=%i)'%(xx.size,yy.size))
        if xyscaler is None and hasattr(self, 'xyscaler'): xyscaler = self.xyscaler
        if xyscaler:
            xx, yy = xyscaler(xx, yy)

        # Params
        kwargs.update(color=color, zorder=zorder)
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        kwsh.setdefault('zorder', zorder)
        kwgl.setdefault('zorder', zorder)
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')

        # Plot
        o = self.axes.plot(xx, yy, **kwargs)
        self.register_obj(o, **kwargs)

        # Effects
        if shadow: self.add_shadow(o, **kwsh)
        if glow: self.add_glow(o, **kwgl)

        return o


    def add_glow(self, objs, gtypes=None, **kwargs):
        """Add glow effect to objects

        :See: :func:`add_glow`
        """
        kwargs['ax'] = self.axes
        anim = kwargs.pop('anim', None)
        return self.register_obj(add_glow(objs, **kwargs), gtypes, anim=anim)

    def add_shadow(self, objs, gtypes=None, **kwargs):
        """Add shadow to objects

        :See: :func:`add_shadow`
        """
        kwargs['ax'] = self.axes
        anim = kwargs.pop('anim', None)
        return self.register_obj(add_shadow(objs, **kwargs), gtypes, anim=anim)

    def savefig(self, figfile, verbose=False, mkdir=True, **kwargs):
        """Save the figure to a file

        :Params:

            - **figfile**: Figure file name. Also accepts
              :class:`~vacumm.misc.remote.OutputWorkFile`.
            - **verbose**, optional: Informs about file name when written.
            - **mkdir**, optional: Make figure directory if it does not exists.
            - Other keywords are passed to :func:`matplotlib.pyplot.savefig`.
        """
        if not figfile:
            return

        # List of files
        if isinstance(figfile, (list, tuple)):
            oo = []
            for ff in figfile:
                oo.append(self.savefig(ff, verbose=verbose, mkdir=mkdir, **kwargs))
            return oo

        # Remote output file
        rem = figfile if isinstance(figfile, OutputWorkFile) else False
        if rem: figfile = figfile.local_file

        # Extension
        backend = P.get_backend().lower()
        ext_standalone = ['ps','pdf','svg']
        ext_others = ['png','gif','jpg','jpeg','bmp']
        basename, ext = os.path.splitext(figfile)
        if ext == '':
            ext = 'png'
            figfile += '.png'
        elif ext == ".py":
            ext = 'png'
            figfile = basename+'.png'
        else:
            ext = ext[1:]
        if not rem and ext.lower() not in ext_others:
            if backend not in ext_standalone and ext.lower() not in ext_standalone:
                figfile += '.png'
            elif backend in ext_standalone and not ext.lower().endwith(backend):
                figfile += '.'+backend

        # Directory
        figdir = os.path.dirname(figfile)
        if mkdir and figdir and not os.path.exists(figdir):
            os.makedirs(figdir)

        # Save
        self.fig.savefig(figfile, **kwargs)
        if verbose: print 'Saved plot to '+figfile
        self._last_figfile = figfile

        # Transfer
        if rem: rem.put()

        return figfile

    def savefigs(self, figfile, **kwargs):
        """Save a figure to png (and optionaly) pdf files using :func:`~vacumm.misc.plot.savefigs`"""
        if figfile is None: return
        from vacumm.misc.plot import savefigs
        savefigs(figfile, fig=self.fig, **kwargs)

    def show(self, **kwargs):
        """Show the current figure

        If the backend does not allow showing the figure using
        :func:`matplotlib.pyplot.show`, it uses an external viewer
        after saving the figure to a temporary file using
        :func:`matplotlib.pyplot.savefig`
        """
        kwanim = kwfilter(kwargs, 'anim_')
        viewers = {
            'pdf':['/usr/bin/kghostview',
                '/usr/bin/evince',
                '/usr/bin/xpdf',
                '/usr/local/bin/acroread'],
            'ps':['/usr/bin/kghostview',
                '/usr/bin/evince',
                '/usr/bin/ghostview'],
            'svg':['/usr/bin/svgdisplay',
                '/usr/bin/konqueror']}
        backend = P.get_backend().lower()
        if backend in viewers.keys():
            if savefig is None:
                tmpfig = mktemp(suffix='.'+backend)
                self.fig.savefig(tmpfig, **kwargs)
            else:
                tmpfig = savefig
            for viewer in viewers[backend]:
                if os.path.exists(viewer):
                    cmd = '%s %s' % (viewer,tmpfig)
                    try:
                        os.system(cmd)
                    except:
                        raise PlotError('Unable to view file with command: '+str(cmd))
                    if savefig is None: os.remove(tmpfig)
                    return
        else:
            # Animation
            if self.anim is not False:
                dict_check_defaults(kwanim, interval=50, repeat_delay=3000, blit=True)
                self.animation = self.animator.make_animation(**kwanim)

            # Show
            P.show()

    def has_data(self):
        """Guess if object has data"""
        return self.data is not None
    def has_valid_data(self):
        """Guess if object has unmasked data"""
        return self.data is not None and not self.masked

    # Some properties

    def get_anim(self):
        """Get the :attr:`anim` attribute`"""
        return getattr(self, '_anim', False)
    def set_anim(self, anim):
        """Set the :attr:`anim` attribute`"""
        self._anim = anim
    anim = property(get_anim, set_anim, doc="""Is each plot saved for final animation?""")

    def get_animator(self):
        """Get the current :class:`Animator` instance or None"""
        if self.anim is False: return
        if not hasattr(self, '_animator'):
            if not hasattr(self.fig, '_vacumm_animator'):
                self.fig._vacumm_animator = Animator(self.fig)
            self._animator = self.fig._vacumm_animator
        return self._animator
    def set_animator(self, animator):
        """Set the current :class:`Animator` instance"""
        self._animator = animator
    animator = property(get_animator, set_animator, doc="""Current :class:`Animator` instance or None""")

    def animator_append(self, obj, anim=None):
        """Append an object to current :class:`Animator`"""
        if self.fig is None: return
        if anim is None: anim = self.anim
        if anim is False: return
        anim = self.animator.append(obj, anim)
        if self.anim is True:
            self.anim = anim
        return anim

    def register_obj(self, obj, gtypes=None, anim=None, **kwargs):
        """Register an object with :meth:`add_obj` and :meth:`animator_append`"""
        if gtypes is None:
            gtypes = []
        elif not isinstance(gtypes, list):
            gtypes = [gtypes]
        if 'plotted' not in gtypes:
            gtypes.append('plotted')
        self.add_obj(gtypes, obj)
        self.animator_append(obj, anim=anim)
        return obj


    def get_uvlat(self, default=45.):
        """Get :attr:`uvlat`

        If a latitude axis is available on X or Y plot axis,
        its mean value is used, else it defaults to ``default``
        """
        lat = self.get_obj('uvlat')
        if self.has_data() and 'y' in self.order:
            return [self.y[:], self.x[:]][self.order.index('y')].mean()
        return lat if lat is not None else default
    def set_uvlat(self, value):
        """Set :attr:`uvlat`"""
        self.set_obj('lat', value)
    def del_uvlat(self):
        """Del :attr:`uvlat`"""
        self.del_obj('uvlat')
    uvlat = property(get_uvlat, set_uvlat, del_uvlat,
        doc="Latitude used for guessing :attr:`uvscaler`")

    def get_uvscaler(self, guess=True, lat=None, raw=False):
        """Get :attr:`uvscaler`

        :Params:

            - **guess**, optional: Guess scaler from axis types
              and data units if not specified.
            - **lat**, optional: Latitude value passed to
              :meth:`get_metric_scale` to guess plot axis metric scale.

        """
        uvscaler = self.get_obj('uvscaler')
        if uvscaler is not None:
            return self._uvscaler_(uvscaler, raw=raw)
        if not guess: return

        # Skip data mode
        if 'd' in self.order: return

        # Guess
        if uvscaler is None:

            # Latitude
            if lat is not None:
                self.uvlat = lat
            else:
                lat = self.uvlat

            # Along X
            xmscale = self.get_metric_scale('x', lat=lat)
            if xmscale is None: return
            uunits = getattr(self.data[1], 'units',  None)
            if uunits is None: return
            umscale = tometric(uunits, munits='m/s')
            if umscale is None: return

            # Along Y
            ymscale = self.get_metric_scale('y', lat=lat)
            if ymscale is None: return
            vunits = getattr(self.data[2], 'units',  None)
            if vunits is None: return
            vmscale = tometric(vunits, munits='m/s')
            if vmscale is None: return

            # Scale
            uvscaler = (vmscale*xmscale*self.x[:].ptp())/(umscale*ymscale*self.y[:].ptp())

        # Parse
        self.uvscaler = uvscaler
        return self._uvscaler_(uvscaler, raw=raw)

    @staticmethod
    def _uvscaler_(uvscaler, raw=False):
        if callable(uvscaler): return uvscaler
        if isinstance(uvscaler, (tuple,list)):
            _ = 1., uvscaler[1]/uvscaler[0]
            if raw: return uvscaler
            uvscaler = _
        if N.isscalar(uvscaler):
            if raw: return uvscaler
            return lambda u, v: (u, uvscaler*v)
        raise TypeError("uvscaler must be a callable, a scalar or a"+
            " 2-element tuple or list")

    def set_uvscaler(self, uvscaler):
        """Set :attr:`uvscaler`"""
        self._uvscaler_(uvscaler)
        self.set_obj('uvscaler', uvscaler)
    def del_uvscaler(self):
        """Del :attr:`uvscaler`"""
    uvscaler = property(get_uvscaler, set_uvscaler,
        del_uvscaler, doc="""Function to rescale U anv V data along X and Y axes.
            ``None`` is returned if no scaling is possible.

            :Example:

                >>> uvscaler = myplot.uvscaler
                >>> if uvscaler is not None:
                    ... u2,v2 = myplot.uvscaler(u, v)

            :Return: A callable function or ``None`` if no scaling is possible
            """)



    def get_metric_scale(self, xy, lat=None):
        """Get units of X or Y plot axis as meters if possible

        Longitude and latitude coordinates are converted using
        :func:`~vacumm.phys.units.deg2m`.
        else :func:`~vacumm.misc.phys.units.tometric` is used
        using axis units.

        :Params:

            - **xy**: Plot axis type (``"x"``, ``"y"``...).
            - **lat**, optional: Latitude for degrees to meter conversion
              of longitude coordinates. It default to :attr:`uvlat`.

        :Return:

            - ``None`` if conversion of possible, else a float value.

        """
        # X or Y plot axis
        if xy not in ['x', 'y']: return

        # From order type
        axtype = getattr(self, xy+'type')
        if axtype=='x':
            if lat is not None:
                self.uvlat = lat
            else:
                lat = self.uvlat
            return deg2m(1., lat=lat)
        if axtype=='y':
            return deg2m(1.)
        units = getattr(self, xy+'units', None)
        if units is None:
            if axtype!='z': return
            units = 'm'
        #if units.startswith('deg'):
            #if getattr(self, xy+'type')=='x':
                #return deg2m(1., lat=lat)
            #if getattr(self, xy+'type')=='y':
                #return deg2m(1.)
        return tometric(units,  munits='m')

    def get_vmin(self, index=0, glob=False):
        """Get :attr:`vmin`"""
        gmin = self.get_obj('vmin')
        if gmin is not None: return gmin
#        if not self.has_data(): return
        if index=='all': # All components (m,u,v)
            if not self.has_valid_data():
                vmins = [None]*3
            else:
                vmins = [var.min() for var in self.data]
            if glob: # Of all plotters
                for b in self.get_brothers(notme=True):
                    for i, v in enumerate(b.get_min(index='all')):
                        if v is not None:
                            vmins[i] = min(vmins[i], v) if vmins[i] is not None else v
        else: # Selected component
            vmins = [self.data[index].min()] if self.has_valid_data() else []
            if glob: # Of all plotters
                vmins = vmins+[b.get_min(index=index) for b in self.get_brothers(notme=True)]
            vmins = [v for v in vmins if v is not None]
            vmins = min(vmins) if len(vmins) else None
        return vmins
    def set_vmin(self, value):
        """Set :attr:`vmin`"""
        self.set_obj('vmin', value)
    def del_vmin(self):
        """Del :attr:`vmin`"""
        self.del_obj('vmin')
    vmin = property(get_vmin, set_vmin, del_vmin, doc="Data min to use for plot")

    def get_vmax(self, index=0, glob=False):
        """Get :attr:`vmax`"""
        gmax = self.get_obj('vmax')
        if gmax is not None: return gmax
#        if not self.has_data(): return
        if index=='all': # All components (m,u,v)
            if not self.has_data() or self.masked:
                maxs = [None]*3
            else:
                maxs = [var.max() for var in self.data]
            if glob: # Of all plotters
                for b in self.get_brothers(notme=True):
                    for i, v in enumerate(b.get_max(index='all')):
                        if v is not None:
                            vmaxs[i] = max(vmaxs[i], v) if vmaxs[i] is not None else v
        else: # Selected component
            vmaxs = [self.data[index].max()] if self.has_valid_data() else []
            if glob: # Of all plotters
                vmaxs = vmaxs+[b.get_max(index=index) for b in self.get_brothers(notme=True)]
            vmaxs = [v for v in vmaxs if v is not None]
            vmaxs = max(vmaxs) if len(vmaxs) else None
        return vmaxs
    def set_vmax(self, value):
        """Set :attr:`vmax`"""
        self.set_obj('vmax', value)
    def del_vmax(self):
        """Del :attr:`vmax`"""
        self.del_obj('vmax')
    vmax = property(get_vmax, set_vmax, del_vmax, doc="Data max to use for plot")




    # X/Y MASKED

    def get_xymasked(self):
        """Get :attr:`xymasked`"""
        gmasked = self.get_obj('xymasked')
        return gmasked is None and True or gmasked
    def set_xymasked(self, value):
        """Set :attr:`xymasked`"""
        self.set_obj('xymasked', value)
    def del_xymasked(self):
        """Del :attr:`xymasked`"""
        self.del_obj('xymasked')
    xymasked = property(get_xymasked, set_xymasked, del_xymasked, doc="Whether X and Y data are not considered if no data at these coordinates")
    def get_xmasked(self):
        """Get :attr:`xmasked`"""
        gmasked = self.get_obj('xmasked')
        if gmasked is None: gmasked = self.xymasked
        return gmasked is None and True or gmasked
    def set_xmasked(self, value):
        """Set :attr:`xmasked`"""
        self.set_obj('xmasked', value)
    def del_xmasked(self):
        """Del :attr:`xmasked`"""
        self.del_obj('xmasked')
    xmasked = property(get_xmasked, set_xmasked, del_xmasked, doc="Whether X data are not considered if no data at these coordinates")
    def get_ymasked(self):
        """Get :attr:`ymasked`"""
        gmasked = self.get_obj('ymasked')
        if gmasked is None: gmasked = self.xymasked
        return gmasked is None and True or gmasked
    def set_ymasked(self, value):
        """Set :attr:`ymasked`"""
        self.set_obj('ymasked', value)
    def del_ymasked(self):
        """Del :attr:`ymasked`"""
        self.del_obj('ymasked')
    ymasked = property(get_ymasked, set_ymasked, del_ymasked, doc="Whether X data are not considered if no data at these coordinates")

    # X MIN/MAX

    def get_xmin(self, glob=True, masked=None):
        """Get :attr:`xmin`"""
        gmin = self.get_axobj('xmin')
        if gmin is not None: return gmin
        if self.masked: masked = False
        if masked is None: masked = self.xmasked
#        if not self.has_data(): return
        if self.order[1]=='d':
            index = max(0, len(self.data)-2)
            xmin = self.get_vmin(index=index)
        else:
            if not self.has_data():
                xmin = None
            elif masked:
                if self.masked:
                    xmin = None
                else:
                    xmin = self.get_xdata(masked=masked).min()
            else:
                xmin = self.get_xdata(masked=False).min()
        if glob:
            xmins = [] if xmin is None else [xmin]
            xmins += [b.get_xmin(glob=False, masked=masked)
                for b in self.get_brothers(notme=True)]
            xmins = [v for v in xmins if v is not None]
            xmin = min(xmins) if len(xmins) else None
        return xmin
    def set_xmin(self, value):
        """Set :attr:`xmin`"""
        self.set_axobj('xmin', value)
    def del_xmin(self):
        """Del :attr:`xmin`"""
        self.del_axobj('xmin')
    xmin = property(get_xmin, set_xmin, del_xmin, doc="Min of X axis to use for plot")

    def get_xmax(self, glob=True, masked=None):
        """Get :attr:`xmax`"""
        gmax = self.get_axobj('xmax')
        if gmax is not None: return gmax
        if self.masked: masked = False
        if masked is None: masked = self.xmasked
#        if not self.has_data(): return
        if self.order[1]=='d':
            index = max(0, len(self.data)-2)
            xmax = self.get_vmax(index=index)
        else:
            if not self.has_data():
                xmax = None
            elif masked:
                if self.masked:
                    xmax = None
                else:
                    xmax = self.get_xdata(masked=masked).max()
            else:
                xmax = self.get_xdata(masked=False).max()
        if glob:
            xmaxs = [] if xmax is None else [xmax]
            xmaxs += [b.get_xmax(glob=False, masked=masked)
                for b in self.get_brothers(notme=True)]
            xmaxs = [v for v in xmaxs if v is not None]
            xmax = max(xmaxs) if len(xmaxs) else None
        return xmax
    def set_xmax(self, value):
        """Set :attr:`xmax`"""
        self.set_axobj('xmax', value)
    def del_xmax(self):
        """Del :attr:`xmax`"""
        self.del_axobj('xmax')
    xmax = property(get_xmax, set_xmax, del_xmax, doc="Max of X axis to use for plot")

    # Y MIN/XMAX

    def get_ymin(self, glob=True, masked=None):
        """Get :attr:`ymin`"""
        gmin = self.get_axobj('ymin')
        if gmin is not None: return gmin
        if self.masked: masked = False
        if masked is None: masked = self.ymasked
#        if not self.has_data(): return
        if self.order[0]=='d':
            index = min(0, len(self.data)-2)
            ymin = self.get_vmin(index=index)
        else:
            if not self.has_data():
                ymin = None
            elif masked:
                if self.masked:
                    ymin = None
                else:
                    ymin = self.get_ydata(masked=masked).min()
            else:
                ymin = self.get_ydata(masked=False).min()
        if glob:
            ymins = [] if ymin is None else [ymin]
            ymins += [b.get_ymin(glob=False,  masked=masked)
                for b in self.get_brothers(notme=True,)]
            ymins = [v for v in ymins if v is not None]
            ymin = min(ymins) if len(ymins) else None
        return ymin
    def set_ymin(self, value):
        """Set :attr:`ymin`"""
        self.set_axobj('ymin', value)
    def del_ymin(self):
        """Del :attr:`ymin`"""
        self.del_axobj('ymin')
    ymin = property(get_ymin, set_ymin, del_ymin, doc="Min of Y axis to use for plot")

    def get_ymax(self, glob=True, masked=None):
        """Get :attr:`ymax`"""
        gmax = self.get_axobj('ymax')
        if gmax is not None: return gmax
        if masked is None: masked = self.ymasked
        if self.masked: masked = False
#        if not self.has_data(): return
        if self.order[0]=='d':
            index = max(0, len(self.data)-2)
            ymax = self.get_vmax(index=index)
        else:
            if not self.has_data():
                ymax = None
            elif masked:
                if self.masked:
                    ymax = None
                else:
                    ymax = self.get_ydata(masked=masked).max()
            else:
                ymax = self.get_ydata(masked=False).max()
        if glob:
            ymaxs = [] if ymax is None else [ymax]
            ymaxs += [b.get_ymax(glob=False)
                for b in self.get_brothers(notme=True)]
            ymaxs = [v for v in ymaxs if v is not None]
            ymax = max(ymaxs) if len(ymaxs) else None
            ymaxs = [ymax]+[b.get_ymax(glob=False, masked=masked)
                for b in self.get_brothers(notme=True)]
            ymaxs = [v for v in ymaxs if v is not None]
            ymax = max(ymaxs) if ymaxs else None
        return ymax
    def set_ymax(self, value):
        """Set :attr:`ymax`"""
        self.set_axobj('ymax', value)
    def del_ymax(self):
        """Del :attr:`ymax`"""
        self.del_axobj('ymax')
    ymax = property(get_ymax, set_ymax, del_ymax, doc="Max of Y axis to use for plot")



    # X/Y TICKS

    def get_xticks(self, raw):
        """Get X ticks"""
        return self.axes.get_xticks()
    def set_xticks(self, ticks):
        """Set X ticks"""
        if ticks is False:
            ticks = []
        elif ticks=='auto':
            if not self.has_data():return
            if self.order[1]=='t':
                ticks = AutoDateLocator2()
            ticks = self.auto_scale(self.xmin, self.xmax)
        if isinstance(ticks, Locator):
            self.axes.xaxis.set_major_locator(ticks)
        elif ticks is not None:
            self.axes.set_xticks(ticks)
    def del_xticks(self):
        """Del X ticks"""
        self.xticks = False
    xticks = property(get_xticks, set_xticks, del_xticks, doc="X ticks to use plot")

    def get_yticks(self, glob=False):
        """Get Y ticks"""
        return self.axes.get_yticks()
    def set_yticks(self, ticks):
        """Set Y ticks"""
        if ticks is False:
            ticks = []
        elif ticks=='auto':
            if not self.has_data():return
            if self.order[1]=='t':
                ticks = AutoDateLocator2()
            ticks = self.auto_scale(self.xmin, self.xmax)
        if isinstance(ticks, Locator):
            self.axes.xaxis.set_major_locator(ticks)
        elif ticks is not None:
            self.axes.set_xticks(ticks)
    def del_yticks(self):
        """Del Y ticks"""
        self.yticks = False
    yticks = property(get_yticks, set_yticks,
        del_yticks, doc="X ticks to use for plot")


    # Generic attribute management
    def _get_xyattr_(self, xy, att, idata=None):
        """Get an attribute of an X/Y axis or current data"""

#        # Nothing
#        if not self.has_data(): return

        # From a variable
        if xy=='d' or getattr(self, xy+'type')=='d':
            if not self.has_data(): return
            if idata is None:
                idata = range(len(self.data))
            elif not isinstance(idata, (list, tuple)):
                idata = [idata]
            for i in idata:
                if i<0:
                    i = len(self.data)+i
                if len(self.data)>i and hasattr(self.data[i], att):
                    return getattr(self.data[i], att)
            return

        # From an axis
        return getattr(getattr(self, xy), att, None)

    def _set_xyattr_(self, xy, att, value, idata=0):
        """Set an attribute to an X/Y axis or current data"""

#        # Nothing
#        if not self.has_data(): return

        # Del
        if value is False:
            self._del_xyattr_(xy, att)

        # Variable or axis?
        if xy=='d' or getattr(self, xy+'type')=='d':
            if not self.has_data():
                return
            target = self.data[idata]
        else:
            target = getattr(self, xy)
#            setattr(self.data[idata], att, value)

        # Set
        setattr(target, att, value)

    def _del_xyattr_(self, xy, att, idata=0):
        """Del an attribute from an X/Y axis or current data"""

        # Nothing
        if not self.has_data(): return

        # To a variable
        if xy=='d' or getattr(self, xy+'type')=='d':
            if not self.has_data():
                return
            target = self.data[idata]
        else:
            target = getattr(self, xy)

        # Del
        if hasattr(target, att):
            delattr(target, att)


    # ID / SHORT NAME

    def get_id(self, idata=None):
        """Get :attr:`id`"""
        return self._get_xyattr_('d', 'id', idata=idata)
    def set_id(self, id=None):
        """Set :attr:`id`"""
        self._set_xyattr_('d', 'id', id)
    def del_id(self):
        """Del :attr:`id`"""
        self._del_xyattr_('d', 'id')
    id = property(get_id, set_id,
        del_id, 'Current data id')

    def get_xid(self):
        """Get :attr:`xid`"""
        return self._get_xyattr_('x', 'id')
    def set_xid(self, id=None):
        """Set :attr:`xid`"""
        self._set_xyattr_('x', 'id', id)
    def del_xid(self):
        """Del :attr:`xid`"""
        self._del_xyattr_('x', 'id')
    xid = property(get_xid, set_xid,
        del_xid, 'Current id of X axis')

    def get_yid(self):
        """Get :attr:`yid`"""
        return self._get_xyattr_('y', 'id')
    def set_yid(self, id=None):
        """Set :attr:`yid`"""
        self._set_xyattr_('y', 'id', id)
    def del_yid(self):
        """Del :attr:`yid`"""
        self._del_xyattr_('y', 'id')
    yid = property(get_yid, set_yid,
        del_yid, 'Current id of Y axis')



    # LONG_NAME

    def get_long_name(self, idata=None):
        """Get :attr:`long_name`"""
        long_name = self._get_xyattr_('d', 'long_name', idata=idata)
        id = self.get_id(idata=idata)
        if long_name is None and id:
            long_name = id.title().replace('_', ' ')
        return long_name
    def set_long_name(self, long_name=None):
        """Set :attr:`long_name`"""
        self._set_xyattr_('d', 'long_name', long_name)
    def del_long_name(self):
        """Del :attr:`long_name`"""
        self._del_xyattr_('d', 'long_name')
    long_name = property(get_long_name, set_long_name,
        del_long_name, 'Current long name')

    def get_xlong_name(self):
        """Get :attr:`xlong_name`"""
        long_name = self._get_xyattr_('x', 'long_name')
        if long_name is None:
            long_name = self.get_xid().title().replace('_', ' ')
        return long_name
    def set_xlong_name(self, long_name=None):
        """Set :attr:`xlong_name`"""
        self._set_xyattr_('x', 'long_name', long_name)
    def del_xlong_name(self):
        """Del :attr:`xlong_name`"""
        self._del_xyattr_('x', 'long_name')
    xlong_name = property(get_xlong_name, set_xlong_name,
        del_xlong_name, 'Current long_name of X axis')

    def get_ylong_name(self):
        """Get :attr:`ylong_name`"""
        long_name = self._get_xyattr_('y', 'long_name')
        if long_name is None:
            long_name = self.get_yid().title().replace('_', ' ')
        return long_name
    def set_ylong_name(self, long_name=None):
        """Set :attr:`ylong_name`"""
        self._set_xyattr_('y', 'long_name', long_name)
    def del_ylong_name(self):
        """Del :attr:`ylong_name`"""
        self._del_xyattr_('y', 'long_name')
    ylong_name = property(get_ylong_name, set_ylong_name,
        del_ylong_name, 'Current long_name of Y axis')


    # UNITS

    def get_units(self, idata=None):
        """Get :attr:`units`"""
        units = self._get_xyattr_('d', 'units', idata)
        if not isinstance(units, basestring):
            return
        if units[0]=='$' and units[-1]=='$':
            units = units[1:-1]
            if not self.isset('latex_units'):
                self.latex_units = True
        return units
    def set_units(self, units=None):
        """Set :attr:`units`"""
        self._set_xyattr_('d', 'units', units)
    def del_units(self):
        """Del :attr:`units`"""
        self._del_xyattr_('d', 'units')
    units = property(get_units, set_units, del_units, 'Current units')

    def get_xunits(self):
        """Get :attr:`xunits`"""
        return self._get_xyattr_('x', 'units')
    def set_xunits(self, units=None):
        """Set :attr:`xunits`"""
        self._set_xyattr_('x', 'units', units)
    def del_xunits(self):
        """Del :attr:`xunits`"""
        self._del_xyattr_('x', 'units')
    xunits = property(get_xunits, set_xunits,
        del_xunits, 'Current units of X axis')

    def get_yunits(self):
        """Get :attr:`yunits`"""
        return self._get_xyattr_('y', 'units')
    def set_yunits(self, units=None):
        """Set :attr:`yunits`"""
        self._set_xyattr_('y', 'units', units)
    def del_yunits(self):
        """Del :attr:`yunits`"""
        self._del_xyattr_('y', 'units')
    yunits = property(get_yunits, set_yunits,
        del_yunits, 'Current units of Y axis')

    # INTERPRET UNITS WITH LATEX?

    def get_latex_units(self):
        """Get :attr:`latex_units`"""
        return self.get_obj('latex_units') or False
    def set_latex_units(self, value):
        """Set :attr:`units`"""
        self.set_obj('latex_units', value)
    def del_latex_units(self):
        """Del :attr:`units`"""
        self.del_obj('latex_units')
    latex_units = property(get_latex_units, set_latex_units, del_latex_units, 'Interpret units with latex')

    re_latex_match = re.compile(r'\$.+\$').match
    def is_latex(self, text):
        """Is this text formatted as latex formula ("$...$")?"""
        return self.re_latex_match(text)


    # TITLE OF THE PLOT

    def _strfill_(self, strpat, name='pattern'):
        """Try to fill strpat with self.sndict() or return strpat with not filling"""
        try:
            return strpat % self.sndict()
        except:
#            warn("Error when filling %(name)s (skipping): %(strpat)s"%locals())
            return strpat


    def get_title(self):
        """Get :attr:`title`"""
        title = self.get_axobj('title')
        if title is False: return False
        if title == '_auto_': title = True
        if isinstance(title, basestring):
            return self._strfill_(title, 'title pattern')
        return self.long_name
    def set_title(self, title=None):
        """Set attr:`title`

        .. note:: If set to ``False``, the title is not plotted.
        """
        self.set_axobj('title', title)
    def del_title(self):
        """Del attr:`title`"""
        self.del_axobj('title')
    title = property(get_title, set_title,
        del_title, 'Preformed title to use for the plot. '
        'It may be formed as a template using other attributes '
        'like :attr:`long_name`, :attr:`units`,  :attr:`xmin`,  etc.')


    # LABELS

    def get_label(self):
        """Get :attr:`label`"""
        label = self.get_axobj('label')
        if label is False: return
        if label is not None:
            return self._strfill_(label, 'label pattern')
        return self.long_name
    def set_label(self, label=None):
        """Set :attr:`label`

        .. note:: If set to ``False``, the label is not plotted.
        """
        self.set_axobj('label', label)
    def del_label(self):
        """Del :attr:`label`"""
        self.del_axobj('label')
    label = property(get_label, set_label,
        del_label, 'Preformed label to use for the plot. '
        'It may be formed as a template using other attributes '
        'like :attr:`long_name`, :attr:`units`,  :attr:`xmin`,  etc.')



    def get_xlabel(self):
        """Get :attr::`xlabel`"""
        label = self.get_axobj('xlabel')
#        label = getattr(self, '_xlabel', None)
        if label is False: return ''
        if label is None or label is True:
            if self.order[1]=='d':
                label = self.get_fmt_lnu(long_name=False)
            elif label is True or self.order[1]=='-':
                label = self.get_fmt_lnu(prefix='x')
            else:
                label = ''
        return self._strfill_(label, 'xlabel pattern')
    def set_xlabel(self, label):
        """set :attr::`xlabel`"""
        self.set_axobj('xlabel', label)
#        self._xlabel = label
    def del_xlabel(self):
        """Del :attr::`xlabel`"""
        self.del_axobj('xlabel')
#        if hasattr(self, '_xlabel'): del self._xlabel
    xlabel = property(get_xlabel, set_xlabel,
        del_xlabel, 'Label of X axis. '
        'It may be formed as a template using other attributes '
        'like :attr:`long_name`, :attr:`units`,  :attr:`xmin`,  etc.')

    def get_ylabel(self):
        """Get :attr::`ylabel`"""
        label = self.get_axobj('ylabel')
#        label = getattr(self, '_ylabel', None)
        if label is False: return ''
        if label is None or label is True:
            if self.order[0]=='d':
                label = self.get_fmt_lnu(long_name=False)
            elif label is True or self.order[0]=='-':
                label = self.get_fmt_lnu(prefix='y')
            else:
                label = ''
        return self._strfill_(label, 'ylabel pattern')
    def set_ylabel(self, label):
        """Set :attr::`ylabel`"""
        self.set_axobj('ylabel', label)
#        self._ylabel = label
    def del_ylabel(self):
        """Del :attr::`ylabel`"""
        self.del_axobj('ylabel')
#        if hasattr(self, '_ylabel'): del self._ylabel
    ylabel = property(get_ylabel, set_ylabel,
        del_ylabel, 'Label of Y axis. '
        'It may be formed as a template using other attributes '
        'like :attr:`long_name`, :attr:`units`,  :attr:`xmin`,  etc.')

    def get_fmt_lnu(self, prefix='', fmtln='%(long_name)s', fmtu='%(units)s',
        fmtlnu='%(long_name)s [%(units)s]', long_name=True, units=True):
        """Format long_name and units as string according to their availability

        :Params:

            - **prefix**, optional: Prefix of the attributes

        :Example:

            >>> myplot.get_fmt_lnu()
            '%(long_name)s [%(units)s]'
            >>> myplot.get_fmt_lnu(prefix='x', fmtlnu='%(long_name)s [%(units)s]')
            '%(xlong_name)s [%(xunits)s]'
            >>> myplot.get_fmt_lnu(long_name=False)
            '%(xunits)s'
        """
        if getattr(self, prefix+'long_name') is None: long_name = False
        if getattr(self, prefix+'units') is None: units = False
        if long_name and units:
            fmt = fmtlnu
        elif long_name:
            fmt = fmtln
        elif units:
            fmt = fmtu
        else:
            fmt=''
        long_name = '%%(%slong_name)s'%prefix
        units = '%%(%sunits)s'%prefix
        return fmt%locals()


    # DICTIONARY OF ATTRIBUTES

    def dict(self, *keys, **items):
        """Get attributes as a dictionary

        .. note::

            :attr:`units` is treated in a special way.
            If :attr:`latex_units` is ``True``, it is formatted
            as ``$<units>$``.

        """
        if len(keys)==0:
            keys = list(self._primary_attributes)
        dd = {}
        for key in keys:
            dd[key] = getattr(self, key)
        items = dict([(key, val) for (key, val) in items.items() if val is not None])
        dd.update(items)
        if self.latex_units and 'units' in dd and not self.is_latex(dd['units']):
            dd['units'] = '$%(units)s$'%dd
        return dd

    def sndict(self, *keys, **items):
        """Get attributes as a dictionary of strings or numbers"""
        dd = self.dict(*keys, **items)
        for key, val in dd.items():
            if val is None or val is False or val is True:
                dd[key] = ''
        return dd


    def sdict(self, *keys, **items):
        """Get attributes as a dictionary of strings"""
        dd = self.sndict(*keys, **items)
        for key, val in dd.items():
            val = str(val)
        return dd


    # ADVANCED GRAPHICAL METHODS


class Plot1D(Plot):

    _order = ['zd', 'd-', '-d']
    rank = 1

    def _check_order_(self, vertical=None, **kwargs):
        """Check order of data

        :Params:

            - **vertical**, optional: Plot vertically.

        :Sea also: :meth:`Plot._check_order_`
        """

        # Force vertical plot
        if vertical is True:
            self._order = [o for o in self._order if o.startswith('d') ]
        elif vertical is False:
            self._order = [o for o in self._order if not o.startswith('d')]

        # Old stuff
        return Plot._check_order_(self, **kwargs)
#    check_order.__doc__ = Plot.check_order.__doc__

    def _set_axes_(self, axis=None, axisatts=None, **kwargs):
        """Get/set used axes

        :Params:

            - **axis**, optional: Change plot axis.
        """

        # Get axis and adjust order
        if axis is not None:
            if not isaxis(axis):
                axis = cdms2.createAxis(axis)
                if self.has_data():
                    cp_atts(self.data[0].getAxis(0), axis, overwrite=False)
            elif self.has_data(): # adjust order
                orders = [var.getOrder() for var in self.data]
                for ivar, var in enumerate(self.data):
                    c = getattr(axis, 'axis', '-').lower()
                    if c!='-' or orders[ivar][0]=='-':
                        set_order(var, c)
        elif self.has_data():
            axis = self.data[0].getAxis(0)
        else:
            raise PlotError('No axis data available for this plot')

        # Change attributes
        if axisatts is not None:
            set_atts(axis, axisatts)

        # Set X axis by default
        self.x = axis
        self.y = None


#    _set_axes_.__doc__ = Plot._set_axes_.__doc__

        # Restore order


    def get_axis_data(self, **kwargs):
        """Return data of the axis"""
        if not self.has_data(): return
        if self.xtype=='d': return self.get_ydata(**kwargs)
        return self.get_xdata(**kwargs)


class ScalarMappable:
    """Abstract class for adding scalar mappable utilities

    :Attribute params:

        - **levels**, optional: Levels to use for contours or colorbar ticks.
          "They can be specified as a single value, a list or array, or "
          "as a tuple used as argument to :func:`numpy.arange`.
          It sets the attribute :attr:`~vacumm.misc.core_plot.ScalarMappable.levels`.
        - **nmax_levels**, optional: Maximal number of levels when
          :attr:`~vacumm.misc.core_plot.ScalarMappable.levels` are computed automatically.
          It sets the attribute :attr:`~vacumm.misc.core_plot.ScalarMappable.nmax_levels`.
        - **nmax**, optional: Same as **nmax_levels**.
          It sets the attribute :attr:`~vacumm.misc.core_plot.ScalarMappable.keepminmax`.
        - **cmap**, optional: Colormap name (see :func:`vacumm.misc.color.get_cmap`).
          If not specified, it is taken from
          config section ``[vacumm.misc.plot]`` and config option ``cmap``, as a string
          that defaults to ``magic``.
        - **levels_mode**, optional: Mode of computing levels if needed.
          It can be specified at initialisation with
          attribute :attr:`~vacumm.misc.core_plot.ScalarMappable.levels_mode`.
          If not specified, it it taken from the config section
          ``[vacumm.misc.plot]`` and config option ``levels_mode``.

            - ``"symetric"``: Min and max are set opposite.
            - ``"positive"``: Min is set to 0.
            - ``"negative"``: Max is set to 0.
            - ``"auto"``: If abs(min) and abs(max) are close,
              ``"symetric"`` is assumed. If min and max are > 0,
              ``"positive"`` is assumed, and the reverse for
              ``"negative"``.
        - **keepminmax**, optional:
          It can be specified at initialisation with
          attribute :attr:`~vacumm.misc.core_plot.ScalarMappable.keepminmax`.
          If not specified, it it taken from the config section
          ``[vacumm.misc.plot]`` and config option ``keepminmax``.
          If False or 0, adjust
          :attr:`~vacumm.misc.core_plot.ScalarMappable.vmin` and
          :attr:`~vacumm.misc.core_plot.ScalarMappable.vmax`
          to first and last values of :attr:`~vacumm.misc.core_plot.ScalarMappable.levels`;
          if 1, do not change :attr:`~vacumm.misc.core_plot.ScalarMappable.vmin`,
          :attr:`~vacumm.misc.core_plot.ScalarMappable.vmax` and :attr:`levels`
          if 2, adjust first and last values of :attr:`levels` or
          :attr:`~vacumm.misc.core_plot.ScalarMappable.vmin`,
          :attr:`~vacumm.misc.core_plot.ScalarMappable.vmax`.
          It sets the attribute :attr:`~vacumm.misc.core_plot.ScalarMappable.cmap`.
        - **cblabel**, optional: Preformed label of the colorbar.
          It may be formed as a template using other attributes
          like :attr:`long_name`, :attr:`~vacumm.misc.core_plot.Plot.units`,
          :attr:`xmin`,  etc. Example: ``"%(long_name)s [%(units)s]"``.


    """
    _primary_attributes = Plot._primary_attributes + ['nmax', 'nmax_levels',
        'levels_mode', 'levels', 'cmap', 'keepminmax']
    _secondary_attributes = Plot._secondary_attributes + ['cblabel']

    def get_nmax_levels(self):
        """Get :attr:`nmax_levels`"""
        nmax = self.get_obj('nmax_levels')
        if nmax is not None: return nmax
        return int(get_config_value('vacumm.misc.plot', 'nmax_levels'))
    def set_nmax_levels(self, value):
        """Set :attr:`nmax_levels`"""
        self.set_obj('nmax_levels', value)
    def del_nmax_levels(self):
        """Del :attr:`nmax_levels`"""
        self.del_obj('nmax_levels')
    nmax_levels = nmax = property(get_nmax_levels, set_nmax_levels,
        del_nmax_levels, doc="Max number of :attr:`levels` for contours and colorbar ticks.")


    def get_levels(self, mode=None, keepminmax=None, nocache=False,
            autoscaling='normal', **kwargs):
        """Get :attr:`levels` for contours and colorbar ticks

        :Params:

            - **mode**, optional: Mode of computing levels if needed.
              It can be specified at initialisation with
              attribute :attr:`levels_mode`.
              If not specified, it it taken from the config section
              ``[vacumm.misc.plot]`` and config option ``levels_mode``.

                - ``"normal"``: Min and max are not preprocessed.
                - ``"symetric"``: Min and max are set opposite.
                - ``"positive"``: Min is set to 0.
                - ``"negative"``: Max is set to 0.
                - ``"auto"``: If abs(min) and abs(max) are close,
                  ``"symetric"``is assumed. If min and max are > 0,
                  ``"positive"`` is assumed, and the reverse for
                  ``"negative"``.

            - **keepminmax**, optional:
              It can be specified at initialisation with
              attribute :attr:`keepminmax`.
              If not specified, it it taken from the config section
              ``[vacumm.misc.plot]`` and config option ``keepminmax``.
              If False or 0, adjust
              :attr:`vmin` and :attr:`vmax` to first and last values of :attr:`levels`;
              if 1, do not change :attr:`vmin`, :attr:`vmax` and :attr:`levels`
              if 2, adjust first and last values of :attr:`levels` or
              :attr:`vmin`, :attr:`vmax`.
            - **nocache**, optional: Once levels are computed, they are stored
              in cache. If ``nocache is True``, first check cache before
              trying to compute levels.
            - **autoscaling**, optional: Autoscaling mode.

                - ``"normal"``: Use :func:`~vacumm.misc.misc.auto_scale`.
                - ``"degrees"``: Use :func:`~vacumm.misc.misc.geo_scale`.
                - `A callable: Use it to auto scale. It should accept
                  the follwing keywords: vmin, vmax, nmax, keepminmax.

        """
        # Cache
        levels = self.get_obj('levels')
        if levels is not None and not isinstance(levels, str):
            if not hasattr(levels, '__len__'):
                levels = N.asarray([levels])
            elif isinstance(levels, tuple):
                levels = N.arange(*levels[:3]).astype('d')
            self._levels = levels
            return levels
        if hasattr(self, '_levels') and not nocache: return self._levels

        # Inits
        if isinstance(levels, str):
            mode = levels
        if mode is not None:
            self.levels_mode = mode
        else:
            mode = self.levels_mode
        if not self.has_data(): return
        if keepminmax is not None:
            self.keepminmax = keepminmax
        keepminmax = self.keepminmax

        # Min and max
        for key in 'positive', 'negative', 'symetric', 'anomaly':
            if kwargs.has_key(key) and kwargs[key]:
                mode = key
        vmin = self.vmin if self.isset('vmin') else None
        vmax = self.vmax if self.isset('vmax') else None
        if mode == 'auto':
            mode = self.minmax2levelsmode(self.vmin, self.vmax)
        self.levels_mode = mode
        if mode=='positive':
            vmin = 0.
        elif mode=='negative':
            vmax = 0.
        elif mode=='symetric':
            vmax = max(N.abs(self.vmin), N.abs(self.vmax))
            vmin = -vmax

        # Compute base levels
        assert autoscaling in ['normal', 'degrees'] or callable(autoscaling), 'Wrong autoscaling parameter'
        if autoscaling=='normal':
            autoscaling = auto_scale
        elif autoscaling=='degrees':
            autoscaling = geo_scale
        levels = autoscaling((self.vmin, self.vmax), vmin=vmin, vmax=vmax,
            nmax=self.nmax_levels, keepminmax=keepminmax==2)

        # Change min and max
        if not keepminmax:
            del self.vmin, self.vmax


        # Cache
        self._levels = levels
        return levels

    def minmax2levelsmode(self, vmin=None, vmax=None):
        """Get auto levels mode from min and max value

        :Return: normal, symetric, positive or negative"""
        if vmin is None:
            vmin = self.vmin
        if vmax is None:
            vmax = self.vmax
        if self.masked:
            mode = 'normal'
        elif (vmin>=0 and vmax > 0) and (vmin<=vmax/3.):
            mode = 'positive'
        elif (vmin<0 and vmax<=0) and (vmax>=vmin/3.):
            mode = 'negative'
        elif (vmin+vmax)< 0.05 * (vmax-vmin):
            mode = 'symetric'
        else:
            mode = 'normal'
        return mode

    def set_levels(self, value):
        """Set :attr:`levels`"""
        self.set_obj('levels', value)
    def del_levels(self):
        """Del :attr:`levels`"""
        self.del_obj('levels')
    levels = property(get_levels, set_levels,
        del_levels, doc="Levels to use for contours or colorbar ticks. "
            "They can be specified as a single value, a list or array, or "
            "as a tuple used as argument to :func:`numpy.arange`.")

    def get_keepminmax(self):
        """Get :attr:`keepminmax`"""
        keepminmax = self.get_obj('keepminmax')
        if keepminmax is None:
            keepminmax = get_config_value('vacumm.misc.plot', 'keepminmax')
        try:
            keepminmax = int(keepminmax)
        except:
            raise PlotError('Error with keepminmax: %s'%keepminmax)
        return keepminmax
    def set_keepminmax(self, value):
        """Set :attr:`keepminmax`"""
        self.set_obj('keepminmax', value)
    def del_keepminmax(self):
        """Del :attr:`keepminmax`"""
        self.del_obj('keepminmax')
    keepminmax = property(get_keepminmax, set_keepminmax,
        del_keepminmax, doc="Do not adjust :attr:`vmin` and :attr:`vmax` when setting :attr:`levels`. If 0, :attr:`vmin` and :attr:`vmax` are set to first and last :attr:`levels`; if 1, they are not adjested to :attr:`levels` ; if 2, first and last :attr:`levels` are adjusted to :attr:`vmin` and :attr:`vmax`.")

    def get_levels_mode(self):
        """Get :attr:`levels_mode`"""
        levels_mode = self.get_obj('levels_mode')
        if levels_mode is None:
            levels_mode = get_config_value('vacumm.misc.plot', 'levels_mode')
        levels_mode = self._check_levels_mode_(levels_mode)
        return levels_mode
    def set_levels_mode(self, mode):
        """Set :attr:`levels_mode`"""
        mode = self._check_levels_mode_(mode)
        self.set_obj('levels_mode', mode)
    def del_levels_mode(self):
        """Del :attr:`levels_mode`"""
        self.del_obj('levels_mode')
    levels_mode = property(get_levels_mode, set_levels_mode,
        del_levels_mode, doc="The way :attr:`levels` are estimated from "
            ":attr:`vmin` and :attr:`vmax`: 'positive'/'negative' means levels "
            "starting/ending from 0, 'anomaly' or 'symetric' means symetric "
            "levels, 'auto' or 'smart' means that mode is estimated "
            "from min and max, and 'normal' means no special treatment.")
    def _check_levels_mode_(self, mode):
        """Check values and aliases, and fallback to config value"""
        mode = str(mode).lower()
        if mode not in ['auto', 'smart', 'normal', 'basic', 'positive',
                'negative', 'symetric', 'anomaly']:
            oldmode = mode
            mode = get_config_value('vacumm.misc.plot', 'levels_mode', user=False)
            warn('Bad value for config value [vacumm.misc.plot] levels_mode: %s. Switched to default: %s'%(oldmode, mode))
        if mode=='anomaly':
            mode = 'symetric'
        elif mode=='smart':
            mode = 'auto'
        elif mode=='basic':
            mode = 'normal'
        return mode

    def get_vmin(self, index=0, glob=False):
        """Get :attr:`vmin`"""
        gmin = self.get_obj('vmin')
        if gmin is not None: return gmin
        if not self.has_data(): return
        if index==0 and hasattr(self, '_levels'):
            mins = self._levels[0]
            if glob: # Of all plotters
                mins = max([mins]+[b.min for b in self.get_brothers(notme=True)])
            return mins
        return Plot.get_vmin(self, index=index, glob=glob)
    vmin = property(get_vmin, Plot.set_vmin,
        Plot.del_vmin, doc="Data min to use for plot")
    def get_vmax(self, index=0, glob=False):
        """Get :attr:`vmax`"""
        gmax = self.get_obj('vmax')
        if gmax is not None: return gmax
        if not self.has_data(): return
        if index==0 and hasattr(self, '_levels'):
            maxs = self._levels[-1]
            if glob: # Of all plotters
                maxs = max([maxs]+[b.max for b in self.get_brothers(notme=True)])
            return maxs
        return Plot.get_vmax(self, index=index, glob=glob)
    vmax = property(get_vmax, Plot.set_vmax,
        Plot.del_vmax, doc="Data max to use for plot")


    @staticmethod
    def _get_config_cmap_(key='cmap'):
            cmap = get_config_value('vacumm.misc.plot', key)
            if cmap is None:
                cmap = get_config_value('vacumm.misc.plot', key, user=False)
            if cmap is None or cmap.lower() in ['none', 'mpl', 'default',
                'normal', 'true', 'false']:
                cmap = None # default from matplotlib
            return cmap

    def get_cmap(self, cmap=None, nocache=False, tint=None, lum=None, sat=None,
            pastel=False, **kwargs):
        """Get :attr:`cmap`

        :Params:

            - **cmap**, optional: Colormap name (see :func:`vacumm.misc.color.get_cmap`).
              It defaults to ``"magic"``.
            - **nocache**, optional: Once cmap is computed, it is stored
              in cache. If ``nocache is True``, first check cache before
              trying to compute cmap.
            - **lum**, optional: Change luminosity, between -1 and 1.
            - **sat**, optional: Change saturation.

        """

        # From cache
        cmap = self.get_obj('cmap')
        if (cmap is None and not nocache and hasattr(self, '_cmap') and
                tint is None and lum is None):
            cmap = self._cmap

        # From config
        if cmap is None or cmap is True :
            cmap = self._get_config_cmap_()

        # Aliases
        if cmap == 'mpl':
            cmap = False
        elif cmap=='mg':
            cmap = 'magic'
        elif cmap=='rb':
            cmap = 'rainbow'
        elif cmap == 'normal' or cmap is True or cmap is False:
            cmap = None

        # Adaptative choices
        cmap_levels = ['positive', 'negative', 'symetric', 'anomaly']
        if cmap in ['auto', 'magic', 'rainbow'] + cmap_levels:

            # Levels mode
            if cmap in cmap_levels:
                levels_mode = self._check_levels_mode_(cmap)
                cmap = 'auto'
            else:
                levels_mode = self.levels_mode
                if levels_mode == 'auto':
                    levels_mode = self.minmax2levelsmode(self.vmin, self.vmax)

            if cmap=='auto':

                if levels_mode == 'normal':
                    cmap = None
                elif levels_mode == 'symetric':
                    cmap = self._get_config_cmap_('cmap_symetric')
                    if cmap is None:
                        cmap = CMAP_SYMETRIC
                elif levels_mode == 'positive':
                    cmap = self._get_config_cmap_('cmap_positive')
                    if cmap is None:
                        cmap = CMAP_POSITIVE
                else:
                    cmap = self._get_config_cmap_('cmap_negative')
                    if cmap is None:
                        cmap = CMAP_NEGATIVE
                cmap = get_cmap(cmap, errmode='warn')

            elif cmap=='magic' or cmap=='rainbow':

                if cmap=='magic':
                    kwargs.setdefault('mode', levels_mode)
                kwargs.setdefault('stretcher', 'reduced_green')
                if getattr(self, 'fill_method', None)=='contourf' or getattr(self, 'fill', '')=='contourf':
                    dict_check_defaults(kwargs, stretch=0, lstretch=0, rstretch=0)
                cmap = getattr(color, 'cmap_'+cmap)(self.levels, **kwargs)


        elif not isinstance(cmap, Colormap):

            cmap = get_cmap(cmap, errmode='warn', **kwargs)

        # Luminosity
        if lum is None:
            lum = .5
        if tint is not None:
            lum = tint*.5 + .5
        if lum != 0.5:
            cmap = change_luminosity(cmap, lum)


        # Saturation
        if sat is not None:
            cmap = change_saturation(cmap, sat)

        # Pastelisation
        if pastel:
            cmap = pastelise(cmap)

        self._cmap = cmap
        return cmap
    def set_cmap(self, cmap):
        """Set :attr:`cmap`"""
        self.set_obj('cmap', cmap)
    def del_cmap(self, cmap):
        """Del :attr:`cmap`"""
        self.del_obj('cmap')
    cmap = property(get_cmap, set_cmap,
        del_cmap, doc="Colomap to use for filled plots and colorbar.")

    def get_cblabel(self):
        """Get :attr:`cblabel`"""
        label = getattr(self, '_cblabel', None)
        if label is None or label is True:
            label = self.get_fmt_lnu(long_name=False)
        if not isinstance(label, basestring): return ''
        return self._strfill_(label, 'cblabel pattern')
    def set_cblabel(self, label):
        """Set :attr:`cblabel`"""
        self._cblabel = label
    def del_cblabel(self):
        """Del :attr:`cblabel`"""
        if hasattr(self, '_cblabel'): del self._cblabel
    cblabel = property(get_cblabel, set_cblabel,
        del_cblabel, doc='Preformed label of the colorbar. '
        'It may be formed as a template using other attributes '
        'like :attr:`long_name`, :attr:`units`,  :attr:`xmin`,  etc.')

    def get_scalar_mappable(self):
        """Get the current scalar mappable or ``None``

        It is useful for instance for :meth:`colorbar`.

        .. note:

            The scalar mappable is search for in current instance only.
        """
        return self.get_obj('scalar_mappable')
    get_sm = get_scalar_mappable

    def get_colorbar(self):
        """Get the current colorbar object        """
        cb = self.get_obj('colorbar')
        if cb is not None:
            return cb

#    def get_colorbars(self):
#        """Get the colorbar of all brothers"""
#        return [b.get_colorbar()
#            for b in self.get_brothers(notme=False, mefirst=True)
#                if cb is not None]

    def colorbar(self, cax=None, fit=False, ticklabels_nmax=12, **kwargs):
        """Add a colorbar

        The colorbar is drawn only if :meth:`get_scalar_mappable` returns a valid
        scalar mappable.

        :Params:

            - **cax**, optional: Axes for the colorbar.
            - **label_<param>**, optional: <param> is passed to
              :meth:`~matplotlib.colorbar.Colorbar.set_label`.
            - **ticklabels_<param>**, optional: <param> is set as a property
              of tick labels.
            - Other keywords are passed to the :meth:`matplotlib.figure.Figure.colorbar`
              method.

        :See also: :meth:`get_colorbar` :meth:`get_scalar_mappable`
        """
        # Get scalar mappable
        sm = self.get_scalar_mappable()
        if sm is None: return

        # Check if already plotted
        cb = self.get_colorbar()
        kwcmap = kwfilter(kwargs, 'cmap')
        if cb is not None: # Update
            cb.set_cmap(sm.get_cmap(**kwcmap))
            cb.set_clim(sm.get_clim())
            cb.update_normal(sm)
            if self.levels is not None:
                cb.set_ticks(self.levels)
            return cb

        # Plot
        kwtl = kwfilter(kwargs, 'ticklabels')
        kwl = kwfilter(kwargs, 'label')
        # - levels
        levels = self.levels
        if levels is not None:
            kwargs.setdefault('ticks', levels)
        # - extend
        if kwargs.get('extend', None) is None:
            kwargs['extend'] = self._get_extend_(sm.get_clim())
        # - axes
        if isinstance(cax, list):
            cax = self.fig.add_axes(cax)
        # - plot it
#        print kwargs
#        xxx
        cb = self.fig.colorbar(sm, ax=self.axes, cax=cax, **kwargs)
        # - fit to axes limits
        if fit:
            axbbox = self.axes.get_position()
            caxbbox = cb.ax.get_position()
            print 'avant', axbbox, caxbbox
            if cb.orientation=='horizontal':
                cb.ax.set_position([axbbox.x0, caxbbox.y0, axbbox.width, caxbbox.height])
            else:
                cb.ax.set_position([caxbbox.x0, axbbox.y0, caxbbox.width, axbbox.height])
            print 'apres', cb.ax.get_position()
        # - ticks
        if levels is not None and 'format' not in kwargs:
            samp = len(levels)/ticklabels_nmax+1
            if samp>1:
                levels = list(levels)
                for i in xrange(len(levels)):
                    if i%samp:
                        levels[i] = ''
            labels = cb.set_ticklabels(levels)

            if kwtl:
                axis = getattr(cb.ax, ('x' if cb.orientation=='horizontal' else 'y')+'axis')
                P.setp(axis.get_ticklabels(), **kwtl)
        # - label
        label = self.cblabel
        if label is not None:
            cb.set_label(label, **kwl)
        return self.set_obj('colorbar', cb)

    def _get_extend_(self, clim):
        """clim: sm.get_clim() or self.levels"""
        cmin = min(clim)
        cmax = max(clim)
        data = self.get_data(scalar=True)
        vmin = data.min()
        vmax = data.max()
        if cmin>vmin and cmax<vmax:
            return 'both'
        if cmin>vmin:
            return 'min'
        if cmax<vmax:
            return 'max'
        return 'neither'

    def post_plot(self, colorbar=True, savefig=None, savefigs=None, show=True, close=False, **kwargs):
        """

        :Params:

            - **colorbar**, optional: Plot the colorbar.
            - **colorbar_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.ScalarMappable.colorbar`.
        """

        # Keywords
        kw = {}
        for kwtype in 'savefig', 'savefigs', 'show', 'colorbar':
            kw[kwtype] = kwfilter(kwargs, kwtype+'_')

        # Base stuff
        kwargs.update(savefig=None, savefigs=None, show=False)
        Plot.post_plot(self, **kwargs)

        # Colorbar
        if 'extend' in kwargs:
            kw['colorbar'].setdefault('extend', kwargs.get('extend'))
        if colorbar: self.colorbar(**kw['colorbar'])

        # Save it
        self.savefig(savefig, **kw['savefig'])
        self.savefigs(savefigs, **kw['savefigs'])

        # Show or close
        if show:
            self.show(**kw['show'])
        elif close:
            P.close(self.fig)



class Curve(Plot1D):
    """Class for plotting simple curve

    :Params:

        - **data**: Data to plot. It may be a single variable or tuple.
          It a tuple is passed, here are the possible forms:

            - ``(M,)``: Simple scalar.
            - ``(U,V)``: Vector coordinates and the modulus is plotted.
            - ``(M,U,V)``: Modulus and vector coordinates and
              only the modulus is used and plotted.

        - **parg**, optional: Argument passed to :func:`~matplotlib.pyplot.plot`
          after the data.
        - **Specific plot params**: See :meth:`plot`.
        - **Other generic params**: See :class:`Plot`.


    :Example:

        >>> c = Curve(mydata, xmin=3, plot=False)
        >>> c.ymax = 5.
        >>> c.plot('-r')
        >>> c.post_plot(savefig='toto.png')
    """

    def load_data(self, *args, **kwargs):
        """Check order of data

        :Params:

            - **data**: A :mod:`MV2` 1D array.

        :See also: :meth:`Plot.load_data`
        """
        return Plot.load_data(self, *args, **kwargs)

    def plot(self, parg=None, nosingle=False, label=None, err=None, fill_between=False,
            shadow=False, glow=False, log=False, **kwargs):
        """Plot of data as a curve

        :Params:

            - **parg**, optional: Argument passed to :func:`~matplotlib.pyplot.plot`
              after the data.
            - **nosingle**, optional: Single point with missing data around
              are not plotted.
            - **fill_between**, optional: Plot curve using
              :func:`~matplotlib.pyplot.fill_between`.
              Reference value defaults to 0. and may be given
              provided by the parameter.
            - **fill_between_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.fill_between`.
            - **err**, optional: Errors as ``(2,nx)`` array to add to the plot
              using :func:`~matplotlib.pyplot.errorbar`.
            - **err_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.errorbar`.
            - **label**, optional: Alternative label for the plot
              (see also :attr:`~vacumm.misc.core_plot.Plot.label`).
            - **shadow**, optional: A shadow is plotted below the line and points.
            - **shadow_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.Plot.add_shadow`.
            - **glow**, optional: A glow effect is plotted below the line and points.
            - **glow_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.Plot.add_glow`.
        """
        if not self.has_valid_data(): return

        # Data
        xx = self.get_xdata(scalar=0)
        yy = self.get_ydata(scalar=0)
        if err is not None:
            if cdms2.isVariable(err):
                err = err.asma()
            else:
                err = N.ma.asarray(err)
            egood = ~N.ma.getmaskarray(err)
            if egood.ndim==2:
                if egood.shape[0]==2:
                    egood = egood[0]&egood[1]
                else:
                    err = err.ravel()
                    egood = egood.ravel()
            if not egood.any():
                err = None
            elif err.shape[-1]!=len(xx):
                raise PlotError('Error array must be 1D and of size %i'%xx.size)

        # Plot keywords
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        kwerr = kwfilter(kwargs, 'err')
        kwfb = kwfilter(kwargs, 'fill_between')
        kwline = {}
        line_keys = ['color','linewidth','linestyle','markerwidth','drawstyle']
        marker_keys = ['markeredgecolor','markeredgesize','markersize','markerfacecolor','marker','zorder','alpha']
        line_keys.extend(marker_keys)
        for key in line_keys:
            if key in kwargs and kwargs[key] is not None:
                kwline[key] = kwargs[key]
        if parg is None:
            parg = []
        elif not isinstance(parg, list):
            parg = [parg]
        kwline['label'] =  self.label if label is None else label

        # Plot func
        if log:
            if self.order[0]=='d':
                pfunc = self.axes.semilogy
            else:
                pfunc = self.axes.semilogx
        else:
            pfunc = self.axes.plot

        # Plot

        # - main
        ll = pfunc(xx, yy, *parg, **kwline)
        self.register_obj(ll, ['curve', 'lines', 'plot'], **kwargs)

        # - fill_between
        if fill_between is not False:
            kwfb.setdefault('zorder', ll[0].get_zorder()-0.01)
            kwfb.setdefault('linewidth', 0)
            kwfb.setdefault('interpolate', True)
            if 'color' in kwline:
                kwfb.setdefault('color', kwline['color'])
            b0 = yy
            if fill_between is True or fill_between is None:
                b1 = 0.
            elif N.asarray(fill_between).ndim == 2:
                b0, b1 = fill_between
            if self.order[0]=='d':
                ff = self.axes.fill_between(xx, b0, b1, **kwfb)
            else:
                ff = self.axes.fill_betweenx(yy, b0, b1, **kwfb)
            self.add_obj('fill_between', ff)
            self.register_obj(ff, 'fill_between', **kwargs)

        # - error
        if err is not None:
            kwerr.setdefault('ecolor', ll[0].get_color())
            kwerr.setdefault('elinewidth', ll[0].get_linewidth())
            if kwerr.get('zorder', None) and kwerr['zorder']<0:
                kwerr['zorder'] = ll[0].get_zorder() - kwerr['zorder']
            ee = self.axes.errorbar(xx.compress(egood), yy.compress(egood),
                err.compress(egood, axis=-1), fmt='none', **kwerr)
            self.register_obj(ee, 'error', **kwargs)

        # - filters
        if shadow:
            self.add_shadow(ll, 'lines_shadow', **kwsh) # 'lines_shadow'
            if err is not None:
                self.add_shadow(ee, 'error_shadow', **kwsh) # 'error_shadow'
        if glow:
            self.add_glow(ll, 'lines_glow', **kwgl) # 'lines_glow'
            if err is not None:
                self.add_glow(ee, 'error_glow', **kwgl) # 'error_glow'
        # - single dots
        if not nosingle and len(xx) > 2:
            mask = self.mask.copy()
            for m in [self.mask[:-2],self.mask[2:]]:
                mask[1:-1] = mask[1:-1] | 1-m
            mask[0] = mask[0] or not self.mask[1]
            mask[-1] = mask[-1] or not self.mask[-2]
            if N.sometrue(1-mask):
                kwmark = {}
                for key in marker_keys:
                    if kwargs.has_key(key):
                        kwmark[key] = kwargs[key]
                if kwargs.has_key('zorder'): kwmark['zorder'] = kwargs['zorder']
                if self.xtype!='d':
                    xx = N.ma.array(xx,mask=mask)
                else:
                    yy = N.ma.array(yy,mask=mask)
                kwmark.update(label='',color=ll[0].get_color())
                dots = self.axes.plot(xx, yy, '.', **kwmark)
                self.register_obj(dots, ['lines_dots', 'plot'], **kwargs)
                if shadow: self.add_shadow(dots, 'lines_dots_shadow', **kwsh) # 'lines_dots_shadow'
                if glow: self.add_glow(dots, 'lines_dots_glow', **kwsh) # 'lines_dots_glow'



class Bar(Plot1D):
    """Class for plotting simple curve

    :Params:

        - **data**: 1D array.
        - **Specific plot params**: See :meth:`plot`.
        - **Other generic params**: See :class:`Plot`.


    :Example:

        >>> c = Bar(mydata, xwidth=.8, plot=False, order='zd')
        >>> c.ymax = 5.
        >>> c.plot()
        >>> c.post_plot(savefig='rain.png')
    """

    def plot(self, width=1., lag=0, align='center', shadow=False, glow=False, offset=None,
        label=None, **kwargs):
        """Plot data as bar plot

        :Params:

            - **width**, optional: Relative width of the bars(``0<width<1``).
              a width of ``1`` means that successive are bars are joined.
            - **lag**, optional: Relative lag to apply to the position
            - **align**, optional: Alignment relative to coordinates.
            - **shadow**, optional: Add a shadow to the bars.
            - **shadow_<param>**, optional: ``<param>`` is passed to :meth:`add_shadow`.
            - **glow**, optional: Add a glow effect to the bars.
            - **glow_<param>**, optional: ``<param>`` is passed to :meth:`add_shadow`.
            - **offset**, optional: Bars start at ``offset``.
            - **label**, optional: Alternative label for the plot
              (see also :attr:`~Plot.label`).
            - **bar_<param>**, optional: ``param`` is passed to :func:`matplotlib.pyplot.bar` (or :func:`matplotlib.pyplot.barh`).

        :Example:

            >>> Bar(rain1).plot(width=.45, align='left', color='cyan')
            >>> Bar(rain2).plot(width=.45, lag=.5, align='left', color='b')
        """
        if not self.has_valid_data(): return

        # Data
        data = self.get_data()[0]
        axis = self.get_axis_data().astype('d')
        bounds = meshbounds(axis)
        widths = N.diff(bounds)
        axis += widths*(lag if lag else 0.)
        widths *= width
        if cdms2.isVariable(offset): offset = offset.asma()

        # Plot keywords
        kwsh = kwfilter(kwargs,'shadow_')
        kwgl = kwfilter(kwargs,'glow_')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        kwbar = kwfilter(kwargs, 'bar_')
        bar_keys = ['color','linewidth','linestyle','markerwidth', 'edgecolor', 'zorder', 'alpha',  'log']
        for key in bar_keys:
            if kwargs.has_key(key):
                kwbar[key] = kwargs[key]
        kwbar['label'] =   self.label  if label is None else label

        # Log data
        if kwbar.get('log', False) and (data<=0.).all(): return

        # Plot
        # - main
        func_name = 'bar'+('h' if self.xtype=='d' else '')
        plot_func = getattr(self.axes, func_name)
        pp = plot_func(axis, data, widths, offset, align=align, **kwbar)
        self.register_obj(pp, ['patches', func_name], **kwargs)
        # - filters
        if shadow: self.add_shadow(pp, **kwsh)
        if glow: self.add_glow(pp, **kwgl)

class QuiverKey:

    def quiverkey(self, qv, value=None, value_mode=80, **kwargs):
        """Add a quiver key to the plot

        See :meth:`Plot.quiverkey` for arguments.

        :Params:

            - **qv**: Results of :func:`~matplotlib.pyplot.quiver`.
            - **pos**, optional: Position of key for arrow .
            - **text**, optional: Text or format with variables 'value' and 'units'.
            - **value**, optional: Numeric value for key (used by text).
            - **units**, optional: Units for key (used by text).
            - **latex_units**, optional: Interpret units using latex.
            - Extra keywords are passed to :func:`~matplotlib.pyplot.quiverkey`.
        """
        # Value
        if value is None:
            m,u,v = self.get_data()
            value = get_quiverkey_value((u, v), mode=value_mode)
            del m, u, v

        return Plot.quiverkey(self, qv, value, **kwargs)


class Stick(QuiverKey, ScalarMappable, Curve):
    """Class for makeing a stick plot (vectors on a line)

    :Params:

        - **udata**: 1D array of intensities along X.
        - **vdata**: 1D array of intensities along V.
        - **Specific data loading params**: See :meth:`load_data`.
        - **Specific plot params**: See :meth:`plot`.
        - **Specific plot finalization params**: See :meth:`~vacumm.misc.core_plot.ScalarMappable.post_plot`.
        - **Generic  params**: See :class:`Plot`.


    :Example:

        >>> c = Stick(u10, v10, savefig='cart.png')
        >>> c = Stick(r, theta, polar=True,width=.8, plot=False, order='zd')
        >>> c.plot()
        >>> c.post_plot(savefig='polar.png')
    """

#    _primary_attributes = Plot._primary_attributes + ['levels', 'nmax', 'nmax_levels', 'cmap']
#    _secondary_attributes = Plot._secondary_attributes + ['cblabel']

    def __init__(self, udata, vdata, polar=False, degrees=True, **kwargs):
        Curve.__init__(self, (udata, vdata), polar=polar, degrees=degrees, **kwargs)


    # Polar case
    def load_data(self, data, **kwargs):
        """Load variables

        :Special params:

            - **udata**: X or radial component of arrows.
            - **vdata**: Y or directional component of arrows.
            - **polar**, optional: Consider polar coordinates: ``(u, v) -> (rho, theta)``
            - **degrees**, optional: If True (default), trat ``theta`` as degrees, else radians.

        :Tasks:

            #. Calls :meth:`Plot.load_data`.
            #. Deals with polar case, if keyword ``polar==True``.
        """
        polar = kwargs.pop('polar')
        degrees = kwargs.pop('degrees')
        Curve.load_data(self, data, **kwargs)

        if polar:
            m, u, v = self.data
            m = u
            angle = (v*N.pi/180.) if degrees else v
            u = m * MV2.cos(angle)
            v = m * MV2.sin(angle)
            del angle
            if hasattr(m, 'units'):
                u.units = v.units = m.units
            if hasattr(m, 'long_name'):
                u.long_name = 'X component of '+m.long_name
                v.long_name = 'Y component of '+m.long_name
            self.data[:] = m, u, v


    # Read vector data and not modulus
    def get_xdata(self, scalar=1, masked=False, bounds=False):
        return Curve.get_xdata(self, scalar=scalar, masked=masked, bounds=bounds)
    get_xdata.__doc__ = Curve.get_xdata.__doc__
    def get_ydata(self, scalar=2, masked=False, bounds=False):
        return Curve.get_ydata(self, scalar=scalar, masked=masked, bounds=bounds)
    get_ydata.__doc__ = Curve.get_ydata.__doc__

    def plot(self, mod=False, pos=None, line=False, color='k', alpha=1, quiverkey=True,
        headwidth=None, headlength=None, headaxislength=None, width=None, scale=None,
        minshaft=None, minlength=None,
        shadow=False, glow=False, cmap=None, levels=None, label=None,
        anomaly=False, **kwargs):
        """Main plot

        :Params:

            - **pos**, optional: Position of the arrow tails.
              It defaults to the middle of the appropriate axis.
            - **mod**, optional: Plot the curve of the modulus.
            - **mod_<param>***, optional: ``<param>`` is passed to
              :meth:`Curve.plot` when plotting the modulus.
            - **color**: Can be either

                - ``"mod"``: The color is function of the modulus.
                - A normal color argument for :func:`~matplotlib.pyplot.quiver`
                  (single or list/array).

            - **line**, optional: Add a transversal line along the arrow tails.
            - **alpha**, optional: Opacity.
            - **scale**, optional: Scale of arrows
              (see :func:`~matplotlib.pyplot.quiver`).
            - **headwidth**, optional: Head width of arrows
              (see :func:`~matplotlib.pyplot.quiver`).
            - **headlength**, optional: Head length of arrows
              (see :func:`~matplotlib.pyplot.quiver`).
            - **headaxislength**, optional: Length of arrow head on axis
              (see :func:`~matplotlib.pyplot.quiver`).
            - **minlength**, optional: See :func:`~matplotlib.pyplot.quiver`.
            - **minshaft**, optional: See :func:`~matplotlib.pyplot.quiver`.
            - **quiverkey**, optional: Add key to scale arrows.
            - **quiverkey_<param>**, optional: ``<param>`` is passed
              to :meth:`~vacumm.misc.core_plot.QuiverKey.quiverkey`.
            - **shadow**, optional: Add a drop shadow.
            - **shadow_<param>**, optional: ``<param>`` is passed to :meth:`add_shadow`
            - **glow**, optional: Add a drop glow effect.
            - **glow_<param>**, optional: ``<param>`` is passed to :meth:`add_glow`
            - **cmap**, optional: Colormap to use when ``color="mod"``.
            - **cmap_<param>**, optional: Passed to
              meth:`~vacumm.misc.core_plot.ScalarMappable.get_cmap` to tune
              colormap.
            - **levels**, optional: Levels of values to use when ``color="mod"``.
            - **levels_<param>**, optional: Passed to
              meth:`~vacumm.misc.core_plot.ScalarMappable.get_levels` to tune levels.

        """
        if not self.has_valid_data(): return

        # Keywords
        kwmod = kwfilter(kwargs, 'mod_')
        kwqv = kwfilter(kwargs, 'quiver_')
        kwqvkey = kwfilter(kwargs, 'quiverkey_')
        kwline = kwfilter(kwargs, 'line_')
        kwcmap = kwfilter(kwargs, 'cmap_')
        kwlevels = kwfilter(kwargs, 'levels_', defaults=dict(anomaly=anomaly))
        kwsh = kwfilter(kwargs, 'shadow_')
        kwgl = kwfilter(kwargs, 'glow_')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        kwqv['label'] =  self.label if label is None else label

        # Cmap and levels
        numcolors = color=='mod' or color is True or \
            (isinstance(color, N.ndarray) and color.dtype.char!='S')
        if numcolors:
            if levels is None:
                levels = self.get_levels(**kwlevels)
            else:
                self.levels = levels
#            color = None
            if cmap is not None:
                self.cmap = cmap
            cmap = self.get_cmap(cmap, **kwcmap)

        # Add modulus
        if mod:
            for att in 'color', 'alpha', 'shadow', 'glow':
                value = eval(att)
                if value is not None:
                    kwmod.setdefault(att, eval(att))
#            if color is not None and isinstance(kwmod['color'], (list, N.ndarray)):
            if kwmod['color'] is not None and isinstance(kwmod['color'], (list, N.ndarray)):
                kwmod['color'] = kwmod['color'][0]

            Curve.plot(self, **kwmod)

        # Data
        mm, uu, vv = self.get_data(scalar=False)

        # Coordinates
        pos = kwargs.pop('ycenter', pos) # compat
        if self.ytype=='d':
            if pos is None: pos = (self.axes.get_ylim()[0]+self.axes.get_ylim()[1])/2.
            xx = self.get_xdata()
            yy = N.zeros(len(xx))+pos
        else:
            if pos is None: pos = (self.axes.get_xlim()[0]+self.axes.get_xlim()[1])/2.
            yy = self.get_ydata()
            xx = N.zeros(len(yy))+pos

        # Plot
        for att in ['color', 'headwidth', 'headlength', 'alpha', 'headwidth',
            'headlength', 'headaxislength', 'minshaft', 'minlength', 'width', 'scale']:
            value = eval(att)
            if value is not None:
                kwqv.setdefault(att, eval(att))
        args = [xx, yy, uu, vv]
        if color is None or color is True:
#            kwqv['vmin'] = self.vmin
#            kwqv['vmax'] = self.vmax
            args.append(mm)
            kwqv.pop('color')
        elif numcolors: # numerical values for colors
            args.append(color)
            kwqv.pop('color')
        elif isinstance(color, (list, N.ndarray)):
            kwqv.setdefault('color', color)
        qv = self.register_obj(self.axes.quiver(*args, **kwqv), ['quiver', 'patches'], **kwargs)
        if numcolors:
            kwqv['cmap'] = cmap
            qv.set_clim(self.vmin, self.vmax)
            self.set_obj('scalar_mappable', qv)
            sel.axes._sci(av)
        qv.set_alpha(alpha)

        # Filters
        if shadow: self.add_shadow([qv], 'quiver_shadow', **kwsh)
        if glow: self.add_glow([qv], 'quiver_glow', **kwgl)

        # Zorder of modulus
        if mod and not kwmod.has_key('zorder'):
            oo = self.get_obj('lines')
            if self.get_obj('lines_shadow') is not None:
                oo.extend(self.get_obj('lines_shadow'))
            if self.get_obj('lines_dot') is not None:
                oo.extend(self.get_obj('lines_dot'))
                if self.get_obj('lines_dot_shadow') is not None:
                    oo.extend(self.get_obj('lines_dot_shadow'))
            dz = qv.get_zorder()-oo[0].get_zorder()
            if dz<=0:
                for o in oo:
                    o.set_zorder(o.get_zorder()+dz-.2)

        # Quiver key
        if quiverkey:
            self.quiverkey(qv, **kwqvkey)

        # Base line
        if line:
            hv, pp = ('h', yy[0]) if self.ytype=='d' else ('v', xx[0])
            funcname = 'ax%sline'%hv
            func = getattr(self.axes, funcname)
            self.register_obj(func(pp, **kwline), [funcname], **kwargs)

        return qv

    def format_axes(self, **kwargs):
        if self.get_obj('quiver') is not None and self.get_obj('curve') is None: #
            xy = 'y' if self.ytype=='d' else 'x'
            kwargs.setdefault(xy+'hide', True)
            kwargs.setdefault(xy+'ticks', [])
#            kwargs.setdefault(xy+'label', '')
        return Curve.format_axes(self, **kwargs)

class Plot2D(ScalarMappable, QuiverKey, Plot):

    _order = ['--']
    _primary_attributes = ScalarMappable._primary_attributes + ['fill']
    _secondary_attributes = ScalarMappable._secondary_attributes + ['cblabel']
    rank = 2
    _plotter = None

    def _set_axes_(self, xaxis=None, yaxis=None, xatts=None, yatts=None,
            x2d=None, y2d=None, x2db=None, y2db=None, **kwargs):
        """Change axes and their attributes of the variables

        The main goal is to deal 2D plots on special grids where the grid
        cannot be properly stored internally in variables (like a section grid).


        :Params:

            - **x/yaxis**, optional: 1D or 2D array or axes to use for axis,
              instead of internal axes of input variables.
            - **x/y2d**, optional: (ny,nx) array of cell center coordinates.
            - **x/y2db**, optional: (ny+1,nx+1) array of cell bounds coordinates
              used in pcolor like plots.
            - **x/yatts**, optional: Dictionnary of alter given ``x/yaxis``.

        """



        # Get axes and adjust order
        if (xaxis is not None or yaxis is not None) and self.has_data():

            # Adjust order
            for ivar, var in enumerate(self.data):
                order = var.getOrder()
                #self.data[ivar] = var2d(var, xaxis, yaxis, xatts=xatts, yatts=yatts)
                for i, axis in enumerate([yaxis, xaxis]):
                    c = getattr(axis, 'axis', '-').lower()
                    if order[i]=='-' and c!='-': order[i] = c
                set_order(self.data[ivar], order)
                self.data[ivar]._nogridorder = True
        if self.has_data():
            xax = get_axis(self.data[0], -1, strict=True)
            yax = get_axis(self.data[0], -2, strict=True)
        if xaxis is None:
            if self.has_data():
                xaxis = xax
            else:
                raise PlotError('No xaxis data available')
        elif not isaxis(xaxis):
            xaxis = MV2.array(xaxis, copy=0)
            if self.has_data():
                cp_atts(xax, xaxis, overwrite=False)
        if yaxis is None:
            if self.has_data():
                yaxis = yax
            else:
                raise PlotError('No yaxis data available')
        elif not isaxis(yaxis):
            yaxis = MV2.array(yaxis, copy=0)
            if self.has_data():
                cp_atts(yax, yaxis, overwrite=False)

        # Change attributes
        if xatts is not None:
            set_atts(xaxis, xatts)
        if yatts is not None:
            set_atts(yaxis, yatts)

        # Set axes
        self.x = xaxis
        self.y = yaxis

        # Mesh form
        xdata = self.get_xdata(masked=False)
        ydata = self.get_ydata(masked=False)
        # - centers
        if x2d is None or y2d is None:
            self.x2d, self.y2d = meshgrid(xdata, ydata)
        if x2d is not None:
            self.x2d = x2d
        if y2d is not None:
            self.y2d = y2d
        self.x2dr,self.y2dr = self.x2d, self.y2d # save raw values
        # - bounds
        if x2db is None or y2db is None:
            self.x2db, self.y2db = meshbounds(xdata, ydata)
        if x2db is not None:
            if x2db.ndim==1:
                self.x2db = N.resize(x2db, (self.x2d.shape[0]+1, x2db.shape[0]))
            elif x2db.shape[0] == self.x2d.shape[0]:
                self.x2db = N.zeros((self.x2d.shape[0]+1, self.x2d.shape[1]+1))
                self.x2db[1:] = shift1d(x2db, 1, axis=0, mode='same')
                self.x2db[:-1] += shift1d(x2db, -1, axis=0, mode='same')
                self.x2db[1:-1] *= 0.5
            else:
                self.x2db = x2db
        if y2db is not None:
            if y2db.ndim==1:
                self.y2db = N.resize(y2db, (self.y2d.shape[1]+1, y2db.shape[0])).T
            elif y2db.shape[1] == self.y2d.shape[1]:
                self.y2db = N.zeros((self.y2d.shape[0]+1, self.y2d.shape[1]+1))
                self.y2db[:, 1:] = shift1d(y2db, 1, axis=1, mode='same')
                self.y2db[:, :-1] += shift1d(y2db, -1, axis=1, mode='same')
                self.y2db[:, 1:-1] *= 0.5
            else:
                self.y2db = y2db
        self.x2dbr,self.y2dbr = self.x2db, self.y2db # save raw values

    def load_data(self, data, **kwargs):
        """Load data and axes

        :Tasks:

            #. Call to :meth:`Plot.load_data`.
            #. Compute :attr:`x2d`, :attr:`y2d`, :attr:`x2db`, :attr:`y2db`

        :Params:

            - **data**: A single or a tuple or 2D :mod:`MV2` arrays.

              - Single variable: Plot 2D scalar field with contours, filled contours, pcolor or image.
              - 2-tuple ``(U,V)``: Plot arrows.
              - 3-tuple ``(M,U,V)``: Plot both scalar field ``M`` and arrows.

        :Specific attributes:

            .. attribute:: x2d

                2D version of the X axis data.

            .. attribute:: y2d

                2D version of the Y axis data.

            .. attribute:: x2db

                2D bounds coordinates of the X axis data.

            .. attribute:: y2db

                2D bounds coordinates of the Y axis data.
        """
#        # Check
#        if data is None: return

        # Normal load
        Plot.load_data(self, data, **kwargs)

    def _transpose_axes_(self):

        # Normal transpose
        Plot._transpose_axes_(self)

        # Mesh form
        # - centers
        self.x2d, self.y2d = self.y2d.T, self.x2d.T
        self.x2dr, self.y2dr = self.y2dr.T, self.x2dr.T
        # - bounds
        self.x2db, self.y2db = self.y2db.T, self.x2db.T
        self.x2dbr, self.y2dbr = self.y2dbr.T, self.x2dbr.T


    def get_xyscaler(self, guess=True):
        """Get :attr:`xyscaler`"""
        return self.get_obj('xyscaler')
    def set_xyscaler(self, scaler):
        """Set :attr:`xyscaler`"""
        self.set_obj('xyscaler', scaler)
    def del_xyscaler(self):
        """Del :attr:`xyscaler`"""
    xyscaler = property(get_xyscaler, set_xyscaler,
        del_xyscaler, doc="""Function to rescale  X and Y coordinates.
            ``None`` is returned if no scaler is available.

            A typical scaler is a :class:`mpl_toolkits.basemap.Basemap` or
            a :class:`mpl_toolkits.pyproj.Proj` instance that converts
            positions from degrees to meters using a geographic
            projection.

            This scaler is used by :meth:`plot_quiver` for undersampling
            based on resolution.

            :Example:

                >>> xyscaler = myplot.xyscaler
                >>> if xyscaler is not None:
                    ... x,y = myplot.xyscaler(x, y)

            :Return: A callable function or ``None`` if no scaling is possible
            """)

    @staticmethod
    def _fill_method_(fill='pcolormesh', pcolor=None, nofill=None, **kwargs):
        if fill is None:
            fill = get_config_value('vacumm.misc.plot', 'fill')
            if fill is not None:
                if fill.isdigit():
                    fill = int(fill)
                elif fill in ['True', 'None', 'False']:
                    fill = eval(fill)
        if fill is True or str(fill).startswith('pcolor'):
            fill = 'pcolormesh'
        if nofill:
            fill = 0
        elif pcolor:
            fill = 'pcolormesh'
        elif pcolor in (False, 0, 2):
            fill = 'contourf'
        elif pcolor == 3:
            fill = 'imshow'
        if fill == 'nofill' or fill is False:
            fill = 'no'
        elif fill == 'contour':
            fill = 'contourf'
        elif isinstance(fill, int):
            fill = N.clip(fill, 0, 4)
            fill = ['no', 'pcolor', 'contourf', 'imshow', 'scatter'][fill]
        return fill


    def get_fill(self, fill=None, pcolor=None, nofill=None, **kwargs):
        """Get the :attr:`fill` attribute"""
        if fill is None:
            fill = self.get_obj('fill')
        return self._fill_method_(fill=fill, pcolor=pcolor, nofill=nofill)
    def set_fill(self, fill):
        """Set the :attr:`fill` attribute"""
        self.set_obj('fill', fill)
    def del_fill(self):
        """Del the :attr:`fill` attribute"""
        self.del_obj('fill')
    fill = property(get_fill, set_fill, del_fill,
        doc="Fill method for 2D plots: 'pcolor', None, True, False, 'contourf', 'imshow', 'pcolormesh'")

    def plot_fill(self, norm=None, shading='flat', alpha=1, extend=None,
        zorder=None, shadow=False, glow=False, **kwargs):
        """Plot filled stuff

        :Params:

            - **fill**, optional: Type filling.
              It defaults to ``fill`` option of the ``[vacumm.misc.plot]`` section
              of the configuration (see :func:`~vacumm.config.edit_config`).

                - ``0`` or ``"no"`` or ``"nofill"`` or ``False``: No filling.
                - ``1`` or ``"pcolor"`` or ``"pcolormesh"``: Fill using :func:`~matplotlib.pyplot.pcolor`
                  or :func:`~matplotlib.pyplot.pcolormesh`.
                - ``2`` or ``"imshow"``: Fill using :func:`~matplotlib.pyplot.imshow`.
                - ``3`` or ``"contourf"``: Fill using :func:`~matplotlib.pyplot.contourf`.

            - **fill_<param>**: ``<param>`` is passed to the corresponding plot function.
            - **nofill**, optional: Implies ``fill=0``.
            - **cmap**, optional: Colormap (defaults to ``"magic"``).
            - **cmap_<param>**, optional: Passed to
              :meth:`~vacumm.misc.core_plot.ScalarMappable.get_cmap` to tune
              colormap.
            - **norm**, optional: :class:`~matplotlib.colors.Normalize` instance.
            - **alpha**, optional: Opacity.
            - **shading**, optional: Shading with :func:`~matplotlib.pyplot.pcolor`.
            - **extend**, optional: Let's :func:`~matplotlib.pyplot.contourf` add
              contours to cover all data range.
            - **zorder**, optional: Plot order.

        :Tasks:

            #. Guess the fill method.
            #. Get the colormap with :meth:`~ScalarMappable.get_cmap`.
            #. Get the scalar data with :meth:`~Plot.get_data`.
            #. Calls :func:`~matplotlib.pyplot.imshow` or
               :func:`~matplotlib.pyplot.pcolor` or
               :func:`~matplotlib.pyplot.pcolormesh` or
               :func:`~matplotlib.pyplot.contourf`.
            #. Calls :func:`~matplotlib.pyplot.clabel`

        .. note::

            :meth:`~ScalarMappable.get_levels` must have been previously called.

        """
        if not self.has_valid_data(): return



        # Keywords
        kwfill = kwfilter(kwargs, 'fill_')
        kwnorm = kwfilter(kwargs, 'norm')
        kwcmap = kwfilter(kwargs, 'cmap', copy=True)
        kwlevels = kwfilter(kwargs, 'levels_')
        if zorder is not None:
            kwfill.setdefault('zorder', zorder)
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')

        # Fill method
        fill = self.fill_method = self.get_fill(**kwargs)
        if fill == 'no': return

        # Levels
        levels = self.get_levels(**kwlevels)

        # Norm
#        from matplotlib.colors import Normalize
        if norm is None:
            from color import StepsNorm
            norm = StepsNorm(levels, **kwnorm)
        elif norm is True:
            norm = Normalize(min(levels), max(levels))
        elif not isinstance(norm, Normalize):
            norm = None

        # Color map
        cmap = self.get_cmap(**kwcmap)

        # Data
        data = self.get_data(scalar=True)

        # Plots
        if fill == 'imshow':
            kwfill.setdefault('origin', 'lower')
            pp = self._plotter.imshow(data, cmap=cmap, alpha=alpha,
                norm=norm, vmin=self.vmin, vmax=self.vmax,
                extent=[self.x2db.min(), self.x2db.max(),
                self.y2db.min(), self.y2db.max()], **kwfill)
            if alpha != 1:
                pp.set_alpha(alpha)
            self.set_obj('imshow', pp)

        elif fill.startswith('pcolor'):
            if fill=='pcolormesh': # self.mask is MV2.nomask or
                func = self._plotter.pcolormesh
                kwfill.setdefault('edgecolors', 'None')
            else:
                kwfill.setdefault('edgecolors', 'None')
                kwfill.setdefault('linewidth', 0)
                func = self._plotter.pcolor
            for att, val in dict(cmap=cmap,alpha=alpha, norm=norm,
                    vmin=self.vmin, vmax=self.vmax, edgecolors='none').items():
                kwfill.setdefault(att, val)
            pp = func(self.x2db, self.y2db, data, **kwfill)
            if alpha != 1:
                pp.set_alpha(alpha)
                pp.set_linewidth(0.)
            self.set_obj('pcolor', pp)


        elif fill == 'contourf':
            if extend is None:
                extend = self._get_extend_(levels)
            dict_check_defaults(kwfill, cmap=cmap,
                levels=levels, alpha=alpha, norm=norm,
                extend=extend, edgecolors='none')
            pp = self._plotter.contourf(self.x2d, self.y2d, data, copy=0, **kwfill)
            self.set_obj('contourf', pp)
            if alpha != 1:
                for col in pp.collections:
                    try:
                        col.set_alpha(alpha)
                    except:
                        pass

        elif fill == 'scatter':
            pp = self._plotter.scatter(self.x2d.ravel(), self.y2d.ravel(),
                c=data.ravel(), cmap=cmap, alpha=alpha,
                norm=norm, vmin=self.vmin, vmax=self.vmax,
                **kwfill)
            if alpha != 1:
                pp.set_alpha(alpha)
            self.set_obj('scatter', pp)
        else:
            raise PlotError('Wrong fill method')


        # Register
        self.set_obj(['fill', 'scalar_mappable'], pp)
        self.register_obj(pp,  'patches', **kwargs)
        self.axes._sci(pp)


        # Filters
        if shadow or glow:
            obj = pp.collections if hasattr(pp, 'collections') else pp
            if shadow: self.add_shadow(obj, 'fill_shadow', **kwsh)
            if glow: self.add_glow(obj, 'fill_glow', **kwgl)

        # Hacks
#        if zorder is not None: pp.set_zorder(zorder)
        for key in 'edgecolors', 'linewidths':
#            if hasatt(pp, 'set_'+key):
#                getattr(pp, 'set_'+key)(kwfill[key])
            if not kwfill.has_key(key): continue
            if hasattr(pp, 'collections'):
                for p in pp.collections:
                    getattr(p, 'set_'+key)(kwfill[key])
#        if kwfill.has_key('edgecolor'):
#            for p in pp.collections:
#                p.set_edgecolor(kwfill['edgecolor'])
#        if kwfill.has_key('linewidth'):
#            for p in pp.collections:
#                p.set_linewidth(kwfill['linewidth'])


    def plot_contour(self, zorder=None, alpha=1, clabel=None, linewidths=None,
            colors='k', shadow=False, glow=False, **kwargs):
        """Plot contour lines

        :Params:

            - **contour_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.contour`.
            - **zorder**, optional: Plot order.
            - **alpha**, optional: Opacity.
            - **linewidths**, optional: Contour linewidths.
            - **clabel**, optional: Add contour labels with
              :func:`~matplotlib.pyplot.clabel`. If not specified, it is taken from
              config section ``[vacumm.misc.plot]`` and config option ``clabel``.
            - **clabel_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.clabel`..
            - **contour_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.contour`.

        :Tasks:

            #. Get the scalar data with :meth:`~Plot.get_data`.
            #. Calls :func:`~matplotlib.pyplot.contour`.
            #. Calls :func:`~matplotlib.pyplot.clabel`

        """
        if not self.has_valid_data():
            return

        # Keywords

        if clabel is None:
            try:
                clabel = eval(get_config_value('vacumm.misc.plot', 'clabel'))
            except:
                clabel = eval(get_config_value('vacumm.misc.plot', 'clabel', user=False))

        kw = kwfilter(kwargs, 'contour', defaults=dict(levels=self.levels))
        kwcl = kwfilter(kwargs, 'clabel')
        kwsh = kwfilter(kwargs, 'shadow')
        kwsh.update(kwfilter(kw, 'shadow_'))
        kwgl = kwfilter(kwargs, 'glow')
        kwgl.update(kwfilter(kw, 'glow_'))
        kwcmap = kwfilter(kwargs, 'cmap', copy=True)
#        kwcmap.update(kwfilter(kw, 'cmap_'))
        shadow = kw.get('shadow', shadow)
        glow = kw.get('glow', glow)
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')

        # Setup
        if kwargs.has_key('fmt'):
            kwcl['fmt'] = kwargs['fmt']
        kwcl.setdefault('fmt', '%g')
        kw.setdefault('alpha', .5 if self.fill_method.startswith('pcolor') or
            self.fill_method.startswith('imshow') else alpha)
        kwcl.setdefault('alpha', alpha)
        if 'linewidth' in kw:
            linewidths = kw.pop('linewidth')
        if linewidths is None and 'linewidth' in kwargs:
            linewidths = kwargs.pop('linewidth')
        if linewidths is not None:
            kw.setdefault('linewidths', linewidths)
        if zorder is not None:
            kw.setdefault('zorder', zorder)

        # Colors
        if not 'cmap' in kw:
            cmap = None
        else:
            cmap = self.get_cmap(kw.get('cmap', None))
        if colors is None and 'color' in kwargs:
            colors = kwargs.pop('color')
        colors = kw.get('colors', colors)
        if cmap is not None:
            colors = None
        else:
            cmap = self.get_cmap(**kwcmap)
            if colors is not None:
                cmap = None
        kw['cmap'] = cmap
        kw['colors'] = colors

        # Data
        data = self.get_data(scalar=True)

        # Contours
        cc = self._plotter.contour(self.x2d, self.y2d, data, **kw)
        self.set_obj('contour', cc)
        self.register_obj(cc, 'lines', **kwargs)
        if shadow: self.add_shadow(cc.collections, 'contour_shadow', **kwsh)
        if glow: self.add_glow(cc.collections, 'contour_glow', **kwgl)


        # Labels
        if len(cc.collections) and not kwcl.pop('hide', not clabel):
            kwcl.setdefault('zorder', cc.collections[0].get_zorder())
            clshadow = kwcl.pop('shadow', shadow)
            kwclsh = kwfilter(kwcl, 'shadow')
            clglow = kwcl.pop('glow', glow)
            kwclgl = kwfilter(kwcl, 'glow')
            dict_copy_items(kwargs, [kwclsh, kwclgl], 'anim')
            cl = self.axes.clabel(cc, **kwcl)
            self.set_obj('clabel', cl)
            self.register_obj(cl, 'text', **kwargs)
            for l in cl:
                l.set_alpha(kwcl['alpha'])
                l.set_zorder(kwcl['zorder'])
            if clshadow: self.add_shadow(cl, **kwclsh)
            if clglow: self.add_glow(cl, **kwclgl)

    def plot_quiver(self, zorder=None, quiverkey=True, barbs=False,
            shadow=False, glow=False, quiver_cmap=None,
            quiver_vmin=None, quiver_vmax=None,
            quiver_samp=None, quiver_xsamp=None, quiver_ysamp=None, quiver_res=None,
            quiver_relres=None, quiver_xres=None, quiver_xrelres=None, quiver_yres=None,
            quiver_yrelres=None, quiver_res_scaler=None, quiver_nauto=None, **kwargs):
        """Plot arrows

        You can undersample arrows using direct undersampling (parameters
        with "samp") or undersampling based on resolution (using
        :func:`~vacumm.misc.grid.masking.resol_mask`).
        Resolution undersampling may be in input (like degrees) or
        transformed (like meters) coordinates. Transformed
        coordinates are deduced from input coordinates using
        ``quiver_res_scaler``, which defaults to :attr:`xyscaler`.

        .. note::

            Direct unsampling is incompatible with resolution undersampling:
            the former prevails against the latter.

        :Params:

            - **quiver_norm**,optional: Normalize/colorize arrows

                - ``0`` or ``None``: No normalization, no colorization (default)
                - ``1``: Normalization, no colorization.
                - ``2``: Normalization, colorization.
                - ``3``: No normalization, colorization.

            - **quiver_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.quiver`.
            - **barbs**, optional: Plot wind barbs instead of arrows
              (see :func:`~matplotlib.pyplot.barbs`).
            - **zorder**, optional: Plot order.
            - **alpha**, optional: Opacity.
            - **shadow**, optional: Add shadow below arrows.
            - **glow**, optional: Add glow effect to arrows.
            - **quiverkey**, optional: Add key to scale arrows.
            - **quiverkey_<param>**, optional: ``<param>`` is passed
              to :meth:`~vacumm.misc.core_plot.QuiverKey.quiverkey`.
            - **cmap**, optional: Colormap to use when ``color="mod"`` (defaults to
              ``"magic"``).
            - **quiver_samp**, optional: Horizontal sampling of arrows (in both directions) [default: 1]
            - **quiver_x/ysamp**, optional: Sampling along X/Y [default: quiver_samp]
            - **quiver_res**, optional: Horizontal resolution of arrows (in both directions)
              for undersampling [default: None]
                If ``'auto'``, resolution is computed so as to have at max ``quiver_nauto``
                arrow in along an axis. If it is a :class:`complex` type, its imaginary part
                set the ``quiver_nauto`` parameter and ``quiver_res`` is set to ``'auto'``.
            - **quiver_x/yres**, optional: Same along X/Y [default: quiver_res]
            - **quiver_relres**, optional: Relative resolution (in both directions).

              - If > 0, = ``mean(res)*relres``.
              - If < -1, = ``min(res)*abs(relres)``.
              - If < 0 and > -1, = ``max(res)*abs(relres)``

            - **quiver_x/yrelres**, optional: Same along X/Y [default: quiver_relres]


        :Tasks:

            #. Get the scalar data with :meth:`~Plot.get_data`.
            #. Calls :func:`~matplotlib.pyplot.quiver`.
            #. Calls :func:`~matplotlib.pyplot.quiverkey`

        """
        if not self.has_valid_data(): return

        # Get data
        if len(self.data)!=3: return
        x2d, y2d = self.x2d, self.y2d
        x2dr, y2dr = self.x2dr, self.y2dr
        mm, uu, vv = self.get_data(scalar=False)
        del mm
        mm = N.ma.sqrt(uu**2+vv**2)

        # Keywords
        kwqv = kwfilter(kwargs, 'quiver')#, defaults=dict(angles='uv'))
        kwqvkey = kwfilter(kwargs, 'quiverkey')
        if zorder is not None:
            kwqv['zorder'] = zorder
#        if kwqv.has_key('width') and kwqv['width'] > .01:
#            kwqv['width'] *= 0.001
        shadow = kwqv.pop('shadow', shadow)
        glow = kwqv.pop('glow', glow)
        kwsh = kwfilter(kwqv, 'shadow')
        kwgl = kwfilter(kwqv, 'glow')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        kwcmap = kwfilter(kwargs, 'cmap', copy=True) if not self.get_scalar_mappable() else {}
        kwcmap.update(kwfilter(kwqv, 'cmap'))

        # Undersampling
        # - setup
        if quiver_samp is not None:
            if quiver_xsamp is None: quiver_xsamp = quiver_samp
            if quiver_ysamp is None: quiver_ysamp = quiver_samp
        if quiver_xsamp is not None:
            quiver_xres = quiver_xrelres = False
        if quiver_ysamp is not None:
            quiver_yres = quiver_yrelres = False
        # - indirect
        rmask = resol_mask((x2dr, y2dr), res=quiver_res, xres=quiver_xres, yres=quiver_yres,
            relres=quiver_relres, xrelres=quiver_xrelres, yrelres=quiver_yrelres,
            scaler = self.xyscaler, compact=True, nauto=quiver_nauto)

        # - apply
        if rmask is not False:
            uu = uu.copy()
            vv = vv.copy()
            mm = mm.copy()
            uu[rmask] = MV2.masked
            vv[rmask] = MV2.masked
            mm[rmask] = MV2.masked
        if quiver_xsamp or quiver_ysamp:
            xslice = slice(None, None, quiver_xsamp)
            yslice = slice(None, None, quiver_ysamp)
            x2d = x2d[yslice, xslice]
            y2d = y2d[yslice, xslice]
            uu = uu[yslice, xslice]
            vv = vv[yslice, xslice]
            mm = mm[yslice, xslice]

        # Norm
        quiver_norm = kwqv.pop('norm', None)
        if quiver_norm in [1, 2] :
            quiverkey = False
            good = mm!=0.
            uu = N.ma.where(good, uu/mm, uu)
            vv = N.ma.where(good, vv/mm, vv)
            del good
            qvargs = [x2d, y2d, uu, vv]
            if quiver_norm == 2:
                qvargs.append(mm)
        else:
            qvargs = [x2d, y2d, uu, vv]
            if quiver_norm == 3:
                qvargs.append(mm)
        if quiver_norm >=2:
            if quiver_cmap is None:
                quiver_cmap = self.get_cmap(**kwcmap)
            kwqv['cmap'] = quiver_cmap
            vmin = quiver_vmin if quiver_vmin is not None else self.vmin
            vmax = quiver_vmax if quiver_vmax is not None else self.vmax
            kwqv['norm'] = Normalize(vmin=vmin, vmax=vmax)
        mask = N.ma.asarray(uu).mask
        if mask is not MV2.nomask:
            good = ~mask.ravel()
            qvargs = [N.compress(good, qa.ravel()) for qa in qvargs] ; del qa

        # Arrows
        if barbs:
            quiver_func = self._plotter.barbs
            quiverkey = False
        else:
            quiver_func = self._plotter.quiver
            angles = 'xy' if self.uvscaler is not None else 'uv'
            kwqv.setdefault('angles', angles)
        qv = quiver_func(*qvargs, **kwqv)
        if isinstance(qv,  tuple): qv = list(qv)
        self.set_obj('quiver', qv)
        self.register_obj(qv, 'patches')
        if quiver_norm>=2 and self.get_scalar_mappable() is None:
            self.set_obj('scalar_mappable', qv)

        # Filters
        if shadow: self.add_shadow(qv, 'quiver_shadow', **kwsh)
        if glow: self.add_glow(qv, 'quiver_glow', **kwgl)

        # Quiver key
        if quiverkey:
            self.quiverkey(qv, **kwqvkey)
        return qv


    def plot_streamplot(self, zorder=None,
        shadow=False, glow=False, streamplot_color=None, streamplot_linewidth=None,
        streamplot_lwmodmin=.5, streamplot_lwmodmax=3,
        streamplot_cmap=None, streamplot_norm=None,
        streamplot_vmin=None, streamplot_vmax=None,
        **kwargs):
        """Plot stream lines with :func:`~matplotlib.pyplot.streamplot`

        :Params:

            - **streamplot_color**,optional:

                - ``None`` or `"default"`: Default colorization.
                - ``"modulus"``: Colorization as a function of the modulus.
                - Else, passed directly to :func:`~matplotlib.pyplot.streamplot`

            - **streamplot_linewidth**,optional:

                - ``None`` or ``"default"``: Default linewidth.
                - ``"modulus"``: Linewidth proportional to the modulus with a maximum
                  linewidth of ``streamplot_lwmod``.
                - Else, passed directly to :func:`~matplotlib.pyplot.streamplot`

            - **streamplot_lwmodmin**, optional: Min linewidth used when ``streamplot_linewidth``
              is set to ``"modulus"``.
            - **streamplot_lwmodmax**, optional: Max linewidth used when ``streamplot_linewidth``
              is set to ``"modulus"``.
            - **streamplot_<param>**, optional: ``<param>`` is passed to
              :func:`~matplotlib.pyplot.streamplot`.
            - **zorder**, optional: Plot order.
            - **alpha**, optional: Opacity.
            - **shadow**, optional: Add shadow below arrows.
            - **glow**, optional: Add glow effect to arrows.
            - **streamplotkey**, optional: Add key to scale arrows.
            - **cmap**, optional: Colormap to use when ``color="modulus"`` (defaults to
              ``"magic"``).


        :Tasks:

            #. Get the scalar data with :meth:`~Plot.get_data`.
            #. Calls :func:`~matplotlib.pyplot.streamplot`.

        """
        if not self.has_valid_data(): return

        # Get data
        if len(self.data)!=3: return
        x2d, y2d = self.x2d, self.y2d
        mm, uu, vv = self.get_data(scalar=False)
        del mm
        mm = N.ma.sqrt(uu**2+vv**2)

        # Keywords
        kwsp = kwfilter(kwargs, 'streamplot')#, defaults=dict(angles='uv'))
        if zorder is not None:
            kwsp['zorder'] = zorder
        if kwsp.has_key('width') and kwsp['width'] > .01:
            kwsp['width'] *= 0.001
        shadow = kwsp.pop('shadow', shadow)
        glow = kwsp.pop('glow', glow)
        kwsh = kwfilter(kwsp, 'shadow')
        kwgl = kwfilter(kwsp, 'glow')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        kwcmap = kwfilter(kwargs, 'cmap', copy=True) if not self.get_scalar_mappable() else {}
        kwcmap.update(kwfilter(kwsp, 'cmap'))

        # Color
        vmax = streamplot_vmax if streamplot_vmax is not None else self.vmax
        vmin = streamplot_vmin if streamplot_vmin is not None else self.vmin
        if streamplot_color in [None, 'default']:
            streamplot_color = None
        elif streamplot_color=='modulus':
            streamplot_color = mm
            if streamplot_cmap is None:
                streamplot_cmap = self.get_cmap(**kwcmap)
            kwsp['cmap'] = streamplot_cmap
            if streamplot_norm is None:
                streamplot_norm  = Normalize(vmin=vmin, vmax=vmax)
            kwsp['norm'] = streamplot_norm
        kwsp['color'] = streamplot_color


        # Linewidth
        streamplot_linewidth = kwsp.pop('lw', streamplot_linewidth)
        streamplot_lwmodmax = kwsp.pop('lwmod', streamplot_lwmodmax)
        if streamplot_linewidth in [None, 'default']:
            streamplot_linewidth = None
        elif streamplot_linewidth=='modulus':
            streamplot_linewidth = (streamplot_lwmodmin+
                (streamplot_lwmodmax-streamplot_lwmodmin)*(mm-vmin)/(vmax-vmin))
        kwsp['linewidth'] = streamplot_linewidth

        # Lines
        sp = self._plotter.streamplot(x2d, y2d, uu, vv, **kwsp)
        self.set_obj('streamplot', sp)
        self.register_obj(sp.arrows, 'patches')
        self.register_obj(sp.lines, 'lines')
        if isinstance(streamplot_color, N.ndarray) and self.get_scalar_mappable() is None:
            sp.arrows.set_array(streamplot_color)
            sp.arrows.set_clim(vmin, vmax)
            sp.arrows.set_norm(kwsp['norm'])
            sp.arrows.set_cmap(kwsp['cmap'])
            self.set_obj('scalar_mappable', sp.arrows)

        return sp

    def plot(self, levels=None, contour=True, anomaly=False, streamplot=False,  **kwargs):
        """Plot filled and/or contoured content

        :Params:

            - **levels**, optional: Levels for contouring and/or colormap.
            - **levels_<param>**, optional: Passed to
              :meth:`~vacumm.misc.core_plot.ScalarMappable.get_levels` to tune levels.
            - **contour**, optional: Plot line contours.
            - **streamplot**, optional: Plot stream lines instead of arrows.
            - Other arguments are passed to :meth:`plot_fill` and :meth:`plot_contour`.

        :Tasks:

            #. Get :attr:`~ScalarMappable.levels`.
            #. Calls :meth:`plot_fill`.
            #. Calls :meth:`plot_contour`.
            #. Calls :meth:`plot_quiver` or :meth:`plot_streamplot`.
        """
        if not self.has_valid_data(): return

        # Get levels
        kwlevels = kwfilter(kwargs, 'levels', defaults=dict(anomaly=anomaly))
        if levels is not None:
            self.levels = levels
        levels = self.get_levels(**kwlevels)

        # Plotter
        if self._plotter is None:
            self._plotter = self.axes

        # Fill
        self.plot_fill(**kwargs)

        # Contours
        if contour: self.plot_contour(**kwargs)

        # Quiver or stream lines
        if not streamplot:
            self.plot_quiver(**kwargs)
        else:
            self.plot_streamplot(**kwargs)

    def add_grid(self, **kwargs):
        """Add the coordinates as a grid with centers and corners

        See :func:`~vacumm.misc.plot.add_grid` for parameters.
        """
        from .plot import add_grid
        kwargs.update(xcorners=self.x2db, ycorners=self.y2db)
        return self.register_obj(add_grid((self.x2d, self.y2d), **kwargs), 'grid')


class Hov(Plot2D):

    _order = ['zT', 'Tx', 'yT', 'T-', '-T']

class Section(Plot2D):

    _order = ['Zx', 'Zy']




class Map(Plot2D):
    """This class is used for plotting a map with or without data

    It uses a :class:`mpl_toolkits.basemap.Basemap` instance, stored in attribute :attr:`map`.
    Some other attributes: :attr:`lon`, :attr:`lat`, :attr:`map_update`.

    .. attribute:: map

        A :class:`~mpl_toolkits.basemap.Basemap` initialized in :meth:`pre_plot`.

    """
    _order = 'YX'
    map = None

    def _set_axes_(self, xaxis=None, yaxis=None, xatts=None, yatts=None, **kwargs):
        Plot2D._set_axes_(self, xaxis=xaxis, yaxis=yaxis, xatts=xatts, yatts=yatts, **kwargs)
        if self.has_data():
            if xaxis is not None:
                for var in self.data:
                    var.getAxis(1).designateLongitude()
            if yaxis is not None:
                for var in self.data:
                    var.getAxis(0).designateLatitude()
    _set_axes_.__doc__ = Plot2D._set_axes_.__doc__

    def format_axes(self, **kwargs):
        pass

    def mapproj(self, x, y, inverse=False):
        """Convert from degrees to meters (or inverse)

        It uses attribute :attr:`map` if set to a
        :class:`~mpl_toolkits.basemap.Basemap` instance with
        a projection different from "cyl".
        If no valid projection is found, et uses functions
        :func:`~vacumm.misc.phys.units.deg2m` and
        :func:`~vacumm.misc.phys.units.m2deg`.

        :Params:

            - **x/y**: Coordinates in degrees or meters
            - **inverse**, optional: Reverse projection (meters->degrees)
            - **current**, optional: If True,
        """
        if getattr(self, 'map', None) is None or self.map.projection=='cyl':
            if inverse:
                return m2deg(x, y), m2deg(x)
            return deg2m(x, y), deg2m(x)
        return self.map(x, y, inverse=inverse)

    def get_xyscaler(self):
        """Get :attr:`xyscaler`"""
        xyscaler = Plot2D.get_xyscaler(self)
        if xyscaler is None:
            if getattr(self, 'map', None) is None:
                return lambda x,y,inverse=True: x,y
            return self.map
        return xyscaler
    xyscaler = property(get_xyscaler, Plot2D.set_xyscaler,
        Plot2D.del_xyscaler, Plot2D.xyscaler.__doc__)

    def get_uvscaler(self, guess=False):
        """Get :attr:`uvscaler`"""
        uvscaler = self.get_obj('uvscaler')
        if uvscaler is not None: return uvscaler
        if guess: return Plot.get_uvscaler(self, guess=True)
        uvscaler = self.uvscaler = lambda u, v: self.map.rotate_vector(u, v, self.x.getValue(), self.y.getValue(), returnxy=False)
        return uvscaler
    uvscaler = property(get_uvscaler, Plot.set_uvscaler,
        Plot.del_uvscaler, Plot.uvscaler.__doc__)

    def load_data(self, data=None, lon=None, lat=None, **kwargs):
        """Data loading


        It performs the following tasks:

            #. Call to the generic :meth:`Plot2D.load_data` method.
            #. Set attributes :attr:`lon` :attr:`lat`.

        :Params:

            - **data**, optional: See :meth:`~vacumm.misc.core_plot.Plot2D.load_data`.
            - **lon**, optional: Longitude interval.
            - **lat**, optional: Latitude interval.

        .. attribute:: lon

            Longitudes as a two-elements tuple or and array.

        .. attribute:: lat

            Latitudes as a two-elements tuple or and array.
        """
        # Default bounds
        default_lon = N.array([-180., 180.])
        default_lat = N.array([-90., 90.])
        if lon is not None:
            lon = N.asarray(lon)
        if lat is not None:
            lat = N.asarray(lat)

        # Get some keys
        self.lon = lon
        self.lat = lat
        self.order = 'yx'

        # Minimal setup
        if data is None:
            if 'xaxis' not in kwargs:
                kwargs['xaxis'] = lon if lon is not None else default_lon
            if 'yaxis' not in kwargs:
                kwargs['yaxis'] = lat if lat is not None else default_lat

        # Basic load
        Plot2D.load_data(self, data, **kwargs)

        # Longitudes and latitudes
        if not self.has_data():
            if self.lon is None:
                self.lon = default_lon
            if self.lat is None:
                self.lat = default_lat
        else:
#            masked = False if self.masked else None
            self.lon = self.get_xdata(masked=True) if self.lon is None else self.lon
            self.lat = self.get_ydata(masked=True) if self.lat is None else self.lat


    def pre_plot(self, map=None, projection='cyl', resolution='auto',
            epsg=None, overlay=False, fullscreen=False, map_update=None,
            lon_min=None, lon_max=None, lat_min=None, lat_max=None,
            lon_center=None, lat_center=None, lat_ts=None,
            nocache=False, cache_dir=None, zoom=None, **kwargs):
        """Plot initialisation.

       :Tasks:

            #. Call to the generic :meth:`Plot.pre_plot` method.
            #. Setup a :class:`~mpl_toolkits.basemap.Basemap` instance
               and store it in :attr:`map`.
            #. Project :attr:`~Plot2D.x2d`, :attr:`~Plot2D.y2d`,
               :attr:`~Plot2D.x2db` and :attr:`~Plot2D.y2db`
               using :attr:`map`.

       :Generic params: See :meth:`Plot.pre_plot`.

       :Special params:

           - **projection**: Map projection, like "merc".
             See :mod:`~mpl_toolkits.basemap.Basemap` for a list of possible projections.
           - **resolution**: GSHHS resolution of shoreline or 's' for Histolitt (SHOM).

            - ``"c"``: Crude.
            - ``"l"``: Low.
            - ``"i"``: Intermediate.
            - ``"h"``: High.
            - ``"f"``: Full.
            - ``"s"``: Histolittfor the french coast
              (from a shapefile file automatically
              downloaded once the license is accepted).

           - **nocache**: Management of cached maps.

             .. only:: html and epub

                   - ``0``: Read and write cached maps.
                   - ``1``: Only write cached maps.
                   - ``2``: Caching is disabled.

             .. only:: latex

                    ``0``: read and write cached maps.;
                    ``1``: Only write cached maps;
                    ``2``: Caching is disabled.

            - **map_update**: Force the update of the map latter by :meth:`post_plot`,
              setting :attr:`map_update`.

            - **zoom**: Zoom on map bounds before creating it.

        .. attribute:: map_update

            An attribute to know if current map must be updated by :meth:`post_plot`.
            It is also a parameter to method.
        """



        # Aliases
        self.map = kwargs.pop('m', map)


        # Common init
        Plot.pre_plot(self, **kwargs)

        # Check if map exists
        if isinstance(self.map, Plot):
            self.map = getattr(self.map, 'map', None)
        if self.map is None and self.get_axobj('map'):
            self.map = self.get_axobj('map')
        if self.map is not None:
            if self.map.ax is not self.axes:
                self.map.ax = self.axes
                map_update = 2
            elif map_update is None:
                map_update = 0
        if map_update is not None:
            self.map_update = int(map_update)
        else:
            self.map_update  = 2


        # Create the map
        if self.map is None:

            kwmap = kwfilter(kwargs, 'basemap')
            kwmap.update(kwfilter(kwargs, 'map_'))
            lon = self.lon if self.xmasked else N.ma.getdata(self.lon)
            lat = self.lat if self.ymasked else N.ma.getdata(self.lat)
            if lon_min is None:
                lon_min = lon.min()
            if lon_max is None:
                lon_max = lon.max()
            if lat_min is None:
                lat_min = lat.min()
            if lat_max is None:
                lat_max = lat.max()
            projection = kwargs.pop('proj', projection)
            from vacumm.misc.grid.basemap import create_map
            resolution = kwargs.pop('res', resolution)

            self.map = create_map(lon_min, lon_max, lat_min, lat_max, projection=projection,
                resolution=resolution, overlay=overlay, fullscreen=fullscreen,
                lon_center=lon_center, lat_center=lat_center, lat_ts=lat_ts, epsg=epsg,
                nocache=nocache, cache_dir=cache_dir, zoom=zoom, ax=self.axes, **kwmap)

        # Projection of coordinates
        if self.has_data():
            self.x2d, self.y2d = self.xyscaler(self.x2d, self.y2d)
            self.x2db, self.y2db = self.xyscaler(self.x2db, self.y2db)

        if not hasattr(self.map, 'res'): self.map.res = self.map.resolution
        self.set_axobj('map', self.map)
        #self._plotter = self.map
        self._plotter = self.axes

    def _check_map_(self):
        if self.map is None: raise PlotError('Map still no set')

    def __call__(self, lon, lat, inverse=False):
        """Convert coordinates with a geographic projection
            using the current :meth:`~mpl_toolkits.basemap.Basemap` instance (:attr:`map`).

        See :meth:`mpl_toolkits.basemap.Basemap.__call__`
        """
        self._check_map_()
        return self.map(lon, lat, inverse=inverse)


    def plot(self, **kwargs):
        """Main plot

        It performs the following tasks:

            #. Call to generic :meth:`~Plot2D.plot` method.
        """



        # Generic 2D plot
        Plot2D.plot(self, **kwargs)

        # Update map limits
        self.map.set_axes_limits(ax=self.axes)


    def plot_quiver(self, quiver_res_scaler=None, quiver_angles='uv', **kwargs):
        if quiver_res_scaler is None:
            quiver_res_scaler = self.mapproj
        return Plot2D.plot_quiver(self, quiver_res_scaler=quiver_res_scaler, **kwargs)
    plot_quiver.__doc__ = Plot2D.plot_quiver.__doc__



    def add_lowhighs(self, size=40, weight='normal', lowtext='L', hightext='H',
        shadow=False, glow=False, smooth=False, va='center', ha='center', **kwargs):
        """Mark position of mins and maxs using letters.

        It is typically used for adding L and H to depressions and anticyclones.
        """

        if not self.has_valid_data(): return

        # Data
        data = self.get_data(scalar=True)
        if data.shape[0]<=2 or data.shape[1]<=2: return

        # Smooth
        if smooth:
            data[:] = generic2d(data, smooth)

        # Gradients
        ndata = data.filled()
        dx = N.sign(N.diff(data[1:-1], axis=1))
        dy = N.sign(N.diff(data[:, 1:-1], axis=0))

        # High and low fields
        lows =  dx[:, :-1]<0
        lows &= dx[:, 1:]>0
        lows &= dy[:-1]<0
        lows &= dy[1:]>0
        highs =  dx[:, :-1]>0
        highs &= dx[:, 1:]<0
        highs &= dy[:-1]>0
        highs &= dy[1:]<0

        # Masking
        if self.mask.any():
            mask = self.mask[:-1] | self.mask[1:]
            mask |= self.mask[:, -1] | self.mask[:, 1:]
            for hl in lows, highs:
                var[mask] = False

        # Positions
        xlows = self.x2d[1:-1, 1:-1][lows]
        ylows = self.y2d[1:-1, 1:-1][lows]
        xhighs = self.x2d[1:-1, 1:-1][highs]
        yhighs = self.y2d[1:-1, 1:-1][highs]

        # Plot text
        kwsh = kwfilter(kwargs, 'shadow')
        kwgl = kwfilter(kwargs, 'glow')
        dict_copy_items(kwargs, [kwsh, kwgl], 'anim')
        for key in 'size', 'va', 'ha', 'weight':
            kwargs.setdefault(key, eval(key))
        tts = []
        for (xx, yy), text in (((xlows, ylows), lowtext), ((xhighs, yhighs), hightext)):
            for x, y in zip(xx, yy):
                tt = self.register_obj(self.axes.text(x, y, text, **kwargs), 'lowhighs', **kwargs)
                tts.append(tt)
                if shadow: self.add_shadow(tt, 'lowhighs_shadow', **kwsh)
                if glow: self.add_glow(tt, 'lowhighs_glow', **kwgl)

        del xlows, ylows, xhighs, yhighs
        return tts


    def add_arcgisimage(self, service, **kwargs):
        """Add an Arcgis image

        Available service aliases (see :attr:`ARCGISIMAGE_ALIASES`): {}
        """

        # Alias
        if not service:
            return
        if service is True:
            service = "esriimagery"
        if isinstance(service, basestring):
            service = ARCGISIMAGE_ALIASES.get(service, service)

        # Pixels
        axext = self.axes.bbox.extents
        fact = 1.
        xpixels = (axext[2]-axext[0])*fact
#        ypixels = (axext[3]-axext[1])*fact
        dpi = self.fig.dpi

        # Call the basemap method
        dict_check_defaults(kwargs, service=service,
            xpixels=xpixels,
#            ypixels=ypixels,
            dpi=dpi)
        try:
            img = self.map.arcgisimage(**kwargs)
            img.set_zorder(kwargs.get('zorder', 1))
            return img
        except Exception, e:
            warn('Error when plotting arcgisimage: {}.\nMessage: {}'.format(
                service, e.message))

    add_arcgisimage.__doc__.format(', '.join(ARCGISIMAGE_ALIASES.keys()))

    def _get_posposref_(self, pos=None, posref=None, xrel=0.1, yrel=0.1, transform='axes',
        xpad=30, ypad=None):
        """Get position of point and reference points in data coordinates (meters)

        :Params:

            - **pos**: Position in data (lon, lat), axes or figure coordinates for placement.
            - **posref**: Position in data (lon, lat) coordinates of reference point.
            - **x/yrel**: Default placement position relative to longitude and latitude ranges.
            - **transform**: Coordinates transform for ``pos`` which defaults to "data".
              Use a valid transform, or "data", "axes" or "figure".
            - **x/ypad**: Padding if dots for placement when given as a string like "upper left".

        Return: ``(x,y),(xref,yref)`` in meters
        """
        # Placement point
        offset = None
        if xrel<0: xrel += 1
        if yrel<0: yrel += 1
        if pos is None:
            lonmin = self.map.lonmin
            lonmax = self.map.lonmax
            latmin = self.map.latmin
            latmax = self.map.latmax
            dlon = lonmax-lonmin
            dlat = latmax-latmin
            lon0 = lonmin+xrel*dlon
            lat0 = latmin+yrel*dlat
            pos = tuple(self(lon0, lat0))
        else:
            if pos in _locations:
                if xpad is not False and ypad is not False:
                    offset = loc2offset(pos, xpad=xpad, ypad=ypad, margin=0.01)
                    pos = loc2tuple(pos, xpad=0, ypad=0)
                else:
                    pos = loc2tuple(pos, xpad=xrel, ypad=yrel)

                transform = self.axes.transAxes
            else:
                transform = self._transform_(transform)
            if transform is self.axes.transAxes or transform is self.fig.transFigure:
                pos = tuple((transform-self.axes.transData).transform_point(pos))
            else:
                pos = tuple(self(*pos))

        # Reference point
        if posref is None:
            posref = pos
        else:
            posref = self(*posref)

        return pos, posref, offset

    def add_mapscale(self, pos='lower left', scale=None, posref=None, barstyle='simple', transform=None,
        xrel=0.1, yrel=0.1, getpos=False, posonly=False, shadow=False, zorder=10, **kwargs):
        """Add a map scale using :meth:`mpl_toolkits.basemap.Basemap.drawmapscale`

        :Params:

            - **pos**, optional: ``(lon,lat)`` where to draw it.
            - **posref**, optional: ``(lon,lat)`` of reference position where the scale
              is estimated.
            - **scale**, optional: Length of the scale in projection coordinates (m).
            - **barstyle**, optional: Bar style: 'simple' or 'fancy'.
            - **transform**, optional: Coordinates transform for ``pos`` which defaults to "data".
              Use a valid transform, or "data", "axes" or "figure".
            - **x/yrel**, optional: Default placement position relative to longitude and latitude ranges.
            - **x/ypad**, optional: Padding if dots for placement when given as a string like "upper left".
            - **shadow**, optional: Add a shadow.
            - **shadow_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.Plot.add_shadow`.
            - Othe keywords are passed to :meth:`~mpl_toolkits.basemap.Basemap.drawmapscale`.

        :Example:

            >>> mymap.add_mapscale((.2,.2), transform='axes')
            >>> mymap.add_mapscale('upper right', posref=(-2,45), fontsize=10)
        """
        if self.map is None: return
        kwsh = kwfilter(kwargs, 'shadow_')

        # Positions
        pos, posref, offset = self._get_posposref_(pos, posref, transform=transform,
            xrel=xrel if xrel<.5 else (1-xrel), yrel=yrel if yrel<.5 else (1-yrel), xpad=False)
        pos = self(pos[0], pos[1], inverse=True)
        if getpos and posonly: return pos
        posref = self(posref[0], posref[1], inverse=True)

        # Scale
        if scale is None:
            scales = basic_auto_scale(self.map.xmin, self.map.xmax, nmax=10)
            scale = (scales[1]-scales[0])*0.001

        # Draw it
        kwargs['barstyle'] = barstyle
        ms = self.map.drawmapscale(pos[0], pos[1], posref[0], posref[1], scale, **kwargs)
        self.add_obj('mapscale', ms)

        # Finalize
        for o in ms:
            o.set_clip_on(False)
            o.set_zorder(zorder)
            if shadow:
                self.add_shadow(o, **kwsh)

        if getpos: return ms, pos
        return ms


    def add_compass(self, pos='lower right', size=40, posref=None, style='simple', transform=None,
        xpad=None, ypad=None, xrel=0.9, yrel=0.9, getpos=False, shadow=False, zorder=10, **kwargs):
        """Add a compass to the map using :func:`~vacumm.misc.plot.add_compass`

        See :func:`~vacumm.misc.plot.add_compass` all options.

        :Params:

            - **pos**, optional: ``(lon,lat)`` where to draw it.
            - **size**, optional: Size in pixels.
            - **posref**, optional: ``(lon,lat)`` of reference position where the north
              is estimated.
            - **style**, optional: Compas style: 'simple' or 'fancy'.
            - **transform**: Coordinates transform for ``pos`` which defaults to "data".
              Use a valid transform, or "data", "axes" or "figure".
            - **x/yrel**: Default placement position relative to longitude and latitude ranges.
            - **x/ypad**: Padding in dots for placement when given as a string like "upper left".
            - **shadow**, optional: Add a shadow.
            - **shadow_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.Plot.add_shadow`.
            - Othe keywords are passed to :func:`~vacumm.misc.plot.add_compass`.

        :Example:

            >>> map2(data, compass=True)
            >>> mymap.add_compass((-4,45), 40, facecolor1='r', text_weight='bold')
            >>> mymap.add_compass(pos="upper right", style="fancy", xpad=100)
        """
        if self.map is None: return
        kwsh = kwfilter(kwargs, 'shadow_')

        # Transform
        transform = self._transform_(transform, 'axes')

        # Positions
        if xpad is None: xpad = int(size/2)+20
        pos, posref, poffset = self._get_posposref_(pos, posref, transform=transform,
            xrel=xrel, yrel=yrel, xpad=xpad, ypad=ypad)

        # Departure angle from 90 for north
        pr = self(posref[0], posref[1], inverse=True)
        uo, vo = self.map.rotate_vector(N.array([0]),N.array([1.]),
            N.array(pr[:1]),N.array(pr[1:]))
        angle = N.degrees(N.arctan2(vo[0], uo[0]))-90.

        # Draw it
        x, y = pos
        kwargs['posref'] = posref
        dict_check_defaults(kwargs, offset=poffset, text_zorder=zorder)
        kwfont = kwfilter(kwargs, 'font', short=True)
        if 'size' in kwfont: kwargs['text_size'] = kwfont['size']
        if 'color' in kwfont: kwargs['text_color'] = kwfont['color']
        cp = add_compass(x, y, size=size, ax=self.axes, angle=angle, anglemode='data',
            transform=self.axes.transData, zorder=zorder, **kwargs)
        self.add_obj('compass', cp)

        # Finalize
        if shadow:
            for o in cp[0]+cp[1]:
                self.add_shadow(o, **kwsh)

        if getpos: return cp, pos
        return cp


    def add_mscp(self, pos=None, posref=None, compass_size=40, transform=None,
        shadow=False, color='k', zorder=10, **kwargs):
        """Plot a mapscale and a compass above


        :Params:

            - **pos**: Position of the mapscale (see :meth:`add_mapscale`).
              The compass is drawn just above.
            - **posref**, optional: Position of reference point for both mapscale and
              compass.
            - **mapscale_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.Map.add_mapscale`.
            - **scompass_<param>**, optional: ``<param>`` is passed to
              :meth:`~vacumm.misc.core_plot.Map.add_compass`.


        Example:

            >>> map2(data, mscp=True, mscp_pos=(-2, 48))
            >>> mymap.add_mscp('lower left')
            >>> mymap.add_mscp((-4,45), mapscale_scale=100, mapscale_barstyle='fancy',
                compass_size=50, compass_style='simple')
            >>> mymap.add_mscp((.9,.1), transform='axes')
        """

        kwms = kwfilter(kwargs, 'mapscale')
        kwms['transform'] = transform
        kwcp = kwfilter(kwargs, 'compass')
        dict_check_defaults(kwcp, facecolor1=color, edgecolor=color, text_color=color,
            zorder=zorder)
        dict_check_defaults(kwms, fontcolor=color, fillcolor2=color,
            zorder=zorder)

        msoffset = kwms.get('offset', 0.02*(self.map.ymax-self.map.ymin))
        cptextpad = kwcp.get('text_pad', 5)
        cpoffset = (0, int(compass_size)/2)

        # Mapscale position at top
        if pos in _locations and 'top' in loc2align(pos)['va']:

            pos = self.add_mapscale(pos=pos, posref=posref, getpos=True, posonly=True, **kwms)
            mpos = self(*pos)
            ppos = self.axes.transData.transform_point(mpos)
            mpos = self.axes.transData.inverted().transform_point(
                (ppos[0], ppos[1]-cpoffset[1]-cptextpad*2))
            pos = self(mpos[0], mpos[1]-3*msoffset, inverse=True)
            del kwms['transform']

        # Mapscale
        kwms.setdefault('pos', pos)
        ms, pos = self.add_mapscale(posref=posref, getpos=True, **kwms)

        # Compass above mapscale
        pos = self(*pos)
        pos = self(pos[0], pos[1]+msoffset*3, inverse=True)
        kwcp.update(transform = 'data', offset=cpoffset)
        kwcp.setdefault('pos', pos)
        cp = self.add_compass(posref=posref, size=compass_size, **kwcp)

        return ms+list(cp)

    def post_plot(self, drawrivers=False, fillcontinents=True,
            meridional_labels=True, zonal_labels=True,
            drawcoastlines=True, drawmapboundary=True,
            meridians=None, parallels=None,
            land_color=None, ticklabel_size=None, refine=0,
            no_seconds=False, fullscreen=False,
            minutes=True, mapscale=False, compass=False, mscp=False,
            bfdeg=None, lowhighs=False, arcgisimage=None,
            **kwargs):
        """Post-processing of the plot

        :Tasks:

            #. Update the map if attribute :attr:`map_update` is True:

                #) Treat continents if attribute :attr:`map_update` = 2 and
                   attribute :attr:`map` has a coastline resolution not ``None``:
                   draw coastlines, fill continents and draw rivers.

                #) Draw parallels and associated labels.

            #. Call to generic post processing of 2D plots (:meth:`Plot2D.post_plot`).

        :Params:

            - **drawrivers**: Draw rivers on the map using method
              :meth:`~mpl_toolkits.basemap.Basemap.drawrivers`.
            - **drawrivers_<param>**: Pass ``<param>`` to the method.
            - **fillcontinents**: Fill continents with color ``land_color``
              using method :meth:`~mpl_toolkits.basemap.Basemap.fillcontinents`.
            - **land_color**: Fill color.
            - **fillcontinents_<param>**: Pass ``<param>`` to the method.
            - **drawparallels**: Display or hide parallels and associated labels using
              method :meth:`~mpl_toolkits.basemap.Basemap.drawparallels`.
            - **drawparallels_<param>**: Pass ``<param>`` to the method.
            - **drawparallels_fmt**: Default to :class:`MinuteLabel`.
            - **drawparallels_gs_<param>**: Passed to :func:`~vacumm.misc.misc.geo_scale`.
            - **parallels**: Parallels to plot.
            - **drawmeridians**: Display or hide medidians and associated labels using
              method :meth:`~mpl_toolkits.basemap.Basemap.drawmeridians`.
            - **drawmeridians_<param>**: Pass ``<param>`` to the method.
            - **drawmeridians_fmt**: Default to :class:`MinuteLabel`.
            - **drawmeridians_gs_<param>**: Passed to :func:`~vacumm.misc.misc.geo_scale`.
            - **meridians**: Meridians to plot.
            - **meridional/zonal_labels**: Display or hide meridional/zonal labels.
              ``meridional/zonal_labels=False`` is equivalent to ``y/xhide=True``.
            - **no_seconds**: Do not display seconds in labels (if applicable).
            - **minutes**: Do not use decimal degrees for labels (if applicable).
            - **bfdeg**: Degrees are in bold (latex text only, if applicable).
            - **x/y/ticklabels_<param>**: Pass ``<param>`` to
              :meth:`~mpl_toolkits.basemap.Basemap.drawmeridians` and
              :meth:`~mpl_toolkits.basemap.Basemap.drawparallels` to change text properties.
            - **fullscreen**: Full screen mode -> no colorbar and no labels.
            - **mapscale**: Add a map scale using
              :meth:`~vacumm.misc.core_plot.Map.add_mapscale`.
            - **mapscale_<param>**: Pass ``<param>`` to the
              :meth:`~vacumm.misc.core_plot.Map.add_mapscale` method.
            - **compass**: Add a compass using
              :meth:`~vacumm.misc.core_plot.Map.add_compass`.
            - **compass_<param>**: Pass ``<param>`` to the
              :meth:`~vacumm.misc.core_plot.Map.add_compass` method.
            - **mscp**: Add a mapscale AND a compass using
              :meth:`~vacumm.misc.core_plot.Map.add_mscp`.
            - **mscp_<param>**: Pass ``<param>`` to the
              :meth:`~vacumm.misc.core_plot.Map.add_mscp` method.


        """



        if self.map_update:

            # Full screen
            if fullscreen:
                kwargs.setdefault('colorbar',  False)
                default_drawparallels = default_drawmeridians = False
            else:
                default_drawparallels = default_drawmeridians = True
            drawparallels = kwargs.pop('drawparallels', default_drawparallels)
            drawmeridians = kwargs.pop('drawmeridians', default_drawmeridians)


            # Map boundaries
            if drawmapboundary and not self.is3d:
                try:
                    self.set_axobj('drawmapboundary', self.map.drawmapboundary(**kwfilter(kwargs,'drawmapboundary')))
                except:
                    print 'Problem with Basemap.drawmapboundary'

            # Continents
            if self.map.res is not None and self.map_update>=2:

                if land_color is None: land_color = land
                kwdrawcoastlines = kwfilter(kwargs, 'drawcoastlines_')
                kwdrawcoastlines.setdefault('linewidth', 0.2)
                kwfillcontinents = kwfilter(kwargs, 'fillcontinents_', defaults={'color':land_color})

                if self.map.resolution is not None:# GSHHS

                    # Draw coastlines
                    if drawcoastlines and not self.get_axobj('drawcoastlines'):
                        o = self.map.drawcoastlines(**kwdrawcoastlines)
                        self.set_axobj('drawcoastlines', o)
                        if self.is3d:
                            self.axes.add_collection3d(o)

                    # Fill continents
                    fcsh = kwfillcontinents.pop('shadow', False)
                    kwfcsh = kwfilter(kwfillcontinents, 'shadow_')
                    if fillcontinents and not self.get_axobj('fillcontinents'):
                        if not self.is3d:
                            try:
                                self.set_axobj('fillcontinents', self.map.fillcontinents(**kwfillcontinents))
                            except:
                                pass
                            if self.get_axobj('fillcontinents') is not None and fcsh:
                                add_shadow(self.get_axobj('fillcontinents'), **kwfcsh)
                        else:
                            polys = []
                            for polygon in self.map.landpolygons:
                                polys.append(polygon.get_coords())
                            kwfillcontinents['facecolor'] = kwfillcontinents.pop('color')
                            kwfillcontinents.pop('lake_color', None)
                            lc = PolyCollection(polys, linewidth=0, closed=False,
                                **kwfillcontinents)
                            del polys
                            self.axes.add_collection3d(lc)


                    # Draw rivers
                    if drawrivers and not self.get_axobj('drawrivers'):
                        o = self.map.drawrivers(**kwfilter(kwargs,'drawrivers'))
                        self.set_axobj('drawrivers', o)
                        if self.is3d:
                            self.axes.add_collection3d(o)

                #elif m.res == 's': # SHOM
                else:
                    from io import Shapes
                    if self.map.res=='s':
                        from vacumm.bathy.shorelines import Histolitt
                        shoreline = Histolitt(m=self.map, proj=True)
                    else:
                        shoreline = self.map.res
                    if isinstance(shoreline, Shapes):
                        try:
                            kwshoreline = {}
                            if not drawcoastlines:
                                kwshoreline['linewidth'] = 0
                            else:
                                kwshoreline.update(kwdrawcoastlines)
                            if not fillcontinents:
                                kwshoreline['fill'] = False
                            else:
                                kwshoreline['fillcolor'] = kwfillcontinents.pop('color')
                                kwshoreline.update(kwfillcontinents)
                            kwshoreline.update(show = False, axes=self.axes)
                            shoreline.plot(**kwshoreline)
                        except Exception, e:
                            import traceback
                            print '*** Error in shaprefile shoreline'
                            print traceback.format_exc()

            # Tick labels and lines
            kwpm_def = dict(linewidth=0.5)
            kwpm_def.update(kwfilter(kwargs, 'ticklabels_'))
            if ticklabel_size is not None: kwpm_def.update(fontsize=ticklabel_size)
            kwp = kwpm_def.copy() ; kwm = kwpm_def.copy()
            if kwargs.get('grid', True) is False:
                kwp['linewidth'] = kwm['linewidth'] = 0
            # - parallels
            if drawparallels:
                if self._xyhide_('y', kwargs.get('yhide', False)) or self.is3d:
                    meridional_labels = 0
                kwp.setdefault('labels', [int(meridional_labels),0,0,0])
#                kwp.update(kwpm_def)
                kwp.update(kwfilter(kwargs, 'yticklabels_'))
                kwp = kwfilter(kwargs,'drawparallels',defaults=kwp)
                kwgs = kwfilter(kwp,'gs',
                    defaults={'vmin':self.map.llcrnrlat, 'vmax':self.map.urcrnrlat,
                        'minutes':minutes})
                if parallels is None:
                    parallels = geo_scale(**kwgs)
                parallels = filter(lambda _:
                    _>= self.map.llcrnrlat and _<= self.map.urcrnrlat, parallels)
                if minutes: kwp.setdefault('fmt',
                    MinuteLabel(parallels, zonal=False, tex=None,
                        no_seconds=no_seconds, bfdeg=bfdeg))
#                self.set_axobj('drawparallels', self.map.drawparallels(parallels,**kwp))
                self.parallels = parallels
            else:
                self.parallels = None
            # - meridians
            if drawmeridians:
                if self._xyhide_('x', kwargs.get('xhide', False)) or self.is3d:
                    zonal_labels = 0
                kwm.setdefault('labels', [0,0,0,int(zonal_labels)])
#                kwm.update(kwpm_def)
                kwm.update(kwfilter(kwargs, 'xticklabels_'))
                kwm = kwfilter(kwargs,'drawmeridians',defaults=kwm)
                kwgs = kwfilter(kwm,'gs',
                    defaults={'vmin':self.map.llcrnrlon, 'vmax':self.map.urcrnrlon,
                        'minutes':minutes})
                if meridians is None:
                    meridians = geo_scale(**kwgs)
                    if (meridians[-1]-meridians[0])>=360.: # to prevent overlaping
                        meridians = meridians[meridians<meridians[0]+360.]
                crnlon = self.map.llcrnrlon,
                meridians = filter(lambda _:
                    _>= self.map.llcrnrlon and _<= self.map.urcrnrlon, meridians)
                if minutes:
                    kwm.setdefault('fmt',
                        MinuteLabel(meridians, zonal=True,  tex=None,
                            no_seconds=no_seconds, bfdeg=bfdeg))
                self.meridians = meridians
            else:
                self.meridians = None

            if self.is3d:
                def addto3d(what, fmt='%g', labstyle=None, **kwargs):
                    """Handle lines ans labels for 3D plots"""
                    coords = []
                    for value, par in self.get_axobj(what).items():
                        for line in par[0]:
                            xy = line.get_xydata()
                            coords.append(xy)
                            x = xy[0, 0]
                            y = xy[0, 1]
                            if what=='drawmeridians':
                                lab = _setlatlab(fmt, value, labstyle)
                                zdir = 'y'
                                ha = 'left'
                            else:
                                lab = _setlonlab(fmt, value, labstyle)
                                zdir = 'x'
                                ha = 'right'
                            self.axes.text(x, y, 0, lab, ha=ha, va='center',
                                           zdir=zdir)
                    if coords:
                        self.axes.add_collection3d(LineCollection(coords,
                            linewidths=[line.get_linewidth()],
                            linestyles=[line.get_linestyle()],
                            colors=[line.get_color()],
                            ))

            # - bfdeg (bold face degrees) homogeneisation
            if minutes:
                if (drawparallels and isinstance(kwp['fmt'], MinuteLabel) and
                    drawmeridians and isinstance(kwm['fmt'], MinuteLabel) and
                    bfdeg is None and
                    (kwp['fmt'].kwargs['bfdeg'] or kwm['fmt'].kwargs['bfdeg'])):
                    kwp['fmt'].kwargs['bfdeg'] = kwm['fmt'].kwargs['bfdeg'] = True
            # - draw
            if drawparallels:
                self.set_axobj('drawparallels', self.map.drawparallels(parallels,**kwp))
                if self.is3d:
                    addto3d('drawparallels', **kwp)
            if drawmeridians:
                # Draw
                lonlabs = self.set_axobj('drawmeridians', self.map.drawmeridians(meridians,**kwm))
                if self.is3d:
                    addto3d('drawmeridians', **kwm)
                # Remove duplicated labels
                llfound = []
                for ll in lonlabs.keys()[:]:
                    llm = ll%360.
                    istart = 1-int(llm in llfound)
                    for lab in lonlabs[ll][1][1:]:
                        lab.set_visible(False)
                    llfound.append(llm)

        # Plot lows and highs
        lowhighs = kwargs.pop('lowhigh', lowhighs)
        if lowhighs:
            self.add_lowhighs(**kwfilter(kwargs, 'lowhighs'))

        # Plot arcgis image
        arcgisimage = kwargs.pop('arcgisimage', arcgisimage)
        if arcgisimage:
            self.add_arcgisimage(arcgisimage, **kwfilter(kwargs, 'arcgisimage'))

        # Map scale and compass
        if mscp:
            kw = kwfilter(kwargs, 'mscp')
            kw.update(kwfilter(kwargs, 'mapscale_', keep=True))
            kw.update(kwfilter(kwargs, 'compass_', keep=True))
            self.add_mscp(**kw)

        else:

            # Map scale
            if mapscale:
                self.add_mapscale(**kwfilter(kwargs, 'mapscale_'))

            # Compass
            if compass:
                self.add_compass(**kwfilter(kwargs, 'compass_'))

        if fullscreen: self.axes.set_frame_on(False)

        # Basic
        Plot2D.post_plot(self, **kwargs)


    def get_best_loc(self, onland=True, **kwargs):
        """Best location on the plot for an object according to land/sea mask"""
        return best_loc_map(self, onland=onland, **kwargs)



############################################################
## Utilities

def get_axis_scale(axis, type=None):
    """Get an axis scale relative standard units

    It searches for the units of the axis,
    and guess conversion factor.

    :Params:

        - **axis**: A :mod:`cdms2` axis.
        - **type**, optional: Its type.
          If ``None`` it is guessed.
          See :func:`~vacumm.misc.axes.axis_type`.

    :Return: A scalar which defaults to ``1.``
    """
    if type is None:
        type = axis_type(axis).lower()
    units = getattr(axis, 'units', '')
    if type=='-' or not units: return 1.
    if type in ['x', 'y']:
        units
        #NOT FINISHED

def get_quiverkey_value(data, mode=80):
    """Get a decent value for a quiver key"""
    # Scalar
    if N.isscalar(data):
        return data

    # From components
    if isinstance(data, tuple):
        data = N.ma.sqrt(data[0]**2+data[1]**2)

    # Reference value
    if mode=='max':
        vref = N.ma.max(data)
    elif mode=='mean':
        vref = N.ma.mean(data)
    else:
        if mode=='median':
            mode = .5
        elif not isinstance(mode, (int, float)):
            raise PlotError('Unkown quiverkey value mode: {}'.format(mode))
        vref = N.percentile(data, mode)
    del data

    # Quiverkey value
    v10 = N.ma.log10(vref)
    if ((v10+1) % 1) > (N.log10(.5) % 1):
        v10 = N.ma.ceil(v10)
    else:
        v10 = N.ma.floor(v10)
    return 10.**v10


############################################################
## Locators
############################################################



class AutoDateLocator2(AutoDateLocator):
    """A clever :class:`matplotlib.dates.AutoDateLocator`"""
    def __init__(self, *args,  **kwargs):
        noweek = kwargs.pop('noweek', 2)
        try:
            AutoDateLocator.__init__(self, *args, **kwargs)
        except:
            kwargs.pop('maxticks', None)
            AutoDateLocator.__init__(self, *args, **kwargs)
        self._oldlocator = None
        self._oldlims = None
        if noweek:
            self.intervald[DAILY] = [1, 2, 3, 7, 14, 21][:2+int(noweek)]
            self.minticks = 3

    def get_locator(self, *args, **kwargs):

        # Check cache
        if self._oldlocator is not None and \
            self._oldlims == self.axis.get_view_interval().tolist()+self.axis.get_data_interval().tolist():
            return self._oldlocator

        # Get normal locator
        locator = AutoDateLocator.get_locator(self, *args, **kwargs)

        # Rectifications
        try: # sampling interval
            sampling = locator.rule._rrule._interval
        except:
            sampling = locator.base._base
        update = False
        if self._freq == DAILY and sampling%7 == 0:
            # Fix daily/7 to weekly/monday locator (start on monday)
            locator.rule.set(byweekday=MO)
            locator.rule.set(interval=1)
            update = True
        elif self._freq == HOURLY and sampling > 1:
            # Fix hourly locator to start at midnight
            good = [2, 3, 4, 6, 8, 12]
            sampling = good[min(N.searchsorted(good, sampling), len(good)-1)]
            for i in xrange(sampling):
                try:
                    locator.rule.set(byhour=range(i, 24, sampling))
                except:
                    pass
            locator.rule.set(interval=1)
            update = True
        if update:
            locator.set_axis(self.axis)
            locator.set_view_interval(*self.axis.get_view_interval())
            locator.set_data_interval(*self.axis.get_data_interval())
        self._oldlocator = locator
        self._oldlims = self.axis.get_view_interval().tolist()+self.axis.get_data_interval().tolist()
        return locator

class AutoDateMinorLocator(AutoDateLocator2):
    """An extension to the :class:`AutoDateLocator2` for minor locators"""
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('minticks', 3)
        kwargs.setdefault('maxticks', 11)
        AutoDateLocator2.__init__(self, *args, **kwargs)
    def viewlim_to_dt(self):
        major_ticks = self.axis.get_majorticklocs()[:2]
        return num2date(major_ticks[0], self.tz), num2date(major_ticks[-1], self.tz)
    def datalim_to_dt(self):
        return self.get_view_interval()


############################################################
## Formatters
############################################################

class AutoDateFormatter2(AutoDateFormatter):
    scaled = {
           365.0  : '%Y',
           30.    : '%b %Y',
           1.0    : '%d %b',
           1./24. : '%Hh',
           1./(60*24.) : '%Hh%M',
           1./(60*24.*60) : "%M'%Ss''",
           }
    def __init__(self, *args, **kwargs):
        AutoDateFormatter.__init__(self, *args, **kwargs)
        self.scaled = AutoDateFormatter2.scaled

class DualDateFormatter(DateFormatter):
    """Special formatter for dates

    Example: ['00h 01/10/2000' '02h' .... '23h'  '00h 01/10/2000'] to mark days and keep hours.
    Here the 'phase' is 0 (00h) and the level is 3 (='day').

    .. todo:: DualDateFormatter: verify if dual_fmt is really used
    """

    def __init__(self, level, fmt=None, dual_fmt=None, phase=None,
            locator=None, **kwargs):
        # Which level?
        slevels = ['year', 'month', 'week', 'day', 'hour', 'minute']
        if not isNumberType(level):
            level = level.strip().lower()
            if level.endswith('s'):
                level = level[:-1]
            if level not in slevels:
                level = -1
            else:
                level = slevels.index(level)
        self.level = level
        self.locator = locator

        # Autoformat
        best_phases = None
        if self.level == 0: # Year
            if fmt is None: fmt = '%b'
            if dual_fmt is None: dual_fmt = fmt+'%n%Y'
        elif self.level == 1: # Month
            if fmt is None: fmt = '%e'
            if dual_fmt is None: dual_fmt = fmt+'%n%b %Y'
            bast_phases = [0,  15]
        elif self.level == 2: # Week
            if fmt is None: fmt = '%a'
            if dual_fmt is None: dual_fmt = fmt+'%n%d/%m/%Y'
            # BAD PHASE FOR WEEKS! should be monday
            self.level = 2
        elif self.level == 3: # Day
            if fmt is None: fmt = '%Hh'
            if dual_fmt is None: dual_fmt = fmt+'%n%d/%m/%Y'
            self.level -=1
            best_phases = [0, 12]
        elif self.level == 4: # Hour
            if fmt is None: fmt = "%M'"
            if dual_fmt is None: dual_fmt = '%Hh'
            self.level -= 1
            best_phases = [0, 30]
        elif self.level == 5: # Minute
            if fmt is None: fmt = "%S''"
            if dual_fmt is None: dual_fmt = "%M'"
            self.level -= 1
            best_phases = [0, 30]
        else:
            if self.level == -2:
                if fmt is None: fmt = '%Y'
                if dual_fmt is None: dual_fmt = fmt
            elif self.level == -3 or self.level>5:
                if fmt is None: fmt = "%M'%Ss''"
                if dual_fmt is None: dual_fmt = fmt
                self.level = -3
            else:
                if fmt is None: fmt = '%Y-%d-%m %H:%M'
                if dual_fmt is None: dual_fmt = fmt
                self.level = -1
        DateFormatter.__init__(self, fmt, **kwargs)
        self.fmt = fmt
        self.dual_fmt = dual_fmt
        self._reqphase = phase
        if self.level<0: return
        if callable(phase):
            self._phase_check = self._phase = phase
        else:
            np = 5-self.level-1
            phases = [1, 1, 0, 0, 0] # m,d,H,M,S
            if isinstance(phase, int):
                phases[self.level] = phase
            phases = tuple(phases)
            if isinstance(phase, tuple):
                if len(phase)<np: # tail is missing
                    phase += tuple(phases[self.level+len(phase):])
                else: # too much head
                    phase = phase[-np:]
            else:
                phase = phases[self.level:]
            self._phase = phase
            self._phase_check_ = lambda dt: (
                dt.timetuple()[self.level+1:5]+(int(round(dt.second)), ) ==
                    tuple(self._phase))
            if locator: # check real case using locator
                dates = num2date(locator())
                if not any([self._phase_check_(date) for date in dates]):
                    self._phase = list(self._phase)
                    for bp in best_phases or []:
                        self._phase[0] = bp
                        if any([self._phase_check_(date) for date in dates]):
                            break
                    else: # fall back to the first date phase
                        self._phase = (dates[0].timetuple()[self.level+1:5]+
                            (int(round(dates[0].second)), ))



    def __call__(self, x, pos=None):
        dt = num2date(x, self.tz)

        # Main label
        if self.level>=0 and self._phase_check_(dt):
            return self.strftime(dt, self.dual_fmt)

        # Intermediate label
        return self.strftime(dt, self.fmt)

    @classmethod
    def factory(cls, locator, **kwargs):
        """Helper to setup an appropriate formatter associated to a specified locator"""
        # Get best arguments for initializing and return instance
        if 'level' in kwargs:
            level = kwargs.pop('level')
        else:
            keys = sorted(AutoDateFormatter2.scaled.keys(), reverse=True)
            tmin, tmax = locator.axis.get_view_interval()
            dtmin, dtmax = locator.viewlim_to_dt()
            if hasattr(locator, 'get_locator'):
                locator = locator.get_locator(dtmin, dtmax)
            if hasattr(locator, '_freq'):
                freq = locator._freq
            else:
                freq = locator.rule._freq
            trange = tmax-tmin
            if freq == YEARLY :
                kwargs.setdefault('fmt', AutoDateFormatter2.scaled[keys[0]])
                level = -2
            elif freq == MONTHLY:
                if trange > 5*365:
                    kwargs.setdefault('fmt', AutoDateFormatter2.scaled[keys[0]])
                    level = -2
                else:
                    level = 'year'
            elif freq == WEEKLY:
                level = 'month'
            elif freq == DAILY:
                if dtmax.month != dtmin.month:
                    level = 'month'
                else:
                    level = 'week'
            elif freq == HOURLY:
                level = 'day'
            elif freq == MINUTELY:
                level = 'hour'
            elif freq == SECONDLY:
                if trange<1./(24.*60):
                    level = -3
                else:
                    level = 'minute'
        kwargs['locator'] = locator
        return cls(level, **kwargs)


class AutoDualDateFormatter(Formatter):
    """Automatic formatter that searches for the best :class:`DualDateFormatter`"""

    def __init__(self, locator, **kwargs):
        self._locator = locator
        self._kwargs = kwargs
        self._formatter = self.get_formatter()

    def get_formatter(self):
        return DualDateFormatter.factory(self._locator, **self._kwargs)

    def __call__(self, x, pos=None):
        self._formatter = self.get_formatter()
        s = self._formatter(x, pos)
        return s

def _get_split_locator_(locator, kwargs):
    ls = locator.split('/')
    lb = locator.split(':')
    if len(ls)>1 and ls[1].isdigit():
        kwargs.setdefault('interval', int(ls[1]))
        locator = ls[0]
    elif len(lb)>1:
        locator = lb[0]
        by = [int(v) for v in lb[1].split(',') if v.isdigit()]
        if by:
            kwargs.setdefault('by'+locator.lower(), by)
    return locator


def setup_time_axis(axis, auto=True, formatter=None, rotation=None,
        locator=None, minor_locator=True, minor_formatter=None, nominor=False,
        maxticks=None, tz=None, interval_multiples=False, **kwargs):
    """Setup tick positions and labels of an time axis

    :Params:

        - **axis**: :class:`~matplotlib.axis.Axis` instance (like ``P.gca().xaxis``)
        - **tz**, optional: Time zone.
        - **auto**, optional: Auto Scaling [default: True]
        - **rotation**, optional: Rotation angle of tick labels.
          If None, automatic [default: None]
        - **formatter**, optional: Date format:

            - 'mpl' or 0 : use the internal Matplotlib auto date locator.
            - "simple" or 1 or "auto": :class:`AutoDateFormatter2`
            - "dual", 2, None, tuple or dict: :class:`AutoDualDateFormatter`
              (default)
            - String: :class:`DateFormatter`
            - else a :class:`matplotlib.ticker.Formatter` instance.

        - **locator**, optional: Major locator. Can be within
          ['year','month','Weekday','day','hour','minute','second'],
          be like :class:`matplotlib.dates.MonthLocator`.
          or have a special value:

            - None or 'auto' or 'vacumm': use the :class:`AutoDateLocator`
              (default)
            - 'mpl': use the internal Matplotlib auto date locator.

        - **minor_locator**, optional: Minor locator.
        - **nominor**, optional: Do not try to add minor ticks [default: False]
        - **locator_<keyword>**, optional: <keyword> is passed to locator if
          locator is a string. If locator = 'month',
          locator = MonthLocator(locator_<keyword>=<value>).
        - **minor_locator_<keyword>**, optional: Same with minor_locator.
        - **maxticks**, optional: Maximal number of major ticks when
          default auto locator is used.
    """
    if isinstance(axis, int):
        axis = 'xy'[axis]
    if isinstance(axis, str):
        axis = getattr(P.gca(), xy+'axis')
    axes = axis.axes
    iaxis = 1-int(axis.axes.xaxis is axis)
    xy = 'xy'[iaxis]
    yx = 'yx'[iaxis]

    # Base
    try:
        axis.axis_date()
    except:
        import datetime
        axis.update_units(datetime.date(2009,1,1))

    # Scale axes
    if auto:
        axes.autoscale_view(**{'scale'+xy:True,'scale'+yx:False})

    # Guess locators
    major_locator = locator or kwargs.get('major_locator', None)
    if nominor: minor_locator = False
    locs = ['year','month','weekday','day','hour','minute','second']
    locs.extend([(loc+'s') for loc in locs])
    kwmjl = kwfilter(kwargs, 'major_locator', defaults=kwfilter(kwargs, 'locator'))
    kwmnl = kwfilter(kwargs, 'minor_locator')
    # - major
    if major_locator is None:
        major_locator = 'auto'
    if isinstance(major_locator, basestring):
        if major_locator.lower()=='mpl':
            major_locator = None
        else:
            major_locator = _get_split_locator_(major_locator, kwmjl)
            if (major_locator.lower().startswith('auto') or
                    major_locator.lower()=='vacumm'):
                maxticks = kwargs.get('nmax_ticks', maxticks) # compat
                kwmjl.setdefault('maxticks', maxticks)
                major_locator = AutoDateLocator2(**kwmjl)
            elif major_locator.lower() in locs:
                major_locator = eval(major_locator.lower().title()+'Locator(**kwmjl)')
    if major_locator:
        axis.set_major_locator(major_locator)
    else:
        major_locator = axis.get_major_locator()
    major_locator.set_axis(axis)
    if hasattr(major_locator, 'interval_multiples'):
        major_locator.interval_multiples = interval_multiples
    # - minor
    if minor_locator is 1 or minor_locator is True:
        minor_locator = 'auto'
    if isinstance(minor_locator, str):
        minor_locator = _get_split_locator_(minor_locator, kwmnl)
        if minor_locator.lower() in locs:
            minor_locator = eval(minor_locator.lower().title()+'Locator(**kwmnl)')
        elif minor_locator.lower().startswith('auto'):
            minor_locator = AutoDateMinorLocator(**kwmnl)
        else:
            minor_locator = None
    if minor_locator:
        axis.set_minor_locator(minor_locator)
    else:
        minor_locator = axis.get_minor_locator()
    minor_locator.set_axis(axis)
    if hasattr(minor_locator, 'interval_multiples'):
        minor_locator.interval_multiples = True

    #print 'minor_locator', minor_locator
    #print 'major_locator', major_locator

    # Tick format
    # - major
    fmt = formatter or kwargs.get('fmt', None) or kwargs.get('major_formatter', None)
    if fmt is None: fmt = 'dual' # default
    if fmt == 'mpl' or fmt==0: # matplotlib
        fmt = None
    elif fmt == 'simple' or fmt == 1 or fmt=="auto":
        fmt = AutoDateFormatter2(major_locator)
    elif fmt == 'dual' or fmt is True or fmt  == 2:
        fmt = AutoDualDateFormatter(major_locator)
    elif isinstance(fmt, basestring):
        fmt = DateFormatter(fmt)
    elif isinstance(fmt, list):
        fmt = AutoDualDateFormatter(major_locator, *fmt)
    elif isinstance(fmt, dict):
        fmt = AutoDualDateFormatter(major_locator, **fmt)
    elif isinstance(fmt, tuple):
        if len(fmt) == 1:
            fmt += ({}, )
        elif not isinstance(fmt[1], dict):
            fmt = (fmt[0], {})
        fmt = DualDateFormatter(fmt[0], **fmt[1])
    if fmt:
        axis.set_major_formatter(fmt)
    if rotation is not None and rotation!=0 \
        and not isinstance(fmt, (DualDateFormatter, AutoDualDateFormatter)):
        P.setp(axis.get_majorticklabels(), "rotation", rotation)
    # - minor
    if not nominor and not isinstance(minor_locator, NullLocator):
        if minor_formatter is True:
            minor_formatter = AutoDateFormatter2()
        elif isinstance(minor_formatter, basestring):
            minor_formatter = DateFormatter(minor_formatter)
        if minor_formatter:
            axis.set_minor_formatter(minor_formatter)



def add_agg_filter(objs, filter, zorder=None, ax=None, add=True):
    """Add a filtered version of objects to plot

    :Params:

        - **objs**: :class:`matplotlib.artist.Artist` instances.
        - **filter**: :class:`vacumm.misc._ext_plot.BaseFilter` instance.
        - **zorder**, optional: zorder (else guess from ``objs``).
        - **ax**, optional: class:`matplotlib.axes.Axes` instance.

    Inspired from http://matplotlib.sourceforge.net/examples/pylab_examples/demo_agg_filter.html .
    """
    # Input
    if not isinstance(objs, (list, tuple)):
        objs = [objs]
    elif len(objs)==0:
        return []

    # Filter
    if ax is None: ax = P.gca()
    shadows = FilteredArtistList(objs, filter)
    if hasattr(add, 'add_artist'):
        add.add_artist(shadows)
    elif add:
        ax.add_artist(shadows)

    # Text
    for t in objs:
        if isinstance(t, Text):
            t.set_path_effects([Normal()])

    # Adjust zorder
    if zorder is None or zorder is True:
        same = zorder is True
        if hasattr(objs, 'get_zorder'):
            zorder = objs.get_zorder()
        else:
            zorder = objs[0].get_zorder()
        if not same:
            zorder -= 0.1
    if zorder is not False:
        shadows.set_zorder(zorder)

    return shadows

def add_shadow(objs, width=3, xoffset=2, yoffset=-2, alpha=0.5, color='k',
        zorder=None, ax=None, add=True):
    """Add a drop-shadow to objects

    :Params:

        - **objs**: :class:`matplotlib.artist.Artist` instances.
        - **width**, optional: Width of the gaussian filter in points.
        - **xoffset**, optional: Shadow offset along X in points.
        - **yoffset**, optional: Shadow offset along Y in points.
        - **color**, optional: Color of the shadow.
        - **zorder**, optional: zorder (else guess from ``objs``).
        - **ax**, optional: class:`matplotlib.axes.Axes` instance.

    Inspired from http://matplotlib.sourceforge.net/examples/pylab_examples/demo_agg_filter.html .
    """
    if color is not None: color = RGB(color)
    try:
        gauss = DropShadowFilter(width, offsets=(xoffset, yoffset), alpha=alpha, color=color)
        return add_agg_filter(objs, gauss, zorder=zorder, ax=ax, add=add)
    except:
        warn('Cannot plot shadows using agg filters')

def add_glow(objs, width=3, zorder=None, color='w', ax=None, alpha=1., add=True):
    """Add a glow effect to text

    :Params:

        - **objs**: Plotted objects.
        - **width**, optional: Width of the gaussian filter in points.
        - **color**, optional: Color of the shadow.
        - **zorder**, optional: zorder (else guess from ``objs``).
        - **ax**, optional: class:`matplotlib.axes.Axes` instance.

    Inspired from http://matplotlib.sourceforge.net/examples/pylab_examples/demo_agg_filter.html .
    """
    if color is not None: color = RGB(color)
    try:
        white_glows = GrowFilter(width, color=color, alpha=alpha)
        return add_agg_filter(objs, white_glows, zorder=zorder, ax=ax, add=add)
    except:
         warn('Cannot add glow effect using agg filters')

def add_lightshading(objs, width=7, fraction=0.5, zorder=None, ax=None, add=True,
        **kwargs):
    """Add a light shading effect to objects

    :Params:

        - **objs**: Plotted objects.
        - **width**, optional: Width of the gaussian filter in points.
        - **fraction**, optional: Unknown.
        - **zorder**, optional: zorder (else guess from ``objs``).
        - **ax**, optional: class:`matplotlib.axes.Axes` instance.
        - Extra keywords are passed to :class:`matplotlib.colors.LightSource`

    Inspired from http://matplotlib.sourceforge.net/examples/pylab_examples/demo_agg_filter.html .
    """
    if zorder is None: zorder = True
    try:
        lf = LightFilter(width, fraction=fraction, **kwargs)
        return add_agg_filter(objs, lf, zorder=zorder, add=add, ax=ax)
    except:
         warn('Cannot add light shading effect using agg filters')


class MinuteLabel:
    def __init__(self, m, zonal=True, **kwargs):
        bfdeg = kwargs.pop('bfdeg', None)
        if not isinstance(m, Basemap) and (bfdeg is None or bfdeg=='auto'):
            bfdeg = False if N.size(m)==0 else (N.asarray(m)%1).ptp()!=0
        if zonal:
            if isinstance(m, Basemap):
                auto_minutes =  int(m.urcrnrlon) != int(m.llcrnrlon)
            else:
                auto_minutes =  int(min(m)) != int(max(m))
            self.func = lonlab
        else:
            if isinstance(m, Basemap):
                auto_minutes =  int(m.urcrnrlat) != int(m.llcrnrlat)
            else:
                auto_minutes =  int(min(m)) != int(max(m))
            self.func = latlab
        kwargs['bfdeg'] = bfdeg
        kwargs['decimal'] = False
        kwargs['no_zeros'] = True
        kwargs['auto_minutes'] = auto_minutes
        self.kwargs = kwargs
    def __call__(self, deg):
        return self.func(deg, **self.kwargs)

class AutoDegreesMinutesLocator(AutoLocator):
    def bin_boundaries(self, vmin, vmax):
        steps = None
        mn = 1/60.
        if (vmax-vmin) > 50.:
            steps = [1, 2, 3, 6, 9, 10]
        elif (vmax-vmin) > 3.:
            steps = [1, 2, 2.5, 3, 5, 10]
        elif vmax//mn != vmin//mn:
            nmn = int(N.ceil((vmax-vmin)*60.))+1
            nmax = 10
            if nmn > nmax:
                steps = [1, 10/6.,15/6.,20/6.,30/6.,10]
        self.set_params(steps=steps)
        return AutoLocator.bin_boundaries(self, vmin, vmax)

class AutoDegreesMinutesFormatter(Formatter):
    """Phase formatter with degrees and minutes

    :Example:

        >>> locator = xaxis.get_major_locator()
        >>> xaxis.set_formatter(AutoDegreesMinutesFormatter(locator, label='lon'))
    """
    def __init__(self, locator=None, **kwargs):
        if not isinstance(locator, Locator):
            data = N.asarray(locator)
            locator = data.min(), data.max()
        self.locator = locator
        self.kwargs = kwargs.copy()
        self.kwargs.setdefault('decimal', False)
    def __call__(self, x, pos=None):
        if isinstance(self.locator, tuple):
            vmin, vmax = self.locator
        else:
            vmin, vmax = self.locator.axis.get_view_interval()
        kwargs = self.kwargs.copy()
        kwargs['auto_minutes'] = vmax//1 != vmax//1
        return phaselab(x, **kwargs)



class DepthFormatter(Formatter):
    def __call__(self, x, pos=None):
        return "{:g}m".format(x)

def twinxy(xy, ax=None, fig=None):
    """Create an :func:`~matplotlib.axes.Axes` instance based on existing one(s)*

    This is an fusion and extension of :func:`matplotlib.pyplot.twinx` and
    :func:`matplotlib.pyplot.twiny`

    :Params:

        - **xy**: A string containing ``"x"`` and/or ``"y"``
        - **ax**, optional: Consider this axes instance instead of the current one.
        - **fig**, optional: Consider this figure instead of the current one.
    """
    if ax is None:
        if fig is None:
            ax = P.gca()
        else:
            ax = fig.gca()
    if not isinstance(xy, str): return ax
    fig = ax.figure
    twinx = 'x' in xy
    twiny = 'y' in xy
    if not twinx and not twiny:
        return ax
    kw = dict(frameon=False)
    if twinx and not twiny:
        kw.update(sharex=ax)
    if twiny and not twinx:
        kw.update(sharey=ax)
    nax = fig.add_axes(ax.get_position(True), **kw)
    if twinx:
        nax.yaxis.tick_right()
        nax.yaxis.set_label_position('right')
        ax.yaxis.tick_left()
        if not twiny: nax.xaxis.set_visible(False)
    if twiny:
        nax.xaxis.tick_top()
        nax.xaxis.set_label_position('top')
        ax.xaxis.tick_bottom()
        if not twinx: nax.yaxis.set_visible(False)
    return nax


def _add_label_(text, ax0, ax1, ay0, ay1, ax, va, ha, rotation, shadow, glow, **kwargs):
    if ax is None: ax = P.gca()
    fig = ax.figure
    b = ax.get_position()
    x = ax0*b.x0 + ax1*b.x1
    y = ay0*b.y0 + ay1*b.y1
    kwsh = kwfilter(kwargs, 'shadow')
    kwgl = kwfilter(kwargs, 'glow')
    o = fig.text(x, y, text, va=va, ha=ha, rotation=rotation, **kwargs)
    if shadow:
        add_shadow(o, ax=ax, **kwsh)
    if glow:
        add_glow(o, ax=ax, **kwgl)
    return o

def add_left_label(text, pos=.1, ax=None, va='center', ha='left',
    rotation=90., shadow=False, glow=False, **kwargs):
    """Add a text label to the left of a plot"""
    return _add_label_(text, 1-pos, pos, .5, .5, ax, va, ha, rotation, shadow, glow, **kwargs)

def add_right_label(text, pos=.1, ax=None, va='center', ha='right',
    rotation=-90., shadow=False, glow=False, **kwargs):
    """Add a text label to the left of a plot"""
    return _add_label_(text, pos, 1-pos, .5, .5, ax, va, ha, rotation, shadow, glow, **kwargs)

def add_top_label(text, pos=.1, ax=None, va='top', ha='center',
    rotation=0., shadow=False, glow=False, **kwargs):
    """Add a text label to the left of a plot"""
    return _add_label_(text, .5, .5, pos, 1-pos, ax, va, ha, rotation, shadow, glow, **kwargs)

def add_bottom_label(text, pos=.1, ax=None, va='top', ha='center',
    rotation=0., shadow=False, glow=False, **kwargs):
    """Add a text label to the left of a plot"""
    return _add_label_(text, .5, .5, 1-pos, pos, ax, va, ha, rotation, shadow, glow, **kwargs)


_locations = ['upper right',
         'upper left',
         'lower left',
         'lower right',
         'center left',
         'center right',
         'lower center',
         'upper center',
         'center center',
         'center']

def loc2tuple(loc, xpad=0.02, ypad=0.02):
    """Location as tuple of (x,y) in relative coordinates"""
    if isinstance(loc, tuple): return loc
    if isinstance(loc, basestring):
        if len(loc.split())==1:
            loc = loc+' '+loc
    assert loc in _locations, 'loc must be a tuple or a string (%s)'% (", ".join(_locations))
    sy, sx = loc.split()
    xpad = abs(xpad)
    if xpad>.5: xpad = 1-xpad
    ypad = abs(ypad)
    if ypad>.5: ypad = 1-ypad
    x = {"left":xpad, 'center':.5, 'right':1-ypad}[sx]
    y = {"lower":xpad, 'center':.5, 'upper':1-ypad}[sy]
    return x, y

def loc2align(loc, ha=None, va=None, margin=1/3.):
    """From location to dict of ha and va

    :Params:

        - **loc**: (x,y) or string like "upper right.
        - **ha/va**, optional: alignments.

    :Return: dict(ha=ha, va=va)
    """
    x, y = loc2tuple(loc, xpad=0, ypad=0)
    if ha is None:
        if x<margin:
            ha = 'left'
        elif x>(1-margin):
            ha = 'right'
        else:
            ha = 'center'
    if va is None:
        if y<margin:
            va = 'bottom'
        elif y>(1-margin):
            va = 'top'
        else:
            va = 'center'
    return dict(ha=ha, va=va)

def loc2offset(loc, xpad, ypad=None, margin=1/3.):
    """From location to (xoffset, ypoffset) where x/yoffset = +/- x/ypad

    :Return: (xoffset,yoffset)
    """
    align = loc2align(loc, margin=margin)
    xoffset = yoffset = 0
    if ypad is None: ypad = xpad
    if align['ha'] == 'left':
        xoffset = xpad
    elif align['ha'] == 'right':
        xoffset = -xpad
    if align['va'] == 'top':
        yoffset = -ypad
    elif align['va'] == 'bottom':
        yoffset = ypad
    return xoffset, yoffset


def best_loc_map(m, onland=True, allowed=_locations):
    """Find the best location on a plot map according to the land/ocean repartition"""
    from collections import OrderedDict
    if isinstance(m, Map): m = m.map
    if not hasattr(m, 'coastpolygons'):
        sloc = allowed[0]
    else:
        from grid.masking import polygons
        from grid import bounds2d
        from _geoslib import Polygon
        fractions = N.zeros((3, 3))
        polys = polygons(m)
        x = m.xmin + (m.xmax-m.xmin)*(0.5+N.arange(3))/3.
        y = m.ymin + (m.ymax-m.ymin)*(0.5+N.arange(3))/3.
        xxb, yyb = bounds2d(x, y)
        for j in xrange(3):
            for i in xrange(3):
                pcell = Polygon(N.array([xxb[j, i], yyb[j, i]]).T)
                for pcoast in polys:
                    if pcell.within(pcoast):
                        fractions[j, i] = 1.
                        break
                    if not pcell.intersects(pcoast):
                        continue
                    for p in pcell.intersection(pcoast):
                        fractions[j, i] += p.area()

        if not onland:
            fractions = 1-fractions
        jmax, imax = N.unravel_index(fractions.argmax())
        sloc = '%s %s'%(['left', 'center', 'right'][imax], ['bottom', 'center', 'top'][jmax])

    return _loc2tuple_(sloc)




def add_param_label(text, loc=(0.01, 0.01), color='0.4', size=8, family='monospace', fig=None, **kwargs):
    """Add parameters description to the bottom/left of the figure

    :Example:

        >>> c = curve2(sst, show=False)
        >>> c.add_param_label(dict(max=.23, kz=0.25))

    :Params:

        - **text**: Either a string or a dictionary.
    """
    # Convert dict to text
    if isinstance(text, dict):
        clsname = text.__class__.__name__
        tt = []
        for item in text.items():
            tt.append('%s=%s'%item)
        text = ', '.join(tt)

    if fig is None: fig = P.gcf()
    loc = loc2tuple(loc, xpad=0.01, ypad=0.01)
    return fig.text(loc[0], loc[1], str(text), color=color, size=size, family=family, **kwargs)



class Animator(object):
    """Animate objects on a figure

    :Example:

        >>> anim = Animator()
        >>> anim.append(P.plot([5,6])
        >>> anim.append(P.plot([4,8])
        >>> anim.make_animation()

    """
    def __init__(self, fig=None):
        self.objs = []
        if fig is None: fig = P.gcf()
        self.fig = fig

    def append(self, obj, frame=None):
        # Set frame index
        nobjs = len(self.objs)
        if frame == 'last':
            frame == -1
        elif frame=='new':
            frame = None
        if nobjs==0:
            frame = None
        if isinstance(frame, int):
            if frame<0:
                frame = max(0, nobjs-frame)
            elif frame>nobjs-1:
                frame = None
            else:
                frame = min(nobjs-1, frame)

        # Get layer
        if frame is None: # New layer
            frameobjs = []
            self.objs.append(frameobjs)
            frame = len(self.objs)-1
        else: # Old layer
            frameobjs = self.objs[frame]

        # Append data
        if isinstance(obj, list):
            if hasattr(obj[0], 'get_figure'):
                if self.fig is None and obj:
                    self.fig = obj[0].get_figure()
                else:
                    assert self.fig == obj[0].get_figure()
            frameobjs.extend(obj)
        else:
            if hasattr(obj, 'get_figure'):
                if self.fig is None:
                    self.fig = obj.get_figure()
                else:
                    assert self.fig == obj.get_figure()
            frameobjs.append(obj)

        return frame

    def make_animation(self, **kwargs):
        """Create the animation object"""
        if self.fig is None: return
        if not self.objs: return
        from matplotlib.animation import ArtistAnimation
        self.animation = ArtistAnimation(self.fig, self.objs, **kwargs)
        return self.animation


def hlitvs(color='.95', axis='x', units='ticks', axes=None, maxticks=10, **kwargs):
    """Highlight intervals between ticks

    :Params:

        - *color*: Background color.
        - *axis*: Matplotlib axis or 'x'/'y'.
        - **units**, optional: Interval types

            - `"ticks"`: Use axis ticks.
            - A list or array: use it as ticks.
            - `"auto"`: Estimate ticks using :class:`~matplotlib.ticker.AutoLocator`.
            - `"date"` or `"time"`: Same but using :class:`AutoDateLocator2`.
            - A time string like `"day"` for using :class:`~matplotlib.dates.DayLocator`.
            - A locator.

        - Other keyparam are passed to :func:`~matplotlib.pyplot.axhspan`
          or :func:`~matplotlib.pyplot.axvspan`.
    """
    # Get axis
    axes = kwargs.pop('ax', axes)
    if isinstance(axis, Axes):
        axes = axis
        axis = 'x'
    if isinstance(axis, basestring):
        if units.endswith('s'): units = units[:-1]
        if axes is None: axes = P.gca()
        axis = getattr(axes, axis+'axis')
    xy = 'xy'[isinstance(axis, YAxis)]
    axes = axis.axes

    # Locators and ticks
    kwloc = kwfilter(kwargs, 'locator_')
    if units=='ticks':
        ticks = axis.get_ticklocs()
    elif isinstance(units, (list, N.ndarray)):
        ticks = units
    else: # locator
        if units == 'auto':
            locator = AutoLocator(**kwloc)
        elif units=='date' or units=='time':
            locator = AutoDateLocator2(**kwloc)
        elif isinstance(units, basestring):
            locator = getattr(matplotlib.dates, units.title()+'Locator')(**kwloc)
        locator.set_axis(axis)
        ticks = locator()
    tmin, tmax = axis.get_view_interval()
    if isinstance(ticks, N.ndarray):
        ticks = ticks.tolist()
    else:
        ticks = list(ticks)
    if tmin!=ticks[0]: ticks.insert(0, tmin)
    if tmax!=ticks[-1]: ticks.append(tmax)

    # Patch
    alpha = RGBA(color)[-1]
    kwargs.setdefault('zorder', 0)
    kwargs.update(facecolor=color)
    kwargs.setdefault('linewidth', 0)
    kwargs.setdefault('label', '_nolegend_')
    kwargs.setdefault('alpha', alpha)
    objs = []
    axspan = getattr(axes, 'ax%sspan'%'vh'[xy=='y'])
    for i in xrange(0, len(ticks)-1, 2):
        objs.append(axspan(ticks[i], ticks[i+1], **kwargs))
#        objs[-1].set_zorder(0)
    return objs


class _PPatch(Patch):
    """
    Base class for patches with size in pixels.

    You must define your coordinates in  attribute _path_coords
    between -5 and 5.
    """
    def __str__(self):
        return "Compass(%s,%s;%s)" % (self.center[0], self.center[1], self.size)

    def __init__(self, xy, size, angle=0, anglemode='pixels', posref=None, offset=None,
        **kwargs):
        """
        *xy*
          center of compass

        *size*
          size of compass


        """
        Patch.__init__(self, **kwargs)

        self.center = xy
        self.size = size
        self.angle = angle
        self.anglemode = anglemode
        self.posref = posref
        self.offset = offset

        path_coords = N.array(self._path_coords, 'f')
        path_coords *= 0.1
        self._path = Path(path_coords, closed=False)

        self._patch_transform = mtransforms.IdentityTransform()


    def _update_transform_(self):
        self._patch_transform = self._get_pixrot_transform_()

    def _get_pixrot_transform_(self, size=None):
        if size is None: size=self.size
        transform = Artist.get_transform(self) if self.is_transform_set() else None
        return _transform_pixrot_(self.center, size, self.angle, self.anglemode,
            ax = self.axes, posref=self.posref, transform=transform)

    def get_patch_transform(self):
        self._update_transform_()
        return self._patch_transform

    def get_path(self):
        return self._path

    def get_transform(self):
        t = Patch.get_transform(self)
        if self.offset is not None:
            t = t+mtransforms.Affine2D().translate(*self.offset)
        return t

    def draw(self, renderer):
        self._update_transform_()
        Patch.draw(self, renderer)

class QuarterRightCompass(_PPatch):
    _path_coords = [[0, 0], [1, 1], [0, 5]]
#    _path_coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
class QuarterLeftCompass(_PPatch):
    _path_coords = [[0, 0], [-1, 1], [0, 5]]

class ArrowCompass(_PPatch):
    _path_coords = [[-.75, -5], [-.25, 3.5], [-1, 3], [0, 5], [1, 3], [.25, 3.5],
        [.75, -5], [0, -4]]

class PPText(Text):
    def __init__(self, *args, **kwargs):
        self._vc_offset = kwargs.pop('offset', 0)
        self._vc_ppatch = kwargs.pop('ppatch', 0)
        Text.__init__(self, *args, **kwargs)

    def get_transform(self):
        t = self._vc_ppatch._get_pixrot_transform_(self._vc_ppatch.size+self._vc_offset)+ \
            Artist.get_transform(self._vc_ppatch)
        if self._vc_ppatch.offset is not None:
            t = t + mtransforms.Affine2D().translate(*self._vc_ppatch.offset)
        return t

def _transform_(transform, default=None, ax=None, fig=None):
    if ax is None: ax = P.gca()
    if default is None: default = ax.transData
    if transform is None: return default
    if isinstance(transform, basestring):
        if transform.lower()=='figure':
            if fig is None:
                fig = P.gcf() if ax is None else ax.fig
            return fig.transFigure
        if hasattr(ax, 'trans'+transform.title()):
            return getattr(ax, 'trans'+transform.title())
    return transform


def _transform_pixrot_(center, size, angle, anglemode, transform=None, posref=None, ax=None):
        """Special transform: size in pixel, rotation, translation in data

        :Params:

            - **trandform**: Transform for the center location.
            - **posref**: Reference point for the angle if anglemode == 'data'.
              It default to the center. Values must be in data coordinates.
        """

        if ax is None: ax = P.gca()
        transform = _transform_(transform, 'data', ax=ax)

        # center
        xp0, yp0 = transform.transform(center)

        # scales
        xp1, yp1 = xp0 + size,  yp0 + size
        xcorner, ycorner = transform.inverted().transform_point((xp1, yp1))
        xscale = xcorner-center[0]
        yscale = ycorner-center[1]

        # angle
        anglemode = str(anglemode)
        if anglemode is None or anglemode.startswith('dot'): anglemode = 'pixels'
        if not anglemode.startswith('pix') and not anglemode.startswith('dat'):
            raise PlotError('Unkown anglemode: %s'%anglemode)
        if anglemode.startswith('dat'): # from data to pixels

            ar = N.radians(angle)

            if posref is None:
                posref = center
                if transform is not ax.transData:
                    posref = (transform-ax.transData).transform_point(posref)
            if posref[0]>1.e20 or posref[1]>1e20:
                raise PlotError('Coordinates of reference seems wrong: %s,%s'%tuple(posref))
            xlen = ylen = min(posref)*0.01
            xr0, yr0 = ax.transData.transform_point(posref)
            xr1, yr1 = ax.transData.transform_point((posref[0]+xlen*N.cos(ar),
                posref[1]+ylen*N.sin(ar)))
            angle = N.degrees(N.arctan2(yr1-yr0, xr1-xr0))

        # Transform
        t = mtransforms.Affine2D() \
            .rotate_deg(angle) \
            .scale(xscale, yscale) \
            .translate(*center)
        return t

def _add_text_(ax, text, **kwtext):
    ax._set_artist_props(text)
    text.update(kwtext)
    ax.texts.append(text)
    return text

def add_compass(x, y, size=40, facecolor1='k', facecolor2='w', edgecolor='k',
    text_color='k', text_size='medium', linewidth=0.5, alpha=1, text_pad=5,
    style='simple',
    angle=0, ax=None, m=None, anglemode='pixels', **kwpatch):
    """Add a compass to the current plot

    :Params:

        - **x/y**: Position in data coordinates.
        - **size**, optional: size in opixels.
        - **style**, optional: "simple", "fancy".
        - **angle**, optional: Angle of the north (trigo).
        - **anglemode**, optional: "pixels" or "data".
        - **facecolor1**, optional: Dark face color.
        - **facecolor2**, optional: Light face color.
        - **edgecolor**, optional: Edges color.
        - **text_color**, optional: Text labels color.
        - **text_size**, optional: Text labels size.
        - **text_<param>**, optional: ``<param>`` is passed to
          :class:`~matplotlib.text.Text` when adding labels.
        - **fontsize/color**, optional: Same as ``text_size/color``.
    """

    kwtext = kwfilter(kwpatch, 'text_')
    kwtext.update(color=text_color, size=text_size, alpha=alpha)
    if ax is None:
        if m:
            ax = m.ax or m._check_ax
        else:
            ax = P.gca()
    if m:
        x, y = m(x, y)
    if style is None: style = 'arrow'
    else: style = str(style)



    # Loop on quarters
    dict_check_defaults(kwpatch, linewidth=linewidth, alpha=alpha, clip_on=False,
        edgecolor=edgecolor)
    patches = []
    texts = []
    if style.startswith('a') or style.startswith('s'):

        # Patch
        kwp = kwpatch.copy()
        kwp['facecolor'] = facecolor1
        patches.append(ax.add_patch(ArrowCompass((x, y), size, angle=angle, anglemode=anglemode, **kwp)))

        # Label (north only)
        text = PPText(0, .5, 'N', offset=text_pad, ppatch=patches[-1], ha='center', va='bottom')
        texts.append(_add_text_(ax, text, **kwtext))

    elif style.startswith('f'):

        for iq, (label, ha, va) in enumerate([
            ('N', 'center', 'bottom'),
            ('W', 'right', 'center'),
            ('S', 'center', 'top'),
            ('E', 'left', 'center')]):

            aa = iq*90.

            # Loop on back/white patches
            for (C, fc) in [
                (QuarterRightCompass, facecolor1),
                (QuarterLeftCompass, facecolor2),
                ][:]:
                kwp = kwpatch.copy()
                kwp['facecolor'] = fc
                patches.append(ax.add_patch(C((x, y), size, angle=aa+angle,
                    anglemode=anglemode, **kwp)))


            # Add text labels
            text = PPText(0, .5, label, offset=text_pad, ppatch=patches[-1], ha=ha, va=va)
            texts.append(_add_text_(ax, text))

    else:

        raise PlotError('Unknown compass style. Please use one of : '+
            ', '.join(['simple', 'fancy']))


    return patches, texts

def _asnum_(xy, atleast_1d=False):
    """Get xy as a number

    If it is a number or a numpy array, it not converted,
    else it is converted with :func:`numtime`.

    xy can also be a list, a list of lists, etc.
    """
    if isinstance(xy, N.ndarray):
        return xy
    single = not isinstance(xy, list)
    xys = xy if not single else [xy]
    out = []
    for xy in xys:
        if is_numtime(xy):
            out.append(xy)
            continue
        if isinstance(xy, list):
            out.extend(xy)
            continue
        xy = numtime(xy)
        out.append(xy)
    if atleast_1d:
        return N.ma.atleast_1d(out)
    if single:
        return out[0]
    return out


docfiller.scan(Plot, Plot.format_axes, Plot.load_data, Plot._check_order_,
    Plot.pre_plot, Plot.post_plot,
    Plot1D._set_axes_, Plot1D._check_order_,
    Curve.plot, Curve.load_data,
    Bar.plot,
    Stick.load_data, Stick.plot,
    ScalarMappable, ScalarMappable.colorbar, ScalarMappable.post_plot,
    ScalarMappable.get_cmap, ScalarMappable.get_levels,
    Plot2D.load_data,
    Plot2D.plot, Plot2D.plot_contour, Plot2D.plot_fill,
    Plot2D.plot_quiver, Plot2D._set_axes_,
    Map.load_data, Map.pre_plot, Map.post_plot,
    quiverkey_=QuiverKey.quiverkey,
    savefig_=Plot.savefig)
