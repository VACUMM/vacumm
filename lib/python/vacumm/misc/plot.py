# -*- coding: utf8 -*-
"""
Generic plots using Matplotlib and taking advantage of CDAT

Tutorial: ":ref:`user.tut.misc.plot`".
"""
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
import cPickle
import re, os, glob
from copy import copy
from itertools import cycle
from operator import isNumberType
from string import Template
from tempfile import mktemp

import cdtime
import numpy as N,MV2, cdms2
import pylab as P
from cdms2.axis import AbstractAxis
from genutil import minmax, statistics
from matplotlib.axes import Subplot
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import ColorConverter,is_color_like, Normalize
from matplotlib.dates import (DateFormatter, MonthLocator, WeekdayLocator, YearLocator,
    DayLocator, HourLocator, MinuteLocator, SecondLocator, MONDAY, WEEKLY, YEARLY, MONTHLY,
    AutoDateLocator, AutoDateFormatter, MO, DAILY, HOURLY, num2date)
from matplotlib.patches import Wedge,Shadow,Circle, Arc
from matplotlib.ticker import FormatStrFormatter, Formatter, FixedLocator
import matplotlib.image as mpimg
from mpl_toolkits.basemap import Basemap

from .misc import (dict_aliases, latlab, lonlab, deplab,
    geo_scale, kwfilter, dict_check_defaults, auto_scale)
from .atime import (mpl, time, axis_add, compress, SpecialDateFormatter,
    )
from .axes import check_axes, isaxis, istime, axis_type
from .color import simple_colors, Scalar2RGB, get_cmap
from .color import get_cmap, Scalar2RGB
from .core_plot import (add_glow, add_shadow, add_agg_filter, hlitvs,
    AutoDateFormatter2,  DualDateFormatter, add_lightshading,
    AutoDateLocator2, AutoDateMinorLocator, AutoDualDateFormatter, add_compass,
    add_right_label, add_left_label, add_top_label, add_bottom_label,
    add_lightshading, add_param_label, get_quiverkey_value)
from .docstrings import docfill
from .grid.misc import (meshgrid, meshbounds, var2d)
from .grid import regridding
from .grid.basemap import gshhs_reslist, gshhs_autores, cache_map, cached_map
from .grid.misc import get_xy
from .misc import is_iterable, broadcast, squarebox, zoombox
from .phys.units import deg2m, m2deg
from ..__init__ import VACUMMError


__all__ = [ 'traj', 'ellipsis',
    'add_colorbar', 'scolorbar','rotate_tick_labels', 'rotate_xlabels', 'rotate_ylabels',
    'scale_xlim', 'scale_ylim', 'xhide', 'yhide', 'xdate', 'ydate',
    'get_locator_from_units', 'set_major_locator', 'get_major_locator',
    'set_minor_locator', 'get_minor_locator', 'add_logo', 'wedge',
    'add_time_mark', 'xi2a', 'yi2a', 'yi2f', 'xi2f', 'gobjs', 'add_key',
    'colorbar', 'hldays',  'xscale', 'yscale', 'xrotate', 'yrotate',
    'add_grid', 'savefigs', 'decorate_axis', 'markers', 'linestyles',
    'make_movie', 'taylor', 'vtaylor', 'get_cls', 'dtaylor', 'rankhist',
    'curve2', 'bar2', 'stick2', 'hov2' , 'map2', 'section2',
    'minimap', 'add_map_point', 'add_map_line', 'add_map_lines', 'add_map_box',  'rshpere_wgs84',
    'add_glow', 'add_shadow', 'add_agg_filter', 'plot2d', 'hlitvs',
    'add_compass', 'add_param_label', 'dtarget', 'add_map_places',
    'add_right_label', 'add_left_label', 'add_top_label', 'add_bottom_label',
    'get_quiverkey_value', 'add_lightshading']
__all__.sort()


MA = N.ma
MV = MV2
cdms = cdms2

#try:
#   from matplotlib.toolkits.basemap import Basemap
#except:



#: WGS radius of earth
rshpere_wgs84 = (6378137.0,6356752.3141)

#: Useful markers list
markers = [
    'o', #  : circl
#   '.', #  : point
#   ',', #  : pixel
    'v', #  : triangle_down
    '^', #  : triangle_up
    '<', #  : triangle_left
    '>', #  : triangle_right
#   '1', #  : tri_down
#   '2', #  : tri_up
#   '3', #  : tri_left
#   '4', #  : tri_right
    '8', #  : 'octogon
    's', #  : square
    'p', #  : pentagon
    'h', #  : hexagon1
#   'H', #  : hexagon2
    '+', #  : plus
    'x', #  : x
#   'D', #  : diamond
    'd', #  : thin_diamond
#   '|', #  : vline
#   '_', #  : hline
]

#: Useful linestyles list
linestyles = [
    '-' , #         : solid',
    '--', #         : dashed',
    '-.', #         : dash_dot',
    ':', #          : dotted',
]



def traj(lons,lats,res=None,dots=True,color=True,m=None,zoom=.9,cmap=None,colorbar=True,**kwargs):

    # Start the plot
    _start_plot_(**kwargs)

    # Time series
    lons = _check_var_(lons,1,'t',**kwargs)
    lats = _check_var_(lats,1,'t',**kwargs)

    # Map
    if m is None:
        lon_min,lon_max = minmax(lons)
        lat_min,lat_max = minmax(lats)
        lon_range = (lon_min+lon_max)*.5 + .5*(lon_max-lon_min)/zoom * N.array([-1,1])
        lat_range = (lat_min+lat_max)*.5 + .5*(lat_max-lat_min)/zoom * N.array([-1,1])
        kwmap = kwfilter(kwargs,'map',defaults=dict(lon=lon_range,lat=lat_range))
        kwmap['show'] = False
        kwmap['res'] = kwargs.pop('resolution',res)
        m = map(**kwmap)

    # Real coordinates
    xx,yy = m(lons,lats)

    #TODO:finish traj
    # Plot
    kwargs.setdefault('linewidth',2)
    if color is not True and is_color_like(color): # Classic one color
        kwdots = kwfilter(kwargs,'dots')
        lines = m.plot(xx,yy,**kwargs)
        kwdots.setdefault('color',lines[0].get_color())
        if dots: m.plot(xx,yy,'o',**kwdots)

    elif color is not False: # Color varying
        if not is_iterable(color) or len(color) != len(xx):
            tt = xx.getAxis(0)
        if istime(tt): tt = mpl(tt)
        mp = ScalarMappable(cmap=cmap)#,norm=no_norm())
        mp.set_array(N.array(tt))
        ttmin = min(tt)
        dtt = max(tt)-ttmin
        color = tuple([tuple(rgba) for rgba in mp.cmap((.5*(tt[:-1]+tt[1:])-ttmin)/dtt)])
        args = []
        for ip in xrange(len(xx)-1):
            args.extend([xx[ip:ip+2],yy[ip:ip+2]])
        lines = []
        for ip,line in enumerate(P.gca()._get_lines(*args,**kwargs)):
            line.set_color(color[ip])
            P.gca().add_line(line)
            lines.append(line)
        if dots:
            kwdots = kwfilter(kwargs,'dots',defaults=dict(linestyle='o',markersize='5'))
            kwdots.update(kwargs)
            for line in P.gca()._get_lines(*args,**kwdots):
                line.set_color(color[ip])
                P.gca().add_line(line)

        if colorbar:
            kwargs.setdefault('colorbar_horizontal',True)
##          kwargs.setdefault('colorbar_aspect',20)
            colorbar(mp,tt, True, **kwargs)

    # End plot
    kwargs.setdefault('title','Trajectory')
    _end_plot_(var=lons,**kwargs)

    return lines


def vtaylor(data, ref=None, ref_color='red', color=None, label=None, title=None, xoffset=5, yoffset=-2, show=False, stdmax=None, savefig=None, **kwargs):
    """Plot a Taylor diagram using :mod:`vcs` module from CDAT

    - **data**: A single or several points where a point is a pair of (std, corr).
    - *ref*: Reference value (standard deviation of observations).
    - *color*: Color of the points. It can be a list to set a different color for each pair.
    - *label*: Label of the points. Obviously, if there are several points, there must have several labels.
    - *x|yoffset*: Offset in font size unit of label position.
    - *title*: A string or ``False``.
    - *savefig*: Save figure to ``savefig``.
    - *savefig_<keyword>*: ``keyword`` is used for saving figure.
    - Other keywords are passed to :meth:`plot`.

    Example:

    >>> taylor([21.,.83], ref=30., title='Single point')
    >>> taylor(([[[21.,.83], [33., .7]], ref=28., color=['blue', 242], label=['Point 1', 'Point2'])


    .. warning::

        This function uses the :mod:`vcs` module for plotting.
        Therefore, yan can't alter the plot with :mod:`matplotlib` functions.
    """
    # Data
    if cdms2.isVariable(data):
        data = data.clone()
    else:
        data = MV2.array(data)
    if data.ndim == 1:
        data = MV2.resize(data, (1, 2))
    if title is False:
        title = None
    elif title is None and hasattr(data, 'long_name') and \
        not data.long_name.startswith('variable'):
        title = data.long_name
    data.id = ' '
    np = len(data)

    # Init canva and diagram
    import vcs
    x=vcs.init(size=1)
    td=x.createtaylordiagram()
    if ref is not None:
        if hasattr(ref, 'std') and callable(ref.std):
            ref = ref.std()
        td.referencevalue = ref
    td.preserveaspectratio='n'
    if color is not None:
        if not isinstance(color, list):
            color = broadcast(list([color]), np)
        td.Marker.id_color = color
    if label is not None:
        if not isinstance(label, list):
            label = list(label)
        td.Marker.id = label
#   td.Marker.symbol=['dot','cross','circle', 'dot']
    if not isinstance(xoffset, list):
        xoffset = broadcast(list([xoffset]), np)
    if not isinstance(yoffset, list):
        yoffset = broadcast(list([yoffset]), np)
    td.Marker.xoffset=xoffset
    td.Marker.yoffset=yoffset

    # Create the template
    t=x.createtemplate(source="deftaylor")
    t.data.x1=.1
    t.data.x2=.9
    t.data.y1=.1
    t.data.y2=.9
    t.xtic1.y1 = t.data.y1
    t.xtic1.y2 = t.xtic1.y1+.01
    t.xmintic1.y1 = t.data.y1
    t.xmintic1.y2 = t.xmintic1.y1+.005
    t.xlabel1.y = t.data.y1-.01
    t.xname.y = t.xtic1.y1 - .05

    # Sets the correlation major line
    l=x.createline()
    l.type='dash'
    l.width=1
    t.ytic2.x1=.1
    t.ytic2.x2=.9
    t.ytic2.line=l
    majorcor=vcs.mklabels([.1,.3,.6,.8,.9,.95,.99])
    td.cticlabels1=majorcor

    # Sets the correlation minor line
    l=x.createline()
    l.type='dot'
    l.width=1
    l.color=[252] # grey
    t.ymintic2.x1=.1
    t.ymintic2.x2=.9
    t.ymintic2.line=l
    minorcor=vcs.mklabels([.2,.4,.5,.7,.85,.91,.92,.93,.94,.96,.97,.98,.998])
    td.cmtics1=minorcor

    # Sets standard dev major tics
    l=x.createline()
    l.type='dash'
    l.width=1
    t.xtic2.y1=.1
    t.xtic2.y2=.9
    t.xtic2.line=l
    t.xtic2.priority=1
    vmax = data[:, 0].max()
    stdticks = auto_scale([0., vmax], vmin=0.)
    mjrstd1=vcs.mklabels(stdticks)
    td.xticlabels1=mjrstd1
    mjrstd2=vcs.mklabels(stdticks)
    td.yticlabels1=mjrstd2

    # Sets the std minor line
    l=x.createline()
    l.type='dot'
    l.width=1
    l.color=[252] # grey
    t.xmintic2.y1=.1
    t.xmintic2.y2=.9
    t.xmintic2.line=l
    t.xmintic2.priority=1
    mnrstd=vcs.mklabels((stdticks[:-1]+stdticks[1:])/2.)
    td.xmtics1=mnrstd
    td.ymtics1=mjrstd1

    # Sets the reference
    l=x.createline()
    l.type='solid'
    l.width=2
#   l.color=[ref_color]
    l.color=[242]
    t.line2.line=l

    # Sets the maximum for arc value
    if stdmax is None: stdmax = stdticks[-1]
    td.max=stdmax

    # Main plot
    kwsf = kwfilter(kwargs, 'savefig')
    x.plot(data,t,td, bg=int(not show), title=title, **kwargs)

    # Save
    if savefig is not None:
        for suf, func in [('eps', 'eps'), ('ps', 'postscript'), ('png', 'png'), ('gif', 'gif')]:
            if savefig.endswith('.'+suf):
                getattr(x, func)(savefig, **kwsf)
                break

def dtaylor(stats, stdref=None, labels=None, colors=None, **kwargs):
    """Directly plot a Taylor diagram

    :Params:

        - **stats**: Statistics is the form ``[std, corr]`` or ``[[std1,corr1],[std2,corr2],...]``
        - **stdref**: Standard deviation of reference
        - **labels**,optional: Label of the points. Obviously, if there are several points, there must have several labels.
          If ``None``, it tries to get it from the `long_name` attribute. If ``False``, it doesn't display it.
        - **reflabel**,optional: Name of the reference label. It defaults to ``"Reference"``.
        $_taylor_
    """
    stats = N.asarray(stats)
    if stats.ndim==1:
        stats.shape = (1,)+stats.shape
    elif stats.shape[1]!=2:
        stats = stats.T
    return _taylor_(stats, labels, colors, stdref=stdref, **kwargs)


def taylor(datasets, ref, labels=False, colors=None, units=None, normalize=True, **kwargs):
    """Plot a Taylor diagram after computing statistics

    :Params:

        - **datasets**: One or several arrays that must be evaluated against the reference.
          Best is to provide :mod:`cdms2` variables with an appropriate `long_name` attribute
          to set labels automatically.
        - **ref**: Reference array of the same shape as all ``datasets`` arrays.
        - **labels**,optional: Label of the points.
          Obviously, if there are several points, there must have several labels.
          If ``None``, it tries to get it from the `long_name` attribute.
          If ``False``, it doesn't display it.
        - **reflabel**,optional: Name of the reference label.
          It defaults to the `long_name` attribute or to ``"Reference"``.
        - **normalize**,optional: Normalize standard deviation and RMSE
        $_taylor_

    :Example:

        Basic usage:

        >>> taylor(var_model, var_obs, title='Single point', reflabel='Observations')

        Two variables with the same reference:

        >>> taylor([var_model1,var_model2], var_ref)

        Two variables without the same reference (like temprature and salinity):

        >>> taylor([sst,sss], [sstref,sssref], labels=True)

        A huge model ensemble of SST with colors varying with averaged SSS:

        >>> taylor([sstm for sstm in sstens], sstsat, colors=cdutil.averager(sss, 'xy'))

    :Tutorial: ":ref:`user.tut.misc.plot.basic.taylor`".
    """
    # Check datasets, units and labels
    # - datasets
    if isinstance(datasets, N.ndarray):
        datasets = [datasets]
    datasets = list(datasets)
    n = len(datasets)
    labels = broadcast(labels, n, mode='last')
    for i, d in enumerate(datasets):
        d = MV2.asarray(d)
        if units in [True, None]:
            units = getattr(d, 'units', None)
        if labels[i] in [True, None]:
            labels[i] = getattr(d, 'long_name', None)
        d = d.asma().ravel()
        datasets[i] = d
    # - references
    if isinstance(ref, N.ndarray):
        ref = [ref]
    ref = list(ref)
    ref = broadcast(ref, n, mode='last')
    for i, r in enumerate(ref):
        r = MV2.asarray(r)
        if units in [True, None]:
            units = getattr(r, 'units', None)
        if labels[i] in [True, None]:
            labels[i] = getattr(r, 'long_name', None)
        ref[i] = r.asma().ravel()

    # Statistics
    stats = []
    stdref = []
    for i,(d, r) in enumerate(zip(datasets, ref)):
        stats.append([d.std(), statistics.__correlation(d, r)])
        stdref.append(r.std())
    stats = N.ma.asarray(stats)

    # Diagram
    return _taylor_(stats, labels, colors, stdref=stdref,
        units=units, normalize=normalize, **kwargs)


def _taylor_(stats, labels=False, colors=None, stdref=None, units=None, refcolor='.5',
    reflabel='Observations', normalize=True, addref=True, legend=True, negative=None,
    minimal=False, stdmax=None, autoresize=True, scatter=False, zorder=12, cmap=None,
    **kwargs):
    """- **stdref**, optional: Standard deviation of the reference used either when
          normalizing statistics or to plot the reference arc and point.
          All statistics can have a different reference standard deviation;
          in this case, statistics are normalized (``normalize`` is set to True).
        - **addref**,optional: Add the reference arc and point to the plot.
        - **normalize**, optional: Force the normalization of statistics (standard
          deviations).
        - **negative**,optional: Force plot of negative correlations.
        - **legend**, optional: Add the legend.
        - **size**,optional: Size(s) of the markers in points.
        - **colors**,optional: Color of the points. It can be a list to set
          a different color for each pair.
        - **cmap**, optional: Colormap use when ``colors`` is a numeric array.
        - **label_colors**,optional: Color of the labels.
          If ``"same"``, color is taken from ``colors``.
        - **label_<keyword>**,optional: ``<keyword>`` is passed
          to :func:`~matplotlib.pyplot.text` for plotting the labels.
        - **symbols**,optional: Symbols for the points (default to ``"o"``).
        - **point_<keyword>**,optional: ``<keyword>`` is passed to
          :func:`~matplotlib.pyplot.plot` (:func:`~matplotlib.pyplot.scatter` if
          ``scatter`` is True) or for plotting the points.
        - **reflabel**,optional: Label of the reference.
        - **refcolor**,optional: Color of the reference.
        - **title**,optional: A string for the general title.
        - **scatter**, optional: Use :func:`matplotlib.pyplot.scatter` to improve
          plot efficiency when the number of points is huge. Labels are not plotted in
          this mode.
        - **savefig**,optional: Save figure to ``savefig``.
        - **savefig_<keyword>**,optional: ``<keyword>`` is used for saving figure.
        - **savefigs**,optional: Save figure to ``savefig`` using :func:`~vacumm.misc.plot.savefigs`.
        - **savefigs_<keyword>**,optional: ``<keyword>`` is used for saving figure using :func:`~vacumm.misc.plot.savefigs`."""
    from matplotlib.transforms import offset_copy

    # Normalization
    np = len(stats)
    scatter = kwargs.get('singleref', scatter)
    if stdref is None:
        normalize=True
        stdref = N.ma.array([1.])
    stdref = N.ma.asarray(stdref)
    if stdref.ndim==0: stdref.shape = 1,
    if stdref.shape[0]>1 and not N.ma.allclose(stdref, stdref.mean()):
        normalize = True
    stats = N.ma.asarray(stats)
    if normalize:
        stats[:, 0] /= stdref
        stdref = N.ma.array([1])

    # Arrays
    if stdmax is None:
        stdmax = max(N.ma.abs(stats[:, 0]).max(), stdref.max())
    corrmin = stats[:, 1].min()
    if negative is None:
        negative = corrmin < 0
    sticks = auto_scale([0, stdmax], vmin=0.)
    stdmax = sticks[-1]
    stdmin = -stdmax if negative else 0
    sticksref = N.arange(0., stdmax*N.sqrt(2.), sticks[1]-sticks[0])
    rticks = N.arange(0, 1., .1).tolist()+[.95, .99, 1.]
    if negative:
        rticks = (-N.asarray(rticks[:0:-1])).tolist()+rticks

    # Check options
    labels = broadcast(labels, np, fillvalue=None)
    cmap = get_cmap(cmap)
    if colors is None:
        colors =  kwargs.pop('color', 'k')
    if not isinstance(colors, N.ndarray):
        colors = broadcast(colors,np)
        if any([isinstance(c,(int,float)) for c in colors]): # list of scalars
            colors = N.array(colors)
    if isinstance(colors, N.ndarray):
        colors = N.resize(colors, (np,))
        if not scatter:
            colors = Scalar2RGB(colors, cmap)(colors)
    for i in xrange(np):
        if colors[i] is None:
            colors[i] = 'k'
    sizes = N.asarray(broadcast(kwargs.pop('size', 8), np))
    alphas = broadcast(kwargs.pop('alpha', 1), np)
    label_size = broadcast(kwargs.pop('label_size', P.rcParams['font.size']), np)
    label_colors = kwargs.pop('label_colors', 'k')
    if label_colors=='same':
        label_colors = colors
    else:
        label_colors = broadcast(colors, np, fillvalue='k')
    symbols = broadcast(kwargs.pop('symbols', 'o'), np)
    if negative and autoresize:
        x0,y0 = P.gcf().get_size_inches()
        x1 = x0*2.
        y1 = y0
        if autoresize != 2:
            x1 /= N.sqrt(2.)
            y1 /= N.sqrt(2.)
        P.gcf().set_size_inches((x1,y1),forward=True)

    # Conversion function
    def xy(std, corr):
        return std*corr, std*N.sin(N.arccos(corr))

    # Start the plot
    if not minimal: _start_plot_(**kwargs)
    fig = P.gcf()
    fig.set_facecolor('w')
    ax = P.gca()

    # Points
    kwpt = kwfilter(kwargs, 'point')
    kwlab = kwfilter(kwargs, 'label')
    if scatter:

        # All points
        x,y = xy(stats[:,0],stats[:,1])
        ppts = P.scatter(x, y, 0.75*sizes**2, colors, marker=symbols[0],
            zorder=zorder, alpha=alphas[0], cmap=cmap, **kwpt)

    else:

        ppts = []
        for i, (std, corr) in enumerate(stats):

            # Single point
            x, y = xy(std, corr)
            ppts.append(P.plot([x], [y], symbols[i], label=labels[i],
                zorder=12, color=colors[i],
                markersize=sizes[i], alpha=alphas[i],**kwpt))

            # Labels
            if isinstance(labels[i], basestring):
                transoffset = offset_copy(ax.transData, fig=fig,
                    x = 0.0, y=-1.2*sizes[i], units='points')
                P.text(x, y, labels[i], va='center', ha='center', transform=transoffset,
                    zorder=zorder-2, size=label_size[i], **kwlab)
    del labels

    # Decorations for normal plots
    if not minimal:


        # Boundary standard deviation arc
        arc = Circle((0., 0.), stdmax, fill=False, zorder=0, facecolor='none')
        ax.add_patch(arc)

        # Intermediates standard deviation arcs
        for radius in sticks[1:-1]:
            pstd = Circle((0., 0.), radius, linestyle='dashed',
                linewidth=P.rcParams['grid.linewidth'], fill=False, edgecolor='b',zorder=0)
            ax.add_patch(pstd)

        # RMS (red) for first ref only
        rmax = N.sqrt(stdmax**2+stdref[0]**2)
        for r in sticksref[1:-1]:
            prms = Circle((stdref[0], 0.), r,
                linestyle='dashed', fill=False, linewidth=P.rcParams['grid.linewidth'],
                edgecolor='r',zorder=0)
            ax.add_patch(prms)
            prms.set_clip_on(True)
            prms.set_clip_path(arc)
            xt = stdref[0]*r/rmax
            yt = stdmax*r/rmax
            P.text(stdref[0]-xt, yt, ('%g'%r).ljust(4), color='r',
                rotation=N.arctan2(xt, yt)*180./N.pi,
                va='center', ha='center', clip_on=True,
                bbox=dict(facecolor='w', linewidth=0, alpha=.5))

    # Reference
    addref = kwargs.get('ref', addref)
    if addref:

        # First Arc
        ax.add_patch(Arc((0., 0.), stdref[0]*2, stdref[0]*2, fill=False,
            edgecolor=refcolor, linewidth=2, zorder=0))

        # All points unless normalized
        imax = 1 if normalize else None
        rslice = slice(0,imax)
        pref = P.plot(stdref[rslice], [0]*len(stdref[rslice]), 'o', color=refcolor,
            markersize=10, clip_on=False)

    if not minimal:

        # Ticks for correlation
        for r in rticks:
            x1, y1 = xy(stdmax, r)
            x0, y0 = xy(stdmax*.98, r)
            dxt = 12*x1/stdmax
            dyt = 12*y1/stdmax
            P.plot([x0, x1], [y0, y1], 'k-', label=False)
            pcor = P.plot([0, x1], [0, y1], 'k--', linewidth=P.rcParams['grid.linewidth'], color='g')
            transoffset = offset_copy(ax.transData, fig=fig, x=dxt, y=dyt, units='points')
            P.text(x1, y1, ('%g'%r).ljust(4), va='center', ha='center',
                transform = transoffset,
                rotation=N.arccos(r)*180./N.pi + (r<0)*180.)
        xt = yt = stdmax/N.sqrt(2.)
        dxt = 29*xt/stdmax
        dyt = 29*yt/stdmax
        transoffset = offset_copy(ax.transData, fig=fig, x=dxt, y=dyt, units='points')
        P.text(xt, yt, 'Correlation', va='center', ha='center',size=P.rcParams['axes.labelsize'],
                transform = transoffset, rotation=-45.)

        # Legend
        if legend:
            legstd = 'Standard dev.'
            legrms = 'RMS error'
            if normalize:
                legstd = 'Rel. '+legstd.lower()
                legrms = 'Rel. '+legrms
            hdls = (pcor[0], pstd, prms)
            labs = ('Correlation', legstd, legrms)
            if addref:
                hdls += pref[0],
                labs += reflabel,
            P.legend(hdls, labs,
                loc='upper right' if not negative else 'upper left',
                #markerscale=10,
                numpoints=1, **kwfilter(kwargs, 'legend_',
                    defaults=dict(shadow=False, fancybox=True))).legendPatch.set_alpha(.5)

        # Final setup
        P.axis([stdmin, stdmax, 0, stdmax])
        ax.set_aspect(1.)
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        if negative:
            ax.spines['left'].set_position('center')
        #ax.set_frame_on(False)
        if negative:
            P.xticks(N.concatenate((-sticks[:0:-1], sticks)))
        else:
            P.xticks(sticks)
        P.yticks(sticks)
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        P.axhline(0, color='k', label=False).set_clip_on(False)
        P.axvline(0, color='k', label=False).set_clip_on(False)
        stdlabel = 'Standard deviation'
        if normalize: stdlabel = 'Relative '+stdlabel.lower()
        if not normalize and units is not None and units is not False:
            stdlabel += ' [%s]'%units
        P.xlabel(stdlabel)
        if not negative:
            P.ylabel(stdlabel)

        # End
        kwargs.update(autoresize=0, units=None, grid=False)
        kwargs.setdefault('title', 'Taylor diagram')
        kwargs.setdefault('title_y', 1.05)
        kwargs.setdefault('savefigs_pdf', True)
        _end_plot_(**kwargs)

    else:

        P.ylim(ymin=0)
        ax.set_aspect(1.)
        ax.set_frame_on(False)
    return ppts



def dtarget(bias, rmsc, stdmod, stdref=None, colors='cyan', sizes=20,
    title='Target diagram', cmap=None, units=None, colorbar=True, **kwargs):
    """Plot a (direct) target diagram from already computed statistics

    :Params:

        - **bias** : Bias.
        - **rmsc**: Centered RMS error.
        - **stdmod**: Standard deviation of the model.
        - **stdref**, optional: Standard deviation of the observations.
        - **colors**, optional: Single or list of colors, or
          array of values, used to color markers.
        - **sizes**, optional: Single or list of sizes of markers.
        - **title**, optional: Title of the plot.
        - **scatter_<keyword>**, optional: ``<keyword>`` is passed
          to :func:`matplotlib.pyplot.scatter`.
        - **circle_<keyword>**, optional: ``<keyword>`` is passed
          change the properties of circles.

    .. todo::

        Make two versions of :func:`target`, as for Taylor
        diagrams: one like with direct statistics, one with
        :mod:`MV2` variables as arguments.
    """

    kwscat = kwfilter(kwargs, 'scatter_')
    kwscat.setdefault('cmap', cmap)
    kwcb =kwfilter(kwargs, 'colorbar_')
    kwcir = kwfilter(kwargs, 'circle_')

    # Convert to masked arrays
    bias = N.ma.asarray(bias)
    if stdref is None: stdref = 1.
    stdref = N.resize(N.ma.asarray(stdref), bias.shape)
    stdmod = N.ma.asarray(stdmod)
    rmsc = N.ma.asarray(rmsc)
    np = len(bias)

    # Make plotted stats
    nrmsc = rmsc / stdref
    nrmsc *= N.sign(stdmod-stdref)
    nbias = bias / stdref

    # Broadcast plot params
    if colors=='rms':
        colors = N.sqrt(rmsc**2+bias**2)
    elif isinstance(colors, N.ndarray):
        colors = N.resize(colors, bias.shape)
    else:
        if not isinstance(colors, list):
            colors = [colors]
        colors = broadcast(colors, np)
    sizes = N.resize(N.asarray(sizes), bias.shape)

    # Start the plot
    _start_plot_(**kwargs)
    ax = P.gca()
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_aspect(1.)

    # Plot
    p = ax.scatter(nrmsc, nbias, sizes, colors, **kwscat)
    xmin, xmax, ymin, ymax = P.axis()
    xmax = max(-xmin, xmax, -ymin, ymax, 1.1)
    P.axis([-xmax, xmax, -xmax, xmax])

    # Arcs
    dict_check_defaults(kwcir, linestyle='dotted', edgecolor='r',
            linewidth=P.rcParams['grid.linewidth'], fill=False, zorder=0)
    for radius in P.xticks()[0]:
        if radius > 0:
            pstd = Circle((0., 0.), radius, **kwcir)
            ax.add_patch(pstd)
    kwcir['linestyle'] = 'solid'
    ax.add_patch(Circle((0., 0.), 1., **kwcir))

    # Colorbar
    if colorbar and isinstance(colors, N.ndarray):
        label = 'RMS'
        if units: label += ' [%s]'%units
        kwcb.setdefault('label', label)
        kwcb.setdefault('pad', 0.15)
        kwcb.setdefault('fraction', 0.04)
        P.colorbar(p, **kwcb)

    # Finish plot
    P.xlabel(r'$RMSC/\sigma_{obs}$', ha='left', va='bottom', position=(1., 0.))
    P.ylabel('$Bias/\sigma_{obs}$', ha='left', va='top', rotation=0, position=(0., 0))
    kwargs['title'] = title
    kwargs.update(autoresize=0, units=None, grid=False)
    _end_plot_(**kwargs)
    return p

target = dtarget


def rankhist(obs, ens, title='Rank histogram', bins=None, **kwargs):
    """Plot rank histogram

    The rank is computed using :func:`~vacumm.misc.stats.ensrank`.
    It may vary from ``0`` to ``nens``.

    - **obs**, (...): Observations or rank
    - **ens**, (nens,...): Ensemble or nens

    """
    from stats import ensrank
    if isinstance(ens, int):
        rank, nens = obs, ens
    else:
        ens = MV2.asarray(ens)
        rank = ensrank(obs, ens)
        nens = ens.shape[0]
    kwhist = kwfilter(kwargs, 'hist_')
    if N.ma.isMA(rank): rank = rank.compressed()
    _start_plot_(**kwargs)
    if bins is None: bins = N.arange(-.5, nens+1.5)
    pdf, bins, patches = P.hist(rank, bins=bins,
        normed=True, align='mid', **kwhist)
    P.axhline(1./(nens+1), color='r', lw=2)
    P.xlim(bins.min(), bins.max())
    P.xlabel('Rank')
    P.ylabel('P(rank)')
    kwargs.setdefault('savefigs_pdf', True)
    _end_plot_(title=title, **kwargs)
    return patches


def ellipsis(xpos,ypos,eAXIS,eaxis=None,rotation=0,
    sign=0,np=100,gobj=None,m=None, scale=1.,fill=False,sign_usearrow=True,
    axes=True,sign_usecolor=False,sign_usewidth=False,units=None,
    key=True,key_align='top',key_orientation='tangential',
    glow=False, shadow=False,
    key_position=1.15,key_format='%(value)g %(units)s',**kwargs):
    """Create a dicretized ellipsis

    - **xpos**: X position of the ellipsis.
    - **ypos**: Y    //
    - **eAXIS**: Length of grand axis
    - *eaxis*: Length of little axis [default: eAXIS]
    - *rotation*: Rotation angle in degrees of grand axis from X axis [default: 0.]
    - *sign*: Indicate a direction of evolution along the ellipsis using a variable line width and/or color [default: 0]. If 0, nothing is done. If positive or negative, it defines the direction relative the trigonometric convention.
    - *sign_usearrow*: Use an arrow to indicate sign [default: True]
    - *sign_usewidth*: Use width variation to indicate sign [default: False]
    - *sign_usecolor*: Use color variation (from grey to default color) to indicate sign [default: False]
    - *scale*: Scale factor to apply to grand and little axes [default: 1.]
    - *np*: Number of discretized segments [default: 100]
    - *gobj*: Graphical object on which to plot the ellipsis (can be a Basemap object) [default: ~matplotlib.pyplot.gca()]. If map object, eAXIS and eaxis are supposed to be in kilometers.
    """

    # Angles for discretization
    pp = N.arange(0,2*N.pi,2*N.pi/np)
    pp = N.concatenate((pp,pp[0:1])) # Closed
    np = len(pp)

    # Axes
    eAXIS = eA_orig = N.absolute(eAXIS)
    if eaxis is None: eaxis = eAXIS
    eaxis = ea_orig = N.absolute(eaxis)
    eaxis *= scale
    eAXIS *= scale

    # Map
    if m: gobj = m
    from core_plot import Map
    if gobj is None:
        gobj = P.gca()
    elif isinstance(gobj, Map):
        gobj = m.map
    ismap = isinstance(gobj,Basemap)
    if ismap:
        iscyl = gobj.projection == 'cyl'
        eaxis *= 1.e3 # km
        eAXIS *= 1.e3 # km

    # Discretization
    xx = N.zeros(np).astype('d')
    yy = N.zeros(np).astype('d')
    rotation *=  N.pi / 180.
    for ip,ap in enumerate(pp):
        # Basic ellipsis
        xxx = eAXIS * N.cos(ap) * 0.5
        yyy = eaxis * N.sin(ap) * 0.5
        # Rotation
        xx[ip] = xxx * N.cos(rotation) - yyy * N.sin(rotation)
        yy[ip] = xxx * N.sin(rotation) + yyy * N.cos(rotation)
        # Cylindrical map: meters to degrees
        if ismap and iscyl:
            xx[ip] = m2deg(xx[ip],ypos)
            yy[ip] = m2deg(yy[ip])

    # Add position
    if ismap: xpos,ypos = gobj(xpos,ypos)
    xx = xx + xpos
    yy = yy + ypos

    # Add axes
    if axes:
        kwaxes = kwfilter(kwargs,'axes',
            defaults=dict(color='k',linestyle=':',alpha=.5))
        for ea,ap in (eAXIS,0.),(eaxis,N.pi/2):
            xaxis = N.array([-1,1]) * N.cos(rotation+ap) * ea/2
            yaxis = N.array([-1,1]) * N.sin(rotation+ap) * ea/2
            if ismap and iscyl:
                xaxis = m2deg(xaxis,ypos)
                yaxis = m2deg(yaxis)
            gobj.plot(xaxis+xpos,yaxis+ypos,**kwaxes)

    # Plot command
    if fill:
        pcmd = gobj.fill
    else:
        pcmd = gobj.plot

    # Plot with sign indicator
    myplot = None
    kwkey = kwfilter(kwargs,'key')
    if int(sign) != 0 :
        sign = float(sign)
        if not sign_usearrow and not fill: # FIXME: color and fill for ellipsis
            linewidth = lw = kwargs.pop('linewidth',P.rcParams['lines.linewidth'])
            color = cl = ColorConverter().to_rgb(kwargs.pop('color',P.rcParams['lines.color']))
            if sign > 0:
                pps = range(0,np-1)
            else:
                pps = range(np-2,-1,-1)
            factor = N.absolute(sign+1.)*3
            for ip in pps:
                if sign_usewidth:
                    lw = linewidth / factor + factor * linewidth * float(ip)/np
                if sign_usecolor:
                    cl = (  (color[0]*ip+(np-ip)*.5)/(np-1),
                            (color[1]*ip+(np-ip)*.5)/(np-1),
                            (color[2]*ip+(np-ip)*.5)/(np-1))
                myplot = pcmd([xx[ip],xx[ip+1]],[yy[ip],yy[ip+1]],
                    linewidth=lw,color=cl,**kwargs)
        else:
            xshaft = N.zeros(3) ; yshaft = N.zeros(3)
            slen = sign*eAXIS/4.
            ap = pp[int(round(np*3./4.))]
            xxx = eAXIS * N.cos(ap) * 0.5
            yyy = eaxis * N.sin(ap) * 0.5
            xshaft[1] = xxx * N.cos(rotation) - yyy * N.sin(rotation)
            yshaft[1] = xxx * N.sin(rotation) + yyy * N.cos(rotation)
            for ia,aa in enumerate([1,-1]):
                aa *= N.pi*.85
                xshaft[ia*2] = xshaft[1] + slen * N.cos(aa+rotation)
                yshaft[ia*2] = yshaft[1] + slen * N.sin(aa+rotation)
            if ismap and iscyl:
                xshaft = m2deg(xshaft,ypos)
                yshaft = m2deg(yshaft)
            kws = kwargs.copy()
            kws['label'] = None
            gobj.plot(xshaft+xpos,yshaft+ypos,**kws)
    if myplot is None:
        myplot = pcmd(xx,yy,**kwargs)
        if glow:
            add_glow(myplot)
        if shadow:
            add_shadow(myplot)

    # Plot key
    if key:

        # Position
        rotation = rotation % N.pi
        key_x = key_position * N.cos(rotation) * eAXIS/2
        key_y = key_position * N.sin(rotation) * eaxis/2
        ha = 'center' ; va = 'center'
        if key_align != 'top':
            key_x *= -1
            key_y *= -1
        if ismap and iscyl:
            key_x = m2deg(key_x,ypos)
            key_y = m2deg(key_y)
        key_x += xpos
        key_y += ypos
        key_rotation = rotation
        if key_orientation == 'tangential':
            key_rotation -= N.pi/2

        # Params
        for key, val in dict(color=myplot[0].get_color(),alpha=myplot[0].get_alpha(),
            rotation=key_rotation*180./N.pi,va=va,ha=ha).items():
            kwkey.setdefault(key, val)

        # Content
        value = eA_orig/2
        if units is None: units = ''
        key_text = key_format % vars()

        # Draw
        kglow = kwkey.pop('glow', glow)
        kwkglow = kwfilter(kwkey, 'glow')
        kshadow = kwkey.pop('shadow', shadow)
        kwkshadow = kwfilter(kwkey, 'shadow')
        tt = P.text(key_x,key_y,key_text,**kwkey)
        if kglow:
            add_glow(tt, **kwkglow)
        if kshadow:
            add_shadow(tt, **kwkshadow)



    # Normal plot
    return myplot


def add_colorbar(var=None,sm=None,horizontal=False,position=None,cmap=None,drawedges=False,axes=None,fig=None,figsize=None,frameon=None,levels=None,barwidth=None,barmargin=.1, norm=None, **kwargs):
    """Standalone colorbar"""

    if fig is None:
        if figsize is None:
            if horizontal:
                figsize = (8,.5)
            else:
                figsize = (.8,6)
        fig = P.figure(figsize = figsize)
    if axes is None:
        #axes = fig.gca()
        if horizontal:
            if barwidth is None: barwidth = .25
            rect = [.05,1-barmargin-barwidth,.9,barwidth]
        else:
            if barwidth is None: barwidth = .18
            rect = [barmargin,.05,barwidth,.92]
        axes = fig.add_axes(rect)
    else:
        axes = fig.add_axes(axes)
    if frameon is not None:
        fig.set_frameon(frameon)
        axes.set_frameon(frameon)

    if levels is not None:
        vmin,vmax = minmax(levels)
    elif var is not None:
        levels,vmin,vmax,tmp = _levels_(var,**kwargs)

    if sm is None and levels is not None:
        sm = ScalarMappable(cmap=cmap)
        sm.set_array(N.array(levels))

    if sm is None:
        raise 'Missing explicit arguments sm'
    for key,val in kwargs.items():
        if key not in ('close','show','savefig','units'):
            kwargs['colorbar_'+key] = kwargs.pop(key)
    cb = colorbar(sm,vars=var,drawedges=drawedges,levels=levels,colorbar_horizontal=horizontal,
        colorbar_position=position,standalone=True,cax=axes,**kwargs)

    kwargs['grid'] = False
    kwargs['title'] = False
    _end_plot_(**kwargs)
    return cb

def scolorbar(*args, **kwargs):
    """Shortcut to :func:`add_colorbar`"""
    return add_colorbar(*args, **kwargs)

def _plot_grid_(ax, xx, yy, samp, alpha, zorder, labels, **kwargs):
        lines = ()
        ny, nx = xx.shape
        for iy in xrange(0, ny, samp):
            lines += zip(xx[iy, ::samp], yy[iy, ::samp]),
        for ix in xrange(0, nx, samp):
            lines += zip(xx[::samp, ix], yy[::samp, ix]),
        oo = LineCollection(lines, **kwargs)
        oo.set_alpha(alpha)
        oo.set_zorder(zorder)
        oo.set_label(labels)
        from mpl_toolkits.mplot3d import Axes3D
        if isinstance(ax, Axes3D):
            ax.add_collection3d(oo)
        else:
            ax.add_collection(oo)


def add_grid(gg, color='k', edges=True, centers=False, m=None, linecolor=None,
        xcorners=None, ycorners=None,
        linewidth=1, alpha=.5, linestyle='solid', zorder=100, marker=',',
        markersize=4, facecolor=None, edgecolor=None, markerlinewidth=0,
        label_edges='Cell edges', label_centers='Cell centers', samp=1,
        ax=None, **kwargs):
    """Add a a 1D or 2D grid to a plot

    Params:
        - **gg**: A cdms grid OR a (lon,lat) tuple of axes OR a cdms variable with a grid
        - **x/ycorners**: Corners coordinates (ny+1,nx+1)
        - *borders*: Display cell borders
        - *centers*: Display cell centers. If == 2, connect centers.
        - *samp*: Undersampling rate
        - *m*: A basemap instance to plot on
        - *color*: Default color (or *c*)
        - *linecolor*: Color of the lines (or *lc*) [defaults to *color*]
        - *linestyle*: Line style (or *ls*)
        - *edgecolor*: Edge color of marker (or *ec*) [defaults to *color*]
        - *facecolor*: Face color of marker (or *fc*) [defaults to *color*]
        - *marker*: Type of marker
        - *markersize*: Size of markers (or *s*)
        - *alpha*: Alpha transparency (or *a*)
        - *zorder*: Order of the grid in the stack (100 should be above all objects)
    """
    # Axes
    xx, yy = get_xy(gg, mesh=True, num=True)

    # Edges
    edges = kwargs.pop('borders', edges)
    if edges:
        if xcorners is None or ycorners is None:
            xx2d_, yy2d_ = meshbounds(xx, yy)
        xx2d = xcorners if xcorners is not None else xx2d_
        yy2d = ycorners if ycorners is not None else yy2d_
    if centers:
        xx, yy = meshgrid(xx, yy)

    # Geographic map
#    if m is None and gobjs('m'): m = gobjs('m')[-1]
    if hasattr(m, 'map'):
        if ax is None: ax = m.axes
        m = m.map
    if isinstance(m, Basemap):
        if centers:
            xx, yy = m(xx, yy)
        if edges:
            xx2d, yy2d = m(xx2d, yy2d)
    if ax is None: ax = P.gca()

    # Style
    color = kwargs.pop('c', color)
    if linecolor is None: linecolor = color
    if edgecolor is None: edgecolor = color
    if facecolor is None: facecolor = color
    linestyle = kwargs.pop('ls', linestyle)
    linecolor = kwargs.pop('lc', linecolor)
    linewidth = kwargs.pop('lw', linewidth)
    alpha = kwargs.pop('a', alpha)
    markersize = kwargs.pop('s', markersize)
    kwlines = {'colors':[linecolor], 'linestyles':linestyle, 'linewidths':linewidth}

    # Centers
    if centers:

        # Markers
        if centers>0:
            kwdef = {}
            kwdef['markerfacecolor'] = facecolor
            kwdef['markeredgecolor'] = edgecolor
            kwdef['linewidth'] = markerlinewidth
            kwdef['label'] = label_centers
            kwdef['alpha'] = alpha
            centers = ax.plot(xx[::samp, ::samp].ravel(), yy[::samp, ::samp].ravel(), marker, markersize=markersize,
                **kwfilter(kwargs, 'center', defaults=kwdef))

        # Lines
        if centers==2 or centers<0:
            kwpg = kwfilter(kwargs, 'center', defaults=kwlines)
            edges = _plot_grid_(ax, xx, yy, samp, alpha, zorder, label_edges, **kwpg)

    else:
        centers=None

    # Create lines
    if edges:
        kwpg = kwfilter(kwargs, 'edge', defaults=kwlines)
        edges = _plot_grid_(ax, xx2d, yy2d, samp, alpha, zorder, label_edges, **kwpg)

#        lines = ()
#        ny, nx = gg.shape
#        for iy in xrange(0, ny+1, samp):
#            lines += zip(xx2d[iy, ::samp], yy2d[iy, ::samp]),
#        for ix in xrange(0, nx+1, samp):
#            lines += zip(xx2d[::samp, ix], yy2d[::samp, ix]),
#        kwdef = {'colors':[linecolor], 'linestyles':linestyle, 'linewidths':linewidth}
#        edges = LineCollection(lines, **kwfilter(kwargs, 'edge', defaults=kwdef))
#        edges.set_alpha(alpha)
#        edges.set_zorder(zorder)
#        edges.set_label(label_edges)
#        ax.add_collection(edges)
#        gobjs(grid_edges=edges)
    else:
        edges=None

    return centers, edges



####################
# Plot alterations #
####################

#def cached_map(*args):
#   """Check if we have a cached map
#
#   @usage:
#   >>> m = cached_map(lon_min, lon_max, lat_min, lat_max, projection, resolution)
#   >>> m = cached_map(m) # Does only caching of map
#   """
#   # We already have one in live memory
#   if isinstance(args[0], Basemap):
#       # Save it
#       cache_map(args[0])
#       # Get it
#       return args[0]
#   # Guess
#   file = _cached_map_file_(*args)
##  print 'Checking', file, os.path.exists(file)
#   if not os.path.exists(file): return None
##  print 'Loadind cached map from '+os.path.basename(file)
#   f = open(file)
#   m = cPickle.load(f)
#   f.close()
#   return m
#
#def cache_map(m):
#   """Cache a map if still not cached"""
#   if m is None: return
#   file = _cached_map_file_(m.llcrnrlon, m.urcrnrlon, m.llcrnrlat, m.urcrnrlat, m.projection, m.resolution)
#   if not os.path.exists(file):
##      print 'Caching map to '+os.path.basename(file)
#       f = open(file, 'wb')
#       cPickle.dump(m, f)
#       f.close()
#
#def _cached_map_file_(lon_min, lon_max, lat_min, lat_max, proj, res):
#   return os.path.join(_cached_map_dir,
#       'cached_map.%(lon_min)g_%(lon_max)g.%(lat_min)g_%(lat_max)g.%(proj)s.%(res)s.pyk' % vars())
#


def rotate_tick_labels(angle,vertical=0,*args,**kwargs):
    """Rotate labels along an axis

    - **angle**: Angle in degrees OR False.
    - *vertical*: On vertical axis if not 0 [default: 0]
    - Other arguments are passed to :func:`~matplotlib.pyplot.setp`
    """
    ax = kwargs.get('ax',P.gca())
    func = eval('ax.get_%sticklabels()' % 'xy'[vertical])
    if angle is False: angle = 0
    P.setp(ax.get_xticklabels(),"rotation",angle,*args)

def xrotate(angle,*args,**kwargs):
    """Rotate xticklabels

    See: :func:`rotate_tick_labels`"""
    rotate_tick_labels(angle,1,*args,**kwargs)
def yrotate(angle,*args,**kwargs):
    """Rotate yticklabels

    See: :func:`rotate_tick_labels`"""
    rotate_tick_labels(angle,0,*args,**kwargs)

def rotate_xlabels(angle,*args,**kwargs):
    """Shortcut to :func:`xrotate`"""
    rotate_tick_labels(angle,1,*args,**kwargs)
def rotate_ylabels(angle,*args,**kwargs):
    """Shortcut to :func:`yrotate`"""
    rotate_tick_labels(angle,0,*args,**kwargs)


def _xyscale_(xylim,factor,keep_min=False,keep_max=False,**kwargs):
    if factor <=0.: factor = 1.
    if factor == 1. or (keep_min and keep_max): xylim()
    mn,mx = xylim()
    delta = (mx-mn)*(factor-1.)
    mn -= delta * 0.5 * (1 + keep_max) * (1 - keep_min)
    mx += delta * 0.5 * (1 + keep_min) * (1 - keep_max)
    return xylim(minmax(auto_scale([mn,mx])))


def xscale(*args,**kwargs):
    """Scale xlim using factor and auto_scale

    - **factor**: Relative factor with 1. meaning no change
    - *keep_min/keep_max*: Kee min/max during scaling [default: False]

    Return: New xlim()
    """
    return _xyscale_(P.xlim,*args,**kwargs)

def scale_xlim(*args,**kwargs):
    """Shortcut to :func:`xscale`"""
    xscale(*args,**kwargs)

def yscale(*args,**kwargs):
    """Scale ylim using factor and auto_scale

    - **factor**: Relative factor with 1. meaning no change
    - *keep_min/keep_max*: Kee min/max during scaling [default: False]

    Return: New ylim()
    """
    return _xyscale_(P.ylim,*args,**kwargs)

def scale_ylim(*args,**kwargs):
    """Shortcut to :func:`yscale`"""
    yscale(*args,**kwargs)


def xhide(choice=True,**kwargs):
    """Hide or not an axis

    - **choice**: What to do [default: True]
    """
    ax = kwargs.get('ax',P.gca())
    P.setp(ax.get_xticklabels(),"visible",not choice)

def yhide(choice=True,**kwargs):
    """Hide or not an axis

    - **choice**: What to do [default: True]
    """
    ax = kwargs.get('ax',P.gca())
    P.setp(ax.get_yticklabels(),"visible",not choice)

def xdate(rotation=45.,**kwargs):
    """Consider x axis as dates

    $_xydate_
    """
    _xydate_('x',rotation=rotation,**kwargs)


def ydate(rotation=0.,**kwargs):
    """Consider y axis as dates

    $_xydate_
    """
    _xydate_('y',rotation=rotation,**kwargs)

def _xydateold_(xy, tz=None, auto=True, fmt=None, rotation=None,
    locator=None, minor_locator=None, nominor=False,
    nmax_ticks=None, intv=None, trange=None, **kwargs):
    """
    - *tz*: Time zone.
    - *auto*: Auto Scaling [default: True]
    - *rotation*: Rotation angle of tick labels. If None, automatic [default: None]
    - *fmt*: Date format.
    - *locator*: Major locator. Can be within ['year','month','Weekday','day','hour','minute','second'] or be like matplotlib.dates.MonthLocator().
    - *minor_locator*: Minor locator.
    - *nominor*: Do not try to add minor ticks [default: False]
    - *nmax_ticks*: Maximal number of ticks
    - *locator_<keyword>*: <keyword> is passed to locator if locator is a string. If locator = 'month', locator = MonthLocator(locator_<keyword>=<value>).
    - *minor_locator_<keyword>*: Same with minor_locator.
    """
    ax = kwargs.get('ax',P.gca())
##  print 'xlimssss',P.gca(),ax
##  if ax.get_xlim() == (0.,1.): raise 'a'
    # Base
    exec "ax.%saxis_date(tz)" % xy

    # Scale axes
    if auto:
        ax.autoscale_view(**{'scale'+xy:True,'scale'+{'x':'y','y':'x'}[xy]:False})

    # Rotates labels
    if rotation is not None:
        eval("rotate_%slabels" % xy)(rotation,**kwargs)

    # String case for locators
    auto_major = locator is None
    major_locator = locator
    locs = ['year','month','weekday','day','hour','minute','second']
    locs.extend([(loc+'s') for loc in locs])
    kwmjl = kwfilter(kwargs,'locator')
    kwmnl = kwfilter(kwargs,'minor_locator')
    if isinstance(major_locator,str):
        if major_locator.lower() in locs:
            major_locator = eval(major_locator.lower().title()+'Locator(**kwmjl)')
        else:
            major_locator = None
    if isinstance(minor_locator,str):
        if minor_locator.lower() in locs:
            minor_locator = eval(minor_locator.lower().title()+'Locator(**kwmnl)')
        else:
            minor_locator = None

    # Tick locations
    guess_major = major_locator is None
    if nmax_ticks is None: nmax_ticks = 24
    if major_locator is None:
        major_locator = get_major_locator(xy,**kwargs)
        nticks = len(eval("ax.get_%sticks()" % xy))
    if nmax_ticks > 0 and not isinstance(major_locator,AutoDateLocator):
        # Major ticks
        locname = str(major_locator).split()[0].split('.')[2]
        nticks = nmax_ticks+1
        i = 1
#       mlr = major_locator.rule._rrule
        while nticks > nmax_ticks:
            try:
                major_locator = set_major_locator(xy,eval("%s(interval=%i,**kwmjl)"%(locname,i)),**kwargs)
            except:
                major_locator = set_major_locator(xy,eval("%s(%i,**kwmjl)"%(locname,i)),**kwargs)
            nticks = len(eval("ax.get_%sticks()" % xy))
            i += 1
    nmax_ticks = abs(nmax_ticks)
    if isinstance(major_locator,AutoDateLocator):
        if N.isfinite(major_locator.axis.get_data_interval()[0]):
            lims = major_locator.datalim_to_dt()
        else:
            lims = major_locator.viewlim_to_dt()
        major_locator = major_locator.get_locator(*lims)

    # Fix phase of some locators
    try:
        sampling = major_locator.rule._rrule._interval # sampling interval
    except:
        sampling = major_locator.base._base
    try:
        freq = major_locator.rule._rrule._freq
    except:
        freq = YEARLY
    if freq == DAILY and sampling == 7:
        # Fix daily/7 to weekly/monday locator (start on monday)
        major_locator = set_major_locator(xy,WeekdayLocator(MO),**kwargs)
    elif freq == HOURLY and sampling > 1:
        # Fix hourly locator to start at mindnight
        good = [2, 3, 4, 6, 8, 12]
        sampling = good[min(N.searchsorted(good, sampling), len(good)-1)]
        byhour = range(0, 24, sampling)
        major_locator = set_major_locator(xy,HourLocator(byhour=byhour),**kwargs)
    major_unit = major_locator._get_unit()

    # Tick format
    if fmt is None:
        # Try a special date formatter
        if trange < 1 and trange > 1/24.: # H/M
            fmt = [3]
        elif trange < 31 and trange > 1: # d/m
            fmt = [2]
        elif trange < 365 and trange > 31: # m/d
            fmt = [1]
        elif trange < 365*5 and trange > 365: # Y/m
                fmt = [0]
    if isinstance(fmt,str):
        fmt = DateFormatter(fmt)
    elif isinstance(fmt,list):
        fmt = SpecialDateFormatter(*fmt)
    elif isinstance(fmt,dict):
        fmt = SpecialDateFormatter(**fmt)
    elif isinstance(fmt,tuple):
        if len(fmt) == 1:
            fmt += ({}, )
        elif not isinstance(fmt[1], dict):
            fmt[1] = {}
        fmt = SpecialDateFormatter(fmt[0], **fmt[1])
    if isinstance(fmt,SpecialDateFormatter) :
        eval("rotate_%slabels" % xy)(0,**kwargs)
    if fmt is not None:
        eval("ax.%saxis.set_major_formatter" % xy)(fmt)

    # Minor ticks
    if not nominor:
#       # Minor formatter for special case
#       if fmt is not None and isinstance(fmt,SpecialDateFormatter):
#           eval("ax.%saxis.set_minor_formatter" % xy)(fmt)
        # Direct
        if minor_locator is not None:
            set_minor_locator(xy,minor_locator,**kwargs)
            return
        # From sampling
        if sampling > 1:
            set_minor_locator(xy,get_locator_from_units(major_locator,**kwargs),**kwargs)
            return
        # From units
        if major_unit == 365: # Year
            for nmaj,itv in [3,1],[5,2],[7,3],[10,4],[nmax_ticks,6]:
                if nticks <= nmaj:
                    set_minor_locator(xy,MonthLocator(interval=itv),**kwargs)
                    break
        elif major_unit == 30: # Month
            if nticks > 10:
                set_minor_locator(xy,DayLocator(bymonthday=15),**kwargs)
            else:
                set_minor_locator(xy,WeekdayLocator(MO),**kwargs)
        elif major_unit == 7: # Week
            if nticks < 15: set_minor_locator(xy,DayLocator(),**kwargs)
        elif major_unit == 1: # Day
            for nmaj,itv in [5,2],[10,4],[nmax_ticks,6]:
                if nticks <= nmaj:
                    set_minor_locator(xy,HourLocator(byhour=range(0, 24, itv)),**kwargs)
                    break
        elif major_unit == 1./24.: # Hour
            for nmaj,itv in [5,10],[nmax_ticks,30]:
                if nticks <= nmaj:
                    set_minor_locator(xy,MinuteLocator(interval=itv),**kwargs)
                    break
        elif major_unit == 1./(24.*60.): # Minutes
            for nmaj,itv in [5,10],[nmax_ticks,30]:
                if nticks <= nmaj:
                    set_minor_locator(xy,SecondLocator(interval=itv),**kwargs)
                    break

def _xydate_(xy, tz=None, auto=True, fmt='dual', rotation=None,
    locator=None, minor_locator=True, minor_formatter=None, nominor=False,
    nmax_ticks=None, intv=None, trange=None, **kwargs):
    """
    - *tz*: Time zone.
    - *auto*: Auto Scaling [default: True]
    - *rotation*: Rotation angle of tick labels. If None, automatic [default: None]
    - *fmt*: Date format.
    - *locator*: Major locator. Can be within ['year','month','Weekday','day','hour','minute','second'] or be like matplotlib.dates.MonthLocator().
    - *minor_locator*: Minor locator.
    - *nominor*: Do not try to add minor ticks [default: False]
    - *nmax_ticks*: Maximal number of ticks
    - *locator_<keyword>*: <keyword> is passed to locator if locator is a string. If locator = 'month', locator = MonthLocator(locator_<keyword>=<value>).
    - *minor_locator_<keyword>*: Same with minor_locator.
    """
    axes = kwargs.get('ax',P.gca())
    axis = getattr(axes, xy+'axis')
##  print 'xlimssss',P.gca(),ax
##  if ax.get_xlim() == (0.,1.): raise 'a'
    # Base
    exec "axes.%saxis_date(tz)" % xy

    # Scale axes
    if auto:
        axes.autoscale_view(**{'scale'+xy:True,'scale'+{'x':'y','y':'x'}[xy]:False})

    # Guess locators
    major_locator = locator
    if nominor: minor_locator = False
    locs = ['year','month','weekday','day','hour','minute','second']
    locs.extend([(loc+'s') for loc in locs])
    kwmjl = kwfilter(kwargs, 'locator')
    kwmnl = kwfilter(kwargs, 'minor_locator')
    # - major
    if isinstance(major_locator, str):
        if major_locator.lower() in locs:
            major_locator = eval(major_locator.lower().title()+'Locator(**kwmjl)')
        else:
            major_locator = None
    else:
        major_locator = AutoDateLocator2()
    if major_locator:
        axis.set_major_locator(major_locator)
#    # - minor
#    if minor_locator is 1 or minor_locator is True:
#        minor_locator = 'auto'
#    if isinstance(minor_locator, str):
#        if minor_locator.lower() in locs:
#            minor_locator = eval(minor_locator.lower().title()+'Locator(**kwmnl)')
#        elif minor_locator.lower().startswith('auto'):
#            minor_locator = AutoDateMinorLocator(**kwmnl)
#        else:
#            minor_locator = None
#    if minor_locator:
#        print 'axis.set_minor_locator'
#        axis.set_minor_locator(minor_locator)



    # Tick format
    # - major
    if major_locator:
        if fmt is None: fmt = 'dual'
        if fmt == 'simple' or fmt == 1:
            fmt = AutoDateFormatter2(major_locator)
        elif fmt == 'dual' or fmt is True or fmt  == 2:
            fmt = AutoDualDateFormatter(major_locator)
        elif isinstance(fmt, str):
            fmt = DateFormatter(fmt)
        elif isinstance(fmt, list):
            fmt = DualDateFormatter(*fmt)
        elif isinstance(fmt, dict):
            fmt = DualDateFormatter(**fmt)
        elif isinstance(fmt, tuple):
            if len(fmt) == 1:
                fmt += ({}, )
            elif not isinstance(fmt[1], dict):
                fmt[1] = {}
            fmt = DualDateFormatter(fmt[0], **fmt[1])
        if fmt:
            axis.set_major_formatter(fmt)
#        if rotation is not None and not isinstance(fmt, (DualDateFormatter, AutoDualDateFormatter)):
#            P.setp(axis.get_majorticklabels(), "rotation", angle, *args)
#
#        # - minor
#        if not nominor and minor_locator:
#            if minor_formatter is True:
#                minor_formatter = AutoDateFormatter2()
#            elif isinstance(minor_formatter, str):
#                minor_formatter = DateFormatter(fmt)
#            if minor_formatter:
#                axis.set_minor_formatter(minor_formatter)


def get_locator_from_units(locator,*args,**kwargs):
    kw = {}
    reps = ['^interval$','^by']
    for key,val in kwargs.items():
        for rep in reps:
            if re.search(rep,key) is not None:
                kw[key] = val
    unit = locator._get_unit()
    units = [365.,30.,7.,1.,1/24.,1/(24*60.),1/(24*3600.)]
    locs = [YearLocator,MonthLocator,WeekdayLocator,
        DayLocator,HourLocator,MinuteLocator,SecondLocator]
    loc = locs[units.index(unit)]
    if loc is HourLocator: # fix hour problem
        kw.setdefault('interval', 1)
        kw['byhour'] = range(0, 24, kw['interval'])
        args = ()
    if loc is WeekdayLocator:
        args = (MO,)+args
    return loc(*args,**kw)

def set_major_locator(xy,locator,**kwargs):
    ax = kwargs.get('ax',P.gca())
    exec("ax.%saxis.set_major_locator(locator)" % xy)
    return get_major_locator(xy,**kwargs)

def get_major_locator(xy,real=False,**kwargs):
    ax = kwargs.get('ax',P.gca())
    locator = eval("ax.%saxis.get_major_locator()" % xy)
    if real and isinstance(locator,AutoDateLocator):
            locator = locator.get_locator(*locator.datalim_to_dt())
    return locator

def set_minor_locator(xy,locator,**kwargs):
    ax = kwargs.get('ax',P.gca())
    exec("ax.%saxis.set_minor_locator(locator)" % xy)
    return get_minor_locator(xy,**kwargs)

def get_minor_locator(xy,**kwargs):
    ax = kwargs.get('ax',P.gca())
    return eval("ax.%saxis.get_minor_locator()" % xy)

def add_logo(logofile, axes=None, fig=None, loc='lower left', scale=None, alpha=1,
    tool=None):
    """Add a logo to the background of the figure

    .. warning:: Except when ``loc='lower left', scale=1`` is passed as arguments,
      the logo is rescaled each time the figure is resized.

    .. note:: It is highly suggested to prefer png to other formats for your logo.

    :Params:

        - *file*: File of the image.
        - *axes*: Axes specs of the image (overwrites loc and scale keywords),
          passed to meth:`~matplotlib.figure.Figure.add_axes`.
        - *loc*: Position of the image (use words in 'lower', 'upper', 'left',
          'right', 'center') .
        - *scale*: Scale the image. By default, it is auto scale so that the logo is
          never greater than 1/4 of the figure width or height.
        - *alpha*: Alpha transparency.
        - *tool*: None, 'mpl' (png only) or 'pil'.
    """
    if not os.path.exists(logofile):
        raise VACUMMError("Logo file not found: %s"%logofile)

    # Read image
    if tool is None:
        if logofile.lower().endswith('png'):
            tool = 'mpl'
        else:
            tool = 'pil'
    else:
        tool = tool.lower()
    if tool not in ['mpl', 'pil']:
        raise VACUMMError('Please choose a valid tool for plotting your logo: '
            '"mpl" (png only), "pil", or None')

    if tool=='mpl':
        try:
            logo  = mpimg.imread(logofile)
            ny, nx = logo.shape[:2]
        except:
            raise VACUMMError("Can't read your logo file. "
                "Please convert it to png format or install PIL and use `tool='pil'`")
    else:
        try:
            import PIL.Image as Image
        except:
            try:
                import Image
            except:
                raise VACUMMError("Can't import PIL to plot your logo. Convert it to png"
                    " and use `tool='mpl'` calling add_logo")
        logo = Image.open(logofile)
        nx, ny = logo.size

    # Dimensions and position
    from matplotlib import rcParams
    dpi = rcParams['figure.dpi']
    if fig is None:
        fig = P.gcf()
    if axes:
        ax = fig.add_axes(axes, aspect=1, frameon=False, autoscale_on='off')
    elif scale==1 and loc=='lower left':
        ax = None
    else:


        dx = nx/(P.gcf().get_figwidth()*dpi)
        dy = ny/(P.gcf().get_figheight()*dpi)
        if scale is not None:
            dx *= scale
            dy *= scale
        else:
            scalex = scaley = 1.
            if dx > .25:
                scalex = .25/dx
            if dy > .25:
                scaley = .25/dy
            if scalex < 1. or scaley < 1.:
                scale = min(scalex,scaley)
                dx *= scale
                dy *= scale
        if 'upper' in loc:
            y = 1-dy
            anchor = 'N'
        elif 'lower' in loc:
            y = 0
            anchor = 'S'
        else:
            y = .5-dy/2
            anchor = 'C'
        if 'right' in loc:
            x = 1-dx
            anchor  = anchor + 'E'
        elif 'left' in loc:
            x = 0
            anchor  = anchor + 'W'
        else:
            x = .5-dx/2
            anchor  = anchor + 'C'
        xc = x+dx/2.
        yc = y+dy/2.
        if anchor=='CC':
            anchor = 'C'
        elif 'C' in anchor:
            anchor = anchor.replace('C', '')
        ax = P.gcf().add_axes([x,y,dx,dy], aspect=1, frameon=False, anchor=anchor)

    # Plot
    if ax is None:
        im = fig.figimage(logo, origin='upper', alpha=alpha)
    else:
        im = ax.imshow(logo, origin='upper', alpha=alpha)
        ax.axis('off')
        ax.set_zorder(-1)
    return im


def wedge(direction,width=18,fig=None,axes=None,figsize=(4,4),shadow=False,
    shadow_alpha=.5,circle=False,strength=None,strength_cmap=None,from_north=False,
    center=True,center_radius=.04,scale=1., **kwargs):
    """Plot a wedge

    - **direction**: Direction of the wedge in degrees from north
    - *width*: Base width of the wedge in degrees [default: 20]
    - *shadow*: Add a shadow [default: False]
    - *shadow_alpha*: Alpha transparency of the shadow [default: .5]
    - *circle*: Add a circle arround the wedge [default: True]
    - *shadow_<keyword>*: <keyword> is passed to :class:`matplotlib.patches.Shadow` class.
    - *circle_<keyword>*: <keyword> is passed to :class:`matplotlib.patches.Circle` class.
    """

    # Figure and axes
    onfig = False
    if fig is None:
        fig = P.figure(figsize=figsize,frameon=False)
    elif fig is True:
        fig = P.gcf()
    else:
        onfig = True
    if axes is None:
        axes = (0,0,1,1)
    axes = fig.add_axes(axes,frameon=False,aspect=1)

    # Relative strength
    # - color
    if strength is not None:
        facecolor = ScalarMappable(cmap=strength_cmap).cmap(strength)
    else:
        facecolor = 'k'
    # - scale
    scale = N.clip(scale, .1, 1.)

    # Plot the wedge
    r = .41*scale
##  direction = (90.-direction+180+360)%360.
    if from_north:
        direction = 90.-direction
    direction = (direction+180+360)%360.
    zorder = kwargs.pop('zorder',100)
    ww = Wedge((-r*1.1*N.cos(direction*N.pi/180.),-r*1.1*N.sin(direction*N.pi/180.)),2*r,
        direction-width/2,direction+width/2,
        **kwfilter(kwargs,['edgecolor','facecolor','linewidth','antialiased'],
            defaults=dict(facecolor=facecolor,linewidth=0)))
    ww.set_zorder(zorder)
    axes.add_patch(ww)

    # Center
    if center:
        cen = Circle((0,0),center_radius*scale,**kwfilter(kwargs,'center',
            defaults=dict(facecolor='#aaaaaa',linewidth=0)))
        cen.set_zorder(zorder+1)
        axes.add_patch(cen)


    # Shadow
    if shadow:
        dx = kwargs.pop('dx',r/10)
        dy = -kwargs.pop('dy',r/10)
        sdw = Shadow(ww,dx,dy,**kwfilter(kwargs,'shadow'))
        sdw.set_alpha(shadow_alpha)
        sdw.set_zorder(zorder-1)
        axes.add_patch(sdw)

    # Circle
    if circle:
        cir = Circle((0,0),.48,**kwfilter(kwargs,'circle'))
        cir.set_zorder(zorder-2)
        axes.add_patch(cir)

    # End
    P.xlim(-.52,.52)
    P.ylim(-.52,.52)
    kwargs['grid'] = False
    axes.set_xticks([])
    axes.set_yticks([])
    _end_plot_(**kwargs)
    return ww

def add_time_mark(date, ymin=0, ymax=1, color='r', line=True,
    marker=True, shadow=False, label=None, zorder=200,
    ax=None, fig=None, **kwargs):#:,alongx=True):
    """Plot a time mark along a time axis using a line and a marker

    :Params:

        - **date**: Date (converted using :func:`~vacumm.misc.atime.mpl`).
        - **color**, optional: Color of the line and marker.
        - **line**, optional: Plot the line?
        - **line_<kw>**, optional: Keyword ``kw`` is passed to
          :func:`~matplotlib.pyplot.axvline`.
        - **marker**, optional: Plot the marker? Or symbol.
        - **line_<kw>**, optional: Keyword ``kw`` is passed to
          :func:`~matplotlib.pyplot.axvline`.

    """
    # Inits
    kwargs.pop('xlim', None)
    kwargs.pop('ylim', None)
    kwline = kwfilter(kwargs, 'line')
    kwmarker = kwfilter(kwargs, 'marker')
    kwshad = kwfilter(kwargs, 'shadow')
    date = mpl(date) # MPL time
    if ax is None:
        if fig is None:
            ax = P.gca()
        else:
            ax = fig.gca()
    res = []
    axis = ax.axis()
    if shadow:
        from core_plot import add_shadow

    # Line
    if line:
        if isinstance(line, basestring):
            kwline.setdefault('linestyle', line)
        o = ax.axvline(date, ymin=ymin, ymax=ymax, color=color, zorder=zorder, **kwline)
        res.append(o)
        if shadow:
            add_shadow(o, ax=ax, **kwshad)

    # Marker
    if marker:
        if isinstance(marker, basestring):
            kwmarker.setdefault('marker', marker)
        kwmarker.setdefault('s', 40)
        o = ax.scatter([date], [axis[2]], c=color,  zorder=zorder, **kwmarker)
        res.append(o)
        if shadow:
            add_shadow(o, ax=ax, **kwshad)
        ax.axis(axis)
    return res

def add_key(key=None, pos=1, fmt='%s)', xmargin=10, ymargin=10, fig=None, axes=None, **kwargs):
    """Add a key to specify the plot number in a figure

    :Params:

        - **key**, optional: A string or a integer. If an integer is given,
          it converted to letter (1->'a').
        - **pos**, optional: Position of the key. It can be an integer or string:

            - ``1 = top left``
            - ``2 = top right``
            - ``3 = bottom right``
            - ``4 = bottom left``

          If negative, push it outside axis box.
          It can also be an explicit position
        - **fmt**, optional: Format of the string.
        - **{x|y}margin**, optional: Margin from axes in points.
        - Other keywords are passed to :func:`~matplotlib.pyplot.text`.

    """
    # Base
    if fig is None: fig = P.gcf()
    if axes is None and kwargs.has_key('ax'): axes = kwargs.pop('ax')
    if axes is None: axes = fig.gca()

    # Text
    if key is False: return
    if key is None or key is True:
        key = max(1, fig.axes.index(axes)+1)
    if isinstance(key, int):
        key = max(key, 1)
        key = 'abcdefghijklmnopqrstuvwxyz'[key-1]
    key = fmt%key

    # Position
    pos = kwargs.pop('loc', pos)
    if isinstance(pos, (list, tuple)):
        xpos, ypos = pos
        gpos = [[1, 2], [3, 4]][xpos < .5]
        pos = gpos[ypos > .5]
    else:
        if isinstance(pos, str):
            pos = pos.lower()
            gpos = [[1, 2], [3, 4]]['bottom' in pos]
            pos = gpos['right' in pos]
        xpos = int(abs(pos) in [2, 3])
        ypos = int(abs(pos) in [1, 2])

    # Alignment
    kwargs.setdefault('va',  ['bottom', 'top'][pos in [1, 2, -3, -4]])
    kwargs.setdefault('ha',  ['right', 'left'][pos in [1, 4, -2, -3]])

    # Margin
    from matplotlib.transforms import ScaledTranslation
    xmargin = 1.*N.abs(xmargin)*(1 - 2*(pos in [-1, 2, -3, 4]))
    ymargin = 1.*N.abs(ymargin)*(1 - 2*(pos in [1, 2, -3, -4]))
    transMargin = ScaledTranslation(xmargin/fig.dpi, ymargin/fig.dpi,
        fig.dpi_scale_trans)

    # Add to plot
    kwargs.setdefault('size', 18)
    return axes.text(xpos, ypos, key, transform=axes.transAxes+transMargin, **kwargs)


def hldays(color='.95', y=False, tmin=None, tmax=None, fig=None, axes=None, **kwargs):
    """Highlight different days with a different background color

    :Params:

        - *color*: Background color.
        - *y*: Work on Y axis.
        - Other keyparam are passed to :func:`~matplotlib.pyplot.axhspan`
          or :func:`~matplotlib.pyplot.axvspan`.
    """
    # Base
    if fig is None: fig = P.gcf()
    if axes is None and kwargs.has_key('ax'): axes = kwargs.pop('ax')
    if axes is None: axes = fig.gca()

    # Bonds
    ttmin, ttmax = getattr(axes, 'get_%slim'%'xy'[y])()
    from vacumm.misc.atime import numtime
    if tmin is None:
        tmin = ttmin
    else:
        tmin = numtime(tmin)
    if tmax is None:
        tmax = ttmax
    else:
        tmax = numtime(tmax)
    if int(tmin)==int(tmax) : return []

    # Patch
    kwargs.setdefault('zorder', 0)
    kwargs.update(facecolor=color)
    kwargs.setdefault('linewidth', 0)
    kwargs.setdefault('label', '_nolegend_')
    objs = []
    axspan = getattr(axes, 'ax%sspan'%'vh'[y])
    for t in xrange(int(tmin), int(tmax)+1, 2):
        objs.append(axspan(t, t+1, **kwargs))
#        objs[-1].set_zorder(0)
    getattr(axes, 'set_%slim'%'xy'[y])(tmin, tmax)
    return objs



def xi2a(inch,data=None):
    """Convert from inch to X axis coordinates"""
    if data is None:
        xmin,xmax = P.gca().xaxis.get_data_interval().get_bounds()
    else:
        xmin,xmax = minmax(data) # data
    trans = (xmax-xmin)/(P.gca().get_position().width*P.gcf().get_figwidth()) # inch to data
    return  trans*inch # data units

def yi2a(inch,data=None):
    """Convert from inch to Y axis coordinates"""
    if data is None:
        ymin,ymax = P.gca().yaxis.get_data_interval().get_bounds()
    else:
        ymin,ymax = minmax(data) # data
    trans = (ymax-ymin)/(P.gca().get_position().height*P.gcf().get_figheight()) # inch to data
    return  trans*inch # data units

def yi2f(inch):
    """Convert from inch to relative height of figure"""
    return inch/P.gcf().get_figheight()

def xi2f(inch):
    """Convert from inch to relative width of figure"""
    return inch/P.gcf().get_figwidth()

def gobjs(*args, **kwargs):
    """Store and retreive graphic objects of the current axes"""
    ax = kwargs.pop('ax', P.gca())
    if not hasattr(ax, 'gobjs'): ax.gobjs = {}
    for arg in args:
        if isinstance(arg, str):
            ax.gobjs.setdefault(arg, [])
            return ax.gobjs[arg]
        if isinstance(arg, dict):
            kwargs.update(arg)
    for key, val in kwargs.items():
        ax.gobjs.setdefault(key, [])
        ax.gobjs[key].append(val)
    if not len(args) and not len(kwargs):
        return ax.gobjs


def savefigs(basename, png=True, pdf=False, verbose=True, dpi=100, nodots=True, fig=None, **kwargs):
    """Save a figure in multiple formats (optimized for sphinx doc generator)

    - **basename**: File name without suffix
    - *png*: If True, save to png format
    - *pdf*: If True, save to pdf format
    - *verbose*: If True, print file names that are created
    - Other keywords are passed to :func:`~matplotlib.pyplot.savefig`

    .. warning::

        Dots ('.') in figure name are converted to dashes  ('-')
        for compatilibity reasons with latex.
        It is optimized for the sphinx doc generator and may
        NOT BE SUITABLE FOR YOU. In this case, use
        :func:`~matplotlib.pyplot.savefig` instead.
    """
    for ext in 'png', 'py':
        if basename.endswith('.'+ext):
            basename = basename[:-len(ext)-1]
            break
    if nodots:
        basename = os.path.join(os.path.dirname(basename), os.path.basename(basename).replace('.', '-'))
    fig = fig if fig is not None else P.gcf()
    figfiles = []
    for ext in ('png', 'pdf'):
        if eval(ext):
            file = '%s.%s'%(basename, ext)
            if ext == 'pdf':
                fig.savefig(file)
            else:
                fig.savefig(file, dpi=dpi)
            figfiles.append(file)
            if verbose: print 'Wrote to '+file
    return figfiles

def make_movie(fig_pattern, outfile, delay=1, clean=False, verbose=False, windows=True):
    """Make a movie from a series of figures

    - **fig_pattern** Unix-like file pattern to select png figures.
    - **outfile**: Output file.
    - *delay**: Delay in seconds between two frames.
    - *windows*: When output is a movie, choose basic mpeg for windows instead of mpeg4 codec.
    - *clean*: If ``True``, remove figures once the movies os created.

    :Usage:

    >>> make_movie('/home/toto/fig*.png', 'movie.mpg', delay=.5, clean=True)
    """
    # Check input files
#   fig_pattern = os.path.abspath(fig_pattern)
#   files = glob.glob(fig_pattern)
#   assert len(files), 'No figure files found'
#   if verbose: print 'Found %i figure files'%len(files)
#   outfile = os.path.abspath(outfile)

    # Command
    if outfile.endswith('gif'):
        delay *= 100.
        if isinstance(fig_pattern, list):
            fig_pattern = ' '.join(fig_pattern)
        cmd = 'convert -delay %(delay)i %(fig_pattern)s %(outfile)s'%vars()
    else:
        delay = 1./delay
        if windows:
            codec = 'msmpeg4v2'
        else:
            codec = 'mpeg4'
        cmd = "mencoder mf://%(fig_pattern)s -lavcopts vcodec=%(codec)s:vbitrate=800 -mf type=png:fps=%(delay)s  -ovc lavc  -oac copy  -o %(outfile)s"%vars()

    # Execution
    import subprocess
    if verbose: print 'Executing command: '+cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err =  p.communicate()
    if verbose:
        print out
        print err
    print 'ret', p.returncode,  err
    if p.returncode:
        raise SystemError


    # Clean files
    if clean:
        if verbose: print 'Cleaning files'
        p = subprocess.Popen('rm -f %(fig_pattern)s'%vars(),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,  err = p.communicate()
        if verbose:
            print out
            print err
        if p.returncode:
            raise SystemError
    return outfile


def get_cls(n, colors=simple_colors, linestyles=linestyles):
    """Get a list of string argument for :func:`~matplotlib.pyplot.plot` to specify the color and the linestyle

    - **n**: Length of the list
    - *colors*: Colors on which to cycle
    - *linestyles*: Linestyles on which to cycle
    """
    cls = []
    cc = []
    for c in cycle(colors):
        cc.append(c)
        if len(cc)>=n: break
    ll = []
    for ls in cycle(linestyles):
        ll.append(ls)
        if len(ll)>=n: break
    cls = []
    for c, ls in zip(cc, ll):
        cls.append(c+ls)
    return cls

#############
# Utilities #
#############

def _check_var_(var,rank,order=None,xaxis=None,yaxis=None,tadd=None,tadd_copy=True,tlocal=False,xatts=None,yatts=None,**kwargs):
    """
    - *xaxis*: Use this X axis.
    - *yaxis*: Use this Y axis.
    - *tadd*: Add value to time (like tadd=1 or tadd=(1,'day'))
    - *tadd*: Make a copy of the axis before any tadd [default: True]
    - *tlocal*: Convert current UTC time axes to local time [default: False]
    """
    if var is None: return

    # Arrows
    istuple = isinstance(var,tuple)
    if not istuple:
        vars = [cdms.createVariable(var),]
    elif len(var) in [1,3]:
        vars = [cdms.createVariable(v) for v in var]
    else:
        # Get vectors and modulus
        vx,vy = [cdms.createVariable(v) for v in var]
        vars = [MV.sqrt(vx**2+vy**2),vx,vy]
        ln = []
        for vv in vx,vy:
            if hasattr(vv,'units'):
                vars[0].units = vv.units
                break
            if hasattr(vv,'long_name'):
                ln.append(vv.long_name)
        ln  = ' and '.join(ln)
        if len(ln) == 0: ln = 'Modulus'
        vars[0].long_name = ln
        vars[0].setAxisList(vars[2].getAxisList())
        vars[0].setGrid(vars[2].getGrid())
        vars[0].id = 'modulus'

# Loop on all set
    for ivar,var in enumerate(vars):
        assert var.rank() >= rank, 'Your variable must have a rank > %i (current rank is %i)' % (rank,var.rank())
        check_axes(var)
        clone = False
        # Reorder axes
        if order is not None:
            var = var.clone() ; clone = True
            order = order.replace('.','')
            order = ('-'*(rank-len(order))+order)[-rank-1:]
            var_order = var.getOrder()
            # Order requested axes if found
            for i,c in enumerate(order[::-1]):
                pos = var_order.rfind(c)
                if c != '-' and pos > -1 and var_order[-i-1] != c:
                    var = var.reorder('-'*(var.rank()-i-1)+c+'-'*i)
            # Move bad axes to front
            i = 0
            for c in 'xyzt':
                pos = var.getOrder()[-rank:].rfind(c)
                if pos > -1 and order[pos-rank] != '-' and order.find(c) == -1:
                    if i >= (var.rank()-rank-order.count('-')):
                        print 'It seems that axis #%i (%s) of var should not be there and cannot be moved' % (i,c.upper())
                    else:
                        var = var.reorder('-'*(var.rank()-i-1)+c+'-'*i)
                    i += 1
        # Reduce rank
        while var.rank() > rank:
            if not clone:
                var = var.clone() ; clone = True
            var = MV.average(var,0)
        # Put variable on new axes if requested
        if rank == 1 and xaxis is not None:
            assert len(xaxis) == len(var), 'Length of new xaxis (%i) is different from length of var (%i)' % (len(xaxis), len(var))
            if type(xaxis) is type([]): xaxis = N.array(xaxis,copy=0)
            if not isinstance(xaxis,AbstractAxis):
                tmpaxis = var.getAxis(0).clone()
                tmpaxis[:] = xaxis
            if not clone:
                var = var.clone() ; clone = True
            var.setAxis(0,xaxis)
        elif xaxis is not None or yaxis is not None:
            if not clone:
                var = var.clone() ; clone = True
            var = var2d(var,xaxis,yaxis,xatts=xatts,yatts=yatts)
        vars[ivar] = var

    # Set to local time
    if tlocal and tadd is None:
        tadd = (-time.altzone,'seconds')

    # Add time
    if tadd not in [None,0]:
        taxes = []
        for var in vars:
            tpos = var.getOrder().find('t')
            if tpos != -1:
                taxis = var.getTime()
                if taxis in taxes: continue # Already treated
                if isinstance(tadd,(list,tuple)):
                    mytaxis = axis_add(taxis,tadd[0],tadd[1],copy=tadd_copy)
                else:
                    mytaxis = axis_add(taxis,tadd,copy=tadd_copy)
                taxes.append(taxis)
                if tadd_copy: # A copy is safer
                    var.setAxis(tpos,mytaxis)

    if len(vars[0].shape) > 1:
        return tuple(vars)
    return vars[0]



def _start_plot_(figure=None,figsize=None,subplot=None,subplots_adjust=None,bgcolor=None,noframe=False, fullscreen=False, **kwargs):
    """
    Misc settings:

    - *figure*: Figure number.
    - *figsize*: Initialize the figure with this size.
    - *subplots_adjust*: Dictionary sent to :func:`~matplotlib.pyplot.subplots_adjust`. You can also use keyparams 'left', 'right', 'top', 'bottom', 'wspace', 'hspace' !
    - *sa*: Alias for subplots_adjust.
    - *bgcolor*: Background axis color.
    - *axes_rect*: [left, bottom, width, height] in normalized (0,1) units to create axes using :func:`~matplotlib.pyplot.axes`.
    - *axes_<keyword>*: <keyword> is passed to :func:`~matplotlib.pyplot.axes`.
    """
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
    axes_rect = kwargs.pop('axes_rect', default_axes_rect)
    if axes_rect is not None:
        axes_rect = [axes_rect]
    else:
        axes_rect = []
    if noframe:
        kwargs['figure_frameon'] = False
        kwargs['axes_frameon'] = False

    # Figure size
    if figsize is not None:
        if figure is None:
            P.gcf().set_size_inches(figsize,forward=True)
        else:
            P.figure(figure,figsize=figsize,**kwfilter(kwargs,'figure'))
    elif figure is not None:
        P.figure(figure,**kwfilter(kwargs,'figure'))
    if noframe:
        P.gcf().set_frameon(False)

    # Axes
    kwaxes = kwfilter(kwargs,'axes')
    if subplot is not None:
        if isinstance(subplot,(list,tuple)):
            P.subplot(*subplot,**kwaxes)
        else:
            P.subplot(subplot,**kwaxes)
    elif axes_rect:
        P.axes(*axes_rect,**kwaxes)
    if subplots_adjust is not None:
        P.subplots_adjust(**subplots_adjust)
    if noframe:
        P.gca().set_frame_on(False)

    # Axis background color
    if bgcolor is not None:
        P.gca().set_axis_bgcolor(bgcolor)





def colorbar(pp=None, vars=None, drawedges=False, levels=None, colorbar_horizontal=False, colorbar=True, colorbar_visible=True, colorbar_position=None,units=None, standalone=False,cax=None,extend='neither', **kwargs):
    """
    Colorbar:

    - *colorbar*: Plot the colorbar [default: True]
    - *colorbar_horizontal*: Colorbar is horizontal [defaults: False]
    - *colorbar_position*: To change the default position. Position is relative to the SPACE LEFT in the form (center,width) with values in [0,1] [defaults: False]
    - *colorbar_visible*: Colorbar is visible [defaults: True]
    - *colorbar_<keyword>*: <keyword> is passed to :func:`~matplotlib.pyplot.colorbar`
    - *units*: Indicate these units on along the colorbar, else it guessed from the variable or suppressed if value is False [default: None]
    """
    if pp is not False and colorbar:

        # Orientation
        if colorbar_horizontal:
            orientation = 'horizontal'
        else:
            orientation = 'vertical'
        xy = 'yx'[int(colorbar_horizontal)]
        # Force the position
        gca = None
        cax = kwargs.pop('colorbar_axes', cax)
        if cax is not None:
            if isinstance(cax, list):
                if not standalone: gca = P.gca()
                cax = P.axes(cax)
        elif colorbar_position is not None:
            if colorbar_position in [True,1,'auto']:
                if colorbar_horizontal:
                    colorbar_position = 0.5
                else:
                    colorbar_position = 0.3
            if not isinstance(colorbar_position,(list,tuple)):
                colorbar_position = (colorbar_position,0.25)
            bbox = P.gca().get_position(original=True)
            cb_start = colorbar_position[0]-colorbar_position[1]
            if colorbar_horizontal:
                cax = [bbox[0],bbox[1]*cb_start,
                    bbox[2],bbox[1]*colorbar_position[1]]
            else:
                bbox_start = bbox[0]+bbox[2]
                cax = [bbox_start+(1-bbox_start)*cb_start,bbox[1],
                    (1-bbox_start)*colorbar_position[1],bbox[3]]
            if not standalone: gca = P.gca()
            cax = P.axes(cax)

        # Plot itself
        kwaxis = kwfilter(kwargs, 'colorbar_'+xy)
        kwlabel = kwfilter(kwargs, 'colorbar_label')
        kwargs.setdefault('colorbar_extend', extend)
        cb = P.colorbar(pp,cax,orientation = orientation,
            **kwfilter(kwargs,'colorbar', defaults=dict(drawedges=drawedges, ticks=levels)))
        if gca is not None: gobjs(ax=gca, colorbar=cb)
        for t in cb.ax.xaxis.get_major_ticks():
            t.tick1On = t.tick2On = False

        # Decorations
        if not isinstance(vars,tuple): vars = (vars,)
        var = vars[0]
        if units is not False:
#           func = getattr(cb.ax,'set_%slabel'%xy)
            if units not in [True,None]:
                cb.set_label(units, **kwlabel)
            elif units is not False and cdms.isVariable(var) and var.attributes.has_key('units'):
                cb.set_label(var.units, **kwlabel)
        if isaxis(var):
            if istime(var):
                var = mpl(var)
            eval('cb.ax.set_%slim'%xy)(minmax(var))
            eval('cb.ax.dataLim.interval%s().set_bounds'%xy)(*minmax(var))
            decorate_axis(var,vertical=not colorbar_horizontal,ax=cb.ax,**kwaxis)
        # Hide it?
        if not colorbar_visible:
            cb.ax.set_visible(False)
        if gca is not None: P.gcf().sca(gca)

        return cb

_colorbar_ = colorbar

def colorbar_new(sm=None, vars=None, drawedges=False, levels=None, horizontal=False, visible=True, position=None,units=None, standalone=False, cax=None, ax=None, extend='neither', **kwargs):
    """
    Plot a colorbar as :func:`~matplotlib.pyplot.colorbar` with some advanced features

    - *horizontal*: Colorbar is horizontal [defaults: False]
    - *position*: To change the default position. Position is relative to the SPACE LEFT in the form (center,width) with values in [0,1] [defaults: False]
    - *visible*: Colorbar is visible [defaults: True]
    - *units*: Indicate these units on along the colorbar, else it guessed from the variable or suppressed if value is False [default: None]
    - Other key are passed to :func:`~matplotlib.pyplot.colorbar`
    """
    if pp is False: return

    # Keywords
    kwcb = kwfilter(kwargs, 'colorbar', defaults=dict(drawedges=drawedges))
    kwargs.update(kwcb)

    # Standalone


    # Orientation
    if horizontal:
        orientation = 'horizontal'
    else:
        orientation = 'vertical'
    xy = 'yx'[int(horizontal)]

    # Force the position
    gca = None
    cax = kwargs.pop('axes', cax)
    if cax is not None:
        if isinstance(cax, list):
            if not standalone: gca = P.gca()
            cax = P.axes(cax)
    elif position is not None:
        if position in [True,1,'auto']:
            if horizontal:
                position = 0.5
            else:
                position = 0.3
        if not isinstance(position,(list,tuple)):
            position = (position,0.25)
        bbox = P.gca().get_position(original=True)
        cb_start = position[0]-position[1]
        if horizontal:
            cax = [bbox[0],bbox[1]*cb_start,
                bbox[2],bbox[1]*position[1]]
        else:
            bbox_start = bbox[0]+bbox[2]
            cax = [bbox_start+(1-bbox_start)*cb_start,bbox[1],
                (1-bbox_start)*position[1],bbox[3]]
        if not standalone: gca = P.gca()
        cax = P.axes(cax)

    # Plot itself
    kwaxis = kwfilter(kwargs, xy)
    kwlabel = kwfilter(kwargs, 'label')
    cb = P.colorbar(pp,cax,ticks=levels,orientation = orientation,extend=extend,**kwargs)
    if gca is not None: gobjs(ax=gca, colorbar=cb)
    for t in cb.ax.xaxis.get_major_ticks():
        t.tick1On = t.tick2On = False

    # Decorations
    if not isinstance(vars,tuple): vars = (vars,)
    var = vars[0]
    if units is not False:
#           func = getattr(cb.ax,'set_%slabel'%xy)
        if units not in [True,None]:
            cb.set_label(units, **kwlabel)
        elif units is not False and cdms.isVariable(var) and var.attributes.has_key('units'):
            cb.set_label(var.units, **kwlabel)
    if isaxis(var):
        if istime(var):
            var = mpl(var)
        eval('cb.ax.set_%slim'%xy)(minmax(var))
        eval('cb.ax.dataLim.interval%s().set_bounds'%xy)(*minmax(var))
        decorate_axis(var,vertical=not horizontal,ax=cb.ax,**kwaxis)

    # Hide it?
    if not visible:
        cb.ax.set_visible(False)
    if gca is not None: P.gcf().sca(gca)

    return cb



def _levels_(var,anomaly=None,levels=None,nmax=8,vmin=None,vmax=None,**kwargs):
    """
    Data levels:

    - *levels*: Force the use of these levels for contours.
    - *anomaly*: Levels must be symetric about zero (useful for anomly plots) and color map is misc.plot.cmap_bwre() [default: None]. If None, is is turned to True if max is near -min (20%)
    - *nmax*: Max number of levels (see misc.auto_scale) [default: 10]
    - *vmin*: Force min value for pcolor.
    - *vmax*: Force max value for pcolor.
    """
    vmin = kwargs.get('min_value', vmin)
    vmax = kwargs.get('max_value', vmax)
    if levels is None:
        if vmin is not None or vmax is not None:
            levels = auto_scale(var,vmin=vmin,vmax=vmax,nmax=nmax)
        else:
            if anomaly is None:
                mm = minmax(var)
                anomaly = N.absolute((mm[1]-mm[0])/(.5*(mm[0]+mm[1]))) < .2
            if anomaly:
                mm = max(N.absolute(minmax(var)))
                levels = auto_scale([-mm,mm],nmax=nmax)
            else:
                levels = auto_scale(var,nmax=nmax)
#    if vmin is None:
    vmin = levels[0]
#    if vmax is None:
    vmax = levels[-1]
    gobjs(levels=levels)
    return levels,vmin,vmax,anomaly

#def _latlab_(deg):
#   return latlab(deg,decimal=False,no_seconds=True,no_zeros=True, )
#
#def _lonlab_(deg):
#   return lonlab(deg,decimal=False,no_seconds=True,no_zeros=True)

def decorate_axis(axis=None,vertical=0,date_rotation=None,date_fmt=None,date_locator=None,date_minor_locator=None,date_nominor=False,nodate=False,values=None,ax=None,**kwargs):
    """
    Axis decoration:

    - *[x/y]title*: Main label of the axis.
    - *[x/y]label*: Same as x/ytitle.
    - *[x/y]hide]*: Hide x/ytick labels.
    - *[x/y]rotation*: Rotation of x/ytick labels (except for time).
    - *[x/y]strict*: Strict axis limit.
    - *[x/y]lim*: Override axis limit.
    - *[x/y]min/max*: Override axis limit by x/ylim.
    - *[x/y]minmax/maxmin*: Set maximal value of x/ymin and minimal value or x/ymax.
    - *[x/y]nmax*: Max number of ticks.
    - *x/yfmt*: Numeric format for x/y axis if not of time type.
    - *date_rotation*: Rotation of time ticklabels.
    - *date_locator*: Major locator for dates.
    - *date_minor_locator*: Minor locator for dates.
    - *date_nominor*: Suppress minor ticks for dates.
    - *date_nmax_ticks*: Max number of ticks for dates
    - *nodate*: Time axis must not be formatted as a time axis [default: False]
    """
    # Current axes
    if ax is not None:
        gca = P.gca()
        P.gcf().sca(ax)
    else:
        gca = None

    # Keywords
    defaults = {}
    props = {}
    vatts = ['min', 'max', 'minmax', 'maxmin', ('label', 'title'), 'units']
    aatts = ['strict', ('fmt', 'format', 'tickfmt', 'tickformat'), ('rotation', 'rotate'),
        'lim', 'hide', 'type', 'locator', 'minor_locator', ('nmax', 'nmax_ticks')]
    defaults = dict()#strict=True)
    for raw_att in aatts+vatts:

        # Merge aliases
        if isinstance(raw_att, tuple):
            dict_aliases(kwargs, raw_att)
            att = raw_att[0]
        else:
            att = raw_att
            raw_att = raw_att,

        # Default values
        # - base
        default = None
        # - get it from cdms variable
        if cdms.isVariable(axis) and raw_att in vatts:
            kwvar = dict_aliases(kwfilter(kwargs, 'v', copy=1), raw_att)
            if not kwargs.has_key(att) and kwvar.has_key(att):
                default = kwvar[att]
        # - special
        if defaults.has_key(att):
            default = defaults[att]

        # Get attribute
        props[att] = kwargs.get(att, default)
    props['label_kwargs'] = kwfilter(kwargs, 'label_', copy=1)
    props['label_kwargs'].update(kwfilter(kwargs, 'title_', copy=1))

    # X or Y axis?
    if props['type'] is None:
        props['type'] = axis_type(axis)
    xy = 'xy'[vertical]

    # Functions
    exec "tick_func = P.%(xy)sticks ; label_func = P.%(xy)slabel ; date_func = %(xy)sdate ; keyaxis = '%(xy)saxis' ; lim_func = P.%(xy)slim" % vars()
    funcs = {}

    # Limits
    if values is None:
        values = axis
    axmin, axmax = None, None
    if props['strict'] is None:
        props['strict'] = props['type'] == 't'
    if props['strict']:
        axmin, axmax = minmax(values)
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
    if props['min'] is not None: axmin = props['min']
    if props['max'] is not None: axmax = props['max']
    if props['minmax'] is not None:
        if axmin is None:
            axmin = min(values)
        axmin = min(axmin, props['minmax'])
    if props['maxmin'] is not None:
        if axmax is None:
            axmax = max(values)
        axmax = max(axmax, props['maxmin'])
    if axmin is not None:
        lim_func(**{xy+'min':axmin})
    if axmax is not None:
        lim_func(**{xy+'max':axmax})
    #axis_save = P.axis()

    # Ticks
    if props['type'] in ['x','y','z']: # Lon, lat or dep axis

        if props['type'] in ['x','y']:
            lab_val = geo_scale(axis, nmax=props['nmax'])
            if props['type'] == 'x':
                lab_func = lonlab
            else:
                lab_func = latlab
        else:
            lab_val = auto_scale(axis)
            lab_func = deplab
        kwtf = {}
        if props['fmt'] is not None: kwtf['fmt'] = props['fmt']
        tick_func(lab_val,lab_func(lab_val,**kwtf))

    elif props['type'] == 't' and not nodate: # Time axis

        if date_fmt is None and hasattr(axis,'units') :
            # Detect climatological plot -> adapt date format
            if len(axis) == 1:
                nmonth = 1
            elif not axis.units.startswith('months'):
                units = 'months '+' '.join(axis.units.split()[1:])
                ct = axis.subAxis(0,len(axis)).asComponentTime(cdtime.DefaultCalendar)
                nmonth = ct[1].torel(units).value - ct[0].torel(units).value + 1
            else:
                nmonth = axis[-1] - axis[0] + 1
            if nmonth == 12 and int(axis.units.split()[2].split('-')[0]) == 1:
                print "DETECTED PLOT OF CLIMATOLOGY"
                if date_rotation is None: date_rotation = 0.
                if date_fmt is None: date_fmt = '%b'
                if date_locator is None:
                    date_locator = MonthLocator()
                if date_minor_locator is None:
                    date_minor_locator = MonthLocator(bymonthday=15)
        # Rotate dates
        if date_rotation is None:
            if vertical:
                date_rotation = 0.
            else:
                date_rotation = 30.
        # Ok let's format
        trange = axis[-1]-axis[0]
        kwdate = kwfilter(kwargs,'date', copy=1)
        kwdate.setdefault('nmax_ticks', props['nmax'])
        kwdate.setdefault('auto', axmin is None or axmax is None)
        date_func(rotation=date_rotation,fmt=date_fmt,locator=date_locator,
            minor_locator=date_minor_locator,nominor=date_nominor,trange=trange,
            nospecial=date_rotation is not None, **kwdate)

    else:
        if props['fmt'] is not None: # Other axes
        # Numeric format
            eval('P.gca().get_%saxis()'%xy).set_major_formatter(FormatStrFormatter(props['fmt']))
    if props['locator']:
        eval('P.gca().get_%saxis()'%xy).set_major_locator(props['locator'])
    if props['minor_locator']:
        eval('P.gca().get_%saxis()'%xy).set_minor_locator(props['minor_locator'])
    if props['type'] != 't' or nodate:
        if props['rotation'] is not None:
            rotate_tick_labels(props['rotation'],vertical=vertical)
    if nodate: props['type'] = '-'

    # Hiding
    if props['hide'] in [True,False] and props['hide']:
        exec xy+'hide()'

    # Label
    elif props['label'] is not False:
        if isinstance(props['label'],(str, unicode)):
            label_func(props['label'], **props['label_kwargs'])
        elif props['type'] not in ['x','y','t','z']:
            ll = getattr(axis,'long_name','')
            if props['units'] is None:
                units = getattr(axis,'units',None)
            slabel = []
            if ll and not cdms.isVariable(axis):
                slabel.append(ll)
            if isinstance(units, (str, unicode)) and units not in ['','-','None']:
                if len(slabel):
                    slabel.append('[%s]'%units)
                else:
                    slabel.append(units)
            slabel = ' '.join(slabel)
            if slabel:
                label_func(slabel, **props['label_kwargs'])


    if axmin is not None:
        lim_func(**{xy+'min':axmin})
    if axmax is not None:
        lim_func(**{xy+'max':axmax})
    #P.axis(axis_save)
    if gca is not None: P.gcf().sca(gca)


#def _strict_lims_(xaxis=None,yaxis=None):
#   if xaxis is not None:
#       P.xlim(minmax(xaxis))
#   if yaxis is not None:
#       P.ylim(minmax(yaxis))


def _end_plot_(var=None,grid=True,figtext=None,show=True,close=False,savefig=None,title=None,logo=False,fullscreen=False,anchor=None, autoresize=2, key=None, dayhl=False, **kwargs):
    """
    - *title*: Title of the figure [defaults to var.long_name or '']
    - *grid*: Plot the grid [default: True]
    - *dayhl*: Add day highlithing [default: False]
    - *figtext*: figtext Add text at a specified position on the figure. Example: figtext=[0,0,'text'] add a 'text' at the lower left corner, or simply figtext='text'.
    - *anchor*: Anchor of the axes (useful when resizing) in ['C', 'SW', 'S', 'SE', 'E', 'NE', 'N', 'NW', 'W'].
    - *logo*: Add a logo to the figure [default: False]. logo can be a file name.
    - *show*: Display the figure [default: True]
    - *savefig*: Save the figure to this file.
    - *savefigs*: Save the figure into multiple formats using :func:`savefigs` and 'savefigs' as the prefix to the files.
    - *autoresize*: Auto resize the figure according axes (1 or True), axes+margins (2). If 0 or False, not resized [default: False=2].
    - *key*: Add a key (like 'a)') to the axes using add_key is different from None [default: None]
    - *close*: Close the figure at the end [default: False]
    - *title_<keyword>*: <keyword> is passed to :func:`~matplotlib.pyplot.title`
    - *key_<keyword>*: <keyword> is passed to :func:`add_key`
    - *logo_<keyword>*: <keyword> is passed to :func:`add_logo`
    - *figtext_<keyword>*: <keyword> is passed to :func:`~matplotlib.pyplot.figtext`
    - *savefig_<keyword>*: <keyword> is passed to :func:`~matplotlib.pyplot.savefig`
    - *savefigs_<keyword>*: <keyword> is passed to :func:`savefigs`
    - *grid_<keyword>*: <keyword> is passed to :func:`~matplotlib.pyplot.grid`
    """

    # Overlay case
    ax = P.gca()
    fig = P.gcf()
    if autoresize and ax.get_aspect() != 'auto' and \
        isinstance(ax, Subplot) and \
        ax.is_first_col() and ax.is_first_row() and \
        ax.is_last_col() and ax.is_last_row():
        r = ax.get_aspect()
        if r=='equal': r=1.
        r *= ax.get_data_ratio()
#       rect = P.gca().get_position(True)
#       if autoresize == 2: r *= rect.width/rect.height
        w,h = ax.get_figure().get_size_inches()
        a = 1.*w*h
        W = N.sqrt(a/r)
        H = r*W
        P.gcf().set_size_inches(W,H)
        if anchor is None: anchor='C'
    if anchor is not None:
        ax.set_anchor(anchor)

    # Grid
    gobjs(grid=P.grid(grid,**kwfilter(kwargs,'grid')))

    # Highlight days
    if dayhl:
        gobjs(hldays=hldays(**kwfilter(kwargs, 'dayhl')))

    # Title
    if title is not False:
        if title not in [True,None]:
            P.title(title,**kwfilter(kwargs,'title'))
        elif var is not None:
            if isinstance(var,tuple): var = var[0]
            if var.attributes.has_key('long_name'):
                gobjs(title=P.title(var.long_name,**kwfilter(kwargs,'title')))

    # Fig text
    if figtext is not None:
        if type(figtext) is type('a'):
            figtext = (figtext,)
        if len(figtext) == 1:
            figtext = [0,0,figtext[0]]
        gobjs(figtext=P.figtext(*figtext,**kwfilter(kwargs,'figtext',defaults={'color':'#888888'})))

    # Key of axes
    if key is not None:
        gobjs(key=add_key(key, **kwfilter(kwargs, 'key')))

    # Logo
    if logo is not False:
        kwlogo = kwfilter(kwargs,'logo')
        if isinstance(logo,str):
            kwlogo['file'] = logo
        gobjs(logo=add_logo(**kwlogo))

    # Save it
    backend = P.get_backend().lower()
    if savefig is not None:
        suff_standalone = ['ps','pdf','svg']
        suff_others = ['png','gif','jpg','jpeg','bmp']
        suff = os.path.splitext(savefig)[1][1:]
        if suff.lower() not in suff_others:
            if backend not in suff_standalone:
                savefig += '.png'
            elif backend in suff_standalone and not suff.lower().endwith(backend):
                savefig += '.'+backend
        P.savefig(savefig,**kwfilter(kwargs,'savefig'))
    if kwargs.has_key('savefigs'):
        if kwargs['savefigs'] is not None:
            savefigs(kwargs['savefigs'], **kwfilter(kwargs, 'savefigs_'))
        else:
            del kwargs['savefigs']

    # Show or close
    if show:
        viewers = {
            'pdf':['/usr/bin/kghostview','/usr/bin/evince','/usr/bin/xpdf','/usr/local/bin/acroread'],
            'ps':['/usr/bin/kghostview','/usr/bin/evince','/usr/bin/ghostview'],
            'svg':['/usr/bin/svgdisplay','/usr/bin/konqueror']}
        if backend in viewers.keys():
            if savefig is None:
                tmpfig = mktemp(suffix='.'+backend)
                P.savefig(tmpfig,**kwfilter(kwargs,'savefig'))
            else:
                tmpfig = savefig
            for viewer in viewers[backend]:
                if os.path.exists(viewer):
                    cmd = '%s %s' % (viewer,tmpfig)
                    try:
                        os.system(cmd)
                    except:
                        print 'Unable to view file with command:',cmd
                    if savefig is None: os.remove(tmpfig)
                    return
        else:
            P.show()
    elif close:#savefig is not None and kwargs.has_key('savefigs') and close:
        P.close()



class CachedBasemap(Basemap):
    def __init__(self, *args, **kwargs):
        from matplotlib import get_home
        cache_dir = kwargs.pop('cache_dir', os.path.join(get_home(), 'basemap', 'cached_maps'))
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)


########
# Docs #
########

def _fill_doc_(*funcs):
    for func in funcs:
        if func.__doc__ is None: continue
        sc = Template(func.__doc__)
        for mm in Template.pattern.findall(func.__doc__):
            func_name = mm[1]
            func_doc = eval(func_name).__doc__.strip(' \t\n')
            func.__doc__ = Template(func.__doc__).safe_substitute(**{func_name:func_doc})

_fill_doc_(xdate, ydate, taylor, dtaylor)



@docfill
def curve2(*args, **kwargs):
    """curve2(data, parg=None, axis=None, title=None, savefig=None, show=True, **kwargs)

    Plot 1D data as a curve and return a :class:`~vacumm.misc.core_plot.Curve` object
    with properties also from :class:`~vacumm.misc.core_plot.Plot` and
    :class:`~vacumm.misc.core_plot.Plot1D`.


    :Examples:

        >>> curve2(sst, 'r-', shadow=True, ymax=25., long_name='SST', subplot=212)

        >>> curve2(xe, label='Full signal', show=False)
        >>> c = curve2(xef, label='Filtered signal', show=False)
        >>> c.legend()
        >>> c.savefig('xe.png')


    :Data params:

        {Curve_load_data[data]}
        {Plot1D__check_order_[vertical]}
        {Plot1D__set_axes_[axis]}
        {Plot[long_name]}
        {Plot[units]}

    :Curve params:

        {Curve_plot[parg]}
        {Curve_plot[nosingle]}
        {Curve_plot[err]}
        {Curve_plot[err_<param>]}
        {Curve_plot[fill_between]}
        {Curve_plot[fill_between_<param>]}
        {Curve_plot[label]}
        {Curve_plot[shadow]}
        {Curve_plot[shadow_<param>]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    :Extra: A few matplotlib plot arguments can be passed (for exemple: ``linewidths``).

    """
    from core_plot import Curve
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    if len(args)>1 and not kwargs.has_key('parg'):
        kwargs['parg'] = args[1]
        args = args[0:1]+args[2:]
    return Curve(*args, **kwargs)

@docfill
def bar2(*args, **kwargs):
    """bar2(data, width=1.,lag=0, align='center', offset=None, title=None, savefig=None, show=True, **kwargs)

    Plot data as a bar plot and return a :class:`~vacumm.misc.core_plot.Bar` object
    with properties also from :class:`~vacumm.misc.core_plot.Plot` and
    :class:`~vacumm.misc.core_plot.Plot1D`.

    :Examples:

        >>> bar2(rain, color='c', align='left', width=0.95, savefig='rain.png')

        >>> bar2(rain, vminmax=5.,show=False)
        >>> bar2(snow, offset=rain, title='Precipitations')

    :Data params:

        {Curve_load_data[data]}
        {Plot1D__check_order_[vertical]}
        {Plot1D__set_axes_[axis]}
        {Plot[long_name]}
        {Plot[units]}

    :Bar params:

        {Bar_plot[width]}
        {Bar_plot[lag]}
        {Bar_plot[align]}
        {Bar_plot[offset]}
        {Bar_plot[label]}
        {Bar_plot[shadow]}
        {Bar_plot[shadow_<param>]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    """
    from core_plot import Bar
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    return Bar(*args, **kwargs)

@docfill
def stick2(*args, **kwargs):
    """stick2(udata, vdata, polar=False, degrees=True, mod=False, pos=None, width=None, scale=None, color='k', line=True, levels=None, cmap=None, shadow=False, **kwargs)

    Plot data as a stick plot and return a :class:`~vacumm.misc.core_plot.Stick` object
    with properties also from :class:`~vacumm.misc.core_plot.Plot`,
    :class:`~vacumm.misc.core_plot.Plot1D`, :class:`~vacumm.misc.core_plot.Curve`,
    :class:`~vacumm.misc.core_plot.ScalarMappable` and
    :class:`~vacumm.misc.core_plot.QuiverKey`.

    :Example:

        >>> stick2(u, v, color='mod', vmax=10., quiver_scale=50., mod=True)
        >>> stick2(r, a, polar=True, degrees=False, quiverkey_value=5)

    :Data params:

        {Stick_load_data[udata]}
        {Stick_load_data[vdata]}
        {Stick_load_data[polar]}
        {Stick_load_data[degrees]}
        {Plot1D__check_order_[vertical]}
        {Plot1D__set_axes_[axis]}
        {Plot[long_name]}
        {Plot[units]}

    :Stick params:

        {Stick_plot[pos]}
        {Stick_plot[mod]}
        {Stick_plot[mod_<param>]}
        {Stick_plot[line]}
        {Stick_plot[alpha]}
        {Stick_plot[headwidth]}
        {Stick_plot[headlength]}
        {Stick_plot[headaxislength]}
        {Stick_plot[minlength]}
        {Stick_plot[minshaft]}
        {Stick_plot[cmap]}
        {Stick_plot[cmap_<param>]}
        {Stick_plot[levels]}
        {Stick_plot[levels_<param>]}
        {Stick_plot[shadow]}
        {Stick_plot[shadow_<param>]}
        {ScalarMappable_post_plot[colorbar]}
        {ScalarMappable_post_plot[colorbar_<param>]}
        {QuiverKey_quiverkey[quiverkey_pos]}
        {QuiverKey_quiverkey[quiverkey_text]}
        {QuiverKey_quiverkey[quiverkey_value]}
        {QuiverKey_quiverkey[quiverkey_units]}
        {QuiverKey_quiverkey[quiverkey_latex_units]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}
        {ScalarMappable_colorbar[cax]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    """
    from core_plot import Stick
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    return Stick(*args, **kwargs)


@docfill
def hov2(*args, **kwargs):
    """hov2(data, order=None, contour=True, fill='pcolor', levels=None, colorbar=True, title=None, xaxis=None, yaxis=None, **kwargs)

    Plot data as a Hovmoller diagram (2D with one axis as time)
    and get a :class:`~vacumm.misc.core_plot.Nov` object.

    :Example:

        >>> hov2(temp[:,-1,:,3], order='ty') # Force time as Y axis
        >>> h = hov2(ssh, show=False)

    :Data params:

        {Plot2D_load_data[data]}
        {Plot__check_order_[order]}
        {Plot2D__set_axes_[x/yaxis]}
        {Plot[long_name]}
        {Plot[units]}

    :2D plot params:

        {Plot2D_plot[contour]}
        {Plot2D_plot_contour[contour_<param>]}
        {Plot2D_plot_fill[fill]}
        {Plot2D_plot_fill[nofill]}
        {ScalarMappable[levels]}
        {Plot2D_plot[levels_<param>]}
        {ScalarMappable[levels_mode]}
        {ScalarMappable[keepminmax]}
        {ScalarMappable[nmax_levels]}
        {ScalarMappable[nmax]}
        {ScalarMappable[cmap]}
        {Plot2D_plot_fill[cmap_<param>]}
        {Plot2D_plot_fill[alpha]}
        {Plot2D_plot_fill[fill_<param>]}
        {Plot2D_plot_fill[shading]}
        {Plot2D_plot_fill[extend]}
        {Plot2D_plot_contour[linewidths]}
        {Plot2D_plot_contour[clabel]}
        {Plot2D_plot_contour[clabel_<param>]}
        {Plot2D_plot_quiver[quiver_<param>]}
        {Plot2D_plot_quiver[quiver_norm]}
        {Plot2D_plot_quiver[quiver_samp]}
        {Plot2D_plot_quiver[quiver_x/ysamp]}
        {Plot2D_plot_quiver[quiver_res]}
        {Plot2D_plot_quiver[quiver_x/yres]}
        {Plot2D_plot_quiver[quiver_relres]}
        {Plot2D_plot_quiver[quiver_x/yrelres]}
        {Plot2D_plot_quiver[quiverkey]}
        {QuiverKey_quiverkey[quiverkey_pos]}
        {QuiverKey_quiverkey[quiverkey_text]}
        {QuiverKey_quiverkey[quiverkey_value]}
        {QuiverKey_quiverkey[quiverkey_units]}
        {QuiverKey_quiverkey[quiverkey_latex_units]}
        {Plot2D_plot_quiver[quiverkey_<param>]}
        {ScalarMappable_post_plot[colorbar]}
        {ScalarMappable_post_plot[colorbar_<param>]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    :More:

        - **Specific params**: see :class:`~vacumm.misc.core_plot.Hov`.
        - **Scalar params**: see :class:`~vacumm.misc.core_plot.ScalarMappable`.
        - **Other generic params**: see :class:`~vacumm.misc.core_plot.Plot`.

    """
    from core_plot import Hov
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    return Hov(*args, **kwargs)

@docfill
def map2(*args, **kwargs):
    """map2(data=None, proj='cyl', res='auto', lon=None, lat=None, contour=True, fill='pcolor', levels=None, colorbar=True, xaxis=None, yaxis=None, title=None, savefig=None, show=True, **kwargs)

    Plot an empty map or data on a map and return a :class:`~vacumm.misc.core_plot.Map`
    object
    with properties also from :class:`~vacumm.misc.core_plot.Plot`,
    :class:`~vacumm.misc.core_plot.Plot2D`,
    :class:`~vacumm.misc.core_plot.ScalarMappable` and
    :class:`~vacumm.misc.core_plot.QuiverKey`.

    :Example:

        >>> map2(xe, resolution='i')
        >>> map2(lon=(-10,0), lat=(45,55), drawrivers=True, drawrivers_color='b')
        >>> m = map(bathy, show=False)
        >>> m.add_place(x, y, 'Brest')
        >>> m.show()

    :Data params:

        {Plot2D_load_data[data]}
        {Map_load_data[lon]}
        {Map_load_data[lat]}
        {Plot2D__set_axes_[x/yaxis]}
        {Plot[long_name]}
        {Plot[units]}
        {Plot[latex_units]}

    :Map params:

        {Map_pre_plot[projection]}
        {Map_pre_plot[resolution]}
        {Map_pre_plot[map_update]}
        {Map_pre_plot[nocache]}
        {Map_pre_plot[zoom]}
        {Map_post_plot[fillcontinents]}
        {Map_post_plot[fillcontinents_<param>]}
        {Map_post_plot[land_color]}
        {Map_post_plot[drawrivers]}
        {Map_post_plot[drawrivers_<param>]}
        {Map_post_plot[meridians]}
        {Map_post_plot[drawmeridians]}
        {Map_post_plot[drawmeridians_<param>]}
        {Map_post_plot[parallels]}
        {Map_post_plot[drawparallels]}
        {Map_post_plot[drawparallels_<param>]}
        {Map_post_plot[meridional/zonal_labels]}
        {Map_post_plot[fullscreen]}
        {Map_post_plot[mapscale]}
        {Map_post_plot[mapscale_<param>]}
        {Map_post_plot[compass]}
        {Map_post_plot[compass_<param>]}
        {Map_post_plot[mscp]}
        {Map_post_plot[mscp_<param>]}

    :2D plot params:

        {Plot2D_plot[contour]}
        {Plot2D_plot_contour[contour_<param>]}
        {Plot2D_plot_fill[fill]}
        {Plot2D_plot_fill[nofill]}
        {ScalarMappable[levels]}
        {Plot2D_plot[levels_<param>]}
        {ScalarMappable[levels_mode]}
        {ScalarMappable[keepminmax]}
        {ScalarMappable[nmax_levels]}
        {ScalarMappable[nmax]}
        {ScalarMappable[cmap]}
        {Plot2D_plot_fill[cmap_<param>]}
        {Plot2D_plot_fill[alpha]}
        {Plot2D_plot_fill[fill_<param>]}
        {Plot2D_plot_fill[shading]}
        {Plot2D_plot_fill[extend]}
        {Plot2D_plot_contour[linewidths]}
        {Plot2D_plot_contour[clabel]}
        {Plot2D_plot_contour[clabel_<param>]}
        {Plot2D_plot_quiver[quiver_<param>]}
        {Plot2D_plot_quiver[quiver_norm]}
        {Plot2D_plot_quiver[quiver_samp]}
        {Plot2D_plot_quiver[quiver_x/ysamp]}
        {Plot2D_plot_quiver[quiver_res]}
        {Plot2D_plot_quiver[quiver_x/yres]}
        {Plot2D_plot_quiver[quiver_relres]}
        {Plot2D_plot_quiver[quiver_x/yrelres]}
        {Plot2D_plot_quiver[quiverkey]}
        {QuiverKey_quiverkey[quiverkey_pos]}
        {QuiverKey_quiverkey[quiverkey_text]}
        {QuiverKey_quiverkey[quiverkey_value]}
        {QuiverKey_quiverkey[quiverkey_units]}
        {QuiverKey_quiverkey[quiverkey_latex_units]}
        {Plot2D_plot_quiver[quiverkey_<param>]}
        {ScalarMappable_post_plot[colorbar]}
        {ScalarMappable_post_plot[colorbar_<param>]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    :Other params:

        - **Specific params**: see :class:`~vacumm.misc.core_plot.Map`.
        - **Specific plot initialization params**: see :class:`~vacumm.misc.core_plot.Map.pre_plot`.
        - **Specific plot params**: see :class:`~vacumm.misc.core_plot.Map.plot`.
        - **Specific plot finalization params**: see :class:`~vacumm.misc.core_plot.Map.post_plot`.
        - **Other generic params**: see :class:`~vacumm.misc.core_plot.Plot`.

    """



    from core_plot import Map
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    if len(args)==0:
        args=[None]
    return Map(*args, **kwargs)

@docfill
def section2(*args, **kwargs):
    """section2(data, contour=True, fill='pcolor', levels=None, colorbar=True, title=None, xaxis=None, yaxis=None, **kwargs)

    Plot geographical data as a vertical-horizontal section
    and return a :class:`~vacumm.misc.core_plot.Section` object
    with properties also from :class:`~vacumm.misc.core_plot.Plot`,
    :class:`~vacumm.misc.core_plot.Plot2D`,
    :class:`~vacumm.misc.core_plot.ScalarMappable` and
    :class:`~vacumm.misc.core_plot.QuiverKey`.

    :Example:

        >>> section2(temp[0,:, :,3])

    :Data params:

        {Plot2D_load_data[data]}
        {Plot2D__set_axes_[x/yaxis]}
        {Plot[long_name]}
        {Plot[units]}

    :2D plot params:

        {Plot2D_plot[contour]}
        {Plot2D_plot_contour[contour_<param>]}
        {Plot2D_plot_fill[fill]}
        {Plot2D_plot_fill[nofill]}
        {ScalarMappable[levels]}
        {Plot2D_plot[levels_<param>]}
        {ScalarMappable[levels_mode]}
        {ScalarMappable[keepminmax]}
        {ScalarMappable[nmax_levels]}
        {ScalarMappable[nmax]}
        {ScalarMappable[cmap]}
        {Plot2D_plot_fill[cmap_<param>]}
        {Plot2D_plot_fill[alpha]}
        {Plot2D_plot_fill[fill_<param>]}
        {Plot2D_plot_fill[shading]}
        {Plot2D_plot_fill[extend]}
        {Plot2D_plot_contour[linewidths]}
        {Plot2D_plot_contour[clabel]}
        {Plot2D_plot_contour[clabel_<param>]}
        {Plot2D_plot_quiver[quiver_<param>]}
        {Plot2D_plot_quiver[quiver_norm]}
        {Plot2D_plot_quiver[quiver_samp]}
        {Plot2D_plot_quiver[quiver_x/ysamp]}
        {Plot2D_plot_quiver[quiver_res]}
        {Plot2D_plot_quiver[quiver_x/yres]}
        {Plot2D_plot_quiver[quiver_relres]}
        {Plot2D_plot_quiver[quiver_x/yrelres]}
        {Plot2D_plot_quiver[quiverkey]}
        {QuiverKey_quiverkey[quiverkey_pos]}
        {QuiverKey_quiverkey[quiverkey_text]}
        {QuiverKey_quiverkey[quiverkey_value]}
        {QuiverKey_quiverkey[quiverkey_units]}
        {QuiverKey_quiverkey[quiverkey_latex_units]}
        {Plot2D_plot_quiver[quiverkey_<param>]}
        {ScalarMappable_post_plot[colorbar]}
        {ScalarMappable_post_plot[colorbar_<param>]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    :More:

        - **Specific params**: see :class:`~vacumm.misc.core_plot.Section`.
        - **Scalar params**: see :class:`~vacumm.misc.core_plot.ScalarMappable`.
        - **Other generic params**: see :class:`~vacumm.misc.core_plot.Plot`.

    """
    from core_plot import Section
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    return Section(*args, **kwargs)

@docfill
def plot2d(*args, **kwargs):
    """plot2d(data, contour=True, fill='pcolor', levels=None, colorbar=True, xaxis=None, yaxis=None, title=None, savefig=None, show=True, **kwargs)

    Generic plot of a 2D variable and return
    a :class:`~vacumm.misc.core_plot.Plot2D` object
    with properties also from :class:`~vacumm.misc.core_plot.Plot`,
    :class:`~vacumm.misc.core_plot.ScalarMappable` and
    :class:`~vacumm.misc.core_plot.QuiverKey`.

    :Example:

        >>> plot2d(xe)

    :Data params:

        {Plot2D_load_data[data]}
        {Map_load_data[lon]}
        {Map_load_data[lat]}
        {Plot2D__set_axes_[x/yaxis]}
        {Plot[long_name]}
        {Plot[units]}
        {Plot[latex_units]}

    :2D plot params:

        {Plot2D_plot[contour]}
        {Plot2D_plot_contour[contour_<param>]}
        {Plot2D_plot_fill[fill]}
        {Plot2D_plot_fill[nofill]}
        {ScalarMappable[levels]}
        {Plot2D_plot[levels_<param>]}
        {ScalarMappable[levels_mode]}
        {ScalarMappable[keepminmax]}
        {ScalarMappable[nmax_levels]}
        {ScalarMappable[nmax]}
        {ScalarMappable[cmap]}
        {Plot2D_plot_fill[cmap_<param>]}
        {Plot2D_plot_fill[alpha]}
        {Plot2D_plot_fill[fill_<param>]}
        {Plot2D_plot_fill[shading]}
        {Plot2D_plot_fill[extend]}
        {Plot2D_plot_contour[linewidths]}
        {Plot2D_plot_contour[clabel]}
        {Plot2D_plot_contour[clabel_<param>]}
        {Plot2D_plot_quiver[quiver_<param>]}
        {Plot2D_plot_quiver[quiver_norm]}
        {Plot2D_plot_quiver[quiver_samp]}
        {Plot2D_plot_quiver[quiver_x/ysamp]}
        {Plot2D_plot_quiver[quiver_res]}
        {Plot2D_plot_quiver[quiver_x/yres]}
        {Plot2D_plot_quiver[quiver_relres]}
        {Plot2D_plot_quiver[quiver_x/yrelres]}
        {Plot2D_plot_quiver[quiverkey]}
        {QuiverKey_quiverkey[quiverkey_pos]}
        {QuiverKey_quiverkey[quiverkey_text]}
        {QuiverKey_quiverkey[quiverkey_value]}
        {QuiverKey_quiverkey[quiverkey_units]}
        {QuiverKey_quiverkey[quiverkey_latex_units]}
        {Plot2D_plot_quiver[quiverkey_<param>]}
        {ScalarMappable_post_plot[colorbar]}
        {ScalarMappable_post_plot[colorbar_<param>]}

    :Plot initialization:

        {Plot_pre_plot[fig]}
        {Plot_pre_plot[figsize]}
        {Plot_pre_plot[subplot]}
        {Plot_pre_plot[subplots_adjust]}
        {Plot_pre_plot[top/bottom/left/right/wspace/hspace]}
        {Plot_pre_plot[axes]}
        {Plot_pre_plot[axes_<param>]}
        {Plot_pre_plot[axes_rect]}
        {Plot_pre_plot[bgcolor]}
        {Plot_pre_plot[twin]}
        {Plot_pre_plot[noframe]}
        {Plot_pre_plot[fullscreen]}

    :Plot finalization:

        {Plot[title]}
        {Plot_post_plot[title_<param>]}
        {Plot[latex_units]}
        {Plot[x/ymin/max]}
        {Plot[x/ymin/max]}
        {Plot_format_axes[x/y/vlim]}
        {Plot_format_axes[x/y/vminmax]}
        {Plot_format_axes[x/y/vmaxmin]}
        {Plot[x/ymasked]}
        {Plot[x/ylong_name]}
        {Plot[x/yunits]}
        {Plot[x/ylabel]}
        {Plot_format_axes[x/yticks]}
        {Plot_format_axes[x/yticklabels]}
        {Plot_format_axes[x/y/vskip]}
        {Plot_format_axes[x/yhide]}
        {Plot_post_plot[grid]}
        {Plot_post_plot[grid_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[figtext]}
        {Plot_post_plot[figtext_<param>]}
        {Plot_post_plot[legend]}
        {Plot_post_plot[legend_<param>]}
        {Plot_post_plot[key]}
        {Plot_post_plot[key_<param>]}
        {Plot_post_plot[anchor]}
        {Plot_post_plot[autoresize]}
        {Plot_post_plot[savefig]}
        {Plot_savefig[savefig_verbose]}
        {Plot_post_plot[savefig_<param>]}
        {Plot_post_plot[savefigs]}
        {Plot_post_plot[savefigs_<param>]}
        {Plot_post_plot[show]}
        {Plot_post_plot[close]}

    :Other params:

        - **Specific params**: see :class:`~vacumm.misc.core_plot.Plot2D`.
        - **Specific plot initialization params**: see :class:`~vacumm.misc.core_plot.Plot.pre_plot`.
        - **Specific plot params**: see :class:`~vacumm.misc.core_plot.Plot2D.plot`.
        - **Specific plot finalization params**: see :class:`~vacumm.misc.core_plot.Plot.post_plot`.
        - **Other generic params**: see :class:`~vacumm.misc.core_plot.Plot`.

    """
    from core_plot import Plot2D
    kwargs.setdefault('plot', True)
    kwargs.setdefault('post_plot', True)
    return Plot2D(*args, **kwargs)


def minimap(gg, bbox= [.85, .85, .14, .14], zoom=1., xmargin=None, ymargin=None,
            lon=None, lat=None, square=False, bgcolor=(0, .8, 1.), fig=None,
            alpha=1, **kwargs):
    """Create a minimap with :func:`map2`

    A minimap is small and generally in a corner of the figure,
    and used to show simple geographic information.

    :Examples:

        >>> minimap(((lonmin,lonmax),(latmin,latmax))).add_point(-5,48)
        >>> minimap(sst.getGrid())
        >>> m = minimap(sst)
        >>> m.add_point(lon, lat)

    :Return: A :class:`~vacumm.misc.core_plot.Map` object.
    """
    from grid import get_xy
    from color import RGB
    data = gg if cdms2.isVariable(gg) else None
    x, y = get_xy(gg)
    if lon is not None:
        x = N.asarray(lon)
    if lat is not None:
        x = N.asarray(lat)
    x = N.asarray(x)
    y = N.asarray(y)
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    xmin, ymin, xmax, ymax = zoombox([xmin, ymin, xmax, ymax], zoom,
        xmargin=xmargin, ymargin=ymargin, square=square)
    kwargs.setdefault('anchor', 'E')
    kwargs.setdefault('colorbar', False)
    kwargs.setdefault('contour', False)
    kwargs.setdefault('fill', 'pcolormesh')
    kwargs.setdefault('cmap', 'jet')
    bgcolor = RGB(bgcolor)
    if alpha:
        bgcolor += alpha,
    oldax = P.gca()
    dict_check_defaults(kwargs, title=False, xhide=True, yhide=True, proj='merc',
        drawparallels_linewidth=.2, drawmeridians_linewidth=.2)
    m = map2(data, lon = (xmin, xmax), lat=(ymin,ymax), show=False,
        axes_rect = bbox, bgcolor=bgcolor, fig=fig, **kwargs)
#    m.axes.set_alpha(alpha)
    fc = m.get_axobj('fillcontinents')
    if alpha!=1 and fc:
        for o in fc: o.set_alpha(alpha)
    m.axes.set_frame_on('off')
    P.sca(oldax)
    return m

def add_map_point(gg, lon, lat, marker='o', color='r', size=40,  m=None, alpha=1,  **kwargs):
    """Add a small map with a point at specified longitude and latitude

    :Params:

        - **gg**: Map limits passed to :func:`minimap`
          or  :meth:`~vacumm.misc.core_plot.Map` instance.
        - **lon/lat**: Coordinates of the point.

    :See also: :func:`minimap` :meth:`~vacumm.misc.core_plot.Plot.add_point`
    """
    kwmap = kwfilter(kwargs, 'map')
    for att in 'bbox', 'bgcolor', 'fig':
        if att in kwargs: kwmap[att] = kwargs.pop(att)
    from core_plot import Map
    if isinstance(gg, Map): m = gg
    if m is None: m = minimap(gg, **kwmap)
    return m.add_point(lon, lat, size=size, color=color, marker=marker, **kwargs)

def add_map_places(gg, lon, lat, txt, marker='o', color='r', size=40,  m=None,
        alpha=1, **kwargs):
    """Add a small map with one or several places at specified longitude and latitude

    :Params:

        - **gg**: Map limits passed to :func:`minimap`
          or  :meth:`~vacumm.misc.core_plot.Map` instance.
        - **lon/lat**: Coordinates arrays of the points.
        - **txt**: Text array associated to the points

    :See also: :func:`minimap` :meth:`~vacumm.misc.core_plot.Plot.add_place`
    """
    kwmap = kwfilter(kwargs, 'map')
    for att in 'bbox', 'bgcolor', 'fig':
        if att in kwargs: kwmap[att] = kwargs.pop(att)
    from core_plot import Map
    if isinstance(gg, Map): m = gg
    if m is None: m = minimap(gg, **kwmap)
    for x,y,text in zip(lon,lat,txt):
      m.add_place(x, y, text, shadow=False, glow=False, **kwargs)
    return m

def add_map_line(gg, extents, color='r', linewidth=1.5, m=None, **kwargs):
    """Add a small map with a line at specified longitudes and latitudes

    :Params:

        - **gg**: Map limits passed to :func:`minimap`
          or  :meth:`~vacumm.misc.core_plot.Map` instance.
        - **extents**: Extents in the forms ``[xmin,ymin,xmax,ymax]``
          ``dict(x=(xmin,xmax),y=xmin,xmax)`` or
              ``dict(lon=(xmin,xmax),lat=xmin,xmax)``.

    :See also: :func:`minimap` :meth:`~vacumm.misc.core_plot.Plot.add_line`
    """
    kwmap = kwfilter(kwargs, 'map')
    for att in 'bbox', 'bgcolor', 'fig':
        if att in kwargs: kwmap[att] = kwargs.pop(att)
    from core_plot import Map
    if isinstance(gg, Map): m = gg
    if m is None: m = minimap(gg, **kwmap)
    return m.add_line(extents, color=color, linewidth=linewidth, **kwargs)

def add_map_lines(gg, xx, yy, color='r', linewidth=1.5, m=None, closed=False, **kwargs):
    """Add a small map with a broken line specified longitudes and latitudes

    :Params:

        - **gg**: Map limits passed to :func:`minimap`
          or  :meth:`~vacumm.misc.core_plot.Map` instance.
        - **xx/yy**: 1D arrays of coordinates.

    :See also: :func:`minimap` :meth:`~vacumm.misc.core_plot.Plot.add_lines`
    """
    kwmap = kwfilter(kwargs, 'map')
    for att in 'bbox', 'bgcolor', 'fig':
        if att in kwargs: kwmap[att] = kwargs.pop(att)
    from core_plot import Map
    if isinstance(gg, Map): m = gg
    if m is None: m = minimap(gg, **kwmap)
    return m.add_lines(xx, yy, color=color, linewidth=linewidth, **kwargs)

def add_map_box(gg, box, color='r', linewidth=1.5, m=None, **kwargs):
    """Add a small map with a box at specified longitudes and latitudes

    :Params:

        - **gg**: Map limits passed to :func:`minimap`
          or  :meth:`~vacumm.misc.core_plot.Map` instance.
        - **box**: Box limits in the forms ``[xmin,ymin,xmax,ymax]``
          ``dict(x=(xmin,xmax),y=xmin,xmax)`` or
          ``dict(lon=(xmin,xmax),lat=xmin,xmax)``.

    :See also: :func:`minimap` :meth:`~vacumm.misc.core_plot.Plot.add_box`
    """
    kwmap = kwfilter(kwargs, 'map')
    for att in 'bbox', 'bgcolor', 'fig':
        if att in kwargs: kwmap[att] = kwargs.pop(att)
    from core_plot import Map
    if isinstance(gg, Map): m = gg
    if m is None: m = minimap(gg, **kwmap)
    return m.add_box(box, color=color, linewidth=linewidth, **kwargs)

#####################################################################
######################################################################
#from grid.misc import meshgrid, meshbounds


