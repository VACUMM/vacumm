# -*- coding: utf8 -*-
"""Variables and utilities about colors and color maps

List of all available colormaps in matplotlib, including VACUMM colormaps
(plotted with :func:`plot_cmaps`)

.. image:: misc-color-cmaps.png

"""
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
import re
import os
import glob
from copy import deepcopy

from matplotlib.cm import ScalarMappable
from matplotlib.colors import (Colormap, Normalize, LinearSegmentedColormap,
    ColorConverter, makeMappingArray, rgb_to_hsv, hsv_to_rgb)
import pylab as P
from matplotlib.cm import cmap_d
import numpy as N
ma = N.ma
import matplotlib.cbook as cbook
from mpl_toolkits.basemap import cm as basemap_cm
try:
    import cmocean.cm as cmoceancm
except:
    cmoceancm = None

import vacumm
from vacumm import vacumm_warn

# Color definitions

#: Land color
bistre = tuple([246./256.,223./256.,180./256.])

#: Alias for land color
land = bistre

#: Alias for land color
jean_pierre = bistre

#: Alias for land color
decimaitre = jean_pierre

#: Ocean color
ocean = (.6, .8, 1)

#: Alias for ocean color
sea = ocean

#: Basic list of colors
simple_colors =  ['k', 'b', 'g', 'r', 'c', 'm', 'y']
#colors = simple_colors

_cpt_dir = os.path.join(os.path.dirname(__file__), 'cpt')

RGB = ColorConverter().to_rgb
RGBA = ColorConverter().to_rgba

__all__ = ['cmap_custom', 'cmap_bwr', 'cmap_bwre', 'cmap_br', 'cmap_wr', 'cmap_wre', 'cmap_bathy',
    'cmap_jete', 'cmap_ajete', 'cmap_jets', 'cmap_wjets', 'cmap_ajets', 'cmap_smoothed_steps',
    'cmap_smoothed_regular_steps', 'cmap_ss', 'cmap_steps', 'cmap_regular_steps', 'cmap_rs',
    'cmap_wjet', 'cmap_pe', 'cmap_grey', 'show_cmap', 'cmaps_mpl', 'cmaps_gmt', 'cmap_gmt',
    'cmaps_vacumm', 'print_cmaps_gmt',  'darken', 'whiten', 'to_shadow', 'StepsNorm',
    'bistre', 'land', 'jean_pierre', 'decimaitre', 'RGB', 'RGBA', 'cmap_srs', 'cmap_rs',
    'cmap_linear', 'ocean', 'sea', 'plot_cmap', 'plot_cmaps', 'get_cmap', 'simple_colors',
    'cmap_land', 'cmap_topo', 'auto_cmap_topo', 'cmap_jet', 'cmap_rainbow', 'rainbow',
    'cmap_magic', 'cmap_mg', 'RangedLinearSegmentedColormap','Scalar2RGB', 'cmap_chla',
    'cmap_previmer', 'cmap_rnb2_hymex','cmap_rainbow_sst_hymex','cmap_dynamic_cmyk_hymex',
    'cmap_white_centered_hymex','cmap_red_tau_hymex', 'cmap_previmer2', 'cmap_ssec',
    'cmap_ncview_rainbow', 'cmap_eke', 'cmap_rb', 'cmap_currents', 'cmaps_registered',
    'cmap_br2','cmap_nice_gfdl', 'anamorph_cmap', 'discretise_cmap', 'to_grey',
    'change_luminosity', 'saturate', 'desaturate', 'change_saturation',
    'change_value', 'pastelise'
    ]
__all__.sort()

# Color maps
# ----------
def rainbow(n=None, mode='auto', first=None, last=None, middle=None, stretcher=None):
    """Get a list of nice rainbow colors

    :Params:

        - **n**, optional: Number of requested colors. It should be > 1. Defaults to ``None``
        - **mode**, optional: Mode for choosing extreme colors.

            - ``"strict"``: Blue and red are first and last colors.
            - ``"first"``: ``(.5, 0, 1)`` is inserted as first color.
            - ``"last"``: ``(1, 0, .5)`` is inserted as last color.
            - ``"auto"`` or None: Switch to ``"extended"`` if ``n>6``.
            - ``"extended"`` or anything else: Like ``"first"`` and ``"last"``.
        - **stretcher**, optional: Function that alter color steps
          from their range [0.->1.] to [0.->1.]. This must be used
          for example to enhance a part of the colormap.
          Its value may also be a predefined function:

            - ``"reduced_green"``: Shrink the green band.

          .. note::

            This function must return a monotonic array of increasing
            values included in the [0,1] interval.


    """
    # Special cases
    if n is None: n=5
    if n==0: return []
    if n==1: return ['.5']

    # Base colors
    rbcols = [(0., 0., 1.), (0., 1., 1.), (0., 1., 0.), (1., 1., 0.), (1., 0., 0.)]
    rbfirst = (.5, 0., 1.)
    rblast = (1., 0., .5)
    rbcyclic = (1., 0, 1.)

    # Mode
    if mode is None: mode = 'auto'
    if mode == 'auto':
        if n>6:
            mode = 'extended'
        else:
            mode = 'strict'
    if mode!='strict':
        if mode=='last':
            rbcols =rbcols+[rblast]
        elif mode=='first':
            rbcols = [rbfirst]+rbcols
        elif mode=='cyclic':
            rbcols = [rbcyclic]+rbcols+[rbcyclic]
        else:
            rbcols = [rbfirst]+rbcols+[rblast]

    # Rainbow colors expanded to n colors
    base_cmap = cmap_custom(_regular_(rbcols, steptype='bounds'), ncol=512)
    positions = N.linspace(0., 1., n)
    if stretcher in ['reduced_green']:
        stretcher = eval('_stretcher_%s_'%stretcher)
    if callable(stretcher):
        positions = N.clip(stretcher(positions), 0, 1)
    colors = [base_cmap(p)[:3] for p in positions]

    # Special values
    if first is not None:
        colors[0] = first
    if last is not None:
        colors[-1] = last
    if middle is not None and n%2:
        colors[int(n/2)] = middle

    return colors

def _stretcher_reduced_green_(x):
    xx = (x-.5)*20
    return .2*P.exp(-N.abs(xx))*xx+x

def _nlev_(n):
    if hasattr(n,  '__len__'):
        lev = n
        n = max(0, len(n)-1)
    else:
        lev = None
    if n is None: n = 5
    elif n==0: n = 1
    return n, lev

def cmap_rainbow(n=None, name='vacumm_rainbow', smoothed=True, mode='auto', **kwargs):
    """Rainbow colormap

    :Params:

        - **n**, optional: Number of colors used (keyword passed to :func:`rainbow`).
        - **stretcher**, **first**, **last**, **middle**, optional: See :func:`rainbow`.
        - **rainbow_<param>**, optional: ``<param>`` is passed to func:`rainbow`
        - Extra keywords are passed to :func:`cmap_rs` or :func:`cmap_srs`

    :Example:

        >>> cmap = cmap_magic(10, first='.5', mode='last')
        >>> cmap = cmap_magic([5.,6.,7,8,9], stretcher='reduced_green')


    :Sample: .. image:: misc-color-vacumm_rainbow.png
    """
    from misc import kwfilter
    kwrb = kwfilter(kwargs, 'rainbow_')
    for att in 'first', 'last', 'middle', 'stretcher':
        kwrb[att] = kwargs.pop(att, None)
    n, pos = _nlev_(n)
    colors = rainbow(n=n, mode=mode, **kwrb)
    if False and pos is not None:
        func = cmap_ss if smoothed else cmap_steps
        norm = kwargs.pop('norm', None)
        if norm is None:
            norm = Normalize(vmin=pos[0], vmax=pos[-1])
        colors = [(c, norm(p)) for (p, c) in zip(pos, colors)]
    else:
        func = cmap_srs if smoothed else cmap_rs
    return func(colors, name=name, **kwargs)
def cmap_rb(*args, **kwargs):
    """Shortcut to :func:`cmap_rainbow`"""
    return cmap_rainbow(*args, **kwargs)

def cmap_magic(n=None, stretch = 0.4, mode='normal', white='.95', name='vacumm_magic', **kwargs):
    """Magic rainbow colormap

    Based upon :func:`cmap_rainbow` with specificities:

        - ``n`` may be for example an array of levels.
        - ``smoothed`` is set to ``False`` by default.
        - ``stretch`` is set negative by default (``-attenuation``).

    :Params:

        - **n**, optional: Number of colors used or
          object with a length (like an array of levels).
        - Extra keywords are passed to :func:`rainbow`

    :Example:

        >>> ssta = MV2.arange(-2.6, 3.1, .1)
        >>> cmap = cmap_magic(ssta, anomaly=True)
        >>> cmap = cmap_magic(len(ssta))

    :Samples:

        - .. image:: misc-color-vacumm_magic.png
        - .. image:: misc-color-vacumm_magic-n10.png
        - .. image:: misc-color-vacumm_magic-anom.png
        - .. image:: misc-color-vacumm_magic-pos.png
        - .. image:: misc-color-vacumm_magic-neg.png
    """
    #print 'cmap_magic n in',n
    n, pos = _nlev_(n)
    kwargs.setdefault('smoothed', False)
    from misc import broadcast
    for att in 'positive', 'negative', 'anomaly',  'symetric': # compat
        if kwargs.has_key(att):
            if kwargs[att]:
                mode = att
            del kwargs[att]
    if mode.startswith('anom'): mode = 'symetric'
    if mode.startswith('norm'):
        kwargs['lstretch'] = -stretch
    elif mode.startswith('pos'):
        kwargs['first'] = white
        kwargs['lstretch'] = [0]+broadcast(-stretch, n-1)
        kwargs['rstretch'] = [-stretch]+[0]*(n-1)
    elif mode.startswith('neg'):
        kwargs['last'] = white
        kwargs['lstretch'] = [0]*(n-1)+[-stretch]
        kwargs['rstretch'] = broadcast(-stretch, n-1)+[0]
    elif mode.startswith('sym'):
        keepc = broadcast(None, n)
        mid = int(n/2)
        if n%2: # odd
            kwargs['lstretch'] = [0]*mid+broadcast(-stretch, mid+1)
            kwargs['rstretch'] = kwargs['lstretch'][::-1]
            kwargs['middle'] = '.8'
            keepc[mid] = False
        else:
            kwargs['lstretch'] = [0]*(mid)+[.98]+broadcast(-stretch, mid-1)
            kwargs['rstretch'] = kwargs['lstretch'][::-1]
            keepc[mid:mid+2] = False, False
    #print 'cmap_magic n out',n
    if pos is not None:
        n = pos
    return cmap_rainbow(n, name=name, **kwargs)

def cmap_mg(*args, **kwargs):
    """Shortcut to :func:`cmap_magic`"""
    return cmap_magic(*args, **kwargs)

def cmap_custom(colors, name='mycmap', ncol=256, ranged=False, register=True, **kwargs):
    """Quick colormap creation

    :Params:

        - **colors**: Like [(color1, position1),(color2, position2), etc...] or
          dict(red=((pos1,r1a,r1b), (pos2,r2a,r2b)),etc...)
        - **ncol/N**: Discretization.
        - **register**: Register the colormap into matplotlib to be accessible with its name?
    """
    ncol = kwargs.pop('N', ncol)
    if isinstance(colors, dict):
        for cname, cvals in colors.items():
            cvals = tuple(cvals)
            if cvals[0][0] != 0:
                colors[cname] = ((0, cvals[0][1], cvals[0][2]), ) + cvals
            if cvals[-1] != 1:
                colors[cname] = cvals + ((1, cvals[-1][1], cvals[-1][2]), )
    else:
        if isinstance(colors, tuple): colors = list(colors)
        tcolors = colors
        colors = dict(red=(), green=(), blue=())
        if tcolors[0][1] != 0:
            tcolors.insert(0, (tcolors[0][0], 0))
        if tcolors[-1][1] != 1:
            tcolors.append((tcolors[-1][0], 1))
        for col, pos in tcolors:
            r, g, b = RGB(col)
            colors['red'] += ((pos, r, r), )
            colors['green'] += ((pos, g, g), )
            colors['blue'] += ((pos, b, b), )
    if ranged:
        cmap = RangedLinearSegmentedColormap(name,colors,ncol, **kwargs)
    else:
        cmap = LinearSegmentedColormap(name,colors,ncol)
    P.register_cmap(name, cmap)
    return cmap

def cmap_linear(c, **kwargs):
    """Basic linear map between colors"""
    kwargs['stretch'] = 0
    return cmap_srs(c, **kwargs)

# blue -> white -> red
def cmap_bwr(wpos=0.5, wcol='w', name='vacumm_bwr'):
    """Blue->white->red colormap

    :Params:

        - *white*: relative position of white color in the map [default: 0.5]

    :Sample: .. image:: misc-color-vacumm_bwr.png
    """
    return cmap_custom((('b', 0), (wcol, wpos), ('r', 1)), name)


# blue -> white -> red for Extremes : violet and pink and ends
def cmap_bwre(wpos=0.5,gap=0.1, wcol='w', name='vacumm_bwre', **kwargs):
    """Returns a violet->blue->white->red colormap->yellow

    - **white**, optional: relative position of white color in the map [default: 0.5]
    - **gap**, optional: Relative width of the pure central white gap [default: 0.1]

    :Sample: .. image:: misc-color-vacumm_bwre.png
    """
    wstart = wpos*(1-gap)
    wstop = wpos + (1.-wpos)*gap
    return  cmap_custom((((1, 0, 1), 0), ('b', .1), (wcol, wstart),
        (wcol, wstop), ('r', .9), ((1, 1, 0), 1)), name, **kwargs)

# blue -> red
def cmap_br(sep=0.5, name='vacumm_br', **kwargs):
    """Blue->red colormap

    :Params:

        - **sep**, optional: relative position of the blue/red transition in the map [default: 0.5]

    :Sample: .. image:: misc-color-vacumm_br.png
    """
    return cmap_custom( (('b', 0), ((.5, 0, .5), sep), ('r', 1)), name, **kwargs)

def cmap_br2(sep=0.5, white=.7, name="vacumm_br2", **kwargs):
    """Blue->light blue|light red-> red

    :Params:

        - **sep**, optional: relative position of the blue/red transition in the map [default: 0.5]
        - **white**, optional: Strenght of the whitening (see :func:`whiten`).

    :Sample: .. image:: misc-color-vacumm_br.png
    """
    bsep = whiten('b', white)
    rsep = whiten('r', white)
    cdict =   {'red':   ((0,1,1),(sep,bsep[0],rsep[0]), (1,0,0)),
               'green': ((0,0,0),(sep,bsep[1],rsep[1]), (1,0,0)),
               'blue':  ((0,0,0),(sep,bsep[2],rsep[2]), (1,1,1)),
               }
    return cmap_custom(cdict, name, **kwargs)

# white -> red
def cmap_wr(name='vacumm_wr', **kwargs):
    """White->red colormap

    :Sample: .. image:: misc-color-vacumm_wr.png
    """
    return cmap_custom((('w', 0), ('r', 1)), name, **kwargs)

# white -> red
def cmap_wre(name='vacumm_wre', **kwargs):
    """White->red->yellow colormap for positive extremes

    :Sample: .. image:: misc-color-vacumm_wre.png
    """
    cdict = {'red':  ((0.,1.,1.),(.8,1.,1.),(1.,1.,1.)),
           'green':((0.,1.,1.),(.8,0.,0.),(1.,1.,1.)),
           'blue': ((0.,1.,1.),(.8,0.,0.),(1.,0.,0.))}
    return cmap_custom(cdict, name, **kwargs)

cy = 0.6 # cyan
ye = 0.95 # yellow
vi = 0.35 # violet
_colors_bathy = (('k',0), ((0,.1,.9),vi), ((0,.8,.8),cy), ((.9,.9,.6),.9), ((.9,.9,.8),1))
def cmap_bathy(start=0., stop=1., name='vacumm_bathy', **kwargs):
    """Colormap for bathymetry maps

    :Sample: .. image:: misc-color-vacumm_bathy.png
    """
    this_cmap = cmap_custom(_colors_bathy, name, ranged=True, start=start, stop=stop, **kwargs)
    this_cmap.set_bad(bistre)
    this_cmap.set_under(bistre)
    this_cmap.set_over(bistre)
    return this_cmap

_colors_land = ((_colors_bathy[-1][0], 0), ('#145D0A',.2),('#62CF60', .4), ('#A1A156',.6),('#6D461D', .85), ('#eeeeff', 1.))
def cmap_land(start=0., stop=1., name='vacumm_land', **kwargs):
    """Colormap for land maps

    :Params:

        - **start/stop**, optional: Positions for a zoom in colormap.

    :Sample: .. image:: misc-color-vacumm_land.png
    """
    this_cmap = cmap_custom(_colors_land, name, ranged=True, start=start, stop=stop, **kwargs)
    this_cmap.set_bad(ocean)
    this_cmap.set_under(ocean)
    this_cmap.set_over(ocean)
    return this_cmap

def cmap_topo(start=0., stop=1., name='vacumm_topo', zero=.5, over=_colors_bathy[0][0],
    under=_colors_land[-1][0], bad='0.5', **kwargs):
    """Colormap for bathy+land maps

    :Params:

        - **start/stop**, optional: Positions for a zoom in colormap.
        - **zero**, optional: Position of 0-depth color.
        - **over**, optional: Color for values over maximal range
          (see :meth:`~matplotlib.colors.Colormap.set_over`).
        - **under**, optional: Color for values under maximal range
          (see :meth:`~matplotlib.colors.Colormap.set_under`).
        - **bad**, optional: Color for bad values
          (see :meth:`~matplotlib.colors.Colormap.set_bad`).

    :Sample: .. image:: misc-color-vacumm_topo.png

    :See also: :func:`cmap_bathy` :func:`cmap_land`
    """
    mycolors = tuple([(c[0], zero*c[1]) for c in _colors_bathy])
    mycolors += tuple([(c[0], zero+(1-zero)*c[1]) for c in _colors_land[1:]])
    this_cmap = cmap_custom(mycolors, name, ranged=True, start=start, stop=stop, **kwargs)
    this_cmap.set_bad(bad)
    this_cmap.set_under(under)
    this_cmap.set_over(over)
    return this_cmap

def auto_cmap_topo(varminmax=(0., 1.), gzoom=1., xzoom=1., **kwargs):
    """Adjusted :func:`cmap_topo` colormap so that altitude 0 fall on center of full colormap by defaut

    :Params:

        - **varminmax**: Define the min and max of topo ; it can be either :

            - a topo variable (array) from which min and max are computed
            - a tuple of (min, max) fo topo

        - **gzoom**, optional: Global zoom in the colormap (>1.)
        - **xzoom**, optional: Zoom out

            - the land colormap if more ocean,
            - the ocean colormap if more land.
    """
    # Data range
    if isinstance(varminmax, tuple):
        topomin, topomax = varminmax
    else:
        topomin = N.asarray(varminmax, 'f').min()
        topomax = N.asarray(varminmax, 'f').max()
    topomin = float(min(topomin, 0))
    topomax = float(max(topomax, 0))
    # More land or ocean?
    xzoom = N.clip(xzoom, 0., 1.)
    if abs(topomin) < topomax:
        # More land
        alpha = -topomin/topomax
        if xzoom<alpha: xzoom=alpha
        zero = xzoom/(1.+xzoom)
        start = zero*(1-alpha)
        stop = 1.
    else:
        # More ocean
        alpha = -topomax/topomin
        if xzoom<alpha: xzoom=alpha
        zero = 1/(1.+xzoom)
        start = 0.
        stop = zero*(1+alpha)
    # Global zoom
    if gzoom<1: gzoom=1.
    start = zero - (zero-start)/gzoom
    stop = zero + (stop-zero)/gzoom
    kwargs.setdefault('register', False)
    return cmap_topo(start, stop, 'vacumm_auto_cmap_topo', zero=zero, **kwargs)



def cmap_jete(name='vacumm_jete', **kwargs):
    """Jet colormap with extremes

    :Sample: .. image:: misc-color-vacumm_jete.png
    """
    cdict =   {'red':((0.,.75,.75),(0.22,0,0),  (0.4,0,0),(0.6,1,1),        (.78,1,1),  (1,1,1)),
               'green':((0.,0.,0.),(0.22,0,0),  (0.375,1,1),(0.64,1,1), (.78,0,0),  (1,0,0)),
               'blue':((0.,1.,1.),(0.22,1,1),   (0.4,1,1),(0.6,0,0),        (.78,0,0),  (1,.75,.75))}
    return cmap_custom(cdict, name, **kwargs)

def cmap_ajete(w=.02, name='vacumm_ajete', **kwargs):
    """Jet colormap with white at center and extremes, for anomaly plots

    :Sample: .. image:: misc-color-vacumm_ajete.png
    """
    cdict =   {'red':((0.,.75,.75),(0.21,0,0),  (0.33,0,0),(.5-w/2,1,1),(.5+w/2,1,1),(0.67,1,1),        (.79,1,1),  (1,1,1)),
               'green':((0.,0.,0.),(0.21,0,0),  (0.3,1,1),(.5-w/2,1,1),(.5+w/2,1,1),(0.71,1,1),     (.79,0,0),  (1,0,0)),
               'blue':((0.,1.,1.),(0.21,1,1),   (0.33,1,1),(.5-w/2,1,1),(.5+w/2,1,1),(0.67,0,0),        (.79,0,0),  (1,.75,.75))}
    return cmap_custom(cdict, name, **kwargs)

def cmap_jet(smoothed=False, name='vacumm_jet', **kwargs):
    """Jet colormap"""
    colors = [(.75,0,1),(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0)]
    kwargs.setdefault('stretch', 0)
    if not smoothed:
        cmap = cmap_regular_steps(colors, name=name, **kwargs)
    else:
        cmap = cmap_srs(colors, name=name, **kwargs)
    return cmap

def cmap_jets(name='vacumm_jets', **kwargs):
    """Jet colormap with smoothed steps

    :Params: Passed to :func:`cmap_smoothed_regular_steps`

    :Samples:

        - .. image:: misc-color-vacumm_jets.png
        - .. image:: misc-color-vacumm_jets+60.png
        - .. image:: misc-color-vacumm_jets-60.png

    :See also: :func:`cmap_smoothed_regular_steps`
    """
    colors = [(.75,0,1),(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0)]
    kwargs.setdefault('stretch', 0)
    return cmap_smoothed_regular_steps(colors,name=name,**kwargs)

def cmap_wjets(wcol=".95", name='vacumm_wjets', **kwargs):
    """White jet colormap with smoothed steps

    :Params:

        - **wcol**, optional: color of "white"
        - other options are passed to :func:`cmap_smoothed_regular_steps`

    :Sample: .. image:: misc-color-vacumm_wjets.png
    """
##  colors = [(1,1,1),(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0),(1,0,.75)]
    wcol = kwargs.get('first_color', wcol)
    colors = [RGB(wcol),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0),(1,0,.75)]
    kwargs.setdefault('stretch', 0)
#   return cmap_linear(colors,name='cmap_wjets',**kwargs)
    return cmap_smoothed_regular_steps(colors, name=name, **kwargs)

def cmap_ajets(wcol="w", name='vacumm_ajets', **kwargs):
    """Jet colormap with smoothed steps and white at center (for anomalies)

    :Params:

        - **wcol**, optional: color of "white"
        - other options are passed to :func:`cmap_smoothed_regular_steps`

    :Sample: .. image:: misc-color-vacumm_ajets.png
    """
##  colors = [(.75,0,1),(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0)]
    colors = [(.75,0,1),(0,0,1),(0,1,1),RGB(wcol),(1,1,0),(1,.5,0),(1,0,0)]
    kwargs.setdefault('stretch', 0)
    return cmap_smoothed_regular_steps(colors,name=name,**kwargs)

def _stretch_(value, stretch):
    """ 0 = no stretch
    -1 = black
    1 = white
    """
    stretch = N.clip(stretch, -1, 1)
    if stretch == 0: return value
    if stretch < 0:
        return (1+stretch) * value
    return value + stretch * (1-value)

def cmap_smoothed_steps(colors, stretch=None, rstretch=0, lstretch=0, name='vacumm_css',
    asdict=False, **kwargs):
    """Smoothed steps

    :Params:

        - **colors**, optional: Central positions for each colors [(col1,pos1),...]

   :See also:

        :func:`cmap_steps` :func:`cmap_regular_steps` :func:`cmap_smoothed_regular_steps`
    """
    if stretch is not None: # compat
        if isinstance(stretch, tuple):
            if len(stretch)>1:
                lstretch, rstretch = stretch[:2]
            else:
                lstretch = stretch[0]
    from misc import broadcast
    ns = len(colors)
    lstretch = broadcast(lstretch, ns)
    rstretch = broadcast(rstretch, ns)
    rr = [] ; gg = [] ; bb = []
    for i,(col,pos) in enumerate(colors):
        r, g, b = RGB(col)
        if i == 0:
            rr = [(0,r,r)]
            gg = [(0,g,g)]
            bb = [(0,b,b)]
        else:
            pcol,ppos = colors[i-1]
            pr,pg,pb = RGB(pcol)
            npos = .5*(pos+ppos)
            rr.append((npos,_stretch_(.5*(r+pr),lstretch[i]),_stretch_(.5*(r+pr),rstretch[i])))
            gg.append((npos,_stretch_(.5*(g+pg),lstretch[i]),_stretch_(.5*(g+pg),rstretch[i])))
            bb.append((npos,_stretch_(.5*(b+pb),lstretch[i]),_stretch_(.5*(b+pb),rstretch[i])))
        if i == (len(colors)-1):
            pos = 1
        rr.append((pos,r,r))
        gg.append((pos,g,g))
        bb.append((pos,b,b))

    cdict = {'red':rr,'green':gg,'blue':bb}
    if asdict: return cdict
    return cmap_custom(cdict, name, **kwargs)

def cmap_ss(*args, **kwargs):
    """Shortcut to :func:`cmap_smoothed_steps`"""
    return cmap_smoothed_steps(*args, **kwargs)



def cmap_smoothed_regular_steps(colors, steptype='center', **kwargs):
    """Smoothed regular steps

    :Params:

        - **colors**: [(r1,g1,b1),...]
        - Other keywords are passed to :func:`cmap_smoothed_steps()`

    :See also:

        :func:`cmap_steps` :func:`cmap_smoothed_steps` :func:`cmap_regular_steps`
        (:func:`cmap_regular_jets` for an example)
    """
    return cmap_smoothed_steps(_regular_(colors,steptype=steptype),**kwargs)

def cmap_srs(*args, **kwargs):
    """Shortcut to :func:`cmap_smoothed_regular_steps`"""
    return cmap_smoothed_regular_steps(*args, **kwargs)


def _regular_(colors, steptype='center', posmin=0., posmax=1.):
    rcolors = []
    dpos = posmax-posmin*1.
    step = dpos/len(colors)
    if steptype == 'center':
        offset = .5
    else:
        if steptype == 'bounds':
            step = dpos/(len(colors)-1)
        offset = 0
    for i,col in enumerate(colors):
        col = RGB(col)
        pos = posmin+(i+offset)*step
#       pos = int(255*pos)/256.
        rcolors.append((col,pos))
    return rcolors

def cmap_steps(cols, stretch=None, lstretch=0., rstretch=0., keepc=None, name='cmap_steps',
    **kwargs):
    """Colormap by steps

    :Params:

        - **cols**: [(col1,pos1),(col2,pos2),...]
        - **lstretch**, optional: Color darkening (<0) or whitening at start of steps (left)
        - **rstretch**, optional: Same but at end of steps (right)
        - **keepc**, optional: If ``lstretch`` and ``rstretch`` are both different from zero,
          it keeps center of steps intact, else it becomes a mean
          of left and right.
        - **ncol/N**: Discretization.

    .See also:

        :func:`cmap_regular_steps` :func:`cmap_smoothed_steps` :func:`cmap_smoothed_regular_steps`
    """
    if stretch is not None: # compat
        if isinstance(stretch, tuple):
            if len(stretch)>1:
                lstretch, rstretch = stretch[:2]
            else:
                lstretch = stretch[0]
    from misc import broadcast
    ns = len(cols)
    lstretch = broadcast(lstretch, ns)
    rstretch = broadcast(rstretch, ns)
    keepc = broadcast(keepc, ns)
    rr = [] ; gg = [] ; bb = []
    pcol, ppos = cols[0]
    pr, pg, pb = RGB(pcol)
    pr = _stretch_(pr, lstretch[0])
    pg = _stretch_(pg, lstretch[0])
    pb = _stretch_(pb, lstretch[0])
    apos = [p for (c, p) in cols]+[1.]
    for i,(col,pos) in enumerate(cols):
        r, g, b = RGB(col)
        if i == 0:
            pos = 0.
        else:
            pr = _stretch_(pr, rstretch[i-1])
            pg = _stretch_(pg, rstretch[i-1])
            pb = _stretch_(pb, rstretch[i-1])
        rr.append((pos, pr, _stretch_(r, lstretch[i])))
        gg.append((pos, pg, _stretch_(g, lstretch[i])))
        bb.append((pos, pb, _stretch_(b, lstretch[i])))
        if keepc[i] is None:
            keepc[i] = lstretch[i] != 0 and rstretch[i] != 0
        if keepc[i]:
            rr.append((.5*(pos+apos[i+1]), r, r))
            gg.append((.5*(pos+apos[i+1]), g, g))
            bb.append((.5*(pos+apos[i+1]), b, b))
        pr, pg, pb = r,g,b
    r = _stretch_(r, rstretch[-1])
    g = _stretch_(g, rstretch[-1])
    b = _stretch_(b, rstretch[-1])
    rr.append((1, r, r))
    gg.append((1, g, g))
    bb.append((1, b, b))
    cdict = {'red':rr,'green':gg,'blue':bb}
    return cmap_custom(cdict, name, **kwargs)

def cmap_regular_steps(colors, steptype='stair', **kwargs):
    """Colormap by regular steps

    :Params:

        - **cols**: [col1,col2,...]
        - **stretch**, optional: Color darkening (<0) or whitening within steps [default: -.6]
        - **steptype**, optional: 'center', 'stair' or 'bounds' [default: 'center']
        - Other keywords are passed to :func:`cmap_steps()`

    :See also:

        :func:`cmap_steps` :func:`cmap_smoothed_steps` :func:`cmap_smoothed_regular_steps`
    """
    return cmap_steps(_regular_(colors, steptype=steptype), **kwargs)

def cmap_rs(*args, **kwargs):
    """Shortcut to :func:`cmap_regular_steps`"""
    return cmap_regular_steps(*args, **kwargs)

def cmap_wjet(wcol='.95', smoothed=True, name='vacumm_wjet', **kwargs):
    """jet colormap with white (or another color) at beginning

    :Sample: .. image:: misc-color-vacumm_wjet.png
    """
    wcol = kwargs.get('first_color', wcol)
    func = cmap_regular_steps if not smoothed else cmap_smoothed_regular_steps
    colors = [wcol, 'b', (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0), (1, 0, .74)]
    return func(colors, name=name, **kwargs)

#    wcol = kwargs.get('first_color', wcol)
#    colors = [RGB(wcol),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0),(1,0,.75)]



def cmap_pe(red=.8, name="vacumm_pe", **kwargs):
    """Colomap for positive extremes (white->grey->red)

    :Sample: .. image:: misc-color-vacumm_pe.png
    """
    cdict = {'red':  ((0.,1.,1.),(red,.7,.7),(1.,1.,1.)),
           'green':((0.,1.,1.),(red,.7,.7),(1.,0.,0.)),
           'blue': ((0.,1.,1.),(red,.7,.7),(1.,0.,0.))}
    return cmap_custom(cdict, name, **kwargs)

def cmap_grey(start=0, end=1., name="vacumm_grey", **kwargs):
    """Grey colormap from position ``start`` to position ``end``

    :Sample: .. image:: misc-color-vacumm_grey.png
    """
    cdict = {'red':  ((0.,start,start),(1.,end,end)),
           'green':((0.,start,start),(1.,end,end)),
           'blue': ((0.,start,start),(1.,end,end))}
    return cmap_custom(cdict, name, **kwargs)

def cmap_chla(name='vacumm_chla', smoothed=True, **kwargs):
    """Colormap for Chlorophyll A

    :Sample: .. image:: misc-color-vacumm_chla.png

    :Source: IDL code from F. Gohin.
    """
    # From?
#    rgb =  [(144, 8, 240), (112, 8, 240), (80, 8, 240), (8, 8, 240), (8, 8, 224), (8, 8, 208), (8, 8, 192), (8, 8, 176), (8, 8, 152), (8, 32, 120), (8, 56, 104), (8, 80, 104), (8, 104, 112), (8, 96, 120), (8, 104, 152), (8, 120, 160), (8, 136, 160), (8, 152, 168), (8, 168, 168), (8, 184, 184), (8, 200, 200), (8, 216, 216), (8, 224, 224), (8, 232, 232), (8, 240, 240), (8, 240, 216), (8, 232, 200), (8, 232, 144), (8, 224, 112), (8, 216, 120), (8, 200, 120), (8, 184, 112), (8, 176, 80), (8, 160, 80), (8, 144, 80), (8, 128, 80), (8, 144, 56), (8, 136, 8), (56, 152, 8), (56, 168, 8), (72, 184, 8), (88, 200, 8), (112, 208, 8), (136, 216, 8), (160, 224, 8), (192, 232, 8), (232, 232, 8), (224, 224, 8), (224, 208, 8), (224, 184, 8), (224, 160, 8), (224, 136, 8), (224, 112, 8), (224, 88, 8), (224, 64, 8), (224, 40, 8), (216, 8, 8), (200, 8, 8), (176, 8, 8), (152, 8, 8), (128, 8, 8), (104, 8, 8)]
#    cc = [(r/255.,  g/255., b/255.) for r, g, b in rgb]
#    colors = _regular_(cc, steptype='bounds', posmax=220./256)
#    colors.append(((200./255, )*3, 1.))

    # From F. gohin (IDL)
    param_red = N.array([0.00,0.00,0.00,0.00,0.00,0.53,0.45,0.35,0.48,
        0.43,0.53,0.67,0.80,0.50,0.60,0.72,0.80,0.87,
        0.43,0.60,0.12,0.00,0.00,0.00,0.00,0.00,0.00,
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.78,0.78,
        0.62,0.62,0.62,0.62,0.62,0.77,0.87,1.00,1.00,
        0.95,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,
        1.00,1.00,1.00,0.87,0.75,0.60,0.53,0.65,0.72,
        0.82,0.92,1.00,1.00,1.00,1.00,0.87,0.80,0.72,
        0.65,0.60,0.57,0.57,0.93,0.80,0.70,0.50,0.40])

    param_green = N.array([0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        0.45,0.55,0.70,0.80,0.80,0.87,0.95,0.97,0.97,
        0.97,0.97,0.97,0.97,0.97,0.97,0.97,1.00,0.88,
        0.75,0.60,0.55,0.65,0.75,0.85,0.93,0.78,0.96,
        0.70,0.70,0.70,0.70,0.00,0.88,0.95,1.00,1.00,
        0.95,1.00,1.00,0.95,0.90,0.83,0.75,0.65,0.58,
        0.50,0.30,0.00,0.00,0.00,0.00,0.38,0.50,0.58,
        0.68,0.78,0.85,0.00,0.65,0.80,0.72,0.65,0.57,
        0.50,0.35,0.28,0.17,0.92,0.80,0.70,0.50,0.40])

    param_blue = N.array([0.55,0.67,0.80,0.90,0.00,0.97,0.87,0.77,0.77,
        0.82,0.92,0.97,0.97,0.97,0.97,0.97,0.97,0.97,
        0.97,0.97,0.90,0.80,0.62,0.48,0.22,0.00,0.00,
        0.00,0.00,0.53,0.57,0.67,0.75,0.82,0.70,0.60,
        0.04,0.16,0.38,0.52,0.62,0.40,0.48,0.62,0.80,
        0.82,0.42,0.00,0.72,0.65,0.62,0.50,0.38,0.20,
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        0.00,0.00,0.00,0.00,0.00,0.00,0.75,0.65,0.70,
        0.60,0.50,0.40,0.18,0.93,0.80,0.70,0.50,0.40])

    indtab = N.array([51,127,124,118,117,116,114,113,
        112,108,109,107,106,94,98,99,100,101,
        102,104,92,86,85,83,84,81,80,79,
        78,77,76,75,74,73,72,70,69,68,
        067,47, 47])
    indtab = -(indtab - 127)

    param_red = param_red[indtab]
    param_green = param_green[indtab]
    param_blue = param_blue[indtab]

    colors = zip(param_red,param_green,param_blue)

    func = cmap_regular_steps if not smoothed else cmap_smoothed_regular_steps
    return func(colors, name=name, **kwargs)

def cmap_previmer(name='vacumm_previmer', **kwargs):
    """Colormap from PREVIMER Website (http://www.previmer.org)

    :Sample: .. image:: misc-color-vacumm_previmer.png

    :Source: http://www.previmer.org (F. Lecornu, G. Charria)
    """
    post = N.array([2,3,4,5,6,7,8,9,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,22,23,24])
    norm = StepsNorm(post)

    r1=N.array([212,191,165,153,97,83,65,48,34,12,8,24,42,56,75,90,105,122,132,134,175,216,250,252,252,252,252,245,231,213,200,185,167,152])
    r1 = r1/255.
    v1=N.array([154,117,95,66,29,46,61,78,94,114,134,150,165,182,199,216,231,244,250,250,250,250,247,231,207,174,141,118,101,87,70,55,40,22])
    v1 = v1/255.
    b1=N.array([229,213,208,209,159,173,119,205,222,242,252,252,252,252,252,252,252,252,228,155,132,132,126,105,82,49,15,4,4,4,4,4,4,4])
    b1 = b1/255.
    colors = []
    for icolor in N.arange(len(r1)):
      colors.append(((r1[icolor],v1[icolor],b1[icolor]),norm.positions[icolor]))
    kwargs.setdefault('ncol', len(colors)-1)
    return cmap_custom(colors, name=name, **kwargs)

def cmap_previmer2(name='vacumm_previmer2', **kwargs):
    """Colormap from PREVIMER Website (http://www.previmer.org)

    Same as :func:`cmap_previmer` but extremes are used for
    :meth:`set_under` and :meth`set_over`.

    :Sample: .. image:: misc-color-vacumm_previmer.png

    :Source: http://www.previmer.org (F. Lecornu, G. Charria)
    """
    r = N.array([212,191,165,153,97,83,65,48,34,12,8,24,42,56,75,90,105,122,132,134,175,216,250,252,252,252,252,245,231,213,200,185,167,152.])
    r /= 255.
    g = N.array([154,117,95,66,29,46,61,78,94,114,134,150,165,182,199,216,231,244,250,250,250,250,247,231,207,174,141,118,101,87,70,55,40,22.])
    g /= 255.
    b = N.array([229,213,208,209,159,173,119,205,222,242,252,252,252,252,252,252,252,252,228,155,132,132,126,105,82,49,15,4,4,4,4,4,4,4.])
    b /= 255.
    colors = zip(r,g,b)
    cmap = cmap_srs(colors[1:-1], name=name, **kwargs)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    return cmap


def cmap_ssec(name='vacumm_ssec', **kwargs):
    """ Colormap from Ncview (http://fossies.org/dox/ncview-2.1.2/colormaps__ssec_8h_source.html)


    :Source: http://fossies.org/dox/ncview-2.1.2/colormaps__ssec_8h_source.html (ncview-2.1.2, G. Charria)
    """


    d = N.array([0,0,45, 0,1,46, 0,2,47, 0,3,48, 0,5,49, 0,6,50, 0,7,51, 0,9,52,
    0,10,53, 0,11,54, 0,13,55, 0,14,56, 0,15,57, 0,17,58, 0,18,59, 0,19,60,
    0,21,62, 0,22,63, 0,23,64, 0,25,65, 0,26,66, 0,27,67, 0,29,68, 0,30,69,
    0,31,70, 0,33,71, 0,34,72, 0,35,73, 0,37,74, 0,38,75, 0,39,76, 0,40,77,
    0,42,79, 0,43,80, 0,44,81, 0,46,82, 0,47,83, 0,48,84, 0,50,85, 0,51,86,
    0,52,87, 0,54,88, 0,55,89, 0,56,90, 0,58,91, 0,59,92, 0,60,93, 0,62,94,
    0,63,96, 0,64,97, 0,66,98, 0,67,99, 0,68,100, 0,70,101, 0,71,102, 0,72,103,
    0,74,104, 0,75,105, 0,76,106, 0,77,107, 0,79,108, 0,80,109, 0,81,110, 0,83,111,
    0,84,113, 0,85,114, 0,87,115, 0,88,116, 0,89,117, 0,91,118, 0,92,119, 0,93,120,
    0,95,121, 0,96,122, 0,97,123, 0,99,124, 0,100,125, 0,101,126, 0,103,127, 0,104,128,
    0,105,130, 0,107,131, 0,108,132, 0,109,133, 0,111,134, 0,112,135, 0,113,136, 0,114,137,
    0,116,138, 0,117,139, 0,118,140, 0,120,141, 0,121,142, 0,122,143, 0,124,144, 0,125,145,
    0,126,147, 0,128,148, 0,129,149, 0,130,150, 0,132,151, 0,133,152, 0,134,153, 0,136,154,
    0,137,155, 0,138,156, 0,140,157, 0,141,158, 0,142,159, 0,144,160, 0,145,161, 0,146,162,
    0,148,164, 0,149,165, 0,150,166, 0,151,167, 0,153,168, 0,154,169, 0,155,170, 0,157,171,
    0,158,172, 0,159,173, 0,161,174, 0,162,175, 0,163,176, 0,165,177, 0,166,178, 0,167,180,
    0,169,181, 0,170,182, 0,171,183, 0,173,184, 0,174,185, 0,175,186, 0,177,187, 0,178,188,
    0,179,189, 0,181,190, 0,182,191, 0,183,192, 0,185,193, 0,186,194, 0,187,195, 0,188,197,
    0,190,198, 0,191,199, 0,192,200, 0,194,201, 0,195,202, 0,196,203, 0,198,204, 0,199,205,
    0,200,206, 0,202,207, 0,203,208, 0,204,209, 0,206,210, 0,207,211, 0,208,212, 0,210,214,
    0,211,215, 0,212,216, 0,214,217, 0,215,218, 0,216,219, 0,218,220, 0,219,221, 0,220,222,
    0,222,223, 0,223,224, 0,224,225, 0,225,226, 0,227,227, 0,228,228, 0,229,229, 8,230,222,
    17,231,214, 26,232,206, 34,233,198, 43,234,190, 52,235,182, 61,236,174, 70,236,166, 78,237,158,
    87,238,150, 96,239,143, 105,240,135, 114,241,127, 122,242,119, 131,242,111, 140,243,103, 149,244,95,
    157,245,87, 166,246,79, 175,247,71, 184,248,63, 193,248,55, 201,249,47, 210,250,39, 219,251,32,
    228,252,24, 237,253,16, 245,254,8, 254,254,0, 255,250,0, 255,245,0, 255,240,0, 255,236,0,
    255,231,0, 255,226,0, 255,221,0, 255,217,0, 255,212,0, 255,207,0, 255,202,0, 255,198,0,
    255,193,0, 255,188,0, 255,183,0, 255,179,0, 255,174,0, 255,169,0, 255,164,0, 255,160,0,
    255,155,0, 255,150,0, 255,145,0, 255,141,0, 255,136,0, 255,131,0, 255,126,0, 255,122,0,
    253,117,0, 249,113,0, 246,109,0, 242,105,0, 239,101,0, 236,97,0, 232,93,0, 229,89,0,
    225,85,0, 222,81,0, 219,77,0, 215,73,0, 212,69,0, 208,65,0, 205,61,0, 202,57,0,
    198,53,0, 195,49,0, 191,45,0, 188,41,0, 185,37,0, 181,33,0, 178,29,0, 175,24,0])
    r = d[0::3]/255.
    g = d[1::3]/255.
    b = d[2::3]/255.

    colors = zip(r,g,b)
    cmap = cmap_srs(colors[1:-1], name=name, **kwargs)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    return cmap

def cmap_ncview_rainbow(name='vacumm_ncview_rainbow', **kwargs):
    """ Colormap Rainbow from Ncview (http://fossies.org/dox/ncview-2.1.2/colormaps__rainbow_8h_source.html)


    :Source: http://fossies.org/dox/ncview-2.1.2/colormaps__rainbow_8h_source.html (ncview-2.1.2, G. Charria)
    """


    d = N.array([38,  0, 50,  38,  0, 51,  38,  0, 51,
         38,  0, 51,  39,  0, 52,  39,  0, 52,  39,  0, 53,  40,  0, 54,  41,  0, 55,
         41,  0, 57,  42,  0, 58,  43,  0, 60,  44,  0, 61,  45,  0, 63,  46,  0, 65,
         47,  0, 67,  48,  0, 70,  50,  0, 72,  51,  0, 75,  52,  0, 77,  53,  0, 80,
         54,  0, 83,  56,  0, 86,  57,  0, 89,  58,  0, 92,  59,  0, 96,  60,  0, 99,
         61,  0,103,  62,  0,106,  63,  0,110,  64,  0,114,  65,  0,117,  65,  0,121,
         66,  0,125,  66,  0,129,  66,  0,133,  66,  0,137,  66,  0,141,  66,  0,145,
         65,  0,150,  65,  0,154,  64,  0,158,  63,  0,162,  62,  0,166,  60,  0,170,
         59,  0,174,  57,  0,179,  55,  0,183,  52,  0,187,  50,  0,191,  47,  0,195,
         44,  0,198,  41,  0,202,  38,  0,206,  34,  0,210,  30,  0,213,  26,  0,217,
         21,  0,220,  17,  0,223,  12,  0,226,   7,  0,229,   2,  0,232,   0,  3,235,
         0,  8,237,   0, 14,240,   0, 20,242,   0, 26,244,   0, 33,246,   0, 39,247,
         0, 46,249,   0, 52,250,   0, 59,251,   0, 66,252,   0, 73,253,   0, 80,253,
         0, 86,253,   0, 93,253,   0,100,253,   0,107,253,   0,114,253,   0,121,253,
         0,128,253,   0,136,252,   0,143,252,   0,150,252,   0,157,252,   0,165,252,
         0,172,252,   0,180,252,   0,187,252,   0,194,251,   0,202,251,   0,210,251,
         0,217,251,   0,225,250,   0,232,250,   0,240,250,   0,247,250,   0,249,244,
         0,249,235,   0,248,227,   0,248,219,   0,247,210,   0,247,202,   0,246,193,
         0,245,185,   0,245,176,   0,244,167,   0,243,159,   0,241,150,   0,240,141,
         0,238,132,   0,237,123,   0,235,114,   0,232,105,   0,229, 96,   0,226, 87,
         0,222, 78,   0,217, 69,   0,212, 61,   0,206, 52,   0,198, 43,   0,190, 35,
         0,181, 27,   0,172, 20,   0,163, 13,   0,157,  7,   0,153,  2,   2,153,  0,
         7,157,  0,  13,163,  0,  20,172,  0,  27,181,  0,  35,190,  0,  43,198,  0,
        52,206,  0,  61,212,  0,  69,217,  0,  78,222,  0,  87,226,  0,  96,229,  0,
       105,232,  0, 114,235,  0, 123,237,  0, 132,238,  0, 141,240,  0, 150,241,  0,
       159,243,  0, 167,244,  0, 176,245,  0, 185,245,  0, 193,246,  0, 202,247,  0,
       210,247,  0, 219,248,  0, 227,248,  0, 235,249,  0, 244,249,  0, 250,247,  0,
       250,240,  0, 250,232,  0, 250,225,  0, 251,217,  0, 251,210,  0, 251,202,  0,
       251,194,  0, 252,187,  0, 252,180,  0, 252,172,  0, 252,165,  0, 252,157,  0,
       252,150,  0, 252,143,  0, 252,136,  0, 253,128,  0, 253,121,  0, 253,114,  0,
       253,107,  0, 253,100,  0, 253, 93,  0, 253, 86,  0, 253, 80,  0, 253, 73,  0,
       252, 66,  0, 251, 59,  0, 250, 52,  0, 249, 46,  0, 247, 39,  0, 246, 33,  0,
       244, 26,  0, 242, 20,  0, 240, 14,  0, 237,  8,  0, 235,  3,  0, 232,  0,  2,
       229,  0,  7, 226,  0, 12, 223,  0, 17, 220,  0, 21, 217,  0, 26, 213,  0, 30,
       210,  0, 34, 206,  0, 38, 202,  0, 41, 198,  0, 44, 195,  0, 47, 191,  0, 50,
       187,  0, 52, 183,  0, 55, 179,  0, 57, 174,  0, 59, 170,  0, 60, 166,  0, 62,
       162,  0, 63, 158,  0, 64, 154,  0, 65, 150,  0, 65, 145,  0, 66, 141,  0, 66,
       137,  0, 66, 133,  0, 66, 129,  0, 66, 125,  0, 66, 121,  0, 65, 117,  0, 65,
       114,  0, 64, 110,  0, 63, 106,  0, 62, 103,  0, 61,  99,  0, 60,  96,  0, 59,
        92,  0, 58,  89,  0, 57,  86,  0, 56,  83,  0, 54,  80,  0, 53,  77,  0, 52,
        75,  0, 51,  72,  0, 50,  70,  0, 48,  67,  0, 47,  65,  0, 46,  63,  0, 45,
        61,  0, 44,  60,  0, 43,  58,  0, 42,  57,  0, 41,  55,  0, 41,  54,  0, 40,
        53,  0, 39,  52,  0, 39,  52,  0, 39,  51,  0, 38,  51,  0, 38,  51,  0, 38,50,  0, 38])
    r = d[0::3]/255.
    g = d[1::3]/255.
    b = d[2::3]/255.

    colors = zip(r,g,b)
    cmap = cmap_srs(colors[1:-1], name=name, **kwargs)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    return cmap

def cmap_nice_gfdl(name='vacumm_nice_gfdl', **kwargs):
    """ GFDL colormap (http://www.gfdl.noaa.gov/visualization)


    :Source: http://www.ncl.ucar.edu/Document/Graphics/ColorTables/nice_gfdl.shtml (included by G. Charria)
    """


    d = N.array([0.996078, 0.984314, 0.964706,
    0.925490, 0.929412, 0.945098,
    0.905882, 0.909804, 0.925490,
    0.862745, 0.882353, 0.901961,
    0.835294, 0.854902, 0.874510,
    0.811765, 0.823529, 0.858824,
    0.784314, 0.796078, 0.831373,
    0.749020, 0.772549, 0.811765,
    0.729412, 0.749020, 0.788235,
    0.694118, 0.717647, 0.768627,
    0.670588, 0.690196, 0.741176,
    0.639216, 0.666667, 0.725490,
    0.611765, 0.639216, 0.698039,
    0.580392, 0.607843, 0.666667,
    0.560784, 0.588235, 0.647059,
    0.517647, 0.560784, 0.623529,
    0.490196, 0.537255, 0.596078,
    0.462745, 0.517647, 0.576471,
    0.435294, 0.490196, 0.545098,
    0.400000, 0.447059, 0.525490,
    0.384314, 0.431373, 0.509804,
    0.352941, 0.407843, 0.486275,
    0.325490, 0.380392, 0.458824,
    0.294118, 0.356863, 0.443137,
    0.270588, 0.329412, 0.415686,
    0.247059, 0.301961, 0.396078,
    0.223529, 0.282353, 0.372549,
    0.196078, 0.254902, 0.360784,
    0.168627, 0.223529, 0.325490,
    0.133333, 0.203922, 0.301961,
    0.113725, 0.180392, 0.274510,
    0.094118, 0.149020, 0.250980,
    0.074510, 0.125490, 0.227451,
    0.050980, 0.109804, 0.203922,
    0.047059, 0.105882, 0.196078,
    0.050980, 0.117647, 0.203922,
    0.062745, 0.129412, 0.219608,
    0.074510, 0.141176, 0.235294,
    0.086275, 0.156863, 0.254902,
    0.094118, 0.176471, 0.258824,
    0.105882, 0.188235, 0.274510,
    0.121569, 0.207843, 0.298039,
    0.133333, 0.219608, 0.309804,
    0.137255, 0.243137, 0.325490,
    0.145098, 0.254902, 0.337255,
    0.160784, 0.270588, 0.356863,
    0.176471, 0.286275, 0.372549,
    0.180392, 0.301961, 0.380392,
    0.196078, 0.313725, 0.396078,
    0.203922, 0.325490, 0.407843,
    0.219608, 0.341176, 0.423529,
    0.223529, 0.360784, 0.427451,
    0.247059, 0.384314, 0.450980,
    0.247059, 0.396078, 0.458824,
    0.262745, 0.415686, 0.478431,
    0.282353, 0.439216, 0.490196,
    0.290196, 0.447059, 0.498039,
    0.298039, 0.462745, 0.513725,
    0.309804, 0.478431, 0.529412,
    0.313725, 0.501961, 0.533333,
    0.329412, 0.517647, 0.549020,
    0.333333, 0.529412, 0.560784,
    0.349020, 0.549020, 0.580392,
    0.356863, 0.564706, 0.592157,
    0.372549, 0.580392, 0.607843,
    0.392157, 0.603922, 0.631373,
    0.403922, 0.615686, 0.643137,
    0.403922, 0.631373, 0.643137,
    0.423529, 0.654902, 0.666667,
    0.431373, 0.662745, 0.674510,
    0.447059, 0.678431, 0.694118,
    0.454902, 0.698039, 0.705882,
    0.474510, 0.717647, 0.725490,
    0.482353, 0.725490, 0.733333,
    0.501961, 0.749020, 0.756863,
    0.505882, 0.772549, 0.752941,
    0.517647, 0.788235, 0.764706,
    0.525490, 0.807843, 0.784314,
    0.541176, 0.819608, 0.800000,
    0.549020, 0.839216, 0.811765,
    0.564706, 0.858824, 0.831373,
    0.580392, 0.874510, 0.847059,
    0.596078, 0.894118, 0.862745,
    0.596078, 0.905882, 0.862745,
    0.596078, 0.905882, 0.862745,
    0.576471, 0.890196, 0.819608,
    0.564706, 0.878431, 0.811765,
    0.549020, 0.866667, 0.760784,
    0.541176, 0.858824, 0.752941,
    0.529412, 0.847059, 0.729412,
    0.517647, 0.835294, 0.713725,
    0.498039, 0.827451, 0.662745,
    0.478431, 0.807843, 0.643137,
    0.470588, 0.803922, 0.607843,
    0.454902, 0.784314, 0.588235,
    0.443137, 0.776471, 0.556863,
    0.431373, 0.764706, 0.545098,
    0.415686, 0.749020, 0.501961,
    0.407843, 0.741176, 0.494118,
    0.392157, 0.729412, 0.458824,
    0.380392, 0.713725, 0.447059,
    0.368627, 0.701961, 0.415686,
    0.352941, 0.682353, 0.400000,
    0.345098, 0.678431, 0.360784,
    0.329412, 0.662745, 0.345098,
    0.317647, 0.647059, 0.325490,
    0.305882, 0.635294, 0.313725,
    0.282353, 0.623529, 0.270588,
    0.274510, 0.615686, 0.262745,
    0.262745, 0.592157, 0.223529,
    0.258824, 0.584314, 0.215686,
    0.247059, 0.576471, 0.180392,
    0.243137, 0.572549, 0.176471,
    0.270588, 0.584314, 0.149020,
    0.282353, 0.600000, 0.160784,
    0.313725, 0.619608, 0.117647,
    0.329412, 0.639216, 0.129412,
    0.372549, 0.654902, 0.098039,
    0.384314, 0.666667, 0.109804,
    0.419608, 0.686275, 0.070588,
    0.435294, 0.701961, 0.086275,
    0.478431, 0.721569, 0.023529,
    0.494118, 0.741176, 0.050980,
    0.529412, 0.756863, 0.000000,
    0.545098, 0.772549, 0.000000,
    0.588235, 0.788235, 0.000000,
    0.603922, 0.807843, 0.000000,
    0.635294, 0.811765, 0.000000,
    0.658824, 0.835294, 0.000000,
    0.698039, 0.850980, 0.000000,
    0.721569, 0.874510, 0.000000,
    0.756863, 0.878431, 0.000000,
    0.780392, 0.905882, 0.000000,
    0.823529, 0.909804, 0.000000,
    0.847059, 0.933333, 0.000000,
    0.878431, 0.945098, 0.000000,
    0.901961, 0.968627, 0.000000,
    0.933333, 0.972549, 0.000000,
    0.960784, 1.000000, 0.000000,
    1.000000, 1.000000, 0.000000,
    1.000000, 1.000000, 0.000000,
    1.000000, 0.984314, 0.000000,
    1.000000, 0.972549, 0.000000,
    1.000000, 0.921569, 0.000000,
    1.000000, 0.905882, 0.000000,
    1.000000, 0.862745, 0.000000,
    1.000000, 0.847059, 0.000000,
    1.000000, 0.803922, 0.000000,
    1.000000, 0.788235, 0.000000,
    1.000000, 0.749020, 0.000000,
    1.000000, 0.733333, 0.000000,
    1.000000, 0.694118, 0.000000,
    1.000000, 0.678431, 0.000000,
    1.000000, 0.631373, 0.000000,
    1.000000, 0.619608, 0.000000,
    1.000000, 0.580392, 0.000000,
    1.000000, 0.568627, 0.000000,
    1.000000, 0.529412, 0.000000,
    1.000000, 0.509804, 0.000000,
    1.000000, 0.466667, 0.000000,
    1.000000, 0.458824, 0.000000,
    1.000000, 0.431373, 0.000000,
    1.000000, 0.407843, 0.000000,
    1.000000, 0.376471, 0.000000,
    0.980392, 0.360784, 0.000000,
    0.952941, 0.333333, 0.000000,
    0.929412, 0.313725, 0.000000,
    0.909804, 0.290196, 0.000000,
    0.886275, 0.270588, 0.000000,
    0.862745, 0.243137, 0.000000,
    0.843137, 0.231373, 0.000000,
    0.819608, 0.203922, 0.000000,
    0.792157, 0.184314, 0.000000,
    0.772549, 0.160784, 0.000000,
    0.749020, 0.145098, 0.000000,
    0.725490, 0.121569, 0.023529,
    0.721569, 0.117647, 0.019608,
    0.686275, 0.125490, 0.023529,
    0.674510, 0.117647, 0.011765,
    0.631373, 0.117647, 0.035294,
    0.627451, 0.117647, 0.031373,
    0.603922, 0.109804, 0.031373,
    0.592157, 0.101961, 0.023529,
    0.549020, 0.105882, 0.035294,
    0.545098, 0.101961, 0.031373,
    0.505882, 0.101961, 0.027451,
    0.501961, 0.098039, 0.023529,
    0.474510, 0.101961, 0.035294,
    0.466667, 0.098039, 0.031373,
    0.431373, 0.094118, 0.039216,
    0.427451, 0.090196, 0.035294,
    0.392157, 0.094118, 0.039216,
    0.388235, 0.090196, 0.035294,
    0.360784, 0.086275, 0.039216,
    0.349020, 0.078431, 0.031373,
    0.313725, 0.086275, 0.047059,
    0.301961, 0.078431, 0.043137,
    0.290196, 0.078431, 0.043137,
    0.278431, 0.070588, 0.039216,
    0.239216, 0.074510, 0.039216,
    0.235294, 0.070588, 0.039216,
    0.215686, 0.066667, 0.043137,
    0.207843, 0.062745, 0.039216,
    0.180392, 0.062745, 0.043137,
    0.160784, 0.050980, 0.031373,
    0.141176, 0.054902, 0.035294,
    0.137255, 0.050980, 0.031373,
    0.113725, 0.050980, 0.035294,
    0.101961, 0.043137, 0.023529,
    0.082353, 0.043137, 0.031373,
    0.070588, 0.031373, 0.019608,
    0.058824, 0.031373, 0.023529,
    0.058824, 0.031373, 0.023529,
    0.054902, 0.031373, 0.019608,
    0.050980, 0.031373, 0.015686,
    0.047059, 0.023529, 0.019608,
    0.050980, 0.027451, 0.023529,
    0.043137, 0.027451, 0.019608,
    0.039216, 0.015686, 0.000000,
    0.035294, 0.019608, 0.015686,
    0.031373, 0.011765, 0.000000,
    0.023529, 0.015686, 0.000000,
    0.023529, 0.015686, 0.000000,
    0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000])
    r = d[0::3]
    g = d[1::3]
    b = d[2::3]

    colors = zip(r,g,b)
    cmap = cmap_srs(colors[1:-1], name=name, **kwargs)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    return cmap


def cmap_eke(name='vacumm_eke', **kwargs):
    """ Colormap for Eddy Kinetic Energy (from Barnier et al., 2006)

    Added by G. Charria
    """


    r = N.array([1,
    0.98792,
    0.97584,
    0.96376,
    0.95169,
    0.93961,
    0.92753,
    0.91545,
    0.90337,
    0.89129,
    0.87922,
    0.86714,
    0.85506,
    0.84298,
    0.8309,
    0.81882,
    0.80675,
    0.79467,
    0.78259,
    0.77051,
    0.75843,
    0.74635,
    0.73427,
    0.7222,
    0.71012,
    0.69804,
    0.68596,
    0.67388,
    0.6618,
    0.64973,
    0.63765,
    0.62557,
    0.61349,
    0.60141,
    0.58933,
    0.57725,
    0.56518,
    0.55122,
    0.53584,
    0.52047,
    0.5051,
    0.48973,
    0.47435,
    0.45898,
    0.44361,
    0.42824,
    0.41286,
    0.39749,
    0.38212,
    0.36675,
    0.35137,
    0.336,
    0.32063,
    0.30525,
    0.28988,
    0.27451,
    0.25914,
    0.24376,
    0.22839,
    0.21302,
    0.19765,
    0.18227,
    0.1669,
    0.15153,
    0.13616,
    0.12078,
    0.10541,
    0.090039,
    0.074667,
    0.059294,
    0.043922,
    0.028549,
    0.013176,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.0035294,
    0.0084706,
    0.013412,
    0.018353,
    0.023294,
    0.028235,
    0.033176,
    0.038118,
    0.043059,
    0.048,
    0.052941,
    0.057882,
    0.062824,
    0.067765,
    0.072706,
    0.077647,
    0.082588,
    0.087529,
    0.092471,
    0.097412,
    0.10235,
    0.10729,
    0.11224,
    0.11718,
    0.12212,
    0.12706,
    0.132,
    0.13694,
    0.14188,
    0.14682,
    0.15176,
    0.15671,
    0.16165,
    0.16659,
    0.17153,
    0.17647,
    0.18643,
    0.20894,
    0.23145,
    0.25396,
    0.27647,
    0.29898,
    0.32149,
    0.344,
    0.36651,
    0.38902,
    0.41153,
    0.43404,
    0.45655,
    0.47906,
    0.50157,
    0.52408,
    0.54659,
    0.5691,
    0.59161,
    0.61412,
    0.63663,
    0.65914,
    0.68165,
    0.70416,
    0.72667,
    0.74918,
    0.77169,
    0.7942,
    0.81671,
    0.83922,
    0.86173,
    0.88424,
    0.90675,
    0.92925,
    0.95176,
    0.97427,
    0.99678,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    0.98824,
    0.96078,
    0.93333,
    0.90588,
    0.87843,
    0.85098,
    0.82353,
    0.79608,
    0.76863,
    0.74118,
    0.71373,
    0.68627,
    0.65882,
    0.63137,
    0.60392,
    0.57647,
    0.54902,
    0.52157,
    0.49412,
    0.46667,
    0.43922,
    0.41176,
    0.38431,
    0.35686,
    0.32941,
    0.30196,
    0.27451,
    0.24706,
    0.21961,
    0.19216,
    0.16471,
    0.13725,
    0.1098,
    0.082353,
    0.054902,
    0.027451,
    0])
    g = N.array([1,
    0.97749,
    0.95498,
    0.93247,
    0.90996,
    0.88745,
    0.86494,
    0.84243,
    0.81992,
    0.79741,
    0.7749,
    0.75239,
    0.72988,
    0.70737,
    0.68486,
    0.66235,
    0.63984,
    0.61733,
    0.59482,
    0.57231,
    0.5498,
    0.52729,
    0.50478,
    0.48227,
    0.45976,
    0.43725,
    0.41475,
    0.39224,
    0.36973,
    0.34722,
    0.32471,
    0.3022,
    0.27969,
    0.25718,
    0.23467,
    0.21216,
    0.18965,
    0.19286,
    0.21537,
    0.23788,
    0.26039,
    0.2829,
    0.30541,
    0.32792,
    0.35043,
    0.37294,
    0.39545,
    0.41796,
    0.44047,
    0.46298,
    0.48549,
    0.508,
    0.53051,
    0.55302,
    0.57553,
    0.59804,
    0.62055,
    0.64306,
    0.66557,
    0.68808,
    0.71059,
    0.7331,
    0.75561,
    0.77812,
    0.80063,
    0.82314,
    0.84565,
    0.86816,
    0.89067,
    0.91318,
    0.93569,
    0.9582,
    0.98071,
    0.99608,
    0.96863,
    0.94118,
    0.91373,
    0.88627,
    0.85882,
    0.83137,
    0.80392,
    0.77647,
    0.74902,
    0.72157,
    0.69412,
    0.66667,
    0.63922,
    0.61176,
    0.58431,
    0.55686,
    0.52941,
    0.50196,
    0.47451,
    0.44706,
    0.41961,
    0.39216,
    0.36471,
    0.33725,
    0.3098,
    0.28235,
    0.2549,
    0.22745,
    0.2,
    0.17255,
    0.1451,
    0.11765,
    0.090196,
    0.062745,
    0.035294,
    0.0078431,
    0.01549,
    0.037176,
    0.058863,
    0.080549,
    0.10224,
    0.12392,
    0.14561,
    0.16729,
    0.18898,
    0.21067,
    0.23235,
    0.25404,
    0.27573,
    0.29741,
    0.3191,
    0.34078,
    0.36247,
    0.38416,
    0.40584,
    0.42753,
    0.44922,
    0.4709,
    0.49259,
    0.51427,
    0.53596,
    0.55765,
    0.57933,
    0.60102,
    0.62271,
    0.64439,
    0.66608,
    0.68776,
    0.70945,
    0.73114,
    0.75282,
    0.77451,
    0.79165,
    0.79741,
    0.80318,
    0.80894,
    0.81471,
    0.82047,
    0.82624,
    0.832,
    0.83776,
    0.84353,
    0.84929,
    0.85506,
    0.86082,
    0.86659,
    0.87235,
    0.87812,
    0.88388,
    0.88965,
    0.89541,
    0.90118,
    0.90694,
    0.91271,
    0.91847,
    0.92424,
    0.93,
    0.93576,
    0.94153,
    0.94729,
    0.95306,
    0.95882,
    0.96459,
    0.97035,
    0.97612,
    0.98188,
    0.98765,
    0.99341,
    0.99918,
    0.97647,
    0.94902,
    0.92157,
    0.89412,
    0.86667,
    0.83922,
    0.81176,
    0.78431,
    0.75686,
    0.72941,
    0.70196,
    0.67451,
    0.64706,
    0.61961,
    0.59216,
    0.56471,
    0.53725,
    0.5098,
    0.48235,
    0.4549,
    0.42745,
    0.4,
    0.37255,
    0.3451,
    0.31765,
    0.2902,
    0.26275,
    0.23529,
    0.20784,
    0.18039,
    0.15294,
    0.12549,
    0.098039,
    0.070588,
    0.043137,
    0.015686,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0])
    b = N.array([1,
    0.99424,
    0.98847,
    0.98271,
    0.97694,
    0.97118,
    0.96541,
    0.95965,
    0.95388,
    0.94812,
    0.94235,
    0.93659,
    0.93082,
    0.92506,
    0.91929,
    0.91353,
    0.90776,
    0.902,
    0.89624,
    0.89047,
    0.88471,
    0.87894,
    0.87318,
    0.86741,
    0.86165,
    0.85588,
    0.85012,
    0.84435,
    0.83859,
    0.83282,
    0.82706,
    0.82129,
    0.81553,
    0.80976,
    0.804,
    0.79824,
    0.79247,
    0.79329,
    0.79906,
    0.80482,
    0.81059,
    0.81635,
    0.82212,
    0.82788,
    0.83365,
    0.83941,
    0.84518,
    0.85094,
    0.85671,
    0.86247,
    0.86824,
    0.874,
    0.87976,
    0.88553,
    0.89129,
    0.89706,
    0.90282,
    0.90859,
    0.91435,
    0.92012,
    0.92588,
    0.93165,
    0.93741,
    0.94318,
    0.94894,
    0.95471,
    0.96047,
    0.96624,
    0.972,
    0.97776,
    0.98353,
    0.98929,
    0.99506,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    0.98039,
    0.95294,
    0.92549,
    0.89804,
    0.87059,
    0.84314,
    0.81569,
    0.78824,
    0.76078,
    0.73333,
    0.70588,
    0.67843,
    0.65098,
    0.62353,
    0.59608,
    0.56863,
    0.54118,
    0.51373,
    0.48627,
    0.45882,
    0.43137,
    0.40392,
    0.37647,
    0.34902,
    0.32157,
    0.29412,
    0.26667,
    0.23922,
    0.21176,
    0.18431,
    0.15686,
    0.12941,
    0.10196,
    0.07451,
    0.047059,
    0.019608,
    0.0021961,
    0.0098824,
    0.017569,
    0.025255,
    0.032941,
    0.040627,
    0.048314,
    0.056,
    0.063686,
    0.071373,
    0.079059,
    0.086745,
    0.094431,
    0.10212,
    0.1098,
    0.11749,
    0.12518,
    0.13286,
    0.14055,
    0.14824,
    0.15592,
    0.16361,
    0.17129,
    0.17898,
    0.18667,
    0.19435,
    0.20204,
    0.20973,
    0.21741,
    0.2251,
    0.23278,
    0.24047,
    0.24816,
    0.25584,
    0.26353,
    0.27122,
    0.2789,
    0.27341,
    0.26573,
    0.25804,
    0.25035,
    0.24267,
    0.23498,
    0.22729,
    0.21961,
    0.21192,
    0.20424,
    0.19655,
    0.18886,
    0.18118,
    0.17349,
    0.1658,
    0.15812,
    0.15043,
    0.14275,
    0.13506,
    0.12737,
    0.11969,
    0.112,
    0.10431,
    0.096627,
    0.088941,
    0.081255,
    0.073569,
    0.065882,
    0.058196,
    0.05051,
    0.042824,
    0.035137,
    0.027451,
    0.019765,
    0.012078,
    0.0043922,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0])


    colors = zip(r,g,b)
    cmap = cmap_srs(colors[1:-1], name=name, **kwargs)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    return cmap

def _local_cmap_cpt_(name, **kwargs):
    """Get a cmap stored in a GMT format in file _cpt_dir/<name>.cpt"""
    if name in cmap_d:
        return cmap_d[name]
    sname = name[7:] if name.startswith('vacumm_') else name
    return cmap_gmt(os.path.join(_cpt_dir, sname+'.cpt'), register='vacumm_'+sname, **kwargs)


def cmap_currents(name='vacumm_currents', **kwargs):
    """A colormap for displaying currents"""
    return _local_cmap_cpt_(name, **kwargs)

def cmap_rnb2_hymex(name="vacumm_rnb2_hymex", **kwargs):
    """RNB2 Colormap for HYMEX

    .. note:: This colormap is registered in matplotlib under the name "vacumm_rnb2_hymex".

    :Sample: .. image:: misc-color-vacumm_rnb2_hymex.png
    """
    data = [
        ((0., 0, .6), 0.),
        ((0, 0.2, 1), .12),
        ((.2, .7, 1), .16),
        ((0, .8, 0), .27),
        ((1, 1, 0), .5),
        ((1, .1, .1), .82),
        ((.8, .1, .1), .9),
        ((.6, .25, .6), 1.)]
    cmap = cmap_smoothed_steps(data, name=name, **kwargs)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    return cmap

def cmap_rainbow_sst_hymex(name="vacumm_rainbow_sst_hymex", **kwargs):
    """RAINBOW_SST Colormap for HYMEX

    .. note:: This colormap is registered in matplotlib under the name "vacumm_rainbow_sst_hymex".

    :Sample: .. image:: misc-color-vacumm_rainbow_sst_hymex.png
    """
    data = [
        ((0., 0, .52), 0.),
        ((0, 0.41, 1), .2),
        ((.09, .98, 0.91), .39),
        ((0.61, 1., 0.38), .53),
        ((1, 0.95, 0), .66),
        ((1, .2, 0), .85),
        ((.61, 0, 0), 1.)]
    cmap = cmap_smoothed_steps(data, name=name, **kwargs)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    return cmap

def cmap_dynamic_cmyk_hymex(name="vacumm_dynamic_cmyk_hymex", **kwargs):
    """DYNAMIC_CMYK Colormap for HYMEX

    .. note:: This colormap is registered in matplotlib under the name "vacumm_dynamic_cmyk_hymex".

    :Sample: .. image:: misc-color-vacumm_dynamic_cmyk_hymex.png
    """
    data = [
        ((0.94902, 0.94510, 0.94902), 0.),
        ((0.94510, 0.88627, 0.94118), .1),
        ((0.69804, 0.62745, 0.80392), .2),
        ((0.50196, 0.65098, 0.80784), .3),
        ((0.40392, 0.75686, 0.60784), .4),
        ((0.33725, 0.69412, 0.20000), .5),
        ((0.80000, 0.90588, 0.04314), .6),
        ((0.98039, 0.96078, 0.09804), .7),
        ((0.97255, 0.50980, 0.06667), .8),
        ((0.98039, 0.07059, 0.04706), .9),
        ((0.59608, 0.05882, 0.03137), 1.)]
    cmap = cmap_smoothed_steps(data, name=name, **kwargs)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    return cmap

def cmap_white_centered_hymex(name="vacumm_white_centered_hymex", **kwargs):
    """RNB2 Colormap for HYMEX

    .. note:: This colormap is registered in matplotlib under the name "vacumm_white_centered_hymex".

    :Sample: .. image:: misc-color-vacumm_white_centered_hymex.png
    """
    data = [
        ((0, 1, 1), 0.),
        ((0.2, 0.2, 1), .2),
        ((1, 1, 1), .35),
        ((1, 1, 1), .65),
        ((1,0.2,0.2), .8),
        ((1, 1, 0), 1.)]
    cmap = cmap_smoothed_steps(data, name=name, **kwargs)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    return cmap

def cmap_red_tau_hymex(name="vacumm_red_tau_hymex", **kwargs):
    """RNB2 Colormap for HYMEX

    .. note:: This colormap is registered in matplotlib under the name "vacumm_red_tau_hymex".

    :Sample: .. image:: misc-color-vacumm_red_tau_hymex.png
    """
    data = [
        ((1.0, 1.0, 1.0), 0.),
        ((1.0, 1.0,   0), .2),
        ((1.0, 0.8,0.05), .4),
        ((1.0, 0.6, 0.1), .6),
        ((1.0, 0.3, 0.1), .9),
        ((1.0, 0.0,   0), 1.)]
    cmap = cmap_smoothed_steps(data, name=name, **kwargs)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    return cmap


def get_cmap(cmap=None, errmode=None, **kwargs):
    """A simple way to get a VACUMM, GMT or MPL colormap

    :Example:

        >>> get_cmap('jet')         # Matplotlib
        >>> pylab.jet()             # Matplotlib
        >>> pylab.get_cmap('jet')   # Matplotlib
        >>> get_cmap('cmap_grey',start=.2)   # VACUMM
        >>> get_cmap('vacumm_grey')   # VACUMM (registered in MPL)
        >>> cmap_grey(start=.2)     # VACUMM
        >>> get_cmap('gmt_gebco')   # GMT (registered in MPL)
        >>> cmap_gmt('gebco')       # GMT (registered in MPL)

    :See also:

        :func:`cmap_gmt` :func:`matplotlib.pyplot.get_cmap`
    """
    if isinstance(cmap, Colormap):
       return cmap
    if isinstance(cmap, basestring):
        if not kwargs: # Try an already registered colormap
            try:
                return P.get_cmap(cmap)
            except:
                pass
        if cmap.startswith('cmap_'):
            cmap = eval(cmap)(**kwargs)
        elif cmap.startswith('gmt_') or cmap.endswith('.cpt'):
            cmap = cmap_gmt(cmap)
        elif  cmap.startswith('vacumm_') and kwargs:
            cmap = eval('cmap_'+cmap[7:])(**kwargs)
        else:
            if errmode=='raise':
                cmap = P.get_cmap(cmap)
            else:
                try:
                    cmap = P.get_cmap(cmap)
                except:
                    if errmode=='warn':
                        vacumm_warn('Invalid colormap: {!s}. Switch to default one'.format(cmap))
                    cmap = P.get_cmap()

    else:
        cmap = P.get_cmap()
    return cmap


def plot_cmap(cmap, ncol=None, smoothed=True,  ax=None, figsize=(5, .25), fig=None,
        show=True, aspect=.05, title=None, title_loc=(.5, .5),
        sa=dict(left=.0, right=1, top=1, bottom=.0),
        savefig=None, savefigs=None, close=True, **kwargs):
    """Display a colormap"""
    from vacumm.misc import kwfilter, dict_check_defaults
    cmap = get_cmap(cmap)
    if ax is None:
        if fig is None:
            fig = P.figure()
        ax = fig.gca()
    else:
        fig = ax.get_figure()
    if figsize:
        fig.set_size_inches(figsize, forward=True)
    ax.set_aspect('auto', anchor='C')
    ax.axis("off")
    if ncol is None:
        ncol = cmap.N if hasattr(cmap, 'N') else 256.
    x = N.arange(0,ncol,1.)
    dx = x[1]-x[0]
    interp = 'bilinear' if smoothed else "nearest"
    p = ax.imshow(N.outer(N.ones(1),x), aspect=aspect*ncol,
        cmap=cmap, origin="lower",
        interpolation=interp, vmin=-.5, vmax=x[-1]+0.5)
    ax.fill([x[0]-dx/2, x[0]-dx/2, x[0]-dx/2-ncol/20.], [-.5, .5, 0],
        facecolor=p.to_rgba(-1), edgecolor='none', linewidth=0)
    ax.fill([x[-1]+dx/2, x[-1]+dx/2, x[-1]+dx/2+ncol/20.], [-.5, .5, 0],
        facecolor=p.to_rgba(ncol+2), edgecolor='none', linewidth=0)
    if title is None: title = cmap.name
    if title is not None and title is not False:
        kwtitle = kwfilter(kwargs, 'title_')
        dict_check_defaults(kwtitle, alpha=.8, ha='center', va='center',
            transform=ax.transAxes, size=9, color='k')
        ax.text(title_loc[0], title_loc[1], title, **kwtitle)
    if sa:
        fig.subplots_adjust(**sa)

    # Save and show
    if savefig is not None:
        fig.savefig(savefig, **kwfilter(kwargs, 'savefig'))
    if savefigs is not None:
        from plot import savefigs as Savefigs
        Savefigs(savefigs, fig=fig, **kwfilter(kwargs, 'savefigs'))
    if show:
        #fig.show()
        P.show()
    if close:
        P.close(fig)


def plot_cmaps(cmaps=None, figsize=None, show=True, savefig=None, ncol=5,
    savefigs=None, aspect=0.05, fig=None, close=True, **kwargs):
    """Display a list of or all colormaps"""

    from vacumm.misc import kwfilter
    kwsf = kwfilter(kwargs, 'savefig')
    kwsfs = kwfilter(kwargs, 'savefigs')
    kwargs.pop('nrow', None)

    # Default colormap list
    from vacumm.misc import kwfilter
    if cmaps is None:
        cmaps = cmaps_mpl(names=False)
    elif isinstance(cmaps, str):
        if cmaps.startswith('vacumm'):
            cmaps = cmaps_vacumm(names=False)
        elif cmaps == 'mpl':
            cmaps = cmaps_mpl(names=False, vacumm=False)
        elif cmaps == 'gmt':
            cmaps = cmaps_gmt()
        else:
            cmaps = [cmaps]

    # Check validity
    goodcmaps = []
    for cmap in cmaps:
        try:
            cmap = get_cmap(cmap)
            if cmap is not None:
                goodcmaps.append(cmap)
        except:
            continue
    cmaps = goodcmaps

    # Uniq
    oldcmaps = cmaps
    cmaps = []
    names = []
    for cmap in oldcmaps:
        if cmap.name not in names:
            cmaps.append(cmap)
            names.append(cmap.name)
    ncmap = len(cmaps)
    assert ncmap, 'No valid cmap'

    # Sort
    cmaps.sort(cmp=lambda a, b: cmp(a.name, b.name))

    # Setup figure
    if fig is None:
        fig = P.figure()
    ncol = min(ncmap, ncol)
    nrow = (ncmap-1)/ncol+1
#    nrow = min(ncmap, nrow)
#    ncol = (ncmap-1)/nrow+1
    if N.isscalar(figsize):
        figwidth = figsize
        figsize = None
    else:
        figwidth = 6.
    if figsize is None:
        onecol = figwidth / ncol
        onerow = figwidth * aspect * .8
        figheight = onerow * nrow
        figsize = (figwidth, figheight)
    if figsize is not False:
        fig.set_size_inches(figsize, forward=True)

    # Loop on colormaps
    ip = N.arange(ncol*nrow).reshape((nrow,ncol)).T.ravel()+1
    for i, cmap in enumerate(cmaps):
        ax = fig.add_subplot(nrow, ncol, ip[i])
        plot_cmap(cmap, show=False, ax=ax, figsize=False,
            aspect=aspect, close=False, title_loc=(.5, 2), title_alpha=1,
            **kwargs)
        ax.set_anchor('C')
#    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    fig.subplots_adjust(wspace=0.0, hspace=0, top=0.98, bottom=0.02, left=0.02, right=0.98)

    # Save and show
    if savefig is not None:
        P.savefig(savefig, fig=fig, **kwsf)
    if savefigs is not None:
        from plot import savefigs as _savefigs
        _savefigs(savefigs, fig=fig, **kwsfs)
    if show: fig.show()
    if close: P.close(fig)

def show_cmap(cmap, *args, **kwargs):
    """Alias for :func:`plot_cmap`"""
    return plot_cmap(cmap, *args, **kwargs)

def cmaps_registered(include=None, exclude=None, names=True):
    """List colormap registered in matplotlib

    :Params:

        - **include**: Include only colormaps that have one of these prefixes.
        - **include**: Exclude colormaps that have one of these prefixes.
        - **names**: Return names OR colormaps.
    """
    cmap_names = cmap_d.keys()
    cmap_names.sort()
    if include is None: include = []
    elif isinstance(include, basestring): include = [include]
    if exclude is None: exclude = []
    elif isinstance(exclude, basestring): exclude = [exclude]
    for inc in include:
        cmap_names = [name for name in cmap_names if name.startswith(inc)]
    for exc in exclude:
        cmap_names = [name for name in cmap_names if not name.startswith(exc)]
    if not names:
        return [cmap_d[name] for name in cmap_names]
    return cmap_names

def cmaps_mpl(vacumm=True, gmt=True, names=True):
    """List available registered colormaps available directly from Matplotlib

    :See also:

        :func:`get_cmap`
    """
    exclude = []
    if not vacumm: exclude.append('vacumm_')
    if not gmt: exclude.append('GMT_')
    return cmaps_registered(exclude=exclude, names=names)

def cmaps_vacumm(names=True):
    """List available VACUMM colormaps

    :See also:

        :func:`get_cmap`
    """
    return cmaps_registered(include='vacumm_', names=names)
cmaps_act = cmaps_vacumm

def cmaps_gmt(names=True):
    """List available GMT colormaps

    :See also:

        :func:`cmap_gmt` :func:`print_cmaps_gmt` :func:`get_cmap`
    """
    return cmaps_registered(include='GMT_', names=names)

def print_cmaps_gmt():
    """List available gmt colormaps

    :See also:

        :func:`cmap_gmt` :func:`cmaps_gmt` :func:`get_cmap`
    """
    print 'List of available GMT colormaps: '+', '.join(cmaps_gmt(names=True))

_re_split_gmt = re.compile(r'[\t/\-]+').split
def cmap_gmt(name, register=True, **kwargs):
    """Get a colormap from GMT

    :Params:

        - **name**: GMT colormap name OR .cpt file name.

    :See also:

        :func:`cmaps_gmt` :func:`print_cmaps_gmt` :func:`get_cmap`
    """
    # Already registered
    if not name.endswith('.cpt'):

#        name = name.lower()
        if name.startswith('gmt_'): name = 'GMT_'+name[4:]
        mname = ('GMT_'+name) if not name.startswith('GMT_') else name
        if mname not in cmaps_gmt():
            raise vacumm.VACUMMError('Invalid GMT colormap: %s. '%name +
                'Please print available name with: print_cmaps_gmt()')
        return P.get_cmap(mname)

    # From file
    import colorsys
    filePath = name
    name = os.path.basename(name[:-4])#.lower()

    if not os.path.exists(filePath):
        raise vacumm.VACUMMError("Wrong GMT colormap file: "+filePath)

    f = open(filePath)
    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = _re_split_gmt(l.strip())
#        ls = l.split()
        if l[0] == "#":
            if ls[-1] == "HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if len(ls)==0: continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        elif len(ls)==8:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])
#        else:
#            return


    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    nTable = len(r)
    x = N.array( x , 'f')
    r = N.array( r , 'f')
    g = N.array( g , 'f')
    b = N.array( b , 'f')
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
        r = r/255.
        g = g/255.
        b = b/255.
    xNorm = (x - x[0])/(x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    if register:
        if isinstance(register, basestring): name = register
    return cmap_custom(colorDict, name, register=register, **kwargs)



def _generic_transform_(c, converter, channels, *args, **kwargs):

    # Colormap case
    if isinstance(c, Colormap):
        c = deepcopy(c)
        if not hasattr(c, '_lut'):
            c._init()
        c._lut[:-3, :3] = _generic_transform_(c._lut[:-3, :3], converter, channels, *args, **kwargs)
        for att in '_under', '_over', '_bad':
            col = getattr(c, '_rgba'+att, None)
            if col is not None:
                fset = getattr(c, 'set'+att)
                fset(_generic_transform_(col, converter, channels, *args, **kwargs)+(col[3], ))
        return c

    # Get RGB
    with_hsv = 'h' in channels or 's' in channels or 'v' in channels
    if isinstance(c, N.ndarray): # array
        r = c[:, 0]
        g = c[:, 1]
        b = c[:, 2]
        if with_hsv:
            hsv = rgb_to_hsv(c[:, :3])
            h = hsv[:, 0]
            s = hsv[:, 1]
            v = hsv[:, 2]
            back_to_rgb = lambda x, y, z: hsv_to_rgb(N.array([x, y, z]).T)[:, :3].T
    else:
        r, g, b, a = RGBA(c)
        if with_hsv:
            h, s, v = rgb_to_hsv((r, g, b))
            back_to_rgb = lambda x, y, z: hsv_to_rgb((x, y, z))

    # Convert channels
    conv = lambda c: N.clip(converter(c, *args, **kwargs), 0, 1)
    if 'r' in channels:
        r = conv(r)
    if 'g' in channels:
        g = conv(g)
    if 'b' in channels:
        b = conv(b)
    if with_hsv:
        if 'h' in channels:
            h = conv(h)
        if 's' in channels:
            s = conv(s)
        if 'v' in channels:
            v = conv(v)
        r, g, b = back_to_rgb(h, s, v)
    if isinstance(c, N.ndarray):
        c = c.copy()
        c[:, 0] = r
        c[:, 1] = g
        c[:, 2] = b
        return c
    return r, g, b


def _pull_toward_value_(x, f, v):
    """Move toward a value

    No effect with f=0
    Max effect with f=1
    """
    return v - (v-x) * (1-f)


def darken(c, f):
    """Darken a color 'c' by a factor 'f' (max when f=1)

    :Params:

        - **c**: Color
        - **f**: Factor between 0 and 1 with null impact at 0

    :Sample:

        >>> darken('r',0)
        (1.0, 0.0, 0.0)
        >>> darken('r',.5)
        (0.5, 0.0, 0.0)

    :See also: :func:`whiten`
    """
    converter = lambda x, f: _pull_toward_value_(x, f, 0)
    return _generic_transform_(c, converter, 'rgb', f)

def whiten(c, f):
    """Whiten a color 'c' by a factor 'f' (max when f=1)

    :Params:

        - **c**: Color
        - **f**: Factor between 0 and 1 with null impact at 0

    :Sample:

        >>> whiten('r',.5)
        (1.0, 0.5, 0.5)
        >>> whiten('r',0)
        (1.0, 0.0, 0.0)

    :See also: :func:`darken`
    """
    converter = lambda x, f: _pull_toward_value_(x, f, 1)
    return _generic_transform_(c, converter, 'rgb', f)

def change_luminosity(c, f):
    """Change the luminosity

    :Params:
        - **c**: color
        - **f**: Factor between 0 and 1 with null impact at 0.5
    """
    if f==0.5:
        return c
    f = 2*f-1
    if f>0:
        return whiten(c, f)
    return darken(c, -f)

def to_shadow(c,att=.3):
    """Return the shadow color of a color

    .. note::

        This is a simple shortcut to :func:`darken`
        with ``f=.3`` be default.

    :Params:

        - **c**: Color.
        - **att**, optional: Attenuation factor
    """
    return darken(c, 1-att)

def to_grey(c, f, g=.5):
    """Pull a color toward a single grey vzlue

    :Params:

        - **c**: Color
        - **f**: Factor between 0 and 1.
          When max, color is converted to medium grey ("0.5").
        - **g**: Value of the grey

    :Sample:

        >>> to_grey('r', 0)
        (1.0, 0.0, 0.0)
        >>> to_grey('r',1)
        (0.5, 0.5, 0.5)

    :See also: :func:`whiten` : :func:`darken`
    """
    converter = lambda x, f, g: _pull_toward_value_(x, f, g)
    return _generic_transform_(c, converter, 'rgb', f, g)

def saturate(c, f):
    """Saturate a color or a colormap"""
    return _generic_transform_(c, _pull_toward_value_, 's', f, 1)
def desaturate(c, f):
    """Desaturate a color or a colormap"""
    return _generic_transform_(c, _pull_toward_value_, 's', f, 0)

def change_saturation(c, f):
    """Change the saturation in HSV mode

    :Params:

        - **c**: Color
        - **f**: Factor between 0 and 1. Null effect at 0.5.
    """
    f = 2*f-1
    if f>0:
        return saturate(c, f)
    return desaturate(c, -f)

def change_value(c, f):
    """Change de value in HSV mode

    :Params:

        - **c**: Color
        - **f**: Factor between 0 and 1. Null effect at 0.5.
    """
    if f>.5:
        return _generic_transform_(c, _pull_toward_value_, 'v', 2*(f-.5), 1)
    return _generic_transform_(c, _pull_toward_value_, 'v', -2*(f-.5), 0)

def pastelise(c, s=.25, v=.9):
    """Make a color more pastel

    Equivalent to::

        >>> c = change_value(c, v)
        >>> c = change_saturation(c, s)

    """
    c = change_value(c, v)
    return change_saturation(c, s)
pastelize = pastelise

class StepsNorm(Normalize):
    """Normalize a given value to the 0-1 range on a stepped linear or log scale

    See tutorial :ref:`user.tut.misc.plot.advanced.stepsnorm`
    """
    def __init__(self, levels, log=False, masked=True, **kwargs):
        levels = N.array(levels)
        Normalize.__init__(self, **kwargs)
        if self.vmin is None:
            self.vmin = levels.min()
        if self.vmax is None:
            self.vmax = levels.max()
        self.levels = N.unique(N.clip(levels, self.vmin, self.vmax)).astype('d')
        self.positions = N.linspace(0, 1., len(self.levels))
        self.log = log
        if self.vmin > self.vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif self.log and self.vmin <= 0.:
            raise ValueError("minvalue must be greater than 0 when using log scale")
        self._masked = masked


    @staticmethod
    def process_value(value):
        """
        Homogenize the input *value* for easy and efficient normalization.

        *value* can be a scalar or sequence.

        Returns *result*, *is_scalar*, where *result* is a
        masked array matching *value*.  Float dtypes are preserved;
        integer types with two bytes or smaller are converted to
        N.float32, and larger types are converted to N.float.
        Preserving float32 when possible, and using in-place operations,
        can greatly improve speed for large arrays.

        Experimental; we may want to add an option to force the
        use of float32.
        """
        if cbook.iterable(value):
            is_scalar = False
            result = ma.asarray(value)
            if result.dtype.kind == 'f':
                if isinstance(value, N.ndarray):
                    result = result.copy()
            elif result.dtype.itemsize > 2:
                result = result.astype(N.float)
            else:
                result = result.astype(N.float32)
        else:
            is_scalar = True
            result = ma.array([value]).astype(N.float)
        return result, is_scalar


    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)
        val = result.copy()
        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax
        mask = ma.getmaskarray(val)

        # if theses 4 lines are not present, color of the 2D scalar field
        # in map2 is black for the values whose range is over the max value in colorbar
#        if cbook.iterable(value):
#            val = N.asarray(value).astype(N.float)
#        else:
#            val = N.array([value]).astype(N.float)
        #

        if vmin==vmax:
            result.fill(0.)
        else:
            if clip:
                result = ma.array(N.clip(val.filled(vmax), vmin, vmax),
                    mask=mask)
            if self.log and N.any(self.levels<=0):
                raise ValueError("All levels must be greater than 0 when using log scale")
            nlev = len(self.levels)

            if self.log:
                val = ma.log(val)

            # Inside
            for ilev in xrange(nlev-1):
                lev0, lev1 = self.levels[ilev:ilev+2]
                p0, p1 = self.positions[ilev:ilev+2]
                if self.log:
                    lev0 = ma.log(lev0)
                    lev1 = ma.log(lev1)
                mincheck = (val>=lev0) if ilev>0 else True
                maxcheck = (val<=lev1) if ilev<nlev-2 else True
                result[:] = ma.where(mincheck&maxcheck,
                    p0+(p1-p0)*(val-lev0)/(lev1-lev0), result)

            # Above (linear extrapolation)
            lev0, lev1 = self.levels[-2:]
            p0, p1 = self.positions[-2:]
            if self.log:
                lev0 = ma.log(lev0)
                lev1 = ma.log(lev1)
            result[:] = ma.where(val>=lev1, p1 + (p1-p0)*(val-lev1)/(lev1-lev0), result)
            result[N.isinf(val)] = N.inf

            # Below (linear extrapolation)
            lev0, lev1 = self.levels[:2]
            p0, p1 = self.positions[:2]
            if self.log:
                lev0 = ma.log(lev0)
                lev1 = ma.log(lev1)
            result[:] = ma.where(val<lev0, p0 + (p0-p1)*(val-lev0)/(lev0-lev1), result)
            result[N.isneginf(val)] = -N.inf

        if self._masked:
            result[mask] = N.ma.masked

        if is_scalar:
            result = result[0]
        return result

    def inverse(self, pos):
        result, is_scalar = self.process_value(pos)
        mask = N.ma.getmaskarray(result)
        pos = result.copy()
        if self.vmin==self.vmax:
            result.fill(0)
        else:
#            result = pos*0.
#            result[pos<0] = self.vmin
#            result[pos>1] = self.vmax

            # Inside
            nlev = len(self.levels)
            for ilev in xrange(len(self.levels)-1):
                lev0, lev1 = self.levels[ilev:ilev+2]
                p0, p1 = self.positions[ilev:ilev+2]
                if self.log:
                    lev0 = ma.log10(lev0)
                    lev1 = ma.log10(lev1)
                mincheck = (pos>=p0) if ilev else True
                maxcheck = (pos<=p1) if ilev<nlev-2 else True
                result[:] = N.where(mincheck&maxcheck, lev0+(pos-p0)*(lev1-lev0)/(p1-p0), result)
                if self.log:
                    result = N.power(result)

            # Above (linear extrapolation)
            lev0, lev1 = self.levels[-2:]
            p0, p1 = self.positions[-2:]
            if self.log:
                lev0 = ma.log(lev0)
                lev1 = ma.log(lev1)
            result[:] = ma.where(pos>=p1, lev1 + (lev1-lev0)*(pos-p1)/(p1-p0), result)
            result[N.isinf(pos)] = N.inf

            # Below (linear extrapolation)
            lev0, lev1 = self.levels[:2]
            p0, p1 = self.positions[:2]
            if self.log:
                lev0 = ma.log(lev0)
                lev1 = ma.log(lev1)
            result[:] = ma.where(pos<p0, lev0 + (lev0-lev1)*(pos-p0)/(p0-p1), result)
            result[N.isneginf(pos)] = -N.inf


        if self._masked:
            result[mask] = N.ma.masked

        if is_scalar:
            result = result[0]
        return result


class RangedLinearSegmentedColormap(Colormap):
    """
    :See also:

        :class:`matplotlib.colors.LinearSegmentedColormap`
    """
    def __init__(self, name, segmentdata, N=256, start=0., stop=1.):
        self.monochrome = False
        Colormap.__init__(self, name, N)
        self._segmentdata = segmentdata
        self._start = start
        self._stop = stop

    def _init(self):
        self._lut = N.ones((self.N + 3, 4), N.float)
        n = int(N.ceil(self.N/(self._stop-self._start)))
        r = makeMappingArray(n, self._segmentdata['red'])
        g = makeMappingArray(n, self._segmentdata['green'])
        b = makeMappingArray(n, self._segmentdata['blue'])
        i0 = int(self._start*n)
        i1 = i0+self.N
        self._lut[:-3, 0] = r[i0:i1]
        self._lut[:-3, 1] = g[i0:i1]
        self._lut[:-3, 2] = b[i0:i1]
        self._isinit = True
        self._set_extremes()

class Scalar2RGB(object):
    """Converter from scalar to colors

    :Params:

        - **vminmax**: Either an array or a tuple of (min,max).
        - **cmap**, optional: A colormap.

    :Example:

        >>> cmap = cmap_srs(['b', 'r'])
        >>> c = Scalar2RGB((1.5, 20.6), cmap)
        >>> print c(1.5, alpha=.5), c(10), c(20.6), c(50)
        (0.0,0.0,1.0,0.5) (0.38,0.0,0.61) (1.0,0.0,0.0) (1.0,0.0,0.0)
        >>> print c([1.5,10])
        [[ 0.          0.          1.        ]
        [ 0.38627451  0.          0.61372549]]

    """

    def __init__(self, vminmax, cmap=None):
        cmap = get_cmap(cmap)
        self.cmap = cmap
        from genutil import minmax
        vmin, vmax = minmax(vminmax)
        self.norm = Normalize(vmin, vmax)
        self.sm = ScalarMappable(cmap=self.cmap, norm=self.norm)

    def __call__(self, value, alpha=None):
        """Convert value to RGB(A)

        :Params:

            - **value**: A scalar or array.
            - **alpha**, optional: Add an alpha value to
              all colors.
        """

        res =  self.sm.to_rgba(value, alpha)
        if alpha is not None: return res
        if isinstance(res, tuple): return res[:3]
        return res[:,:3]

def anamorph_cmap(cmap, transform, name=None):
    """Tranform a colormap with anamorphim

    :Params:

        - **cmap**: Colormap.
        - **transform**: Sorted array of float between 0 and 1.
        - **name**: name of the colormap
    """
    #  Input cmap
    cmap = P.get_cmap(cmap)

    # Transform
    transform = N.clip(N.atleast_1d(transform).astype('f'), 0, 1).tolist()
    if transform[0]!=0:
        transform = [0.] + transform
    if transform[-1]!=1:
        transform.append(1.)

    # First colors
    segmentdata = {'red':[], 'green':[], 'blue':[]}

    # Input segmentdata
    input_segmentdata = cmap._segmentdata.copy()
    if callable(input_segmentdata['red']):
        for cname, cfunc in input_segmentdata.items():
            xind = N.linspace(0, 1, cmap.N) ** cmap._gamma
            lut = N.clip(N.array(cfunc(xind), dtype=N.float), 0, 1)
            input_segmentdata[cname] = []
            for i, xi in enumerate(xind):
                input_segmentdata[cname].append((xi, lut[i], lut[i]))

    # Loop on intervals
    dxo = 1./(len(transform)-1)
    for ii, xi0 in enumerate(transform[:-1]):
        xi1 = transform[ii+1]
        dxi = xi1-xi0
        xo0 = ii*dxo
        xo1 = (ii+1)*dxo
        for cname, cvals in input_segmentdata.items():
            for cval in cvals:
                xi = cval[0]
                if xi<xi0:
                    continue
                if xi>=xi1:
                    break
                xo = (xi-xi0)*dxo/dxi+xo0
                segmentdata[cname].append((xo,)+cval[1:])
            else:
                col = cmap(xi1)[['red', 'green', 'blue'].index(cname)]
                segmentdata[cname].append((xo1, col, col))
    for cname, cval in input_segmentdata.items():
        segmentdata[cname].append(cval[-1])

    # Name
    if name is None:
        name = cmap.name
        if name is None:
            name = 'vacumm'
        name += '_anamorph'

    # Create it
    cmapo = LinearSegmentedColormap(name, segmentdata, cmap.N)

    # Register it
    P.register_cmap(name, cmapo)

    return cmapo

def discretise_cmap(cmap, bounds, name=None, **kwargs):
    """Make discret an existing colormap

    :Examples:

        >>> discretise_cmap('jet', [.25, .5, .9]) # two not evenly spaced colors
        >>> discretise_cmap('jet', 10) # ten evenly spaced colors

    :Params:

        - **cmap**: Colormap.
        - **bounds**: An array of limits that will normalized.
          If a scalar, it is converted into an array of 'scalar' values ranging
          from 0 to 1.
        - **name**, optional: Name of the colormap.
        - Other params are passed to :func:`cmap_custom`.
    """
    from vcmq import P, cmap_custom, N, plot_cmap
    old_cmap = P.get_cmap(cmap)

    if N.isscalar(bounds):
        bounds = N.linspace(0, 1, bounds)
    else:
        bounds = N.asarray(bounds)
        bounds = (bounds-bounds[0])/(bounds[-1]-bounds[0])

    centers = 0.5*(bounds[:-1]+bounds[1:])
    colors = [old_cmap(c) for c in centers]

    data = []
    for ic, color in enumerate(colors):
        data.extend([(color, bounds[ic]), (color, bounds[ic+1])])

    if name is None:
        name = old_cmap.name+'_discrete'
    new_cmap = cmap_custom(data, name=name, **kwargs)
    return new_cmap

discretize_cmap = discretise_cmap

# Register colormaps
# - vacumm
cmap_bwr()
cmap_bwre()
cmap_br()
cmap_br2()
cmap_wr()
cmap_wre()
cmap_bathy()
cmap_land()
cmap_topo()
cmap_jete()
cmap_ajete()
cmap_jet()
cmap_jets()
cmap_wjets()
cmap_ajets()
cmap_wjet()
cmap_pe()
cmap_grey()
cmap_chla()
cmap_previmer()
cmap_previmer2()
cmap_rnb2_hymex()
cmap_rainbow_sst_hymex()
cmap_dynamic_cmyk_hymex()
cmap_white_centered_hymex()
cmap_red_tau_hymex()
cmap_ncview_rainbow()
cmap_eke()
cmap_currents()
cmap_nice_gfdl()
cmap_ssec()
# - basemap
for _name in basemap_cm.datad.keys():
#    sname = _name.lower() if _name.startswith('GMT_') else _name
    _cmap = getattr(basemap_cm, _name)
    P.register_cmap(_name, _cmap)
del _cmap, _name
# - cmocean
if cmoceancm is not None:
    for cmname in cmoceancm.cmapnames:
        for suffix in '', '_r':
            cmapname = cmname + suffix
            cmap = getattr(cmoceancm, cmapname)
            P.register_cmap('cmocean_' + cmapname, cmap)
            P.register_cmap(cmapname, cmap)

#: Colormap for posivite data
CMAP_POSITIVE = 'speed'

#: Colormap for negative data
CMAP_NEGATIVE = 'tempo_r'

#: Colomap for anomalies / symetric data
CMAP_SYMETRIC = 'balance'
CMAP_ANOMALY = CMAP_SYMETRIC
