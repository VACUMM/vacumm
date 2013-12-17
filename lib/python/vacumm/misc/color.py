# -*- coding: utf8 -*-
"""Variables and utilities about colors and color maps

List of all available colormaps in matplotlib, including VACUMM colormaps 
(plotted with :func:`plot_cmaps`)

.. image:: misc-color-cmaps.png

"""
# Copyright or © or Copr. Actimar (contributor(s) : Stephane Raynaud) (2010)
# 
# raynaud@actimar.fr
# 
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
import os,glob
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize, LinearSegmentedColormap, ColorConverter, makeMappingArray
import pylab as P
from matplotlib.cm import cmap_d
import numpy as N
ma = N.ma
import matplotlib.cbook as cbook

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

_gmt_cpt = '/soft/gmt/GMT4.1.1/share/cpt/' # !!!

RGB = ColorConverter().to_rgb
RGBA = ColorConverter().to_rgba

__all__ = ['cmap_custom', 'cmap_bwr', 'cmap_bwre', 'cmap_br', 'cmap_wr', 'cmap_wre', 'cmap_bathy', 
    'cmap_jete', 'cmap_ajete', 'cmap_jets', 'cmap_wjets', 'cmap_ajets', 'cmap_smoothed_steps', 
    'cmap_smoothed_regular_steps', 'cmap_ss', 'cmap_steps', 'cmap_regular_steps', 'cmap_rs', 
    'cmap_wjet', 'cmap_pe', 'cmap_grey', 'show_cmap', 'cmaps_mpl', 'cmaps_gmt', 'cmap_gmt',
    'cmaps_act', 'print_cmaps_gmt',  'darken', 'whiten', 'to_shadow', 'StepsNorm',
    'bistre', 'land', 'jean_pierre', 'decimaitre', 'RGB', 'RGBA', 'cmap_srs', 'cmap_rs', 
    'cmap_linear', 'ocean', 'sea', 'plot_cmap', 'plot_cmaps', 'get_cmap', 'simple_colors', 
    'cmap_land', 'cmap_topo', 'auto_cmap_topo', 'cmap_jet', 'cmap_rainbow', 'rainbow', 
    'cmap_magic', 'cmap_mg', 'RangedLinearSegmentedColormap','Scalar2RGB', 'cmap_chla', 
    'cmap_previmer', 'cmap_rnb2_hymex','cmap_rainbow_sst_hymex','cmap_dynamic_cmyk_hymex',
    'cmap_white_centered_hymex','cmap_red_tau_hymex', 'cmap_previmer2', 'cmap_ssec', 'cmap_ncview_rainbow']
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

def cmap_rainbow(n=None, name='rainbow', smoothed=True, mode='auto', **kwargs):
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
    
def cmap_magic(n=None, stretch = 0.4, mode='normal', white='.95', name='magic', **kwargs):
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

def cmap_custom(colors, name='mycmap', ncol=256, ranged=False, **kwargs):
    """Quick colormap creation
    
    :Params:
    
        - **colors**: Like [(color1, position1),(color2, position2), etc...] or 
          dict(red=((pos1,r1a,r1b), (pos2,r2a,r2b)),etc...)
    """
    if isinstance(colors, dict):
        for cname, cvals in colors.items():
            if cvals[0][0] != 0:
                colors[cname] = ((0, cvals[0][1], cvals[0][2]), ) + cvals
            if cvals[-1] != 1:
                colors[cname] = colors[cname] + ((1, cvals[-1][1], cvals[-1][2]), )
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
        return RangedLinearSegmentedColormap(name,colors,ncol, **kwargs)
    return LinearSegmentedColormap(name,colors,ncol)

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
    cmap = cmap_custom((('b', 0), (wcol, wpos), ('r', 1)), 'vacumm_bwr')
    P.register_cmap(name, cmap)
    return cmap


# blue -> white -> red for Extremes : violet and pink and ends
def cmap_bwre(wpos=0.5,gap=0.1, wcol='w', name='vacumm_bwre'):
    """Returns a violet->blue->white->red colormap->yellow
    
    - **white**, optional: relative position of white color in the map [default: 0.5]
    - **gap**, optional: Relative width of the pure central white gap [default: 0.1]
    
    :Sample: .. image:: misc-color-vacumm_bwre.png
    """
    wstart = wpos*(1-gap)
    wstop = wpos + (1.-wpos)*gap
    cmap =  cmap_custom((((1, 0, 1), 0), ('b', .1), (wcol, wstart), 
        (wcol, wstop), ('r', .9), ((1, 1, 0), 1)), 'vacumm_bwre')
    P.register_cmap(name, cmap)
    return cmap

# blue -> red
def cmap_br(sep=0.5, name='vacumm_br'):
    """Blue->red colormap
    
    :Params:
    
        - **sep**, optional: relative position of the blue/red transition in the map [default: 0.5]
    
    :Sample: .. image:: misc-color-vacumm_br.png
    """
    cmap = cmap_custom( (('r', 0), ((.5, 0, .5), sep), ('b', 1)), 'vacumm_br')
    P.register_cmap(name, cmap)
    return cmap
    
def cmap_br2(sep=0.5, white=.7, name="vacumm_br2"):
    """Blue->light blue|light red-> red
    
    :Params:
    
        - **sep**, optional: relative position of the blue/red transition in the map [default: 0.5]
        - **white**, optional: Strenght of the whitening (see :func:`whiten`).
    
    :Sample: .. image:: misc-color-vacumm_br.png
    """
    bsep = whiten('b', white)
    rsep = whiten('r', white)
    cdict =   {'red':   ((0,0,0),(sep,bsep[0],rsep[0]), (1,1,1)),
               'green': ((0,0,0),(sep,bsep[1],rsep[1]), (1,0,0)),
               'blue':  ((0,1,1),(sep,bsep[2],rsep[2]), (1,0,0)),
               }
    cmap =  LinearSegmentedColormap(name, cdict, 256)
    P.register_cmap(name, cmap)
    return cmap

# white -> red
def cmap_wr(name='vacumm_wr'):
    """White->red colormap
    
    :Sample: .. image:: misc-color-vacumm_wr.png
    """
    cmap = cmap_custom( (('w', 0), ('r', 1)), name)
    P.register_cmap(name, cmap)
    return cmap

# white -> red
def cmap_wre(name='vacumm_wre'):
    """White->red->yellow colormap for positive extremes
    
    :Sample: .. image:: misc-color-vacumm_wre.png
    """
    cdict = {'red':  ((0.,1.,1.),(.8,1.,1.),(1.,1.,1.)),
           'green':((0.,1.,1.),(.8,0.,0.),(1.,1.,1.)),
           'blue': ((0.,1.,1.),(.8,0.,0.),(1.,0.,0.))}
    cmap = LinearSegmentedColormap(name,cdict,256)
    P.register_cmap(name, cmap)
    return cmap

cy = 0.6 # cyan
ye = 0.95 # yellow
vi = 0.35 # violet
_colors_bathy = (('k',0), ((0,.1,.9),vi), ((0,.8,.8),cy), ((.9,.9,.6),.9), ((.9,.9,.8),1))
def cmap_bathy(start=0., stop=1., name='vacumm_bathy'):
    """Colormap for bathymetry maps
    
    :Sample: .. image:: misc-color-vacumm_bathy.png
    """
    this_cmap = cmap_custom(_colors_bathy, name, ranged=True, start=start, stop=stop)
    this_cmap.set_bad(bistre)
    this_cmap.set_under(bistre)
    this_cmap.set_over(bistre)
    P.register_cmap(name, this_cmap)
    return this_cmap
    
_colors_land = ((_colors_bathy[-1][0], 0), ('#145D0A',.2),('#62CF60', .4), ('#A1A156',.6),('#6D461D', .85), ('#eeeeff', 1.))
def cmap_land(start=0., stop=1., name='vacumm_land'):
    """Colormap for land maps
    
    :Params:
    
        - **start/stop**, optional: Positions for a zoom in colormap.
    
    :Sample: .. image:: misc-color-vacumm_land.png
    """
    this_cmap = cmap_custom(_colors_land, 'vacumm_land', ranged=True, start=start, stop=stop)
    this_cmap.set_bad(ocean)
    this_cmap.set_under(ocean)
    this_cmap.set_over(ocean)
    P.register_cmap(name, this_cmap)
    return this_cmap

def cmap_topo(start=0., stop=1., name='vacumm_topo', zero=.5, over=_colors_bathy[0][0], under=_colors_land[-1][0], bad='0.5'):
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
    this_cmap = cmap_custom(mycolors, name, ranged=True, start=start, stop=stop)
    this_cmap.set_bad(bad)
    this_cmap.set_under(under)
    this_cmap.set_over(over)
    P.register_cmap(name, this_cmap)
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
    return cmap_topo(start, stop, 'vacumm_auto_cmap_topo', zero=zero, **kwargs)
    
    

def cmap_jete(name='vacumm_jete'):
    """Jet colormap with extremes
    
    :Sample: .. image:: misc-color-vacumm_jete.png
    """
    cdict =   {'red':((0.,.75,.75),(0.22,0,0),  (0.4,0,0),(0.6,1,1),        (.78,1,1),  (1,1,1)),
               'green':((0.,0.,0.),(0.22,0,0),  (0.375,1,1),(0.64,1,1), (.78,0,0),  (1,0,0)),
               'blue':((0.,1.,1.),(0.22,1,1),   (0.4,1,1),(0.6,0,0),        (.78,0,0),  (1,.75,.75))}
    cmap = LinearSegmentedColormap(name,cdict,256)
    P.register_cmap(name, cmap)
    return cmap

def cmap_ajete(w=.02, name='vacumm_ajete'):
    """Jet colormap with white at center and extremes, for anomaly plots
    
    :Sample: .. image:: misc-color-vacumm_ajete.png
    """
    cdict =   {'red':((0.,.75,.75),(0.21,0,0),  (0.33,0,0),(.5-w/2,1,1),(.5+w/2,1,1),(0.67,1,1),        (.79,1,1),  (1,1,1)),
               'green':((0.,0.,0.),(0.21,0,0),  (0.3,1,1),(.5-w/2,1,1),(.5+w/2,1,1),(0.71,1,1),     (.79,0,0),  (1,0,0)),
               'blue':((0.,1.,1.),(0.21,1,1),   (0.33,1,1),(.5-w/2,1,1),(.5+w/2,1,1),(0.67,0,0),        (.79,0,0),  (1,.75,.75))}
    cmap = LinearSegmentedColormap(name,cdict,256)
    P.register_cmap(name, cmap)
    return cmap

def cmap_jet(smoothed=False, name='vacumm_jet', **kwargs):
    """Jet colormap"""
    colors = [(.75,0,1),(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0)]
    kwargs.setdefault('stretch', 0)
    if not smoothed:
        cmap = cmap_regular_steps(colors, name=name, **kwargs)
    else:
        cmap = cmap_jets(colors, name=name, **kwargs)
    P.register_cmap(name, cmap)
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
    cmap = cmap_smoothed_regular_steps(colors,name=name,**kwargs)
    P.register_cmap(name, cmap)
    return cmap

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
    cmap = cmap_smoothed_regular_steps(colors, name=name,**kwargs)
    P.register_cmap(name, cmap)
    return cmap

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
    cmap = cmap_smoothed_regular_steps(colors,name=name,**kwargs)
    P.register_cmap(name, cmap)
    return cmap

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
    
def cmap_smoothed_steps(colors, stretch=None, rstretch=0, lstretch=0, name='cmap_css', ncol=256, asdict=False):
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
    return LinearSegmentedColormap(name, cdict, ncol)
    
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

def cmap_steps(cols, stretch=None, lstretch=0., rstretch=0., keepc=None, name='cmap_steps'):
    """Colormap by steps
    
    :Params:
    
        - **cols**: [(col1,pos1),(col2,pos2),...]
        - **lstretch**, optional: Color darkening (<0) or whitening at start of steps (left)
        - **rstretch**, optional: Same but at end of steps (right)
        - **keepc**, optional: If ``lstretch`` and ``rstretch`` are both different from zero,
          it keeps center of steps intact, else it becomes a mean
          of left and right.
    
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
    return LinearSegmentedColormap(name, cdict, 256)

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
    cmap = func(colors, **kwargs)
    P.register_cmap(name, cmap)
    return cmap

#    wcol = kwargs.get('first_color', wcol)
#    colors = [RGB(wcol),(0,1,1),(0,1,0),(1,1,0),(1,.5,0),(1,0,0),(1,0,.75)]



def cmap_pe(red=.8, name="vacumm_pe"):
    """Colomap for positive extremes (white->grey->red)
    
    :Sample: .. image:: misc-color-vacumm_pe.png
    """
    cdict = {'red':  ((0.,1.,1.),(red,.7,.7),(1.,1.,1.)),
           'green':((0.,1.,1.),(red,.7,.7),(1.,0.,0.)),
           'blue': ((0.,1.,1.),(red,.7,.7),(1.,0.,0.))}
    cmap = LinearSegmentedColormap(name,cdict,256)
    P.register_cmap(name, cmap)
    return cmap

def cmap_grey(start=0, end=1., name="vacumm_grey"):
    """Grey colormap from position ``start`` to position ``end``
    
    :Sample: .. image:: misc-color-vacumm_grey.png
    """
    cdict = {'red':  ((0.,start,start),(1.,end,end)),
           'green':((0.,start,start),(1.,end,end)),
           'blue': ((0.,start,start),(1.,end,end))}
    cmap = LinearSegmentedColormap(name, cdict, 256)
    P.register_cmap(name, cmap)
    return cmap

def cmap_chla(name='vacumm_chla', smoothed=True):
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
    cmap = func(colors, name=name)
    P.register_cmap(name, cmap)
    return cmap

def cmap_previmer(name='vacumm_previmer'):
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
    cmap = cmap_custom(colors,name=name,ncol=len(colors)-1)
    P.register_cmap(name, cmap)
    return cmap

def cmap_previmer2(name='vacumm_previmer2'):
    """Colormap from PREVIMER Website (http://www.previmer.org)
    
    Same as :func:`cmap_previmer` but extremes are used for 
    methd:`set_under` and :meth`set_over`.
    
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
    cmap = cmap_srs(colors[1:-1], name=name)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    P.register_cmap(name, cmap)
    return cmap    

def cmap_ssec(name='vacumm_ssec'):
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
    cmap = cmap_srs(colors[1:-1], name=name)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    P.register_cmap(name, cmap)
    return cmap

def cmap_ncview_rainbow(name='vacumm_ncview_rainbow'):
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
    cmap = cmap_srs(colors[1:-1], name=name)
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    P.register_cmap(name, cmap)
    return cmap

def cmap_rnb2_hymex(name="vacumm_rnb2_hymex"):
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
    cmap = cmap_smoothed_steps(data, name=name)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    P.register_cmap(name, cmap)
    return cmap

def cmap_rainbow_sst_hymex(name="vacumm_rainbow_sst_hymex"):
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
    cmap = cmap_smoothed_steps(data, name=name)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    P.register_cmap(name, cmap)
    return cmap

def cmap_dynamic_cmyk_hymex(name="vacumm_dynamic_cmyk_hymex"):
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
    cmap = cmap_smoothed_steps(data, name=name)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    P.register_cmap(name, cmap)
    return cmap

def cmap_white_centered_hymex(name="vacumm_white_centered_hymex"):
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
    cmap = cmap_smoothed_steps(data, name=name)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    P.register_cmap(name, cmap)
    return cmap

def cmap_red_tau_hymex(name="vacumm_red_tau_hymex"):
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
    cmap = cmap_smoothed_steps(data, name=name)
    cmap.set_under(data[0][0])
    cmap.set_over(data[-1][0])
    P.register_cmap(name, cmap)
    return cmap

def get_cmap(cmap=None, **kwargs):
    """A simple way to get a VACUMM, GMT or MPL colormap
    
    :Example:
    
        >>> get_cmap('jet')         # Matplotlib
        >>> pylab.jet()             # Matplotlib
        >>> pylab.get_cmap('jet')   # Matplotlib
        >>> get_cmap('cmap_grey',start=.2)   # VACUMM
        >>> get_cmap('vacumm_grey')   # VACUMM (registered in MPL)
        >>> cmap_grey(start=.2)     # VACUMM
        >>> get_cmap('gmt_gebco')   # GMT
        >>> cmap_gmt('gebco')       # GMT
    
    :See also:

        :func:`cmap_gmt` :func:`matplotlib.pyplot.get_cmap`
    """
    if isinstance(cmap, Colormap):
       return cmap 
    if isinstance(cmap, str):
        if cmap.startswith('cmap_'):
            cmap = eval(cmap)(**kwargs)
        elif cmap.startswith('gmt_'):
            cmap = cmap_gmt(cmap)
        elif  cmap.startswith('vacumm_') and kwargs:
            cmap = eval('cmap_'+cmap[7:])(**kwargs)
        else:
            cmap = P.get_cmap(cmap)
    else:
        cmap = P.get_cmap()
    return cmap


def plot_cmap(cmap, ncol=None, smoothed=True,  ax=None, figsize=(5, .25), show=True, aspect=.05, title=None, 
    sa=dict(left=.0, right=1, top=1, bottom=.0), savefig=None, savefigs=None, **kwargs):
    """Display a colormap"""
    from vacumm.misc import kwfilter
    cmap = get_cmap(cmap)
    if ax is None:
        if figsize is not False:
            P.figure(figsize=figsize)
        ax = P.gca()
    P.subplots_adjust(**sa)
    ax.set_aspect('auto', anchor='C')
    P.axis("off")
    if ncol is None:
        ncol = cmap.N if hasattr(cmap, 'N') else 256.
    x = N.arange(0,ncol,1.)
    dx = x[1]-x[0]
    interp = 'bilinear' if smoothed else "nearest"
    p = P.imshow(P.outer(N.ones(1),x), aspect=aspect*ncol, cmap=cmap, origin="lower", 
        interpolation=interp, vmin=-.5, vmax=x[-1]+0.5)
    P.fill([x[0]-dx/2, x[0]-dx/2, x[0]-dx/2-ncol/20.], [-.5, .5, 0], facecolor=p.to_rgba(-1), 
        edgecolor='none', linewidth=0)
    P.fill([x[-1]+dx/2, x[-1]+dx/2, x[-1]+dx/2+ncol/20.], [-.5, .5, 0], facecolor=p.to_rgba(ncol+2), 
        edgecolor='none', linewidth=0)
    if title is None: title = cmap.name
    if title is not None and title is not False:
        from core_plot import add_glow
        add_glow(P.text(ncol/2., 0., title, color='k', ha='center', va='center'), alpha=0.5)
        
    # Save and show
    if savefig is not None:
        P.savefig(savefig, **kwfilter(kwargs, 'savefig'))
    if savefigs is not None:
        from plot import savefigs as Savefigs
        Savefigs(savefigs, **kwfilter(kwargs, 'savefigs'))
    if show:
        P.show()


def plot_cmaps(cmaps=None, figsize=None, show=True, savefig=None, nrow=50, savefigs=None, 
    aspect=0.05, **kwargs):
    """Display a list of or all colormaps"""
    
    from vacumm.misc import kwfilter
    kwsf = kwfilter(kwargs, 'savefig')
    kwsfs = kwfilter(kwargs, 'savefigs')
    
    # Default colormap list
    from vacumm.misc import kwfilter
    if cmaps is None:
        cmaps = cmaps_mpl(names=False)+cmaps_gmt()
    elif isinstance(cmaps, str):
        if cmaps.startswith('act'):
            cmaps = cmaps_act(names=False)
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
    nrow = min(ncmap, nrow)
    ncol = (ncmap-1)/nrow+1
    if figsize is None:
        figsize = (5*ncol, nrow*5*aspect*1.1)
    if figsize is not False:
        P.figure(figsize=figsize)
        
    # Loop on colormaps
    ip = N.arange(ncol*nrow).reshape((nrow,ncol)).T.ravel()+1
    for i, cmap in enumerate(cmaps):
        ax = P.subplot(nrow, ncol, ip[i])
        plot_cmap(cmap, show=False, ax=P.gca(), sa=dict(bottom=.01, top=.99, left=.01, right=.99), 
            aspect=aspect, **kwargs)
        ax.set_anchor('C')
    P.subplots_adjust(wspace=0.05, top=0.98, bottom=0.02, left=0.02, right=0.98)
    
    # Save and show
    if savefig is not None:
        P.savefig(savefig, **kwsf)
    if savefigs is not None:
        from plot import savefigs as Savefigs
        Savefigs(savefigs, **kwsfs)
    if show: P.show()
    
def show_cmap(cmap, *args, **kwargs):
    """Alias for :func:`plot_cmap`"""
    return plot_cmap(cmap, *args, **kwargs)
    
def cmaps_mpl(vacumm=True, names=True):
    """List available Matplotlib colormaps
    
    :See also:
    
        :func:`get_cmap`
    """
    cmap_names = cmap_d.keys()
    cmap_names.sort()
    if not vacumm:
        cmap_names = [name for name in cmap_names if not name.startswith('vacumm_')]
    if names:
        return cmap_names
    return [cmap_d[name] for name in cmap_names]

def cmaps_act(names=True):
    """List available VACUMM colormaps
    
    :See also:
    
        :func:`get_cmap`
    """
    cmap_names = cmap_d.keys()
    cmap_names.sort()
    cmap_names = [name for name in cmap_names if name.startswith('vacumm_')]
    if names:
        return cmap_names
    return [cmap_d[name] for name in cmap_names]
    
def cmaps_gmt():
    """List available GMT colormaps
    
    :See also:
        
        :func:`cmap_gmt` :func:`print_cmaps_gmt` :func:`get_cmap`
    """
    cmaps = []
    for file in glob.glob(_gmt_cpt+'GMT_*.cpt'):
        f = open(file)
        lines = f.readlines()
        f.close()
        cmaps.append(os.path.basename(file)[4:-4])
    return cmaps
        
def print_cmaps_gmt():
    """List available gmt colormaps
    
    :See also:
        
        :func:`cmap_gmt` :func:`cmaps_gmt` :func:`get_cmap`
    """
    print 'List of available GMT colormaps :'
    for file in glob.glob(_gmt_cpt+'GMT_*.cpt'):
        f = open(file)
        lines = f.readlines()
        f.close()
        name = os.path.basename(file)[4:-4]
        description = lines[2].strip('# \t\n')
        print '%s %s' % ((name+':').ljust(20),description)
        

def cmap_gmt(name):
    """Get a colormap from GMT
    
    :See also:
        
        :func:`cmaps_gmt` :func:`print_cmaps_gmt` :func:`get_cmap`
    """
    import colorsys
    filePath = _gmt_cpt+'GMT_'+name+'.cpt'
    if not os.path.exists(filePath):
        print "Bad color map name (%s). Please use cmap_gmt_list() to get the list of available GMT colormaps." % name
        return None
    f = open(filePath)
    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if l[0] == "#":
            if ls[-1] == "HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

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
    
    return LinearSegmentedColormap('gmt_'+name,colorDict)

def darken(c, f):
    """Darken a color 'c' by a factor 'f' (max when f=1)
    
    :Params:
    
        - **c**: Color
        - **f**: Factor between 0 and 1
        
    :Sample:
    
        >>> darken('r',0)
        (1.0, 0.0, 0.0)
        >>> darken('r',.5)
        (0.5, 0.0, 0.0)

    :See also: :func:`whiten`
    """
    if isinstance(c, Colormap):
        c._init()
        c._lut[:-3, :-1] = [darken(rgb,f) for rgb in c._lut[:-3, :-1]]
        for att in '_under', '_over', '_bad':
            col = getattr(c, '_rgba'+att, None)
            if col is not None:
                func = getattr(c, 'set'+att)
                func(darken(col,f)+(col[3], ))             
        return c
    r,g,b,a = RGBA(c)
    f = 1-f
    r = f*r
    g = f*g
    b = f*b
    return r,g,b
    
def whiten(c, f):
    """Whiten a color 'c' by a factor 'f' (max when f=1)
    
    :Params:
    
        - **c**: Color
        - **f**: Factor between 0 and 1
        
    :Sample:
    
        >>> whiten('r',.5)
        (1.0, 0.5, 0.5)
        >>> whiten('r',0)
        (1.0, 0.0, 0.0)

    :See also: :func:`whiten`
    """
    if isinstance(c, Colormap):
        c._init()
        c._lut[:-3, :-1] = [whiten(rgb,f) for rgb in c._lut[:-3, :-1]]
        for att in '_under', '_over', '_bad':
            col = getattr(c, '_rgba'+att, None)
            if col is not None:
                func = getattr(c, 'set'+att)
                func(whiten(col,f)+(col[3], ))             
        return c
    r,g,b,a = RGBA(c)
    f = 1-f
    r = 1.-f*(1-r)
    g = 1.-f*(1-g)
    b = 1.-f*(1-b)
    return r,g,b
    

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



class StepsNorm(Normalize):
    """Normalize a given value to the 0-1 range on a stepped linear or log scale
    
    See tutorial :ref:`user.tut.misc.plot.advanced.stepsnorm`
    """
    def __init__(self, levels, log=False, **kwargs):
        levels = N.array(levels)
        Normalize.__init__(self, **kwargs)
        if self.vmin is None:
            self.vmin = levels.min()
        if self.vmax is None:
            self.vmax = levels.max()
        self.levels = N.unique(N.clip(levels, self.vmin, self.vmax))
        self.positions = N.linspace(0, 1., len(self.levels))
        self.log = log
        if self.vmin > self.vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif self.log and self.vmin <= 0.:
            raise ValueError("minvalue must be greater than 0 when using log scale")
        
        
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
        
        if cbook.iterable(value):
            val = N.asarray(value).astype(N.float)
        else:
            val = N.array([value]).astype(N.float)
        
        if vmin==vmax:
            result.fill(0.)
        else:
            if clip:
                mask = ma.getmask(val)
                result = ma.array(N.clip(val.filled(vmax), vmin, vmax),
                    mask=mask)
#            else:
#                result = ma.array(val)*ma.nomask
#                print result
#                result[val>self.vmax] = 1.1
#                result[val<self.vmin] = -0.1
            if self.log and N.any(self.levels<=0):
                raise ValueError("All levels must be greater than 0 when using log scale")
            nlev = len(self.levels)
            for ilev in xrange(nlev-1):
                lev0, lev1 = self.levels[ilev:ilev+2]
                p0, p1 = self.positions[ilev:ilev+2]
                if self.log:
                    val = ma.log(val)
                    lev0 = ma.log(lev0)
                    lev1 = ma.log(lev1)
                mincheck = (val>=lev0) if ilev>0 else True
                maxcheck = (val<=lev1) if ilev<nlev-2 else True
                result[:] = ma.where(mincheck&maxcheck, p0+(p1-p0)*(val-lev0)/(lev1-lev0), result)
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, pos):
        result, is_scalar = self.process_value(pos)
        pos = result.copy()
        if self.vmin==self.vmax:
            result.fill(0)
        else:
#            result = pos*0.
#            result[pos<0] = self.vmin
#            result[pos>1] = self.vmax
            nlev = len(self.levels)
            for ilev in xrange(len(self.levels)-1):
                lev0, lev1 = self.levels[ilev:ilev+2]
                p0, p1 = self.positions[ilev:ilev+2]
                if self.log:
                    val = ma.log10(val)
                    lev0 = ma.log10(lev0)
                    lev1 = ma.log10(lev1)
                mincheck = (pos>=p0) if ilev else True
                maxcheck = (pos<=p1) if ilev<nlev-2 else True
                result[:] = N.where(mincheck&maxcheck, lev0+(pos-p0)*(lev1-lev0)/(p1-p0), result)
                if self.log:
                    result = N.power(result)
            
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
        if cmap is None: cmap = P.get_cmap()
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
        

# Register colormaps
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


#print 'importing colors 1'
