# -*- coding: utf8 -*-
"""Common utilities"""
# Copyright or © or Copr. Actimar/IFREMER (2013-2015)
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

from vcmq import *
from ConfigParser import SafeConfigParser
from vacumm.misc.color import cmap_srs
from matplotlib.transforms import offset_copy, blended_transform_factory
import matplotlib.image as mpimg
import locale

lang = 'fr_FR.UTF-8'
os.environ['LANG'] = lang
locale.setlocale(locale.LC_ALL,lang)
os.environ['LC_NUMERIC'] = 'en_US.UTF-8'
locale.setlocale(locale.LC_NUMERIC, 'en_US.UTF-8')


class PREVIMERError(Exception):
    pass

# Load configuration here
cfgfiles = [
    os.path.join(os.path.dirname(__file__[:-2]), 'ifroco.cfg'),
    'previmer.cfg', 'ifroco.cfg'
]
cfg = SafeConfigParser()

def load_cfg(cfiles=None):
    """Load configuration

    :Params:

        - **cfiles**: single or list of alternative config files.
    """
    if not cfiles: cfiles = cfgfiles
    return cfg.read(cfiles)

cfgfiles = load_cfg()


def cfgget(option, section='DEFAULT', alt=None, reload=False):
    """Get an evaluated config option

    .. note:: Convert strings to unicode
    """
    if reload: load_cfg(reload)
    if not cfg.has_option(section, option): return alt
    value = eval(cfg.get(section, option))
    if isinstance(value, str):
        value = value.decode('utf8')
    return value

def set_mpl_defaults():
    """Set some Mtpotlib plot defaults"""
    grid_color = cfgget('grid_color')
    rc('grid', color=grid_color, linestyle='dashed')
    rc('font', size=9)

def add_background(fig=None):
    """Add a linear background to the figure"""
    if fig is None: fig = P.gcf()
    colorb = cfgget('bgcolor_bottom')
    colort = cfgget('bgcolor_top')
    cmap = cmap_srs([colorb, colort])
    bgaxes = fig.add_axes([0, 0, 1, 1])
    bgaxes.set_axis_off()
    bgaxes.imshow([[0, 0], [1, 1]], origin='lower', interpolation='bilinear',
        cmap=cmap, aspect='auto', vmin=0, vmax=1)

def add_title(texts, fig=None, addlegtime=1):
    """Add a multiline title"""

    # Some inits
    if fig is None: fig = P.gcf()
    y = 1-0.01
    color = cfgget('title_color')
    family = cfgget('title_font')
    title_size = cfgget('title_size')
    title_weight = cfgget('title_weight')
    title_legsize = cfgget('title_legsize')
    ren = fig.canvas.get_renderer()
    heights = []

    # Form inputs
    if not isinstance(texts, list):
        texts = [[texts]]
    if addlegtime: # add legal time
        legtext = strftime("(heure légale) mise à jour du %m/%d/%Y %Hh%M",now())
        legtext = legtext.decode('utf8')
        legtext = ' '+legtext, title_legsize, None
        if addlegtime==1: # on last line
            last = texts[-1]
            if not isinstance(last, list):
                last = [last]
            last.append(legtext)
            texts = texts[:-1]
            texts.append(last)
        else: # new line
            texts.append(legtext)

    # Loop on lines
    for i, texts in enumerate(texts):
        if not isinstance(texts, list):
            texts = [texts]

        # Several text on sameline?
        multi = len(texts)!=1
        sizes = []
        if multi:
            ha = ['right', 'left']
            widths = []
            tt = []
        else:
            ha = ['center']
        heights = []

        # Loop on text specs
        for it, text in enumerate(texts):

            # Specs
            if isinstance(text, tuple):
                text, size, weight = text
            else:
                size = title_size
                weight = title_weight

            # Alignment and xoffset
            if ha[it]=='left':
                tr = offset_copy(fig.transFigure, x=5, y=0, fig=fig, units='points')
            elif ha[it]=='right':
                tr = offset_copy(fig.transFigure, x=-5, y=0, fig=fig, units='points')
            else:
                tr = fig.transFigure

            # Plot it
            t = fig.text(.5, y, text, transform=tr, family=family, ha=ha[it],
                color=color, size=size, va='top', weight=weight)

            # Adjust
            tw = t.get_window_extent(renderer=ren).inverse_transformed(fig.transFigure)
            if multi:
                widths.append(tw.width)
                tt.append(t)
                if it==1:
                    x0 = 0.5*(1+widths[0]-widths[1])
                    tt[0].set_x(x0)
                    tt[1].set_x(x0)

            heights.append(tw.height)

        # Next position
        y -= max(heights)*1.2
    return y



def add_logos(logofiles=None, fig=None):
    """Add logos using its config name"""
    if logofiles is None: logofiles = cfgget('logo_def')
    lastx = 0.005
    if not isinstance(logofiles, list):
        logofiles = [logofiles]
    logodir = cfgget('logo_dir')
    for logofile in logofiles:
        if os.path.exists(logofile):
            imgfile = logofile
        else:
            imgfile = os.path.join(logodir, logofile)
            if not os.path.exists(imgfile):
                raise PREVIMERError("Can't find logo: "+logofile)
        imgm = mpimg.imread(imgfile)
        ny, nx = imgm.shape[:2]
        if fig is None: fig = P.gcf()
        xsize, ysize = fig.get_size_inches()
        dpi = fig.dpi*1.
        ax = fig.add_axes([lastx, 0.005, nx/dpi/xsize, ny/dpi/ysize], aspect=1,frameon=False)
        ax.set_axis_off()
        ax.imshow(imgm, origin='upper')
        lastx += nx/dpi/xsize

def add_copyright(text, fig=None):
    """Add a copyright info at the bottom right"""
    if fig is None: fig = P.gcf()
    size = cfgget('copyright_size')
    trans = offset_copy(fig.transFigure, x=-3, y=6, fig=fig, units='points')
    fig.text(.99, 0.005, text, size=size, ha='right', transform=trans, va='bottom')


def set_xlabel(text, ax=None, **kwargs):
    """Set the X label with appropriate default font properies"""
    if ax is None: ax = P.gca()
    dict_check_defaults(kwargs,
        color = cfgget('title_color'),
        family = cfgget('title_font'),
        size = cfgget('label_size'),
    )
    ax.set_xlabel(text, **kwargs)

def set_ylabel(text, ax=None, **kwargs):
    """Set the Y labe with appropriate font properiesl"""
    if ax is None: ax = P.gca()
    dict_check_defaults(kwargs,
        color = cfgget('title_color'),
        family = cfgget('title_font'),
        size = cfgget('label_size'),
    )
    ax.set_ylabel(text, **kwargs)

def add_axlabel(x, y, text, ax=None, dx=0, dy=0, ha='center', va='center', **kwargs):
    """A label to an axis with appropriate default font properies"""
    if ax is None: ax = P.gca()
    dict_check_defaults(kwargs,
        color = cfgget('title_color'),
        family = cfgget('title_font'),
        size = cfgget('label_size'),
    )
    trans = offset_copy(ax.transAxes, x=dx, y=dy, fig=ax.get_figure(), units='points')
    ax.text(x, y, text, transform=trans, ha=ha, va=va, **kwargs)

def format_label(labelobj, **kwargs):
    """For a label object with appropriate default font properies"""
    dict_check_defaults(kwargs,
        color = cfgget('title_color'),
        family = cfgget('title_font'),
        size = cfgget('label_size'),
    )
    P.setp(labelobj, **kwargs)
