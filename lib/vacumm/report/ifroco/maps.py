# -*- coding: utf8 -*-
"""Plotting maps"""
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

from common import *


def nested_axes(rect1=None, margins=None, fig=None):
    """Setup two nested axes

    :Params:

        - **rect1**, optional: Rectangle specs for the inner axes
          [left,bottom,width,height] in figure coordinates.
          Default from config ``[maps] cax``.
        - **margin**, option: Margins in pixel around the inner axes
          [left,bottom,right,top].
           Default from config ``[maps] cax_margins``.

    :Return: `ax0, ax1`
    """
    if fig is None: fig = P.gcf()

    # As points
    if rect1 is None: rect1 = cfgget('cax', 'maps')
    if margins is None: margins = cfgget('cax_margins', 'maps')
    p0 = rect1[:2]
    p1 = (p0[0]+rect1[2], p0[1]+rect1[3])
    trans = fig.transFigure.transform_point

    # Convert to pixels and add margin
    p0 = trans(p0) - margins[:2]
    p1 = trans(p1) + margins[2:]

    # Convert back to figure coordinates
    trans = fig.transFigure.inverted().transform_point
    p0 = trans(p0)
    p1 = trans(p1)
    rect0 = [p0[0], p0[1], p1[0]-p0[0], p1[1]-p0[1]]

    # create axes and add them to the figure
    ax0 = fig.add_axes(rect0)
    ax1 = fig.add_axes(rect1)

    return ax0, ax1


def get_plot_specs(lon, lat, xfigsize=None, axes_rect=None):
    """Get plot area specs

    :Params:

        - **lon/lat**: Data coordinates to plot.
    """
    # X-figsize
    xfigsize = xfigsize or cfgget('figsize')[0]

    # Mercator extensions of the area
    lon = N.asarray([lon.min(),lon.max()])
    lat = N.asarray([lat.min(),lat.max()])
    proj = get_proj((lon, lat), proj='merc')
    p0 = proj(lon[0], lat[0])
    p1 = proj(lon[1], lat[1])
    dx = p1[0]-p0[0]
    dy = p1[1]-p0[1]
    hv = 'hor' if dx>dy else 'ver'

    # Axes rectangle
    if axes_rect is None:
        left = cfgget('axes_rect_left')
        top = cfgget('axes_rect_top')
        right = cfgget('axes_rect_right_'+hv, 'maps', cfgget('axes_rect_right'))
        bottom = cfgget('axes_rect_bottom_'+hv, 'maps', cfgget('axes_rect_bottom'))
    axes_rect = [left, bottom, right-left, top-bottom]

    # Figure size
    W,H = axes_rect[2:]
    alpha = W/H*dy/dx
    figsize = xfigsize, xfigsize*alpha


    return figsize, axes_rect, hv


def plot_map(data, levels=None, units=None, long_name=None,
    title=None, axes_rect=None, xfigsize=None, extend=None,
    savefig=None, copyright=None, logos=None, close=True,  **kwargs):
    """Plot data on a map

    :Params:

        - **data**: Single 2D variable or tuple of two variables for quiver.
        - **units**, optional: Change data units displayed in the colorar label.
        - **long_name**, optional: Name dispay in the default title.
        - **title**, optional: Title as a string or list of strings.
          It is dsipayed using :func:`~vacumm.report.common.add_title`.
        - **copyright**, optional: Add a copyright (multi-)line to the plot
          using :func:`~vacumm.report.common.add_copyright`.
        - **logos**, optional: Add logos to the plot
          using :func:`~vacumm.report.common.add_logos`.
        - **savefig**, optional: Save the figure to  this file.
        - Extra keywords are passed to the :meth:`~vacumm.misc.plot.map2`
          method.

    """

    # Domain
    quiver = isinstance(data, tuple)
    if quiver:
        scalar =  data[0].clone()
        scalar[:] = MV2.sqrt(data[0]**2+data[1]**2)
    else:
        scalar = data
    lon = data.getLongitude().getValue()
    lat = data.getLatitude().getValue()

    # Create figure
    set_mpl_defaults()
    figsize, axes_rect, hv = get_plot_specs(lon, lat,
        axes_rect=axes_rect, xfigsize=xfigsize)
    cb_horiz = hv == 'hor'
    xy = 'x' if cb_horiz else 'y'
    P.figure(figsize=figsize)

    # Add background
    add_background()

    # Add title
    time = scalar.getTime()
    if title is None:
        title = long_name or scalar.long_name
        title = [title[0].upper()+title[1:]]
        ctime = time.asComponentTime()[0]
        title.append(strftime('le %d/%m/%Y %Hh%m', ctime))
    add_title(title, addlegtime=1)

    # Axes for colorbar
    cb_shrink = cfgget('cb_shrink', 'maps')
    if cb_horiz: # Horizontal
        bt = cfgget('cax_rect_bottomtop_hor', 'maps')
        base_width = axes_rect[2]
        width = cb_shrink * base_width
        left = axes_rect[0] + 0.5*(1-cb_shrink)*base_width
        cax_rect = [left, bt[0], width, bt[1]-bt[0]]
    else: # Vertical
        lr = cfgget('cax_rect_leftright_ver', 'maps')
        base_height = axes_rect[3]
        height = cb_shrink * base_height
        bottom = axes_rect[1] + 0.5*(1-cb_shrink)*base_height
        cax_rect = [lr[0], bottom, lr[1]-lr[0], height]
    cax_margins = cfgget('cax_margins_'+hv,'maps')
    if quiver: # More space for quiverkey
        cax_margins[1] += cfgget('cax_margins_quiverkey', 'maps')
    caxb, cax = nested_axes(cax_rect, cax_margins)
    caxb.set_xticks([])
    caxb.set_yticks([])
    cb_orientation = 'horizontal' if cb_horiz else "vertical"

    # No time
    if time:
        if quiver:
            data = data[0][0], data[1][0], data[2][0]
        else:
            data = data[0]

    # Main plot
    grid_color = cfgget('grid_color')
    grid_lw = cfgget('grid_linewidth')
    grid_dashes = [5,5]
    m = map2(data, axes_rect=axes_rect, title=False, show=False,
        cmap='vacumm_previmer2', colorbar_cax=cax, autoresize=0, proj='merc',
        extend=extend, levels=levels,  basemap_anchor='N', land_color='w',
        drawcoastlines_linewidth=1, bgcolor=cfgget('bgcolor', 'maps'),
        colorbar_format="%g", linewidths=0.5,
        drawmeridians_dashes=grid_dashes, drawparallels_dashes=grid_dashes,
        drawmeridians_color=grid_color, drawparallels_color=grid_color,
        drawmeridians_linewidth=grid_lw, drawparallels_linewidth=grid_lw,
        **kwargs)

    # Tune colorbar
    cb = m.get_colorbar()
    cb.set_label('')
    pad = 5
    cax.tick_params(xy, labelsize=8, length=0, pad=pad)
    add_shadow(cb.patch, ax=cax, width=3, xoffset=3, yoffset=-3)
    cb.ax.artists.remove(cb.outline)
    if extend:
        args = []
        if hv=="hor":
            if cb._extend_lower():
                args.append([0,0,dict(dx=0, dy=-pad, ha='left', va='top')])
            if cb._extend_upper():
                args.append([1,0,dict(dx=0, dy=-pad, ha='right', va='top')])
        else:
            if cb._extend_lower():
                args.append([1,0,dict(dx=pad, dy=0, ha='left', va='bottom')])
            if cb._extend_upper():
                args.append([1,1,dict(dx=pad, dy=0, ha='left', va='top')])
        for arg in args:
            add_axlabel(arg[0], arg[1], "$\infty$", ax=cax, size=8,
                color='k', family=None, **arg[2])
    units = units or scalar.units
    add_axlabel(1, 0, u"Unités : %s"%units, dy=-2, va='top', ha='center',
        ax=cax, size=cfgget('units_size', 'maps'))

    # Statistics
    if hv=='hor':
        axref = cax
        dy = -8
    else:
        axref = m.axes
        dy = -18
    vmean, vmax, vmin = scalar.mean(), scalar.max(), scalar.min()
    units = scalar.units
    text = "Moy. : %(vmean).2f %(units)s - Min. : %(vmax).2f %(units)s / Max. : %(vmin).2f %(units)s"%locals()
    set_xlabel(text, ax=axref, labelpad=abs(dy))

    # Logos
    if logos: add_logos(logos)

    # Copyright
    if copyright: add_copyright(copyright)

    # Save
    if savefig:
        m.savefig(savefig)
    if close:
        m.close()

    return m




#TOTO:
#    MAPS:
#        COLORBAR:
#            CENTRER COLORBAR SUR AXE
#            UNITS
#            (QUIVER)
#            IFINITY MAX
#        ESPACE TITRE
#
#    CURVES:
#
#        TOUT
