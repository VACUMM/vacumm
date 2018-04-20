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
from vacumm.misc.core_plot import hlitvs

class Day12hFormatter(DateFormatter):
    """Day formatter at 12h only"""
    def __call__(self, x, pos=0):
        if x==0:
            raise ValueError('DateFormatter found a value of x=0, which is an illegal date.  This usually occurs because you have not informed the axis that it is plotting dates, eg with ax.xaxis_date()')
        dt = P.num2date(x, self.tz)
        if dt.hour!=12: return ''
        return self.strftime(dt, self.fmt)



def plot_curves(data,
    units=None, long_name=None, vmin=None, vmax=None,
    title=None, axes_rect=None, figsize=None,
    date_locator='auto', date_minor_locator='auto',
    date_formatter=None, hlitvs_units='ticks',
    savefig=None, copyright=None, logos=None, close=True, **kwargs):
    """Plot one or several curves

    :Params:

        - **data**: Single 1D variable or a list of variables or tuples.
          Each tuple has its first item as a variable, and second item
          as plot line specification such as 'r-' for red solid line.

          .. warning:: All variables must have a proper long_name attribute
            that will be displayed in the legend if there are more
            than one variable.

        - **units**, optional: Change data units displayed in the left label.
        - **long_name**, optional: Name dispay in the left label.
        - **title**, optional: Title as a string or list of strings.
          It is dsipayed using :func:`~vacumm.report.common.add_title`.
        - **copyright**, optional: Add a copyright (multi-)line to the plot
          using :func:`~vacumm.report.common.add_copyright`.
        - **logos**, optional: Add logos to the plot
          using :func:`~vacumm.report.common.add_logos`.
        - **savefig**, optional: Save the figure to  this file.
        - Extra keywords are passed to the :meth:`~vacumm.misc.core_plot.Plot.post_plot`
          method.

    :Examples:

        >>> plot_curve(sst)
        >>> plot_curve([sst_insitu, (sst_model, 'ob')])
    """
    # Convert to list of tuples
    if not isinstance(data, list):
        data = [data]
    lsc = []
    lines = cfgget('lines', 'curves')
    for i in xrange(len(data)):
        if not isinstance(data[i], tuple):
            for sc in lines:
                if sc not in lsc:
                    data[i] = data[i], sc
                    break
            else:
                raise PREVIMERError("Can't find a suitable line color or symbol")
        lsc.append(data[i][1])

    # Create the figure
    set_mpl_defaults()
    P.figure(figsize=figsize or cfgget('figsize'))

    # Add background
    add_background()

    # Add title
    if title is None:
        title = scalar.long_name
        title = [title[0].upper()+title[1:]]
    add_title(title, addlegtime=2)

    # Plot curves
    axes_rect = [cfgget('axes_rect_'+pos, 'curves') for pos in ['left', 'bottom', 'right', 'top']]
    axes_rect[2] -= axes_rect[0]
    axes_rect[3] -= axes_rect[1]
    ms = cfgget('markersize', 'curves')
    for var, lsc in data:
        p = curve2(var, lsc, axes_rect=axes_rect, title=False, markersize=ms,
            shadow=True, post_plot=False)

    # Logos
    if logos: add_logos(logos)

    # Copyright
    if copyright: add_copyright(copyright)

    # Finalize
    units = units or data[0][0].units
    long_name = long_name or data[0][0].long_name
    p.post_plot(show=False, xlabel='Date', ylabel='%(long_name)s (%(units)s)'%locals(),
        date_locator=date_locator, date_formatter = date_formatter,
        minor_locator=date_minor_locator,
        savefig=savefig, vmin=vmin, vmax=vmax,
        legend=len(data)>1, legend_loc='lower right', legend_bbox_to_anchor=[1, 1],
        legend_frameon=False, legend_ncol=len(data),
        hlitvs=True, hlitvs_color=cfgget('hlitvs_color', 'curves'),
        hlitvs_units=hlitvs_units, hlitvs_axis=p.axes.xaxis,
        grid_linewidth=cfgget('grid_linewidth'),
        close=close,
        **kwargs
    )
    p.axes.tick_params('y', labelright=True)
    format_label(p.axes.xaxis.label)
    format_label(p.axes.yaxis.label)


    return p
