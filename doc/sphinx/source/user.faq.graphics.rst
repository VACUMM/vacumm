.. _user.faq.graphics:

Graphiques
==========

How to modify the plot parameters?
----------------------------------

Edit the Matplotlib configuration file.
It is typically located here: :file:`~/.matplotlib/matplotlibrc`.

For more information: http://matplotlib.sourceforge.net/users/customizing.html


My map takes time to draw
-------------------------

The map layouts that take time are those of a reduced area
which use a high resolution coastline ("f" from GSHHS or "s" for Histolitt from SHOM).
In the case of GSHHS traîts, a card caching system
plotted with mod:`~mpl_toolkits.basemap` is used
(see :mod:`vacumm.misc.grid.basemap`):
if a map using the same zoom and resolution
must be traced, it is loaded from the cache.
For your information, cached maps are stored in the directory given by
the function :func:`~vacumm.misc.grid.basemap.get_map_dir` (a priori
:file:`~/.matplotlib/basemap/cached_maps`).

.. note::

    This cache system can generate errors during use
    of a saved card. See:`user.faq.graphics.err_basemap`.

In addition, a map may want to appear on the screen
while you're doing a remote plot, which can take time.
If you want to disable this display, read :
:ref:`user.faq.graphics.window`.



How to tuner a graph to plot?
-----------------------------

VACUMM's graphical routines transfer keywords to
most of the matplotlib functions used.
We specify a keyword to transfer a prefixing this one of the name of
the function suffixed with a"_" (use: func:`~vacumm.misc.misc.kwfilter` to
filter keywords).

For example, it is possible to change the size of the title of a graph,
and the thickness of the contours as follows::

    >>> map2(sst, title_size=20, contour_linewidths=3)

In this specific case, to have the list of the other possible options
we refer to the documentation of functions :func:`matplotlib.pyplot.title`
and func:`matplotlib.pyplot.contour`.

In addition, Matplotlib plotted objects are often accessible using
the :meth:`~vacumm.misc.core_plot.Plot.get_obj` method that returns a list::

    >>> m = map2((u, v), show=False)
    >>> myquiver = m.get_obj('quiver')[0]
    >>> mycontour = m.get_obj('contour')[0]
    >>> myquiver.set_zorder(10)
    >>> m.show()

Note you must set the ``show`` parameter to ``False`` if your
want to tune your plot (write mode), which is not necessary if you just
want to gather come parameters.

How to place a colorbar yourself?
---------------------------------

Turn off the automatic colorbar plot and figure display,
then draw it by specifying its position (cf. : func:`matplotlib.pyplot.axes`) :

    m = map2(data, colorbar=False, right=.9, show=False)
    m.colorbar(cax=[.91, .2, .3, .6])
    m.show()

or directly

    map2(data, colorbar_cax==[.91, .2, .3, .6], right=.9, show=False)

See: func:`matplotlib.pyplot.colorbar` for plotting options.


.. _user.faq.graphics.window:

My figure is not displayed / turn off the display of figures like in batch mode
-------------------------------------------------------------------------------

The creation of figures requires the use of a graphical software layer.

These software layers are called backends. The minimum for the
saving the figures is to use the backend ``Agg``, this allows
also to bypass errors when using batch scripts
when the display is not available.

The other backends require a display, the ``TkAgg`` backend
is almost systematically available.
The qt4agg backend is also available with UVCDAT.

If you want to make a display of your figures, use ``qt4agg``,
if you use ``Agg``, no display will be made.

There are two ways to modify the backend:

    1) The permanent one consists in modifying the corresponding line
       in the Matplotlib configuration file::

           backend : Agg

    2) The second way is to do it on the fly,
       at the beginning of a script (before any loading of :mod:`matplotlib`)
       and :mod:`vacumm`)::

           from matplotlib import use
           use('Agg')

How to disable automatic zoom?
------------------------------

For example when you draw a map with func:`~vacumm.misc.plot.map2`,
if the unmasked data does not cover the entire domain,
automatic zoom will avoid drawing areas that are not covered.

You can avoid this by using the ``xmasked=False`` options,
symasked=False`, see ``xymasked=False`` to disable zoom on both axes
(in the case of a 2D plot)::

    map(sst, xymasked=False)

How to add shadow or glow effects?
--------------------------------------------------

Most of the time, keywords are provided for this purpose for the
plotting functions. For example, if you draw contours, you
can do this to add a shadow to these and an effect
glow to labels::

    map2(data, contour_shadow=True, clabel_glow=True)

You can pass keywords to modify effects (see
:func:`~vacumm.misc.core_plot.add_shadow` and :func:`~vacumm.misc.core_plot.add_glow`
to know the options)::

    map2(data, contour_shadow_width=4, contour_shadow_xoffset=3)

You can also proceed manually for an external trace::

    m = map2(data, show=False)
    m.add_shadow(m.axes.plot(x,y)) # m.axes.plot or pylab.plot

or::

    from vacumm.misc.core_plot import add_shadow
    add_shadow(P.plot(x,y))




.. _user.faq.graphics.err_basemap:

I have an error with class :class:`~mpl_toolkits.basemap.Basemap`
-----------------------------------------------------------------

It may be that the use of a map (typically when plotting with
:func:`~vacumm.misc.plot.map2`) produces an error (for example a missing attribute,
such as :attr:`celestial` or :attr:`_mapboundarydrawn`).
This is probably related to the automatic caching of already drawn maps.
This process saves time when plotting a map with exactly the same
characteristics (domain, projection, coastline).
However, the cached card may not be compatible with the current version of
:mod:`~mpl_toolkits.basemap`, following a python package update.
To solve the problem radically, you can delete all cached cards.
They are found by default in the directory given by the function
:func:`~vacumm.misc.grid.basemap.get_map_dir`.

There is a function to perform this operation::

    >>> from vacumm.misc.grid.basemap import reset_cache
    >>> reset_cache()
