.. _user.scripts.showbatsec:

:program:`showbatsec.py` -- Plot section of bathymetry
======================================================

.. include:: bin/showbatsec.help.txt

:Examples:   
    
    .. code-block:: bash
    
        $> showbatsec.py  -o showbatsec1.png
        Please set longitude of first point: -7
        Please set latitude of first point: 43
        Please set longitude of second point: -1
        Please set latitude of second point: 47

    .. figure:: showbatsec1.png

        Figure created using the first example.
            
    .. code-block:: bash
    
        $> showbatsec.py  --x0=-2 --y0=44 --x1=-4 --y1=48 \
        --along=m --title="Section across shelf" -o showbatsec2.png
        
    .. figure:: showbatsec2.png

        Figure created using the second example.
        
:See also: :class:`~vacumm.bathy.bathy.NcGriddedBathy`, :func:`~vacumm.misc.plot.curve2`, :func:`~vacumm.misc.plot.minimap`, :func:`~vacumm.misc.grid.misc.transect_specs`.
