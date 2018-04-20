.. _user.scripts.showres:

:program:`showres.py` -- Display the resolution of a netcdf grid
================================================================

.. include:: bin/showres.help.txt

:Examples:
    
    These examples are executed from the :file:`bin` directory of the distribution.
    
    .. code-block:: bash
    
        $> showres.py ../data/mars3d.xy.nc
        dx=0.0108562°, dy=0.0071907°
        $> showres.py --meters ../data/mars3d.xy.nc
        dx=0.801031km, dy=0.799714km

:See also: :ref:`user.scripts.showgrid`, :func:`~vacumm.misc.grid.misc.resol`, 
           :func:`~vacumm.misc.grid.misc.get_xy`.
