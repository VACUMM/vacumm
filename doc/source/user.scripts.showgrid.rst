.. _user.scripts.showgrid:

:program:`showgrid.py` -- Show and display a netcdf grid
========================================================

.. include:: bin/showgrid.help.txt

:Example:   
    
    .. code-block:: bash
    
        $> showgrid.py --plot --figsize=12 -s f --xmin=-5 ../data/mars3d.xy.nc
        Dimensions           : nx=227 ny=261
        Axes                 : lon="longitude" lat="latitude"
        Zonal extent         : -6.56181 -> -4.10829  [06°33'43''W -> 04°06'30''W]
        Meridional extent    : 47.4916 -> 49.3615  [47°29'30''N -> 49°21'41''N]
        Zonal resolution     : 0.0108562° / 0.801031km
        Meridional resolution: 0.0071907° / 0.799714km
        
:See also: :ref:`user.scripts.showres`, :func:`~vacumm.misc.grid.misc.resol`, 
           :func:`~vacumm.misc.grid.misc.get_xy`, :func:`~vacumm.misc.plot.map2`,
           :func:`~vacumm.misc.plot.add_grid`.