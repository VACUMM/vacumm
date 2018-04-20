.. _user.scripts.showtime:

:program:`showtime.py` -- Display information of a netcdf time axis
===================================================================

.. include:: bin/showtime.help.txt


:Examples:
    
    These examples are executed from the :file:`bin` directory of the distribution.
    
    .. code-block:: bash
    
        $> showtime.py -v high_tides --ncol=3 ../data/tide.sealevel.BREST.mars.nc
        276 values:
        2006-07-15 06:15:15, 2006-07-15 18:45:15, 2006-07-16 07:15:16
        2006-07-16 19:30:16, 2006-07-17 08:15:17, 2006-07-17 20:30:17
        2006-07-18 09:15:18, 2006-07-18 21:45:18, 2006-07-19 10:45:19
        ...
        $> showtime.py --format='%b %m %Y' --minmax ../data/mars3d.t.nc
        min: 07 Jan 2008
        max: 09 Jan 2008

