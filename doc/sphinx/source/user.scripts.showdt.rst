.. _user.scripts.showdt:

:program:`showdt.py` -- Display de time step of netcdf time axis
================================================================

.. include:: bin/showdt.help.txt


:Examples:
    
    These examples are executed from the :file:`bin` directory of the distribution.
    
    .. code-block:: bash
    
        $> showdt.py ../data/mars3d.t.nc
        3600 seconds
        $> showdt.py -u h --format="Time step in hours: %(value).1f" ../data/mars3d.t.nc
        Time step in hours: 1.0

