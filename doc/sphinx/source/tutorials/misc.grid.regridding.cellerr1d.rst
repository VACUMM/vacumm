.. _user.tut.misc.grid.regridding.cellerr1d:

Weigthed interpolation 1D based on errors
=========================================

Voir : :func:`~vacumm.misc.grid.regridding.regrid1d`  :func:`~vacumm.misc.grid.regridding.cellerr1d`.

.. _fig.misc.grid.regridding.cellerr1d:
.. figure:: ../../../../scripts/tutorials/misc-grid-regridding-cellerr1d.png
    
    Interpolation of radial speed to an hourly time axis using weights that
    are based on measurement error and lag errors: a) lag error model, 
    b) result of interpolation with a time zoom.

.. literalinclude:: ../../../../scripts/tutorials/misc.grid.regridding.cellerr1d.py
