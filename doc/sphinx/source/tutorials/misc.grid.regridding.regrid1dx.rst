.. _user.tut.misc.grid.regridding.regrid1dx:

Regrillage 1D étendu
====================

Voir : :func:`~vacumm.misc.grid.regridding.regrid1d`  :func:`~vacumm.misc.grid.regridding.interp1d` :func:`~vacumm.misc.grid.regridding.nearest1d' :func:`~vacumm.misc.grid.regridding.cubic1d`, et le tutoriel ":ref:`user.tut.misc.grid.regridding.interp1d`.

.. _fig.misc.grid.regridding.regrid1dx:
.. figure:: ../../../../scripts/tutorials/misc-grid-regridding-regrid1dx.png
    
    Interpolation linéaire sur selon un axe (vertical) dont les coordonnées varient selon 1 ou plusieurs des autres axes (zonal).
    On fait alors appel aux mots clés ``xmap`` désignant les axes sur lesquels les coordonnées varient, 
    et ``xmapper`` pour spécifier ces coordonnées.

.. literalinclude:: ../../../../scripts/tutorials/misc.grid.regridding.regrid1dx.py
    :encoding: utf-8
