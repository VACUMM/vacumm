.. _user.tut.misc.grid.regridding.remap_vs_interp_1d:

Remapping versus interpolation - 1D
===================================

Voir : :func:`~vacumm.misc.grid.regridding.regrid1d`  :func:`~vacumm.misc.grid.regridding.interp1d` :func:`~vacumm.misc.grid.regridding.remap1d`.

.. _fig.misc.grid.regridding.remap_vs_interp_1d:
.. figure:: ../../../../scripts/tutorials/misc-grid-regridding-remap_vs_interp_1d.png
    
    Un champ avec des valeurs manquantes sur un axes de temps basse résolution est interpolé vers un axe à plus haute résolution, par *remapping* et par interpolation linéaire.

.. literalinclude:: ../../../../scripts/tutorials/misc.grid.regridding.remap_vs_interp_1d.py
