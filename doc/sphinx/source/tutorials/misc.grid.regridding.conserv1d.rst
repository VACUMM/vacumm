.. _user.tut.misc.grid.regridding.conserv1d:


Remapping conservatif 1D
========================

Le remapping conservatif est très similaire au remapping simple (``method="remap"``), 
et a pour but de conserver la somme du champ d'entrée vers le champ de sortie.
Dans le cas conservatif, on n'effectue pas des moyennes mais des sommes.

Cette méthode est utile pour regriller par exemple des précipitations dans le temps.

Voir : :func:`~vacumm.misc.grid.regridding.regrid1d`.

.. _fig.misc.grid.regridding.conserv1d:
.. figure:: python/misc-grid-regridding-conserv1d.*

    Des précipitations horaires sont réparties sur un axe de temps bi-horaire.

.. literalinclude:: python/misc.grid.regridding.conserv1d.py
