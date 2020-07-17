.. _tut.misc.grid.masking.coast:

Masque d'après traît de côte
============================

Voir : :func:`~vacumm.misc.grid.masking.polygon_mask` :func:`~vacumm.misc.plot.map` :func:`~vacumm.misc.plot.add_key`

.. _fig.misc.grid.masking.coast:
.. figure:: ../../../../scripts/tutorials/misc-grid-masking-coast.png

    En a), une cellule est masquée si son centre est sur la terre. En b), la cellule n'est masquée que si plus de 50% se sa surface est sur de la terre. En c), comme en b) mais seulement 30% de la surface suffit à masquer la cellule. En d), comme en c) mais si la cellule a plusieurs fois de terre, il faut que la terre couvrant 70% de la cellule pour que celle-ci soit masquée.

.. literalinclude:: ../../../../scripts/tutorials/misc.grid.masking.coast.py

