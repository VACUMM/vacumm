.. _user.tut.bathy.bathy.gridded.base:

Les bases d'une bathy grillée
=============================

Voir : :class:`~vacumm.bathy.bathy.GriddedBathy` :func:`~vacumm.bathy.bathy.plot_bathy` :class:`~vacumm.bathy.shorelines.GSHHS` :func:`~vacumm.misc.plot.map` :func:`~vacumm.misc.grid.regridding.regrid2d` :func:`~vacumm.misc.grid.masking.polygon_mask`

.. _fig.bathy.bathy.gridded.base:
.. figure:: ../../../../scripts/tutorials/bathy-bathy-gridded-base.png

    En **a)**, topographie original ; 
    en **b)**, la même topographie, mais masquée sur
    la terre grâce au trait de côte GSHHS fin ;
    en **c)**, topographie masquée puis interpolée
    sur une grille plus fine ;
    en **d)**, la topographie non masquée est d'abord
    interpolée, puis elle est masquée à l'aide du
    traît de côte.
    Cette dernière méthode évite d'interpoler un
    masque grossier et est donc préférable.
    
.. literalinclude:: ../../../../scripts/tutorials/bathy.bathy.gridded.base.py



