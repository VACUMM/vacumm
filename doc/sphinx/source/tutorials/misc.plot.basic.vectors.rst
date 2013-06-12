.. _user.tut.misc.plot.basic.vectors:

Tracé de vecteurs sur carte
===========================

Voir : :func:`~vacumm.misc.plot.map`.

.. _fig.misc.plot.basic.vectors:
.. figure:: ../../../../scripts/tutorials/misc-plot-basic-vectors.png
    
    Un tracé de vecteurs vitesse sur une carte. 
    Le module des vecteurs est automatiquement tracé en fond mais 
    peut être désactivé avec ``nofill=True,contour=False``. 
    Il est aussi possible de tracer un autre champ 
    en fond avec par exemple ``map((u,v,sst))``

.. literalinclude:: ../../../../scripts/tutorials/misc.plot.basic.vectors.py
