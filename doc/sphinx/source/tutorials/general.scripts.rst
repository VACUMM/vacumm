.. _user.tut.generic.scripts:


How to write a script
=====================

This tutorial shows how to write the same script with different level of coding.

In this example, we want to compute the mixed layer (MLD) depth from MARS3D outputs,
and plot it on a map.

All these scripts produce a figure close to the one below.

.. figure:: ../../../../scripts/tutorials/general-scripts-level1.png

    Figure created with the first level script.
    It shows the MARS3D MLD in the Gulf of Lion.

.. note::

    We wont detail here how the algorithm for computing the MLD.
    Please refer to the :func:`vacumm.diag.thermdyn.mixed_layer_depth`
    for more information.


First level: simple case
-------------------------

In this example, we import all the needed functions for reading file,
computing the MLD and plotting it.
All parameterd are defined from within the script.


.. literalinclude:: ../../../../scripts/tutorials/general.scripts.level1.py



Second level: using a child of :class:`~vacumm.data.misc.dataset.Dataset` class
-------------------------------------------------------------------------------

Most of the operations of the first level script are directly integrated 
in the children of the :class:`~vacumm.data.misc.dataset.Dataset` class,
and in particular in the :class:`~vacumm.data.model.mars.MARS3D` class.
This new formulation of the script takes advantage of this.


.. literalinclude:: ../../../../scripts/tutorials/general.scripts.level2a.py


In this the follwing formulation of the script, the mixed layer depth is also
directly computed by the :class:`~vacumm.data.model.mars.MARS3D`: 
it automatically read the temperature and the depths before computing the MLD.

.. literalinclude:: ../../../../scripts/tutorials/general.scripts.level2b.py


.. _user.tut.generic.scripts.cfgarg:
    
Third level: add commandline arguments and options
--------------------------------------------------

We know want the script to be generic and use it without editing it.
Therefore we add options and arguments management using the :mod:`argparse` module.

.. literalinclude:: ../../../../scripts/tutorials/general.scripts.level3.py



