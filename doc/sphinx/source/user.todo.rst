.. _user.todo:


What remains to be done
***********************


.. _user.todo.gen:

Generic list
============

.. todolist::

Stabilisation of the code
=========================

You can refer to :ref:`appendix.risks`  for a list of routines that require special attention, and therefore such improvements.


.. _user.todo.more:

Supplements
===========

This section lists important improvements needed to the library, and possibly the key to achieve or possible ways.



Interpolations and reggridding
------------------------------

- The :mod:`~vacumm.misc.grid.kriging` module must 
  be used by the :mod:`~vacumm.misc.grid.regridding` module, and especially the 
  :class:`~vacumm.misc.grid.regridding.GridData` class.
- The :mod:`~vacumm.misc.grid.regridding.regrid2d` must better handle masks, especially
  in the case of curcilinear grids handled by 
  the :class:`~vacumm.misc.grid.regridding.CDATRegridder` class.
- The :class:`~vacumm.misc.grid.regridding.GriddedMerger` must be checked.


:class:`~vacumm.data.misc.dataset.Dataset` classes
--------------------------------------------------

- A special :class:`GriddedDataset` class must be created.
- Vertical levels treatments must be generalized to be usable for
  :class:`~vacumm.data.misc.dataset.AtmosDataset`, 
  and a :class:`ZDataset` must be created for that.
- A :class:`ZGridded` like class must be implemented to serve as base class
  for :class:`~vacumm.data.misc.dataset.AtmosDataset` and
  :class:`~vacumm.data.misc.dataset.OceanDataset`.
  
:mod:`vacumm.diag` and grid staggering
--------------------------------------

All diagnostics on grids must properly handle grid staggering.
Different schemes can be proposed depending on the grid type and 
the grid point location of variables.
One way is to handle a limited number of configurations,
like for instance for kinetic energy :

- The case of a C grid with u and v and U and V locations recpectively.
- The case of a C grid with u and v and V and U locations recpectively (geostrophic currents).
- The generic case of an A grid when there is no staggering.

If the configuration is not explicitly specified, it must be guessed
using :meth:`vacumm.data.misc.arakawa.ArakawaGrid.factory` for the grid type,
and :func:`vacumm.data.cf.get_loc` for the location.
Another option is to estimate grid staggering from coordinates.


Tidal tools
-----------

These must use the Tidal toolbox thanks to f2py, 
and integrate analysis and prediction capabilities.




