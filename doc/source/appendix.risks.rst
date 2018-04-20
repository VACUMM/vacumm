.. _appendix.risks:

Unstable routines
*****************

This section lists the modules, functions, classes and methods 
whose stability is little or not proven.
The user is cautioned as to their uses.

.. warning::

    This list is probably not exhaustive or updated because 
    there is no direct link between the technical development of a code, 
    and the update of this part of the documentation. 
    It is likely that you can glean some information about the stability 
    of the code in its own documentation.


.. list-table:: *List of instabilities*
    :widths: 20 25 15 10
    :header-rows: 1
    :stub-columns: 1

    * - Routine
      - Statut
      - TODO
      - Priorité (1>2>...)
    * - :mod:`~vacumm.data.cf`
      - Beta version
      - Check post-processing of specifications
      - 1
    * - :class:`~vacumm.data.misc.dataset.Dataset`
      - Beta version
      - Test the behavior of getting variables at different grid horizontal positions.
      - 1
    * - :meth:`~vacumm.data.misc.dataset.Dataset.get_dx` et al
      - Beta version
      - Must be tester
      - 1
    * - :meth:`~vacumm.data.misc.dataset.OceanDataset.get_depth`
      - Beta version
      - Test the behavior of getting depth at different grid vertical positions
      - 1
    * - :meth:`~vacumm.misc.grid.regridding.regrid2d` and :meth:`~vacumm.misc.grid.regridding.CDATRegridder`
      - Beta version
      - Test it on all kind of grids (values and mask)
      - 1
    * - :func:`~vacumm.misc.grid.regridding.regular`
      - Pas utilisée depuis longtemps
      - Needs more tests
      - 2
    * - :func:`~vacumm.misc.grid.regridding.regular_fill1d`
      - Pas utilisée depuis longtemps
      - Needs more tests
      - 2
    * - :func:`~vacumm.misc.grid.masking.check_poly_islands`
      - Not really used 
      - Need tests or rewriting
      - 3
    * - :func:`~vacumm.misc.grid.masking.check_poly_straits`
      - Not really used 
      - Need tests or rewriting
      - 3
    * - :func:`~vacumm.misc.filters.running_average`
      - Obsolete
      - À virer ?
      - 3
    * - :func:`~vacumm.misc.plot.add_colorbar`
      - Little used
      - To rewrite
      - 3
