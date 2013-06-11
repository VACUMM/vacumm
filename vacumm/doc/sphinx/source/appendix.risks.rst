.. _appendix.risks:

Unstable routines
*****************

Cette section liste les modules, fonctions, classes et méthodes (*routines*) 
dont la stabilité
n'est peu ou pas éprouvée.
L'utilisateur est donc mise en garde quant à leur utilisations.

.. warning::

    Cette liste n'est probablement pas exhaustive ni à jour car il n'existe pas
    de lien technique direct entre le développement d'un code, et la mise
    à jour de cette partie de la documentation.
    Il est probable que vous puissiez
    glaner quelques informations sur la stabilité du code 
    dans sa propre documentation.


.. list-table:: *Liste des routines ou modules instables*
    :widths: 20 25 15 10
    :header-rows: 1
    :stub-columns: 1

    * - Routine
      - Statut
      - TODO
      - Priorité (1>2>3)
    * - :func:`~vacumm.misc.grid.regridding.regular`
      - Pas utilisée depuis longtemps
      - Besoin de plus de tests
      - 2
    * - :func:`~vacumm.misc.grid.regridding.regular_fill1d`
      - Pas utilisée depuis longtemps
      - Besoin de plus de tests
      - 2
    * - :func:`~vacumm.misc.grid.masking.check_poly_islands`
      - Pas vraiment utilisée 
      - Besoin de tests voire réécriture
      - 3
    * - :func:`~vacumm.misc.grid.masking.check_poly_straits`
      - Pas vraiment utilisée 
      - Besoin de tests voire réécriture
      - 3
    * - :func:`~vacumm.misc.filters.running_average`
      - Obsolète
      - À virer ?
      - 3
    * - :func:`~vacumm.misc.plot.add_colorbar`
      - Peu utilisée
      - À réécrire
      - 3
    * - :mod:`vacumm.misc.archive`
      - Utilisé une seule fois !
      - Test + doc
      - 4
