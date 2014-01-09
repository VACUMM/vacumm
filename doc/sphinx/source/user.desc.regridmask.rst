Regridding and masking
**********************

     
.. _user.desc.regridding:

Regridding
==========


1D and 2D regridding
--------------------

Generalities
^^^^^^^^^^^^

See also: :mod:`~vacumm.misc.grid.regridding`

Regridding classes
~~~~~~~~~~~~~~~~~~

There are two major classes of regridding:

    **The interpolation**
        Elle peut être typiquement par plus proche voisin, linéaire ou cubique, ou faire appel à des
        techniques plus sophistiquées en 2D.
        L'interpolation est à utiliser dans deux cas de figure :
        
            Pour passer **d'une grille (ou axe en 1D) basse résolution vers une grille haute résolution**.
            Si elle est utilisée dans le sens opposé, elle peut être biaisée par de l'aliasing.
            
    **Le remapping**
        Cette technique considère les points de grille comme des cellules : 
        pour chaque cellule cible, les cellules source recouvrant cette dernière
        sont moyénées avec une pondération proportionnelle aux aires de recouvrement.
        Elle est à utiliser pour passer **d'une grille (ou axe en 1D) 
        haute résolution vers une grille basse résolution**.

Le tutoriel ":ref:`user.tut.misc.grid.regridding.remap_vs_interp_1d`" donne un aperçu du comportement des deux classes de regridding.

Le choix de la méthode peut être fait automatiquement.
Dans ce cas, la fonction :func:`~vacumm.misc.grid.regridding.regrid_method` va calculer les résolutions
moyennes en entrée et en sortie, et choisir en fonction du rapport entre les deux.

Le problème des valeurs manquantes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ce problème est résolu en interpolant les masques eux-mêmes et en utilisant des méthodes d'ordre moins élévé.
La méthode d'ordre le moins élévée, et donc celle de référence, est celle avec pondération par la distance (ou par plus proche voisin).
Par exemple, près d'une valeur manquante, l'interpolation linéaire est approximée par l'interpolation avec pondération par la distance dans la cellule.
De même, celle cubique est approximée par celle linéaire, puis celle avec pondération par la distance.

Dans le cas du remapping, les valeurs non masquées sont celles où le masque regrillé à une valeur 
par exemple inférieure à .5.


Regrillage 1D
^^^^^^^^^^^^^

Le regridding 1D se fait sur l'un des axes d'une variable qui peut être multi-dimensionnelle.
Par exemple, on peut chercher à regriller sur la verticale une variable 'tzyx'.

Le regridding 1D supporte en outre ce que l'on pourrait appeler les "axes étendus" : 
un regridding toujours 1D, mais variable suivant certains des autres axes.
On utilise pour cela les mots clés ``xmap`` et ``xmapper``.

Pour plus d'information, voir :func:`~vacumm.misc.grid.regridding.regrid1d`.


.. _user.desc.regridding.regrid2d:

Regrillage 2D
^^^^^^^^^^^^^

Le regridding 2D se fait sur les deux derniers axes, qui sont généralement
associés aux axes de longitude et latitude en degrés.

L'utilisation du proche voisin à proximité immédiate de la côte 
(valeur manquante pour la cellule voisine)
pose problème si l'on désire un rendu propre : des "marches d'escalier" apparaissent 
sous le trait de côte.
Il est possible de prévenir cet inconvénient en utilisant une interpolation proportionnelle
à la distance à la place du proche voisin.
Les valeurs et le masque résultants sont en outre plus réalistes et naturels.

Pour plus d'information, voir :func:`~vacumm.misc.grid.regridding.regrid2d`.


Ressources
^^^^^^^^^^

Le module :mod:`~vacumm.misc.grid.regridding` met à disposition des fonctions et classes
afin de faciliter le regridding avec les méthodes évoquées ci-dessus.
Néanmoins, le coeur du regridding est délégué à d'autre modules.

    - Le regridding sur grilles rectangulaires est effectué notamment par les routines
      Fortran du fichier :file:`interp.f90`. Ce fichier est converti en module
      python grâce à `f2py <http://docs.scipy.org/doc/numpy/user/c-info.python-as-glue.html#f2py>`_.
    - Le remapping sur grille rectangulaire est effectué par le module CDAT :mod:`regrid2`
    - Le regridding sur des grilles curvilinéaires fait appel au regrilleur
      `SCRIP <http://climate.lanl.gov/Software/SCRIP>`_ et au module CDAT :mod:`regrid2`.

Interpolation de nuages de points
---------------------------------

Les seules techniques utilisée sont dérivées de la triangulation
et procèdent généralement à une interpolation dans un plan défini par trois points.
les algorithmes utilisés sont les suivants :
    
    - **Cargen** : Développé par l'IFREMER (plus d'infos?), adapté et inclu dans le fichier
      :file:`interp.f90`.
    - **Natgrid** : Dérivé de la librairie `NgMath <http://www.dkrz.de/ngdoc/ng/ngmath/index.html>`_
      et inclu dans CDAT. Elle utilise un algorithme basé sur les voisin naturels.
    - **Csagrid** : Même origine que Natgrid mais par fit de splines.
    
Si des valeurs masquées sont contenues dans les données d'entrée, le masque est lui aussi interpolé.

Sachant que ces méthodes peuvent demander des ressouces importante,
le domaine est découpé par blocs de même taille.

Pour plus d'information, voir :func:`~vacumm.misc.grid.regridding.griddata`.


.. _user.desc.masking:

Masking
=======

Le masquage des nuages de points et des données grillées est généralement effectué
soit par polygones, soit par interpolation de masque.

Polygones
---------

Les polygones sont définis par une séries de coordonnées X/Y,
et prennent la forme d'objet de type :class:`_geoslib.Polygon`.
Le module :mod:`_geoslib` est fourni avec :mod:`~mpl_toolkits.basemap`,
et permet aussi de définir des objets de type :class:`~_geoslib.Line` ou :class:`~_geoslib.Point`.
Il est alors possible d'effectuer des tests d'inclusion (par exemple un point dans un polygone),
et de calculer des intersections (entre lignes et/ou polygones).

Si on définit par exemple un ensemble de polygones à partir d'un trait de côte,
nous pouvons masquer des points (:func:`~vacumm.misc.grid.masking.polygon_select`)
ou des cellules d'une variable grillée (:func:`~vacumm.misc.grid.masking.masked_polygon`).
Dans le cas des cellules, il est possible de masquer un point si le centre est sur la terre,
ou par exemple si plus de 50% de sa surface est de la terre.

Interpolation de masque
-----------------------

Si les données d'entrée possèdent des valeurs manquantes, le masque de 
sortie peut être estimé en interpolant celui d'entrée converti en réels (de 0. et 1.),
et à l'aide d'un seuil (typiquement 0.5).
Cette approche est valable aussi bien pour les nuages de points que les données grillées.
Il convient néanmoins pour les données grillées d'adapter l'interpolation
près des cellules masquées (voir :ref:`user.desc.regridding.regrid2d`)
pour que le masque issu de l'interpolation des données soit moins restrictif
celui issu de l'interpolation du masque d'origine.

    
 