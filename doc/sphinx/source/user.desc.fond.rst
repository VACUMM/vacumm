.. _user.desc.fond:

Fondements de la librairie
**************************

    
Architecture générale
=====================

.. _fig.arch:
    
.. figure:: static/modules.*

    Vue schématique des modules de la librairies (en vert), ainsi que de certains
    paquets dont elle dépend.



La librairie générique :mod:`vacumm.misc` est organisée de la manière suivante :

    - :mod:`~vacumm.misc.misc` : Contient ce qui n'est pas inclus dans les
      autres modules.
    - :mod:`~vacumm.misc.axes` : Utilitaires concernant les axes au sens de CDAT.
    - :mod:`~vacumm.misc.atime` : Gestion des différents types de temps.
    - :mod:`~vacumm.misc.color` : Couleurs et palettes de couleurs au sens de matplotlib.
    - :mod:`~vacumm.misc.config` : Gestion de configurations avancées.
    - :mod:`~vacumm.misc.core_plot` : Classes définissant le nouveau cœur des routines de plots de :mod:`~vacumm.misc.plot`.
    - :mod:`~vacumm.misc.filters` : Filtres numériques 1D et 2D.
    - :mod:`~vacumm.misc.plot` : Fonction graphiques.
    - :mod:`~vacumm.misc.stats` : Quelques fonctions statistiques spéciales.
    - :mod:`~vacumm.misc.grid` : À propos de données grillées.

        - :mod:`~vacumm.misc.grid.misc` : Routines générales sur les grilles.
        - :mod:`~vacumm.misc.grid.regridding` : Interpolation et regrillage.
        - :mod:`~vacumm.misc.grid.masking` : Masquage des données.
        - :mod:`~vacumm.misc.grid.basemap` : Objets et projections cartographiques.

    - :mod:`~vacumm.misc.phys` : À propos de données physiques.

        - :mod:`~vacumm.misc.phys.constants` : Quelques constantes utiles.
        - :mod:`~vacumm.misc.phys.units` : Conversions d'unités.
        
        

  .. _user.desc.cdat:
    
CDAT comme base numérique
=========================

Le choix a été fait de prendre `CDAT <http://www2-pcmdi.llnl.gov/cdat>`_ comme base pour l'ensemble
des développements de la librairie.

Les tableaux [:mod:`MV2`]
-------------------------

Les tableaux numériques utilisés sont générés par le module :mod:`MV2` (:mod:`cdms2`) de CDAT,
s'ils représentent une quantité localisée dans l'espace et/ou le temps.
Ces tableaux ne sont pas purement numériques et contiennent notamment :

    - les valeurs numériques,
    - le masque associé,
    - les axes (1D ou 2D) qui localisent la variable (longitude, temps...),
    - des attributs (nom, unités...).

L'avantage est de ne pas avoir à passer explicitement et séparément 
aux fonctions toutes ces informations en plus du tableau numérique.
En outre, sont associées à ces tableaux "étendus" un grand nombre de fonctionnalités 
développées par CDAT, adaptées à un cadre océanographique et météorologique.
Pour finir, notons que CDAT est une excellente interface pour la lecture et 
l'écriture de fichiers netcdf. 

En conséquence, lire et tracer une carte de SST se fera 
par exemple de la manière suivante : 

    >>> import cdms2
    >>> f = cdms2.open('file.nc')
    >>> sst = f('sst') # Lecture
    >>> f.close()
    >>> from vacumm.misc.plot import map
    >>> map(sst) # Plot
    
    
Le temps [:mod:`cdtime`]
------------------------

Le temps au sens du module :mod:`cdtime` de CDAT permet de gérer
plusieurs types de calendrier, de considérer des temps absolus (année, mois, etc)
ou relatifs (nombre et unités depuis une date donnée),
et de s'intégrer naturellement dans les variables :mod:`MV2` à 
travers les axes de temps.

Il s'agit ainsi du type de temps utilisé comme base dans la librairie.
Néanmoins, de nombreuses fonctions (voir :mod:`~vacumm.misc.atime`)
permettent de faire des conversions vers d'autres
types de temps (chaîne de caractères, numérique, :class:`datetime.datetime`, etc).



.. _user.desc.mpl:
    
Matplotlib comme base graphique
===============================

La librairie se base presque exclusivement sur `Matplotlib <http://matplotlib.sourceforge.net>`_ et
`Basemap <http://matplotlib.sourceforge.net/basemap/doc/html/>`_ pour l'ensemble des graphiques.
Matplotlib est une librairie permettant d'effectuer des représentations
graphiques de données numériques, dans un environnement proche de celui Matlab.
Les graphiques générés sont d'excellente qualité, 
et sont disponibles potentiellement en de nombreux formats.

Le module :mod:`vacumm.misc.plot` opère une fusion entre CDAT et matplotlib.
Il contient  un ensemble de fonctions prenant comme argument des variables CDAT.

.. note::
    
    Certaines fonctions sont en cours de portage pour utiliser 
    le module :mod:`vacumm.misc.core_plot`.
    Elles ont le même nom que leur consœur mais avec le suffixe "2" 
    (exemple :func:`~vacumm.misc.plot.curve2` est la nouvelle version de 
    :func:`~vacumm.misc.plot.curve`),
    et leur arguments et paramètres restent inchangés.
    :mod:`~vacumm.misc.core_plot` contient des classes spécialisées héritées directement ou indirectement
    de la même classe :class:`~vacumm.misc.core_plot.Plot`.
    Le but est d'éviter la duplication de code,
    d'organiser plus facilement le développement de nouvelles fonctionnalités,
    et de pouvoir intégrer aisément des fonctionnalités graphiques
    dans des infrastructures opérationnelles.


 
