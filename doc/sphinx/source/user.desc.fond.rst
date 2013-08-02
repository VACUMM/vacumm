.. _user.desc.basis:

Basis of the librairie
**********************


_user.desc.basis.arch:
    
General architecture
====================

.. _fig.arch:
    
.. figure:: static/modules.*

    Schematic view of (most of) the library modules (in green) 
    and of some required external modules.



La librairie générique :mod:`vacumm.misc` est organisée de la manière suivante :

    - :mod:`~vacumm.misc.misc` : Generic place where to put what cannot be placed elsewhere.
    - :mod:`~vacumm.misc.axes` : Utilities about CDAT axes.
    - :mod:`~vacumm.misc.atime` : Working with time.
    - :mod:`~vacumm.misc.color` : Working with colors and colorbars (:mod:`matplotlib`).
    - :mod:`~vacumm.misc.config` : Managing advanced configuration (:mod:`configobj`).
    - :mod:`~vacumm.misc.core_plot` : Core classes used by :mod:`~vacumm.misc.plot`.
    - :mod:`~vacumm.misc.filters` : 1D et 2D filters.
    - :mod:`~vacumm.misc.plot` : Plot functions.
    - :mod:`~vacumm.misc.stats` : Special statistical functions.
    - :mod:`~vacumm.misc.grid` : About CDAT grids.

        - :mod:`~vacumm.misc.grid.misc` : General functions.
        - :mod:`~vacumm.misc.grid.regridding` : Interpolation and regidding.
        - :mod:`~vacumm.misc.grid.kriging` : Pure kriging.
        - :mod:`~vacumm.misc.grid.masking` : Masking data.
        - :mod:`~vacumm.misc.grid.basemap` : Utilities related to geographic maps 
          (:mod:`~mpl_toolkits.basemap`).

    - :mod:`~vacumm.misc.phys` : About physics.

        - :mod:`~vacumm.misc.phys.constants` : Useful constants.
        - :mod:`~vacumm.misc.phys.units` : Units conversions.
        
    - :mod:`~vacumm.data` : Advanced data use.
    - :mod:`~vacumm.diag` : Advanced diagnostics.
        

.. _user.desc.cdat:
    
CDAT as numeric basis
=====================

Le choix a été fait de prendre `CDAT <http://www2-pcmdi.llnl.gov/cdat>`_ 
comme base pour l'ensemble des développements de la librairie.

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
    
Matplotlib as graphic basis
===========================

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


 
