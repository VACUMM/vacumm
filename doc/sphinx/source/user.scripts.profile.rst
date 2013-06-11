.. _user.scripts.profile:

:program:`profile.py`
=====================

Ce script permet l'affichage d'informations sur un fichier de profils (histogramme, carte, distribution, profil)

Ce script permet également de convertir un fichier de profil dans le format attendu, voir plus bas.

Préalables
----------

Les profils pris en charge doivent respecter les conditions suivantes:
    
    - avoir les dimensions ou variables représentant le code plateforme, le temps, la position et la profondeur.
    - les variables à traiter doivent être représentées sur deux dimensions, la première
      correspondant au profil, la seconde à la profondeur.

Fonctionnement
--------------

Ce script permet d'afficher:
    
    - les informations générales sur une liste de profils.
    - le tracé d'un ou plusieurs profils.
    - le tracé de l'histogramme des profils.
    - le tracé d'une carte affichant les positions des profils.
    - le tracé d'une carte affichant la distribution des profils (nombre par polygone).

Utilisation
-----------

Les options --time et --variables permettent respectivement la sélection d'une période de temps et les variables (par alias) à considérer.

L'attribut bbox est spécifique au tracé de carte, il définit l'emprise de la zone visualisée et n'agit pas sur la sélection des donnée.

La restriction de zone peut être réalisée par un ou plusieurs polygones.
Dans le cas du tracé de carte:
    
    - si aucun polygone n'est fourni, la carte affiche simplement la position des profils
    - si un ou plusieurs polygones sont fourni, la carte affiche ceux ci avec la distribution des profils de chacun d'entre eux représenté par une couleur.

Définition des polygones
~~~~~~~~~~~~~~~~~~~~~~~~

Un polygone est défini par un ensemble de cordonnées x,y (longitude,latitude).

L'option -p permet de spécifier un polygone sur ligne de coommande
Exemple:

.. code-block:: none
    
    -p x1,y1,x2,y2,x3,y3,...


L'option --polyfile permet de spécifier un fichier de polygone

.. code-block:: none
    
    -3 50 -6 49
    -6 48.5 -3 48.5

L'option --polysfile permet de spécifier un fichier de polygones avec un polygone par ligne

Exemple:

.. code-block:: none
    
    #-3 50 -6 49 -6 48 -3 47
    -3 50 -6 49 -6 48.5 -3 48.5
    -3 48.5 -6 48.5 -6 48 -3 47

Les lignes vides ou commencant par le caractère '#' sont ignorées.

Les coordonées de polygone sont délimitées par des espaces ou virgules ','

La partie décimale des coordonées est séparée de la partie entière par un point '.'

Toutes ces options de polygones peuvent être répétées.

Conversion d'un fichier de profils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

profiles.py permet aussi de mettre en forme un fichier de profil dans un format attendu pour les tracés.

Les formats d'entrée et de sortie sont les même que ceux de merge_profiles.py.

Cependant, les fonctionnalités de cette conversion sont moindres que celles de merge_profiles.py.

L'exemple suivant montre comment convertir un fichier, en specifiant:
    - que les variables en sortie sont temperature et salinity, à partir de TEMP et PSAL
    - force a utiliser la variable DEPH comme etant la variable de profondeur

.. code-block:: none
    
    profile.py input.nc -v temperature,TEMP -v salinity,PSAL -v depth,DEPH --convert output.nc
    

TODO: permettre de forcer l'utilisation d'une variable de pression pour le calcul des profondeurs

Usage
~~~~~

.. code-block:: none
    
    Usage: profile.py [options] [file1] [file2] [...]
    
    Show profiles informations
    
    Examples:
      This will produce a map with the profiles position for the given variable:
        profile.py profiles.nc -v temperature --dist 
      This will produce a map with the profiles distribution for the given variable and set of polygons:
        profile.py profiles.nc -v temperature -t 2000,2005 --hist --dist -P data/polygon_manga/polygon_manga.txt
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      --cfgfile=CFGFILE     Configuration file
      --hist=n,units        Plot profiles histogram with interval n,units (ex:
                            "1,month" , "2,weeks" , ...)
      --dist                Plot profiles distribution (position and, if polygons
                            are specified, the colored count of profiles per
                            polygon)
      --drawpro             Plot profiles position in distribution plot
      --pro=PRO             Plot a profile at specified comma separated indexes.
                            Ranges may also be specified using start:stop[:step]
                            notation, resulting in [start,...,stop[ indexes.
                            Example: 0,1,2:5,5:11:1 would plot profiles 0 to 10
      -t TIME, --time=TIME  Restrict timerange to be processed, format is
                            "datemin,datemax" with date format: "YYYY-mm-
                            ddTHH:MM:SS"
      -b BBOX, --bbox=BBOX  Set the map bounding box, format is
                            lon_min,lat_min,lon_max,lat_max
      -v VARIABLES, --variables=VARIABLES
                            Variables to be processed. Ex: -v temperature.
                            Aliases may be specified using comma separated names
                            Ex: -v temperature,TEMP
                            Repeat this option for each required variable.
      -p POLY, --poly=POLY  Add a polygon for profile restriction and accounting,
                            polygon format is x1,y1,x2,y2,x3,y3,x4,y4... At least
                            4 couples must be set, otherwise if 4 values are set,
                            it is considered as lonmin,latmin,lonmax,latmax
      --polyfile=POLYFILE   Specify a polygon file containing couples coordinates
                            (in lat, lon order), with values separated by a
                            whitespace or coma
      -P POLYSFILE, --polysfile=POLYSFILE
                            Specify a polygons file (one per line)
      -o pattern, --output=pattern
                            Output files pattern  (default: %(plot)s.png)
      --convert=converted_profiles.nc
                            Convert given profile file(s) into a standard netcdf
                            output file. Variables must be specified using the -v
                            option.
      --show                Show figures
      --nonetdcf3           Do not force saving converted file in NetCDF version 3
      --debug               Debug mode
    

**Exemples d'utilisation:**


    * Production d'une carte avec le positionnement des profils:
      
    .. code-block: bash
        
        profile.py profiles.nc -v temperature --dist 
    
    * Production d'une carte avec la distribution des profils dans des polygones:
    
    .. code-block: bash
        
        profile.py profiles.nc -v temperature -t 2000,2005 --hist --dist -P data/polygon_manga/polygon_manga.txt
    

Aperçu des sorties
~~~~~~~~~~~~~~~~~~

**Profiles:**

.. image:: profiles.png
    :width: 90%

**Histogramme:**

.. image:: profiles_hist.png
    :width: 90%

**Distribution (position):**

.. image:: profiles_dist.png
    :width: 90%

**Distribution (polygones):**

.. image:: profiles_dist_polygons.png
    :width: 90%


