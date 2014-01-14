.. _user.scripts.section:

.. note:: This script is sligthly obsolete and must be updated to be able to work with other model outputs, and renamed.


:program:`section.py`
=====================

Il s'agit du script à utiliser pour le tracé du "bourrelet d'eau froide".

Ce script permet la visualisation d'une section 2D (verticale) pas nécessairement méridionnale ou zonale.

La section 2D (profondeur,position) peut être faite sur des données de modèle ou de climatologie,
plus généralement sur des données 4D (tzyx).

Les tracés résultants ont en abscisse la position le long de la section, en ordonnée la profondeur et
les valeurs de la variable en couleur de remplissage.

Une mini carte est également affichée permettant de localiser la section.

La description des données à représenter est définie par le fichier de configuration.
Les options en ligne de commande (voir usage), sont aussi paramétrables dans le fichier de
configuration.

Le script produit un tracé par intervalle de temps et par variable.

Usage
~~~~~

.. code-block:: none
    
    Usage: section.py [options]
    
    Produce vertical section plots
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      --cfgfile=CFGFILE     Configuration file
      -v varname[,alias]*, --variables=varname[,alias]*
                            Variables to be processed. Use aliases when varname
                            differ between datasets. This option may be repeated
                            to produce figures for each variable definition.
      -t min,max,[bb],step,unit, --time=min,max,[bb],step,unit
                            Time selection: - min,max: specify the time range to
                            operate. - bb: optionnal, time bounds open/closed
                            selection. - step,unit: period covered by each plot.
                            Ex:  "2001-01,2001-01-15T00,7,days"
                            "2001-06,2001-09,co,1,month"
      -c x1,y1,x2,y2, --coords=x1,y1,x2,y2
                            Section coordinates
      -b lonmin,latmin,lonmax,latmax, --bbox=lonmin,latmin,lonmax,latmax
                            Restrict processed zone to the specified bounding box
      -o pattern, --output=pattern
                            Output files pattern (default:
                            section-%(var)s-%(tmin)s-%(tmax)s.png)
      --show                Show figures


Aperçu des sorties
~~~~~~~~~~~~~~~~~~

.. image:: section.png
    :width: 90%


