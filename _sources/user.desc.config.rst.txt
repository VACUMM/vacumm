.. _user.desc.config:

Using advanced configurations
*****************************

Initial configurations
======================

Certains modules de la librairie utilisent des fichiers de configuration
pour définir leur comportement par défaut.

Il s'agit de fichiers proches du format **".ini"**.
Les fichiers .ini définissent des sections et variables associées mais
aussi le type de chacune de ces variables.
(cf :class:`ConfigObj`).

Ces fichiers se situent par défaut dans le même répertoire que le module
concerné et portent le même nom, avec l'extension ".py" remplacée par ".ini".

La classe :class:`vacumm.misc.bases.Object` a, entre autre, pour but d'automatiser
le chargement des fichiers de configuration par défaut ainsi que des fichiers
utilisateur (voir plus bas).

Chaque classe d'un module qui hérite de la classe :class:`~vacumm.misc.bases.Object` chargera automatiquement
sa configuration intiale si le fichier existe, et qu'une section correspondante
(par défaut de même nom que la classe) est présente.

Le chargement de cette configuration initiale est réalisé une seule fois pour
chaque classe, **lors du premier import du module concerné**.

La configuration initiale, aussi nommée configuration par défaut, est appliquée
à la classe concernée, et est recopiée par ses objets (instances de cette classe)
lors de leur création.


User configuration
==================

Les configurations utlisateur **".cfg"** permettent d'écraser les valeurs
de la configuration initiale.

La classe :class:`~vacumm.misc.bases.Object` permet d'automatiser leur chargement.

Les script de la librairie proposent généralement l'option :option:`--cfgfile`
afin de charger un fichier de configuration.

Le fichier de configuration peut alors contenir différentes sections
correspondantes au composants (classes) utilisées par le script.


Classes of configuration
========================

.. note::

    Les classes concernée sont celles héritant de Object.

    L'héritage est respécté lors des chargements configurations:
    une classe MyObj héritant de Object pourra contenir les mêmes
    configurations que la classe Object, une classe MySubObj pourra
    contenir les mêmes configurations que MyObj et Object, ...

    Ainsi, la configuration de la classe Object sera aussi appliquée
    aux classes qui la spécialisent (qui en héritent).


La classe :class:`~vacumm.misc.bases.Object`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exemple:

.. code-block:: none
    
    [Object]
        cfg_debug = False
        log_obj_stats = True
        [[Logger]]
            level = 'verbose'

Variables spécifiques au objets:

  - **cfg_debug**: Si True, des logs concernant le chargement des configurations seront émis.
  - **log_obj_stats**: Si True, les valeur min, max, average et count (nombre de valeurs valides)
                        seront incluses lors des appel à Object.describe(myvar). Une attention
                        particulière doit être apportée à cette variable car elle peut avoir un
                        impact sur la durée des traitements.

Variables spécifiques au logging:

    - **level**: Niveau de log (valeurs possibles: debug,verbose,info,notice,warning,error,critical)


La classe :class:`~vacumm.data.misc.dataset.Catalog`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cette classe est utilisée pour peupler les objets Dataset.

Elle permet de spécifier la vue logique des fichiers de données afin de pouvoir traiter
des aggrégations de fichiers, actuellement:
    
    - jeux de données en best estimate (aka best time serie)
    - jeux de données de profiles

Ces deux types de jeux de données sont actuellement configurable de la même manière
via les variables:
    
    - **files**: Une liste explicite de fichiers.
    
    - **filepattern**: Le motif des fichiers à prendre en compte, pouvant contenir
                       un ou plusieurs marqueurs de date (%Y: année sur 4 chiffre, %m: mois, ...).
    - **time**: Une liste de deux chaines de caractères, spécifiant la date
                minimale et maximale **des fichiers (et non de la couverture temporelle des données
                qu'ils contiennent)** à prendre en compte.
                Les dates doivent être dans un format compatible avec le module
                :mod:`vacumm.misc.atime`

Le jeu de données est alors constitué de la liste des fichier **files** à laquelle s'ajoute
la liste des fichiers établie par la fonction :func:`~vacumm.misc.io.list_forecast_files`
prennant en argument (**filepattern**, **time**)

.. todo:: documentation du format des dates de la librairie vacumm

Exemple
-------

.. code-block:: none

    [[Catalog]]
        files = 'file_2004-01.nc','file_2004-02.nc'
        filepattern = 'data/file_%Y-%m.nc'
        time = '2004-01-01','2004-03-01'

La classe :class:`~vacumm.data.misc.dataset.Dataset`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

La classe Dataset peut contenir une section Catalog afin de décrire
l'accès aux données.

Dataset exploite la liste des fichier établie par Catalog en utilisant :class:`vacumm.misc.io.NcReadBestEstimate`.

Cette classe fourni les fonctionnalités génériques d'exploitation d'un jeu de données:
    - exploitation d'une liste de fichiers de donneés
    - récupération de grille
    - détermination de la résolution temporelle
    - lecture des variables en best estimate
    - lecture de l'axe de temps en best estimate
    - lecture transparente des données de profondeur (directe ou coordonnées sigma)
    - chargement de section, hovmoller, couche à une certaine profondeur, ...
    

Bien qu'orientée vers l'exploitation de données grillées, Dataset doit pouvoir constituer une base de lecteur
de données. La classe :class:`vacumm.data.misc.profiles.ProfilesDataset` redéfinie certaines méthodes afin
de traiter non plus des données grillées mais de profiles.


Exemple
-------

.. code-block:: none
    
    [MARS3D]
        log_obj_stats = True
        [[Logger]]
            level = 'debug'
        [[Catalog]]
            filepattern = 'data/champs_%m_BOBI.nc'
            time = '2004-01-01','2004-03-01'
    

La classe :class:`~vacumm.data.misc.profile.ProfilesDataset`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


La classe :class:`~vacumm.data.misc.coloc.Colocator`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

La classe Colocator permet de colocaliser des données de différents types.
    - Données grillées de modèle
    - Données de profiles in situ
    
Colocator fourni actuellement une méthode de lecture des données de modèle colocalisées
en fonction des données de profiles avec deux méthodes:
    - Lecture des données modèle les plus proche en temps et position
    - Interpolation spatiale et temporelle


Example of a configuration file
===============================

Voici un exemple de fichier de configuration minimale premettant de régler le niveau
de logging et l'accès aux données.

Ce type de fichier peut être utilisés avec les différents scripts de diagnostiques
du modèle Mars.

.. todo::
    - Etendre cette configuration avec des sections spécifiques à chacun des scripts
      pour affiner les paramétrage (paramètres actuellement obligatoire en arguments des scripts,
      critères pour le calcul de MLD, paramétrage des tracés, ...)


.. code-block:: none
    
    [Object]
        cfg_debug = False
        log_obj_stats = False
        [[Logger]]
            level = 'info'
    
    [MARS3D]
        log_obj_stats = True
        [[Logger]]
            level = 'debug'
        [[Catalog]]
            filepattern = 'data/mars_%m_%Y.nc'
            time = '2007-01-01','2007-03-01'
    
    [Satellite]
        [[Catalog]]
            files = 'data/satellite.nc'
    
    [ProfilesDataset]
        [[Catalog]]
            files = 'data/merged_profiles_manga_coriolis.nc','data/merged_profiles_manga_sismer.nc'
    

