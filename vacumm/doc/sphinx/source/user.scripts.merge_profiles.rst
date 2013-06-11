.. _user.scripts.merge_profiles:


:program:`merge_profiles.py`
============================

Ce script est destiné à agréger des fichiers de profils verticaux plus ou moins homogènes dans
l'organisation des fichiers NetCDF (convention, nom des variables et attributs).

Initialement les profils exploités sont issus des bases CORIOLIS et SISMER.

Préalables
----------

Les profils pris en charge doivent respecter les conditions suivantes:
    
    - avoir les dimensions ou variables représentant le code plateforme, le temps, la position et la profondeur.
    - les variables à traiter doivent être représentées sur deux dimensions, la première correspondant au profil, la seconde à la profondeur.
    - en option, des variables de code qualité peuvent être utilisées et servir de filtres.

Fonctionnement
--------------

Le script effectue les actions suivantes:
    
    - chargement des fichiers de profils.
    - exclusion des profils ne répondant pas aux critères de qualité des données.
    - exclusion des données masquées.
    - détection des dates partielles (heure manquante).
    - détection des duplicats.
    - homogénéisation du format des profils.
    - enregistrement d'un fichier de profils. uniques et un fichier de profils dupliqués.

Les types de profils sont déterminés automatiquement à partir des noms de fichiers.

Le script fonctionne avec un type générique de profil, qui peut être spécialisé pour les cas particuliers.

Les codes plateforme sont déterminés soit par une variable associée, un attribut global, ou le nom de fichier.

La sélection de variables à traiter est réalisée par alias; une table de correspondance est définie dans le
code source.

Ex: en spécifiant comme variable "temperature", la première recherche dans les fichiers sera sur une variable nommée TEMP.

**Filtrage qualité**

Si une liste de codes qualité valides est renseignée, les profils pour lesquels les informations de
temps, position ou profondeurs ne correspondent pas à l'un des codes sont exclus.

**Détection de dates partielles**

Certains profils peuvent contenir une information partielle en temps, la partie horaire étant indéterminée.
Dans ce cas la partie horaire est fixée à 12H et le script cherchera à remplacer ce profil par un profil avec 
une date complète si disponible, sinon le profil est conservé avec cette date ajustée.

**Détection des duplicats**

Les critères de duplicité de deux profils sont:
    
- codes plateforme égaux.
- dates et heures égales.
- positions égales.

Si tous ces critères sont vérifiés, tous les profils correspondant sont ajoutés au fichier de duplicats, le premier profil est quand à lui ajouté au fichier profils uniques.

**Homogénéisation des profils**

Le script adapte les données des profils pour un format commun en sortie.

Les informations temporelles, issus de variables de type chaîne de caractères, ou numérique avec unité sont converties
en valeurs numérique avec unité (par défaut "seconds since 1900-01-01 00:00:00").

Les informations de profondeurs, issues de variables de profondeur ou de pression sont converties si nécessaire en profondeur.

Utilisation
-----------

Le paramétrage peut être réalisé soit par la ligne de commande, par un fichier de configuration ou les deux.
Par exemple, la configuration :file:`profiles.cfg`:

.. code-block:: ini

    [ProfilesMerger]
        input-files = profile1.nc,profile2.nc
        [[load]]
            variables = temperature,salinity

Ce fichier de configuration :file:`profiles.cfg`: serait utilisé de la manière suivante:

.. code-block: bash

    merge_profiles.py --cfgfile profiles.cfg

Cela équivaut à la ligne de commande:

.. code-block:: sh

    merge_profiles.py --ProfilesMerger-load-variables temperature --ProfilesMerger-load-variables salinity --ProfilesMerger-input-files profile1.nc --ProfilesMerger-input-files profile2.nc

Ou encore:

.. code-block:: sh

    merge_profiles.py --ProfilesMerger-load-variables temperature --ProfilesMerger-load-variables salinity profile1.nc profile2.nc

Pour une utilisation en ligne de commande:

Le choix des variables à traiter peut être spécifié avec l'option --ProfilesMerger-load-variables (option à répeter pour chaque variable)

Le choix du filtre des codes qualité (des profils à conserver) peut être spécifié avec l'option --ProfilesMerger-load-qualities (option à répeter pour chaque code qualité)

Exemple si on veut garder les profils ayant les codes qualité 1 ou 2:

.. code-block:: sh

    merge_profiles.py --cfgfile profiles.cfg --ProfilesMerger-load-qualities 1 --ProfilesMerger-load-qualities 2

Le nom du fichier de profils uniques en sortie est déterminé par --ProfilesMerger-merge-filter-file

Le noom du fichier de profils dupliqués en sortie est déterminé par --ProfilesMerger-merge-reject-file

Le choix des fichiers de profils en entrée à agréger peut être réalisé de trois manières différentes, pouvant être combinées:
    - un fichier par argument de la ligne de commande
    - un fichier par option "--ProfilesMerger-input-files" de la ligne de commande ou de la configuration
    - par recherche de fichier en utilisant les options --ProfilesMerger-find-* de la ligne de commande ou de la configuration

Si la méthode de recherche de fichier est utilisée, deux modes sont disponibles selon l'utilisation de l'une ou l'autre des options/configurations:
    - --ProfilesMerger-find-regex / regex: utilisation d'une expression régulière (ex: '.*PR_(BA|CT).*\.nc')
    - --ProfilesMerger-find-pattern / pattern: utilisation d'une expression wildcard (ex: ('*PR_BA*.nc', '*PR_CT*.nc'))

Usage
~~~~~

.. code-block:: none
    
    Usage: merge_profiles.py [options] [profilesA.nc] [profilesB.nc] [...] 
    
    Merge various vertical profiles into a netcdf file
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      --nonetcdf3           Do not force saving files in NetCDF version 3
      --nonumpywarnings     Hide all numpy warnings
      --debug               Set logging level to DEBUG
      --verbose             Set logging level to VERBOSE
      --loglevel=level      Set logging level, available levels are: NOTSET,
                            DEBUG, VERBOSE, INFO, NOTICE, WARNING, ERROR,
                            CRITICAL. You may restrict this level for specific
                            class with filters:
                            --loglevel=MyBuggedClass=debug,MyVerbose Class=warning
      --logformat=format    Set logging format. Can also use filters.
      --logdateformat=format
                            Set logging date format. Can also use filters.
      --lognocolor          Disable logging colors
      --cfgfile=CFGFILE     Configuration file [default: "none"]
    
      Merger configuration:
        --ProfilesMerger-find-path=PROFILESMERGER_FIND_PATH
                             Root directory. [default: "None"]
        --ProfilesMerger-find-regex=PROFILESMERGER_FIND_REGEX
                             File inclusion regular expression. This is mutually
                            exclusive with the pattern option. [default: "None"]
        --ProfilesMerger-find-pattern=PROFILESMERGER_FIND_PATTERN
                             File inclusion wildcard pattern(s). This is mutually
                            exclusive with the regex option. [default: ""]
        --ProfilesMerger-find-matchall=PROFILESMERGER_FIND_MATCHALL
                             Evaluate whole files paths (True) or only files names
                            (False). [default: "False"]
        --ProfilesMerger-find-depth=PROFILESMERGER_FIND_DEPTH
                             Maximum recursive search depth (None: no limit, 0:
                            direct children, N: N subdirectories). [default: "0"]
        --ProfilesMerger-load-safe=PROFILESMERGER_LOAD_SAFE
                             Use generic profile on unknown profile type.
                            [default: "False"]
        --ProfilesMerger-load-variables=PROFILESMERGER_LOAD_VARIABLES
                             Variable filter. [default: ""]
        --ProfilesMerger-load-qualities=PROFILESMERGER_LOAD_QUALITIES
                             Quality code filter applied to time, position and
                            depth. [default: "1 2"]
        --ProfilesMerger-merge-filter=PROFILESMERGER_MERGE_FILTER
                             Filter class name. [default:
                            "ProfilesDuplicatesFilter"]
        --ProfilesMerger-merge-filter-file=PROFILESMERGER_MERGE_FILTER_FILE
                             Filtered outputs file. [default:
                            "merged_profiles.nc"]
        --ProfilesMerger-merge-filter-sort=PROFILESMERGER_MERGE_FILTER_SORT
                             Sort filtered (currently 'time' or None). [default:
                            "time"]
        --ProfilesMerger-merge-reject-file=PROFILESMERGER_MERGE_REJECT_FILE
                             Rejected outputs file. [default:
                            "duplicated_profiles.nc"]
        --ProfilesMerger-merge-reject-sort=PROFILESMERGER_MERGE_REJECT_SORT
                             Sort filtered (currently 'time' or None). [default:
                            "time"]

Format des sorties:
-------------------

Les sorties sont au format NetCDF.

Pour chaque profil, la variable filename donne le nom de fichier d'origine du profil et la dimension
profile donne l'indice du profil dans ce fichier d'origine.

La structure des fichiers générés suit celle de l'exemple suivant:

.. code-block:: none
    
    netcdf merged_profiles {
    dimensions:
        profile = 826 ;
        filename_string = 18 ;
        platform_code_string = 4 ;
        level = 511 ;
    variables:
        int profile(profile) ;
            profile:description = "Index of profile in original file" ;
        double time(profile) ;
            time:long_name = "Time" ;
            time:standard_name = "time" ;
            time:units = "seconds since 1900-01-01 00:00:00" ;
            time:missing_value = 1.e+20 ;
            time:axis = "T" ;
        int filename_string(filename_string) ;
        char filename(profile, filename_string) ;
            filename:missing_value = "N/A" ;
            filename:description = "File name of the original file" ;
        double latitude(profile) ;
            latitude:standard_name = "latitude" ;
            latitude:long_name = "Latitude" ;
            latitude:valid_min = -90. ;
            latitude:units = "degree_north" ;
            latitude:missing_value = 1.e+20 ;
            latitude:valid_max = 90. ;
            latitude:axis = "Y" ;
        double longitude(profile) ;
            longitude:standard_name = "longitude" ;
            longitude:long_name = "Longitude" ;
            longitude:valid_min = -180. ;
            longitude:units = "degree_east" ;
            longitude:missing_value = 1.e+20 ;
            longitude:valid_max = 180. ;
            longitude:axis = "X" ;
        int platform_code_string(platform_code_string) ;
        char platform_code(profile, platform_code_string) ;
            platform_code:missing_value = "N/A" ;
        int level(level) ;
            level:long_name = "Level" ;
            level:standard_name = "level" ;
            level:axis = "Z" ;
        double depth(profile, level) ;
            depth:standard_name = "depth" ;
            depth:long_name = "Depth" ;
            depth:valid_min = -4453.39918618524 ;
            depth:positive = "down" ;
            depth:units = "m" ;
            depth:missing_value = 1.e+20 ;
            depth:valid_max = 0. ;
            depth:axis = "Z" ;
        double temperature(profile, level) ;
            temperature:long_name = "sea water temperature" ;
            temperature:standard_name = "sea_water_temperature" ;
            temperature:units = "degree_Celsius" ;
            temperature:missing_value = 1.e+20 ;
        double salinity(profile, level) ;
            salinity:long_name = "sea water salinity" ;
            salinity:standard_name = "sea_water_salinity" ;
            salinity:units = "PSU" ;
            salinity:missing_value = 1.e+20 ;

    // global attributes:
            :Conventions = "CF-1.0" ;
            :software_version = "merge_profiles.py 2012-11-22" ;
            :creation_date = "2012-12-07T17:33:56Z" ;
            :history = "2012-12-07T17:33:56Z : Creation" ;
            :title = "Profiles" ;
            :southernmost_latitude = 43. ;
            :northernmost_latitude = 50. ;
            :westernmost_longitude = -15. ;
            :easternmost_longitude = -1. ;

Exemple
-------

Les profils actuels sont:

Emprise manche-gascogne
~~~~~~~~~~~~~~~~~~~~~~~

La racine de ces données se trouve ici:

/home11/caparmor/mars/VALID/DATA/MANGA

+-----------------------------------------------------------+---------------------------------------------+
| Description                                               | Localisation                                |
+===========================================================+=============================================+
| Profils BATHY provenant du GTS                            | HYDRO_CORIOLIS/MANGA_CORIOLIS_PR_BA         |
+-----------------------------------------------------------+---------------------------------------------+
| Profils CTD                                               | HYDRO_CORIOLIS/MANGA_CORIOLIS_PR_CT         |
+-----------------------------------------------------------+---------------------------------------------+
| Profils flotteurs Argo                                    | HYDRO_CORIOLIS/MANGA_CORIOLIS_PR_PF         |
+-----------------------------------------------------------+---------------------------------------------+
| Recopesca                                                 | HYDRO_CORIOLIS/MANGA_CORIOLIS_PR_RE         |
+-----------------------------------------------------------+---------------------------------------------+
| Profils TESAC provenant du GTS                            | HYDRO_CORIOLIS/MANGA_CORIOLIS_PR_TE         |
+-----------------------------------------------------------+---------------------------------------------+
| Profils XBT ou XCTD                                       | HYDRO_CORIOLIS/MANGA_CORIOLIS_PR_XB         |
+-----------------------------------------------------------+---------------------------------------------+
| Profils possiblement en double avec les données SISMER    | HYDRO_CORIOLIS_DOUBLES_POTENTIELS           |
+-----------------------------------------------------------+---------------------------------------------+
| Profils Bouteilles                                        | HYDRO_SISMER/MANGA_SISMER_PR_BO.nc          |
+-----------------------------------------------------------+---------------------------------------------+
| Profils CTD                                               | HYDRO_SISMER/MANGA_SISMER_PR_CT.nc          |
|                                                           | HYDRO_SISMER/MANGA_SISMER_PR_CT_DEPTH.nc    |
+-----------------------------------------------------------+---------------------------------------------+
| Profils XBT                                               | HYDRO_SISMER/MANGA_SISMER_PR_XB.nc          |
+-----------------------------------------------------------+---------------------------------------------+


Emprise méditéranée nord occidentale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

La racine de ces données se trouve ici:

/home11/caparmor/mars/VALID/DATA/MANGA

+--------------------------------------------------------+----------------------------------------+
| Description                                            | Localisation                           |
+========================================================+========================================+
| Profils BATHY provenant du GTS                         | HYDRO_CORIOLIS/MENOR_CORIOLIS_PR_BA    |
+--------------------------------------------------------+----------------------------------------+
| Profils CTD                                            | HYDRO_CORIOLIS/MENOR_CORIOLIS_PR_CT    |
+--------------------------------------------------------+----------------------------------------+
| Profils flotteurs Argo                                 | HYDRO_CORIOLIS/MENOR_CORIOLIS_PR_PF    |
+--------------------------------------------------------+----------------------------------------+
| Profils TESAC provenant du GTS                         | HYDRO_CORIOLIS/MENOR_CORIOLIS_PR_TE    |
+--------------------------------------------------------+----------------------------------------+
| Profils XBT ou XCTD                                    | HYDRO_CORIOLIS/MENOR_CORIOLIS_PR_XB    |
+--------------------------------------------------------+----------------------------------------+
| Profils possiblement en double avec les données SISMER | HYDRO_CORIOLIS_DOUBLES_POTENTIELS      |
+--------------------------------------------------------+----------------------------------------+
| Profils Bouteilles                                     | HYDRO_SISMER/MENOR_SISMER_PR_BO.nc     |
+--------------------------------------------------------+----------------------------------------+
| Profils CTD                                            | HYDRO_SISMER/MENOR_SISMER_PR_CT.nc     |
+--------------------------------------------------------+----------------------------------------+













