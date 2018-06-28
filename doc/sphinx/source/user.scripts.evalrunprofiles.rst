.. _user.scripts.evalrunprofiles:

:program:`eval_run_profiles2.py`
==================================

Ce script permet de compiler des statistiques concernant la température et la salinité dans la colonne d'eau sur une période de temps donnée choisie ou en mode opérationnel (temps réel).

Ce script permet également de produire des figures (format png) ainsi que des fichiers NetCDF (un par statistique), entre autre à destination du site previmer.org.

Préalables
----------

Les données d'entrée nécessaires au script sont:

- les résultats de simulation numériques disponibles au format NetCDF.
- les observations RECOPESCA disponibles au format .csv.

Fonctionnement
--------------

Ce script permet de calculer et d'afficher:

- les différences de température dans la colonne d'eau entre le modèle et les observations RECOPESCA,
- les différence en salinité dans la colonne d'eau entre le modèle et les observations RECOPESCA,
- les moyennes temporelles et spatiales,
- les sections verticales d'une suite temporelle dans une sous-région donnée.




Utilisation
-----------

Les réglages s'effectuent au moyen du fichier config_profiles.cfg. Les attributs suivants de ce fichiers sont à vérifier avant de lancer une évaluation: ::
    [Time Period]
	Les dates (année, mois, jours, heures) de début et de fin de la période considérée.
    [Action]
	Définit les statistiques et les figures à génerer.
    [Output]
	Définit le répertoire cible où seront enregistrés les fichiers binaires .nc

    [Histogram]
	Produit ou non les histogrammes de séries temporelles des profils observés.
    [Env]
	Prise en compte des sorties modèles.

    [Domain]
	Définit le domaine étudié

    [SubArea]
	Définit le sous-domaine étudié

.. warning::
 Pour le bon fonctionnement de ce script il faut que soient présents dans le même répertoire les fichiers suivants:

         cprofile.py;

         config_profiles.cfg ;

         config_profiles.ini;



Exemple d'utilisation:
----------------------

Configuration du Config.cfg pour générer toutes les stats possibles globalement et régionalement.

.. code-block:: bash

    [Time Period]

    mdeb=6
    jdeb=01
    andeb=2012
    mfin=6
    jfin=10
    anfin=2012

    [Action]

    single_profiles=False
    map_profiles=True
    hist_profiles=True
    timeserie_profile=True
    Valid_Stats= True
    Control_Stats= False
    Vertical_Section=True
    [Output]
    title= mars_manga
    rep_ecriture= /work/sskrypni/netcdf/RECOPESCA

    [Histogram]

    tstep=1
    tstep_unit=day

    [Env]

    use_model=True

    [Domain]

    Lamax = 52
    Lomax = 4.25
    Lomin = -8
    Lamin = 43.2

    [SubArea]

    Lamax = 47
    Lomax = -2
    Lomin = -6
    Lamin = 46

Le code doit être exécuté de la manière suivante sur caparmor:

.. code-block:: bash

        python eval_run_profiles2.py cfgfile="config_profiles.cfg"

.. note::

	cfgfile="config_profiles.cfg" définit le fichier de configuration spécifique à chaque utilisateur qui implémentera le fichier de configuration par défaut config_profiles.ini.
	Chaque utilisateur peut ainsi définir ses propres options de configuration.


Aperçu des sorties
~~~~~~~~~~~~~~~~~~

.. image:: images/Map_RECOPESCA_profile_20120609172349_1-72099995613W_43-4146995544N.png
    :width: 80%

.. image:: images/Vertical_Section_temp2012060100000020120610230000.png
    :width: 80%

.. image:: images/Control_stats2012060100000020120610230000.png
    :width: 80%

