.. _user.scripts.evalrunsst:

:program:`init_eval_run_sst_op.py`
==================================

Ce script permet de compiler des statistiques concernant la température de surface de la mer (SST) à partir de données NetCDF (utilisation de la librairie cdms2) observéeset simulées sur une période de temps donnée choisie ou en mode opérationnel (temps réel).

Ce script permet également de produire des figures (format png) ainsi que des fichiers NetCDF (un par statistique), entre autre à destination du site previmer.org.

Préalables
----------

Les données d'entrée nécessaires au script sont:
    
- les résultats de simulation numériques disponibles au format NetCDF.
- les observations satellites disponibles au format NetCDF. Le script actuel est capable de récupérer les données satellites SEVIRI et NAR sur site ftp.

Fonctionnement
--------------

Ce script permet de calculer et d'afficher:
    
- les différences de SST entre le modèle et les observations satellites,
- les moyennes temporelles et spatiales,
- l'écart type pour le modèle et les observations(et la différence des deux) pour chaque pas de temps,
- l'evolution du biais moyen entre le modèle et les observations,
- le biais moyen, minimum et maximum et systématique entre le modèle et les observations,
- le taux de couverture des observations en chaque pas de temps et l'évolution temporelle de ce taux de couverture,
- la RMS et la RMS centrée,
- la corrélation temporelle,
- un rapport de Validation au format html.

Ces analyses peuvent être réalisées sur tout le domaine ou par sous régions(exemple pour Manche Gascogne: Vendée, Bretagne Sud, Manche Ouest, Manche Est, pays Basque).

Utilisation
-----------

Les réglages s'effectuent au moyen d'un fichier config.cfg. Les attributs suivants de ce fichiers sont à vérifier avant de lancer une évaluation: ::

        [Statistics]
          Définit les statistiques à compiler.
        [Domain]
          Définit le domaine étudié
        [Model Description] 
          Définit le type de modèle considéré et son mode de rapatriement.
        [MARS F1]
          Répertoire où récupérer les sorties du modèle.
        [Observations]
          Définit le type d'observations considérées et leur mode de rapatriement
        [Report]
          Générer ou non un rapport html avec les fichiers image (.png) produits. 
        [Output]
          Générer ou non un fichier NetCDF par statistique et le répertoire où ils seront écrits. 
        [Time Period]
          Les dates (année, mois, jours, heures) de début et de fin de la période considérée.
        [Env]
          Répertoire de travail, où sont écrits et chargés les fichiers NetCDF. 

.. warning::
 Pour le bon fonctionnement de ce script il faut que soient présents dans le même répertoire les fichiers suivants: 

         init_eval_run_sst_op_previmer.py; 

         eval_run_sst_op_previmer.py ; 

         global_stat.py ; 

         config.cfg ; 

         congig_r.cfg


Exemple d'utilisation:
----------------------

Configuration du Config.cfg pour générer toutes les stats possibles globalement et régionalement.

.. code-block:: bash

    [Control]
    interp_map = True
    tseries = True

    [Domain]
    lamax = 52
    lomax = 4.25
    lomin = -8
    lamin = 43.2

    [Statistics]
    allstat = True
    regionalstat = True
    extreme_bias = True
    evol_biais = True
    mean_biais = True
    spatial_std = True
    spatial_rmsc = True
    to_do = True
    obs_coverage = True
    temporal_rmsc = True
    temporal_rms = True
    spatial_mean = True
    obs_spatialcoverage = True
    temporal_mean = True
    temporal_corr = True
    extrema = True
    temporal_std = True
    spatial_rms = True
    systematic_bias = True

    [Model Description]
    download = nfs
    name = mars_manga

    [RECOPESCA]
    pwd = ifre3005
    data_dir = co052404/monthly/PR_RE
    user = c1f1a3
    url_cdoco = eftp.ifremer.fr
    data_dir = products/gridded/experimental-cms/netcdf/msg
    url_cersat = ftp.ifremer.fr/pub/ifremer/cersat
    compress = bz2

    [MARS F1]
    url_f1 = eftp.ifremer.fr/f1_4000/best_estimate/
    fic_prefix = PREVIMER_F1-MARS3D-MANGA4000_
    pwd = Ef0XmnZ4
    time_res = 1
    user = c1e975
    dir = /home/coriolis_exp/spool/co01/co0123/co012302/co01230207_v2/f1_4000/best_estimate/

    [MARS F2]
    url_f2 = eftp.ifremer.fr/f2/best_estimate/
    fic_prefix = PREVIMER_F2-MARS3D-MENOR_
    pwd = c1sisisi
    time_res = 3
    user = c1ef32

    [Observations]
    download = True
    product = SEVIRI SST
    download_detail = nfs
    clean_data = False

    [Env]
    workdir = /work/sskrypni/data_test/

    [Report]
    rep_html = True

    [Output]
    net_cdf = True
    title = Mars_Manga4000
    rep_ecriture = /work/sskrypni/netcdf

    [Time Period]
    mdeb = 02
    andeb = 2012
    jdeb = 21
    jfin = 28
    hfin = 23
    mfin = 02
    anfin = 2012

Le code peut être exécuté de la manière suivante sur caparmor:

.. code-block:: bash

        python init_eval_run_sst_op_previmer.py

.. note::

        Le programme ``init_eval_run_sst_op_previmer.py`` est codé pour utilisation sur la machine caparmor (service4). Il inclut une soumission en batch des calculs (qsub).



Aperçu des sorties
~~~~~~~~~~~~~~~~~~

.. image:: result_mean_bias_all_20120221_20120228.png 
    :width: 80%

.. image:: result_spatial_statmean_manche_est_20120221_20120228.png 
    :width: 80%
