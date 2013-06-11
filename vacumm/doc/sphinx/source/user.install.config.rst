.. _user.install.config:

Configuration de la librairie
=============================

.. sidebar:: Quelle est ma configuration ?

    Si votre librairie n'est pas trop ancienne,
    tapez:
        
        >>> from vacumm.config import print_config
        >>> print_config()
    
    Ou directement dans le shell :
        
    .. code-block:: bash
    
        vacumm_print_config.py
        
    Voir :ref:`user.scripts.vacumm_print_config`.

But
---

Certains modules ont par exemple besoin d'avoir accès à des
fichiers externes (comme un fichier bathymétrique) qui ne sont
pas inclus directement dans la librairie car trop volumineux.
Afin de rendre l'installation générique et portable sur 
différents réseaux, il est nécessaire de ne pas mettre ces 
chemin d'accès en dur dans le code : celui-ci fait donc appel
à un système de configuration. 
Il existe ainsi une configuration par défaut (celle du réseau IFREMER),
et il est possible d'altérer cette configuration pour l'adapter par exemple
à son réseau ou ses besoins.
Par exemple, Actimar possède sa configuration propre.
Et chaque utilisateur peut lui aussi avoir sa configuration qui
n'altère qu'une partie des paramètres.

La notion de configuration intervient donc à deux niveaux :

- Au niveau du développeur qui va faire appel à l'API (:mod:`vacumm.config`)
  pour définir et utiliser la configuration dans son code.
- Au niveau de l'utilisateur qui peut modifier son fichier de
  configuration personnel.


Fonctionnement
--------------

Un certain nombre de modules possèdent une configuration par défaut
prenant la forme d'un fichier :file:`config.cfg` se trouvant
dans le même répertoire que le module, ou dans le répertoire racine de la librairie.
La configuration du module est alors accessible à la section
portant son nom. Par exemple, si le module :mod:`vacumm.toto.tutu`
est configurable, son fichier de configuration peut être situé dans 
:file:`vacumm/toto/config.cfg`, avec le contenu:

.. code-block:: ini

    [vacumm.toto.tutu]
    masection = mavaleur
    
    
En outre, les fichiers de configuration supporte l'expansion de variable
sous deux formes :
    
    #. Une option présente dans une valeur est automatiquement remplacée par
       sa valeur (elle est "étendue") quand elle est présente sous la forme
       classique ``%(option)s``. Par exemple :
           
       .. code-block:: ini
       
           [vacumm.toto.tutu]
           option1 = foo
           option2 = %(option1)bar
           
       ``option2`` prend alors la valeur ``foobar``.
       Il existe plusieurs noms de répertoires qui sont
       actuellement systématiquement étendus :
           
           - ``lib_dir`` : Répertoire racine de la librairie 
             (voir :func:`~vacumm.config.get_lib_dir`).
           - ``data_dir`` : Répertoire où sont stockées les données utilisées
             par la librairie et les tutoriels 
             (voir :func:`~vacumm.config.get_data_dir`).
           - ``tut_dir`` : Répertoire contenant les scripts de tutoriel 
             (voir :func:`~vacumm.config.get_tut_dir`).
           - ``dist_dir`` : Répertoire du paquet de distribution dans le cas
             des développeurs 
             (voir :func:`~vacumm.config.get_dist_dir`).
           - ``conf_dir`` : Répertoire général contenant les fichiers
             de configuration 
             (voir :func:`~vacumm.config.get_conf_dir`).
           - ``user_conf_dir`` : Répertoire utilisateur contenant les fichiers
             de configuration 
             (voir :func:`~vacumm.config.get_user_conf_dir`).
           - ``mod_dir`` : Répertoire du module concerné 
             (voir :func:`~vacumm.config.get_mod_dir`).
           
    #. Le chemin d'accès ``~`` est automatiquement remplacé par le 
       home de l'utilisateur, et les variables d'environnement
       présentent sous la forme ``$HOME`` ou ``${HOME}``
       sont remplacées par leur valeur.
       
       .. note:: Les développeurs amenés à gérer des configuration
           peuvent se référer à la documentation du module :mod:`vacumm.config`.

On accède dans le code aux valeurs de la configuration grâce à la fonction
:func:`vacumm.config.get_config_value`.


Configuration par défaut
------------------------

La configuration générale par défaut (chargée grâce à :func:`~vacumm.config.get_default_config`) 
est la suivante :

.. literalinclude:: vacumm.cfg
    :language: ini


Configuration utilisateur
-------------------------

Il est possible d'altérer la configuration par défaut d'un ou plusieurs modules.
Le meilleur moyen est d'éditer le fichier :file:`$HOME/.config/vacumm/vacumm.cfg`.
Pour **éditer** rapidement le fichier : :ref:`user.scripts.vacumm_edit_config`.

.. obsol
    Il existe en outre d'autres possibilités, mais qui doivent être considérées comme obsolètes.
    Les voici par ordre de priorité (voir :func:`~vacumm.config.get_config_files`) :
        
        #. Créer un fichier :file:`vacumm.cfg` dans votre **répertoire courant** d'utilisation.
        #. Installer un fichier :file:`site.cfg` ou :file:`vacumm.cfg` dans le **répertoire
        principal de la librairie téléchargée** si vous l'utilisez directement sans l'installer
        (vous êtes donc un développeur !).
        #. Installer un fichier :file:`site.cfg` dans le **répertoire du module** que vous ciblez.
        #. Installer un fichier :file:`site.cfg` dans le **répertoire racine de la librairie**.
    
Si vous installez la librairie, vous pouvez la reconfigurer pour les
utilisateur en spécifiant l'option *--cfgfiles* à l'execution de :file:`setup.py`
(voir :ref:`user.install.install.config`).

Vous pouvez enfin altérer cette configuration en ligne grâce à la fonction
:func:`vacumm.config.set_config_value`.

