.. _user.install.modenv:

Using environment modules
=========================

Introduction
------------

Comme cela a été évoqué dans les deux sections précédentes,
l'accès au python fourni par CDAT nécessite la modification
de la variable d'environnement :envvar:`PATH`, et l'accès à la librairie
nécessite la modification
de la variable d'environnement :envvar:`PYTHONPATH` si elle a été installée dans un répertoire tiers.

Sachant qu'il est possible de faire coéxister plusieurs versions de CDAT et de la librairie,
en fonction des besoins (dans un cadre opérationnel, pour un projet donné, 
pour les utilisateurs standards ou ceux voulant les dernières fonctionnalités),
il peut devenir difficile à un administrateur de faire en sorte que
les utilisateurs puissent accéder facilement aux bonnes versions ;
de même qu'il serait désagréable aux utilisateurs de devoir changer eux même
toutes les variables d'environnement (c'est à dire de savoir quelles variables
et connaître tous les chemins d'accès).
Un moyen de gérer facilement ce problème est l'utilisation de `modules d'environnement <http://modules.sourceforge.net>`_.
Cette approche permet à l'administrateur de gérer de nombreuses versions et leur dépendances
de façon transparente pour l'utilisateur qui n'a alors qu'a insérer deux lignes
dans son :file:`.cshrc` (ou :file:`.bashrc`) et le fonctionnement 
de la commande :command:`module` charger l'environnement python qu'il souhaite.

Les documentations de l'utilitaire se trouve ici:

    - http://modules.sourceforge.net/
    - http://modules.sourceforge.net/man/module.html
    - http://modules.sourceforge.net/man/modulefile.html
    - http://sourceforge.net/p/modules/wiki/FAQ/

Module par défaut
-----------------

L'utilitaire de module permet de faire cohabiter plusieurs version de module pour un même logiciel.

Si on a par exemple :

.. code-block:: bash

    shell> ls /my/modules/vacumm
    0.1
    1.0
    trunk

Pour positionner une version comme étant celle par défaut  il faut alors créer dans ce répertoire un fichier :file:`.version` contenant :

.. code-block:: tcl

    #%Module
    set ModulesVersion "1.0"

Ce fichier :file:`.version` indique alors quel fichier de module par défaut il faut charger quand seul le nom de
répertoire est spécifié. La commande ``module load vacumm`` chargera alors le fichier de module :file:`/my/modules/vacumm/1.0`.

Le côté administrateur
----------------------

Il est d'abord nécessaire que l'administrateur configure l'utilisation
des modules d'environnement.

La première étape consiste à créer un réportoire dédié aux fichiers de modules, 
puis des sous-répertoires dédiés aux versions de CDAT et de vacumm et aux dépendances :
    
.. code-block:: bash

    shell> mkdir /my/modules
    shell> mkdir /my/modules/cdat
    shell> mkdir /my/modules/vacumm
   
Il convient ensuite de créer les fichiers de module pour CDAT et ses dépendances
(voir :ref:`user.install.prereq`). 
Par exemple, SCRIP est une des dépendances. 
S'il est installé dans le répertoire :file:`/my/soft/scrip`, 
il faut créer le fichier de module :file:`/my/modules/scrip` avec le contenu suivant:
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "SCRIP regridder"
    }
    
    prepend-path PATH /my/soft/scrip
    
Le module python :mod:`~mpl_toolkits.basemap` utilise la librairie GEOS.
On crée alors le fichier de module :file:`/my/modules/geos/3.2` :
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "GEOS library"
    }
    
    set version 3.2
    setenv GEOS_LIB /my/soft/geos/$version/lib
    prepend-path LD_LIBRARY_PATH /my/soft/geos/$version/lib
    

**Note**: geos est désormais inclus avec cdat >= 6.0, il n'est donc pas nécessaire
d'avoir une installation et un module de geos à part.

Le fichier de configuration de CDAT :file:`/my/modules/cdat/5.2` peut alors être créé 
en tenant compte des dépendances et des conflits avec d'autres versions :

.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "SCRIP regridder"
    }
    
    conlict cdat
    module load scrip
    module load geos/3.2
    
    set version 5.2
    set base_path /my/soft/cdat
    prepend-path PATH $base_path/$version/bin
    prepend-path LD_LIBRARY_PATH  $base_path/$version/lib
    prepend-path LD_LIBRARY_PATH  $base_path/$version/Externals/lib
    prepend-path C_INCLUDE_PATH $base_path/$version/include
    
Si la version par défaut de vacumm est installée dans le python de CDAT par défaut,
le fichier de module correspondant :file:`/my/modules/vacumm/1.0` contient alors :
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "VACUMM python library"
    }
    
    conlict vacumm
    module load cdat

Ici on ne fait que charger CDAT étant donné que vacumm est installé dans le répertoire par défaut
de l'installation python de CDAT.

Si une autre version de VACUMM existe (installée avec un --prefix égal à :file:`/my/soft/vacumm-recent`), 
et utilisant une autre version de CDAT, nous pourrions avoir le fichier :file:`/my/modules/vacumm/recent` :
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "VACUMM python library"
    }
    
    conlict vacumm
    module load cdat/6.0
    prepend-path PATH /my/soft/vacumm-recent/bin
    prepend-path PYTHONPATH /my/soft/vacumm-recent/lib/python2.7/site-packages
    
Avec un module cdat/6.0.0 configuré sur cet autre CDAT.

Le côté utilisateur
-------------------
    
Le fichier :file:`$HOME/.cshrc` (ou :file:`$HOME/.bashrc`) va typiquement contenir les lignes suivantes :
    
.. code-block:: bash

    # Modules
    # - initialisation
    source /usr/share/modules/init/csh # ou bash
    # - ajout du repertoire des modules
    module use /my/modules

.. warning::
    
    Il faut cependant faire attention si l'on charge le module cdat dans un des fichiers de profile tel
    que :file:`$HOME/.bashrc` ou :file:`$HOME/.cshrc`. Ces fichiers étant chargés à l'ouverture de session
    (graphique ou non), des effets de bord peuvent se produire car la version de python, et les packages
    qui y sont installés ne sont plus ceux du système et certains logiciels voir même la session graphique
    complète pourrait être inutilisable.
    Si tel est le cas, il faudra alors se résoudre à charger le module cdat manuellement dans votre terminal.

L'utilsateur a alors accès à tous les modules d'environnement définis
dans la section précédente.
On vérifie les modules disponibles :
    
.. code-block:: bash

    shell> module avail
     
Le chargement de la librairie python vacumm se fait de la façon suivante :
    
.. code-block:: bash

    shell> module load vacumm
    
.. note::

    Il est possible de charger automatiquement la librairie en plaçant ``module load vacumm``
    directement dans le fichier :file:`$HOME/.cshrc` après la ligne ``module use ...``
    
Puis on vérifie :
    
.. code-block:: bash

    shell> module list
    Currently Loaded Modulefiles:
      1) scrip      2) geos/3.2    3) cdat/5.2     4) vacumm/recent
    shell> echo $PATH
    /my/soft/cdat/5.2/bin:/my/soft/scrip:... #etc
    shell> which python
    /my/soft/cdat/5.2/bin/python
    shell> python -c "import vacumm"
    
    
Pour changer de version de vacumm pour par exemple une installation du trunk qui serait régulièrement mise à jour :
    
.. code-block:: bash

    shell> module switch vacumm vacumm/trunk
    
.. _user.install.modenv.dev:

Le côté développeur
-------------------

Le checkout de vacumm contient un module :file:`etc/modulefiles/vacumm`, celui ci permet d'exploiter directement le checkout
de vacumm :

.. code-block:: bash

    shell> module use /path/to/my/vacumm/etc/modulefiles
    shell> module load vacumm
    
Pour ne pas rentrer en conflit avec (masquer) un autre jeu de modules, vous pouvez créer votre espace de modules
personnels et lier le module de votre checkout dans cet espace : 

.. code-block:: bash

    shell> mkdir -p ~/etc/modulefiles/vacumm
    shell> ln -s /path/to/my/vacumm/etc/modulefiles ~/etc/modulefiles/vacumm/dev

Ce module requiert un module nommé cdat, une erreur sera affichée si celui ci est introuvable.

.. warning::

    Si vous avez besoin de spécifier un module cdat particulier pour votre version de développement, ne modifier pas
    :file:`etc/modulefiles/vacumm` (pour ne pas le commiter !), créez plutôt un répertoire cdat dans votre espace de
    modules en créant un module qui chargera le cdat nécessaire et éventuellement en utilisant la méthode de module
    par défaut décrite ci dessus.

.. code-block:: bash

    shell> mkdir -p ~/etc/modulefiles/cdat
    shell> vi ~/etc/modulefiles/cdat/dev

Renseignez ce nouveau module :file:`cdat/dev` avec le contenu suivant :

.. code-block:: tcl

    #%Module
    module load cdat/6.0


Vous pouvez maintenant charger votre version de développement :

.. code-block:: bash

    shell> module load vacumm/dev




Et voilà, that's all folks !
    
