.. _user.dev:

For developpers
***************

Use your own version
====================

To use the library without installing it, see doc sections :ref:`user.install.install.dev` and 
:ref:`user.install.modenv.dev`.



Coding rules
============

Please follow the :ref:`appendix.conventions`.



Architecture
============

L'architecture doit de la librairie ne doit pas s'éloigner
de celle mentionnée à la section :ref:`user.algo.arch`.
Si un nouveau module est ajouté, celui-ci doit être intégré
à :mod:`vacumm.misc` s'il est de nature générique.
Si ce module est très spécialisé, il peut intégré au niveau supérieur, 
c'est-à-dire directement dans :mod:`vacumm`.

.. _user.dev.doc:

Documentation
=============

Compilation
-----------

La documentation est écrite en langage  `rst <http://docutils.sourceforge.net/rst.html>`_
et est compilée par `Sphinx <http://sphinx.pocoo.org>`_ .
Les sources se trouvent dans le sous-répertoire :file:`doc/sphinx/source`.
Pour la compiler : 

.. code-block:: bash

    shell> cd doc/sphinx
    shell> make      # html + pdf
    shell> make html # html seul
    shell> make pdf  # pdf seul via latex

La documentation compilée se trouve alors dans les répertoires 
:file:`doc/sphinx/build/html` et :file:`doc/sphinx/build/latex`.

Les exemples de code présentés dans les :ref:`tutoriels <user.tut>` de la documentation 
sont des scripts situés dans le sous répertoire :file:`doc/sphinx/source/tutorials/python`.
Certains d'entre eux utilisent des fichiers de données situés dans le 
sous répertoire :file:`data`.


Dépendances
-----------
    
Certains des scripts de tutoriel génèrent des figures ou fichiers inclus ensuite dans
la documentation. Il est nécessaire de les exécuter afin de cette dernière puisse
être compilée correctement.

En outre, la documentation des modules peut aussi utiliser des figures
(notamment :mod:`~vacumm.misc.color`).
Celles-ci sont créées grâce au script :program:`doc/sphinx/source/library/make.py`,
qui doit donc être exécuté avant la compilation de la documentation 
(:file:`doc/sphinx/Makefile` s'en charge automatiquement).

.. code-block:: bash

    shell> cd doc/sphinx/library
    shell> python make.py
    
Enfin, certaines figures de la documentation (comme le logo) sont créées par :command:`pdflatex`.  

.. code-block:: bash

    shell> cd doc/sphinx/source
    shell> make

Il est possible que certaines pages incluent des formules mathématiques.
Afin qu'elle soit compilées et affichées correctement par l'extension sphinx
:class:`sphinx.ext.pngmath`, vous devez avoir le programme :program:dvipng
accessible.


Génération des figures TikZ
---------------------------

Cette documentation contient plusieurs figures 
dessinées grâce à :program:`pdflatex` et 
`PGF/TikZ <http://pgf.sourceforge.net>`_ (logo, architexture 
de la librairie).
L'avantage est de pouvoir placer les sources des figures
sur le dépôt svn, et que tout le monde puisse re-générer ces dernières.

Il est nécessaire d'avoir une version récente de PGF/TikZ, 
que vous pouvez par exemple obtenir ici : http://www.texample.net/tikz/builds/  (voici `une version <http://media.texample.net/pgf/builds/pgfCVS2010-09-28_TDS.zip>`_).
Pour l'installation, procédez ainsi :
    
.. code-block:: bash

    shell> mkdir -p ~/texmf
    shell> cd ~/texmf
    shell> wget http://media.texample.net/pgf/builds/pgfCVS2010-09-28_TDS.zip
    shell> unzip pgfCVS2010-09-28_TDS.zip
    shell> rm pgfCVS2010-09-28_TDS.zip

.. sidebar:: Qu'est-ce que PGF/TikZ ?

    Il s'agit d'une librairie permettant de créer des figures
    de haute qualités à partir d'un code source intégré dans du TeX.
    Le meilleur aperçu est fourni par le site qui recense les
    exemples (tutoriels) :  http://www.texample.net/tikz/examples .
    La plupart d'entre eux se base sur une version plus récente que celle
    installée par défaut sur un système.
    C'est donc sur une *build* CVS que les figures de cette documentation
    sont basées.

Les figures sont ensuite générées avec la commande suivante :  

.. code-block:: bash

    shell> cd doc/sphinx/sources
    shell> make

Le code latex est alors compilé, générant un pdf qui est ensuite converti
au format ppm, puis au format png.



Tutoriels et vérifications
==========================

Il fortement suggéré aux développeurs de la librairie de créer des tutoriels
sur les fonctionnalités importantes qu'ils développent.
Ces tutoriels ont deux intérêts :
    
    - Ils compléte la documentation.
    - Il permettent d'effectuer des tests de vérification,
      grâce au programme :program:`check.py`.
      
Le programme :program:`check.py` est situé dans le répertoire de tutoriels
(:file:`doc/sphinx/source/tutorials/python`).
Son utilisation est la suivante :

.. code-block:: bash

    shell> check.py [options] [pattern1 [pattern2] ...]
    
Celui-ci prend comme argument un (ou plusieurs) *global pattern*
permettant de lister les scripts à tester.
La valeur par défaut est ``"*.py"``.
Il est ensuite possible d'exclure des scripts de cette liste
grâce à l'option :option:`-e`.

Ce script affiche en console tout ou partie des informations sur les tests, 
et logue toutes les informations dans le fichier :file:`check.log`.
Le niveau d'affichage des information est console est modifiable
avec l'option :option:`-l`.

Exemples d'utilisation :
    

.. code-block:: bash

    shell> check.py -e "misc.color.py" -e misc.grid.masking.* misc.*.py
    shell> check.py --loglevel=debug

    
.. rubric:: Options de :program:`check.py`
    
.. program:: check.py :

.. option:: -h, --help

    Affiche l'aide.
    
.. option:: -e, --exclude

    Ajoute un *global pattern* listant des scripts à exclure des tests.
    
.. option:: -l, --loglevel

    Définit le niveau de journalisation dans la console.
    Celui-ci peut avoir les valeurs suivantes :
        
        - ``"debug"``: Affiche des sorties standard et d'erreur.
        - ``"info"``: Affiche le nom des scripts ayant passé le test (choix par défaut).
        - ``"error"``: Affiche le nom des scripts n'ayant pas passé le test.


Distribution d'un paquet
========================

Il est possible de créer des paquets correspondant typiquement à des versions précises (par exemple stables).
La procédure est la suivante :
    
.. code-block:: bash

    shell> python setup.py bdist
    
Cette commande va alors créer un fichier distribuable, dont le nom est proche de 
:file:`vacumm-0.9-svn128.linux-x86_64.tar.gz`.
Ce fichier peut alors être déposé dans la section fichiers du site gforge du projet
(`à cette adresse <https://forge.ifremer.fr/frs/admin/qrs.php?package=&group_id=93>`_).


