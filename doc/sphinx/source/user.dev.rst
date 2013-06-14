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

The architecture of the library must follow the spirit
of what is explained at section :ref:`user.desc.basis.arch`.
For instance, generic stuff must go in the :mod:`vacumm.misc`
module.
When a new module is specialized, it should find a place elsewhere.

.. _user.dev.doc:

Generation of the sphinx documentation
======================================

Generation of the figures of scripts
------------------------------------

Some of the scripts located in the :file:`scripts` directory 
create figure that are included in the documentation.
You must run them before compiling the documentation to make
sure that all needed figures are available.
For that, just run a :program:`make` in the :file:`scripts` directory:
    
    
.. code-block:: bash

    $ cd scripts
    $ make


.. note:: It is not necessary to re-generate these figure each time you compile!


Compilation
-----------

The documentation is written in `rst <http://docutils.sourceforge.net/rst.html>`_ language,
and compilated with `Sphinx <http://sphinx.pocoo.org>`_ .
Source files are located in the  :file:`doc/sphinx/source` directory.
To compile it: 

.. code-block:: bash

    $ cd doc/sphinx
    $ make      # html + pdf
    $ make html # html only
    $ make pdf  # pdf only using pdflatex

The documentation is generated in directories 
:file:`doc/sphinx/build/html` and :file:`doc/sphinx/build/latex`.


Regeneration of TikZ figures 
----------------------------

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

    $ mkdir -p ~/texmf
    $ cd ~/texmf
    $ wget http://media.texample.net/pgf/builds/pgfCVS2010-09-28_TDS.zip
    $ unzip pgfCVS2010-09-28_TDS.zip
    $ rm pgfCVS2010-09-28_TDS.zip

.. sidebar:: Qu'est-ce que PGF/TikZ ?

    Il s'agit d'une librairie permettant de créer des figures
    de haute qualités à partir d'un code source intégré dans du TeX.
    Le meilleur aperçu est fourni par le site qui recense les
    exemples (tutoriels) :  http://www.texample.net/tikz/examples .
    La plupart d'entre eux se base sur une version plus récente que celle
    installée par défaut sur un système.
    C'est donc sur une *build* CVS que les figures de cette documentation
    sont basées.

TikZ figure can now be generated with:  

.. code-block:: bash

    $ cd doc/sphinx/sources
    $ make

Le code latex est alors compilé, générant un pdf qui est ensuite converti
au format ppm, puis au format png.


Writing tutorials
=================

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

    $ check.py [options] [pattern1 [pattern2] ...]
    
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

    $ check.py -e "misc.color.py" -e misc.grid.masking.* misc.*.py
    $ check.py --loglevel=debug

    
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


Distributing the library as a package
=====================================

Il est possible de créer des paquets correspondant typiquement à des versions précises (par exemple stables).
La procédure est la suivante :
    
.. code-block:: bash

    $ python setup.py bdist
    
Cette commande va alors créer un fichier distribuable, dont le nom est proche de 
:file:`vacumm-0.9-svn128.linux-x86_64.tar.gz`.
Ce fichier peut alors être déposé dans la section fichiers du site gforge du projet
(`à cette adresse <https://forge.ifremer.fr/frs/admin/qrs.php?package=&group_id=93>`_).


