Leçon sur la génération de documentation via sphinx
===================================================

Introduction
------------

`Sphinx <http://sphinx.pocoo.org>`_ est un générateur de documentation qui donne envie de documenter.

- Langage ascii -> subversion -> partage.
- Syntaxe rapidement comprise.
- Auto-documentations (python, voire fortran avec :mod:`~vacumm.sphinxext.fortran_autodoc`).
- Intégration d'exemples.
- Sorties propres en PDF et HTML (et +).
- Assistants pour la création.
- Personnalisation du rendu.
- Extensions.
- De plus en plus répandu.


Une syntaxe de base très claire :

.. code-block:: rst

    Section
    =======
    Ref. : :func:`plot.map`
    

Exemple
-------

Exemple des différentes étapes pour créer sa documentation :

1. Lancez l'assistant de création puis répondez aux questions.

   ::
   
        $ sphinx-quickstart

   Pour l'utilisation des extensions Sphinx : repondez ``y``.
2. Éditez le fichier :file:`conf.py` si vous le souhaitez.
3. Créez ou éditez vos fichiers de documentation.
   Par exemple :
   
   - :file:`index.rst` : Fichier de base.
   - :file:`usersguide.rst` : Manuel de l'utilisateur.
   - :file:`library.rst` : Guide du développeur.
4. Générer la documentation ::
   
   Pour le html ::
   
        $ make html
        
   Pour le pdf ::
   
        $ make latexpdf

5. Visualisez les documentations dans le répertoire :file:`build`.


Exemple pour le fichier :file:`index.rst` :

.. code-block:: rst

    .. doctree::

        usersguide
        library

Exemple pour le fichier :file:`library.rst` :

.. code-block:: rst

    Guide du developpeur
    ####################

    :mod:`mymodule` --- Mon module
    ******************************

    .. automodule:: mymodule


