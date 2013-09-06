Leçon sur les arguments de scripts
==================================


Introduction
------------

Le but est de filtrer la ligne de commande python pour récupérer des options
et paramètres. La base des arguments :


>>> import sys
>>> print sys.argv

Écrivez un script avec ce code puis testez le contenu de ``sys.argv``.


Utilisation de :mod:`argpase`
-----------------------------

On va maintenant utilisez le module :mod:`argpase` qui permet de gérer facilement les arguments
de la ligne de commande et avoir un comportement à la linux.
Pour comprendre l'intérêt et le fonctionnement, 
**allez voir** le tutoriels :ref:`user.tut.generic.scripts`.



Mélanger ligne de commande et configuration
-------------------------------------------

Le but est d'utiliser un fichier de spécification pour générer des options en ligne de commande.
Les paramètres sont ainsi définis dans l'ordre de priorité suivant :

    1. Argument en ligne de commande: ``--logger-level=debug``.
    2. Option dans le fichier de config utilisateur : ``[logger] level=debug``
    3. Valeur par défaut du fichier de specs : ``[logger] level=choice(debug,info,warning,error,default=info)``

**Allez voir** le tutoriel :ref:`user.tut.misc.config.argparse` pour un exemple simple.

