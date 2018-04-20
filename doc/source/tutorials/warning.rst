.. _tut.warn:
    
Avertissement
=============

Chacun des tutoriels existe sous la forme d'un script python 
pouvant être excécuté tel quel dans son répertoire.
Certains d'entre eux aboutissent à la création de figures et/ou de fichiers.

.. rubric:: Figures et ``__file__``

Les scripts qui créent des figures utilisent systématiquement
la variable ``__file__`` (example : ``map(data, savefigs=__file__``).
Cette variable n'est pas accessible directement en ligne de commande python,
mais uniquement lors de l'execution de scripts en tant que fichier.
Ces scripts ne peuvent donc être excécutés entièrement depuis
l'interpéteur python si la variable ``__file__`` est utilisée.
Vous pouvez contourner le problème en ramplaçant l'option
``savefigs=__file__`` par ``show=True`` 
(ou ``savefigs(__file__)`` par ``P.show()``).

.. note::
    
    L'option ``savefigs`` et la fonction :func:`~vacumm.misc.plot.savefigs`
    sont une version particulière de l'option ``savefig`` 
    et de la fonction :func:`~matplotlib.pyplot.savefig` : 
    elles permettent la sauvegarde de figure à la fois au format png et pdf,
    et avec une résolution adaptée pour la documentation.
    Vous n'avez généralement pas à utiliser cette version spécialisée
    dans vos développements.


.. rubric:: Fichiers de données

Certains scripts lisent des fichiers de données
situés dans un répertoire
relatif à leur emplacement : vous ne pouvez donc les exécuter
correctement en dehors de ce répertoire si le chemin d'accès 
à ces fichiers n'est pas corrigé.
    

