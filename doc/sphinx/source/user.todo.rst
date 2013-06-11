.. _user.todo:


Ce qu'il reste à faire
**********************


.. _user.todo.gen:

Liste générique
===============

.. todolist::

Stabilisation du code
=====================

Vous pouvez vous référer à :ref:`appendix.risks` pour avoir la liste des routines 
nécessitant une attention particulière, et donc par exemple des améliorations.


.. _user.todo.more:

Compléments
===========

Cette section liste les améliorations importantes à apporter à la librairie,
et éventuellement les clés pour y parvenir ou les voies possibles.



Interpolations
--------------

Quatree
~~~~~~~

Il serait souhaitable d'optimiser les interpolation de nuages de points 
vers grilles rectangulaire en utilisant une méthode de décomposition
du domaine avec un algorithme de type quadtree.
Actuellement, les blocs sont tous de taille similaire et peuvent
contenir un nombre de points très inhomogène.

Krigeage
~~~~~~~~

Il serait intéressant d'avoir la possibilité d'interpoler par krigeage.

