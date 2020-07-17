Leçon d'intégration du fortran dans python via :program:`f2py`
==============================================================

La commande :program:`f2py` et le module associé permettent de convertir
un code fortran en modules et fonctions python.
Cela peut être utile pour deux raisons principales :

- Vous avec un code fortran complexe que vous ne voulez pas convertir en python.
- Vous voulez augmenter les performances d'un algorithme, notamment avec de la parallélisation.

Avec cette approche, vous pouvez donc utiliser du code fortran
sans avoir à faire des appels système et passer par
des fichiers intermédiaires. 

.. note:: Notez cependand qu'une bonne utilisation de :mod:`numpy` et :mod:`scipy` permet
    souvent d'obtenir des performances très élevées.

Voici un exemple minimal.

Fichier :file:`myfile.f90` :

.. literalinclude:: ../../../../scripts/courses/myfile.f90


Compilation et création du fichier :file:`mymodule.so` ::

    $ f2py -c -m mymodule myfile.f90


Utilisation (fixhier :file:`courses_f2py.py`) :

.. literalinclude:: ../../../../scripts/courses/courses_f2py.py


    