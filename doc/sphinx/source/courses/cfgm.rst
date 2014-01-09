.. _courses.cfgm:

Leçon sur la gestion des configurations avancées
================================================

Voir :mod:`vacumm.misc.config`, :mod:`configObj` (http://www.voidspace.org.uk/python/configobj.html),
:ref:`user.tut.generic.scripts`.

Le but est de pouvoir gérer des fichiers de configuration plus évolués qu'avec :mod:`ConfigParser` :

- Sections imbriquées.
- Valeurs par défaut via le fichier de spécification.
- Vérification des valeurs.
- Accès comme un dictionnaire.

Le module :mod:`vacumm.misc.config` offre une interface pratique pour gérer ces configurations,
et des outils pour les exploiter encore plus loin (ex: utilisation avec :mod:`argparse`).

Fichier :file:`courses_advanced_cfgm.py`

.. literalinclude:: ../../../../scripts/courses/courses_advanced_cfgm.py


Fichier de spécifications :file:`courses_advanced_cfgm.ini`

.. literalinclude:: ../../../../scripts/courses/courses_advanced_cfgm.ini


Fichier de configuration :file:`courses_advanced_cfgm.cfg`

.. literalinclude:: ../../../../scripts/courses/courses_advanced_cfgm.cfg

