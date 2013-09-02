Introduction à VACUMM
=====================

Organisation
------------

.. figure:: ../static/modules.*

    Architecture (incomplète) de VACUMM

Pourquoi VACUMM ?
-----------------

- Origine : :file:`srtools.py`
- Une librairie et des scripts orientés océano/météo qui regroupent le plus important.
- Quand l'utiliser :

    - Vous utilisez à la fois UVCDAT et Matplotlib.
    - Vous y trouver des fonctions et classes qui font déjà tout le travail.
    - Vous ne voulez pas réinventer la roue.
    - vous cherchez à partager des outils pouvant s'intégrer dans la librairie. 


Configuration
-------------

Voir :ref:`user.install.config`.

La configuration définit les paramètres de la librairie.
Un paramètre peut être par exemple le chemin d'accès au trait de côte du SHOM
(pas inclus dans la librairie),
ou le type de remplissage par défaut pour les plots 2D.

Il existe des valeurs par défaut définies dans des fichiers 
situés au même endroit que les sources.
D'autres fichiers peuvent être introduits à l'installation.
Enfin, chaque utilisateur peut definir sa propre configuration
dans :file:`~/.config/vacumm/vacumm.cfg`.

Vous pouvez afficher et éditer votre configuration notamment depuis le shell :

.. literalinclude:: ../../../../scripts/courses/courses_vacumm_config.sh



Où trouver de l'aide
--------------------

- **Site web :** http://www.ifremer.fr/vacumm ou http://relay.actimar.fr/vacumm/
- **Listes de diffusion :** https://forge.ifremer.fr/mail/?group_id=93
- **Mantis :** https://forge.ifremer.fr/mantis/my_view_page.php
- **En ligne :**

    >>> from vcmq import map2
    >>> help(map2)
    >>> import vacumm
    >>> vacumm.help('map2') # ouvre le navigateur
    
- **En direct :** Guillaume Charria  ou Stéphane Raynaud.



Utilisation
-----------


Par sa librairie :

>>> from vacumm.misc.plot import map2   # import complet
>>> from vcmq import map2               # raccourci


Par ses scripts :

.. code-block:: shell

    $ vacumm_print_config.py
    $ showtime.py -m myfile.nc
    
