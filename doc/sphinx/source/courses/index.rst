Cours Python/CDAT/VACUMM
========================

Ce cours est consitué principalement de TP, 
avec quelques généralités en guise d'introductions.

Pour charger l'environnement Python/CDAT/VACUMM sous Caparmor:

.. code-block:: bash

    $ module purge
    $ module use /home11/caparmor/mars/PYTHON/modulefiles
    $ module load uvcdat121

Vous pouvez créer un alias dans votre .cshrc pour vos futures utilisations:

.. code-block:: bash

    alias uvcdat121 'module purge ; module use /home11/caparmor/mars/PYTHON/modulefiles ; module load uvcdat121'

(Vous pouvez retrouver l'environnement Python/CDAT/VACUMM sur ops en exécutant la commande: source /export/home/logiciels/modulefiles/uvcdat120.csh )

Récupérez les scripts en exécutant le script :program:`vacumm_courses.py` :

.. code-block:: bash

    $ vacumm_courses.py init

Pour les installer dans un autre répertoire que le répertoire courant, tapez par exemple :

.. code-block:: bash

    $ vacumm_courses.py init workdir

Et pour revenir à cette aide :

.. code-block:: bash

    $ vacumm_courses.py show


Pour l'aide du script :

.. code-block:: bash

    $ vacumm_courses.py -h

.. toctree::

    python
    numpy
    f2py
    cdat
    vacumm_intro
    time
    io
    grids
    plot
    filters
    masking
    regrid
    interp
    bigfiles
    bathy
    shorelines
    phys
    advanced
    sphinx
    
    
