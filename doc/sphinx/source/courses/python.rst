Leçon d'introduction à python
=============================


Pourquoi python ?
-----------------

Diverses raisons :

- Pas envie de payer une licence.
- Language interpété.
- Facilement interfaçable avec du fortran (pas d'appel système), et du C.
- Librairies : système, réseau, calcul numérique, graphique, etc.
- Pour être efficace ou faire comme les autres.


Environnement
-------------

Choisissez celui qui vous convient.

- **vim :** notamment pour les anciens.
- **emacs :** possibilité d'exécution directe + édition à distance.
- **ipython :** console pratique avec interface vim.
- **eric4 :** environnement complet avec gestion de projet, exécution et deboguage.
- Plein d'autres IDE.


Exécution
---------

Soit en direct ::

    $ python
    Python 2.7.3 (default, Jan 25 2013, 15:38:17) 
    [GCC 4.5.1 20100924 (Red Hat 4.5.1-4)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> 
    
Soit en exécution de script via python ::

    $ python myscript.py
    

Soit en direct exécution de script (si exécutable avec entête) ::

    $ myscript.py

Soit en soumission de job ::

    PBS -q sequentiel
    PBS -N mypython
    module load mypython
    python myscript.py
    
.. warning:: Si vous passez par un job, placez la ligne suivante en début de script en cas de tracés

    >>> from matplotlib import use ; use('Agg')

Le language
-----------

Fichier :file:`courses_python.py`

.. literalinclude:: ../../../../scripts/courses/courses_python.py


Attention aux règles de codage : :ref:`appendix.conventions`.


Modules importants
------------------

Librairie standard :

- Outils système : :mod:`os`, :mod:`sys`, :mod:`shutil`, :mod:`glob`, :mod:`stat`, :mod:`subprocess`, :mod:`popen2`.
- Entrées sorties : :mod:`ConfigParser`, :mod:`argparse`, :mod:`csv`, :mod:`logging`.
- Divers : :mod:`re`, :mod:`collections`, :mod:`pickle`, :mod:`datetime`, :mod:`time`.

Librairies scientifiques :

- Numérique : :mod:`numpy`.
- Scientifique : :mod:`scipy`.
- Graphiques : :mod:`matplotlib`, :mod:`basemap`.
- Océano/Météo (UVCDAT) : :mod:`cdms2`, :mod:`cdutil`, :mod:`genutil`.

Divers :

- Configurations avancées : :mod:`configobj`.


La librairie :mod:`scipy` regroupe beaucoup d'algorithmes scientifiques intéressants.
Il est important de parcourir son contenu avant d'aller voir ailleur.

-  Clusters (:mod:`~scipy.cluster`)
-  Constantes (:mod:`~scipy.constants`)
-  Transformées de Fourier (:mod:`~scipy.fftpack`)
-  Integration et ODEs (:mod:`~scipy.integrate`)
-  Interpolation (:mod:`~scipy.interpolate`)
-  Entrés/sorties (:mod:`~scipy.io`)
-  Algèbre linéaire (:mod:`~scipy.linalg`)
-  Modèle d'entropie max (:mod:`~scipy.maxentropy`)
-  Traitement d'images (:mod:`~scipy.ndimage`)
-  *Orthogonal distance regression* (:mod:`~scipy.odr`)
-  Optimisation  (:mod:`~scipy.optimize`)
-  Traitement du signal (:mod:`~scipy.signal`)
-  Matrices creuses (:mod:`~scipy.sparse`)
-  Algèbre linéaire pour matrices creuses (:mod:`~scipy.sparse.linalg`)
-  Algorithmes spatiaux (:mod:`~scipy.spatial`)
-  Statistiques (:mod:`~scipy.stats`)
-  Manipulation de tableaux d'image (:mod:`~scipy.stsci`)
-  C/C++ integration (:mod:`~scipy.weave`)



Avoir sa propre librairie
-------------------------

Créez le répertoire :file:`myscripts`::

    mkdir myscripts
    
    
Créez un module ::

    $ echo "earth_radius = 6371006." > myscripts/mydata.py
 

Ajustez votre variable :envvar:`PYTHONPATH` ::

    $ setenv PYTHONPATH $PWD/myscripts

   
Testez ::

    $ python -c 'import mydata; print mydata.earth_radius`






Où trouver de l'aide ?
----------------------

Sources génériques :

- Google.
- Les documentations (générée avec sphinx) des paquets (modules).
- Les didactitiels.

Documentation importantes :

- Librairies de base : http://docs.python.org/2.7/library/index.html
- Manuel de référence de :mod:`numpy` : http://docs.scipy.org/doc/numpy/reference
- Manuel de référence de :mod:`scipy` : http://docs.scipy.org/doc/scipy/reference
- Galerie de :mod:`matplotlib` : http://matplotlib.org/gallery.html
- CDAT (5) : http://www2-pcmdi.llnl.gov/cdat/manuals/cdms5.pdf
