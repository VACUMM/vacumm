.. _user.faq.install:

Installations
=============


.. _user.faq.install.intel:

"This Intel <math.h> is for use with only the Intel compilers!"
---------------------------------------------------------------


Si à l'installation (compilation) de la librairie, vous rencontrez un message d'erreur du type :


.. code-block:: bash

    /appli/intel/Compiler/11.1/073/include/math.h:27:3: error: #error "This Intel <math.h> is for use with only the Intel compilers!"
    In file included from /home1/caparmor/sraynaud/soft/cdat-vacumm/include/python2.5/pyport.h:231,
                     from /home1/caparmor/sraynaud/soft/cdat-vacumm/include/python2.5/Python.h:57,
                     from build/src.linux-x86_64-2.5/fortranobject.h:7,
                     from build/src.linux-x86_64-2.5/fortranobject.c:2:
    /appli/intel/Compiler/11.1/073/include/math.h:27:3: error: 
        #error "This Intel <math.h> is for use with only the Intel compilers!"
    error: Command "gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O3 -Wall 
        -Wstrict-prototypes -fPIC -Ibuild/src.linux-x86_64-2.5 
        -I/home1/caparmor/sraynaud/soft/cdat-vacumm/lib/python2.5/site-packages/numpy/core/include 
        -I/home1/caparmor/sraynaud/soft/cdat-vacumm/include/python2.5 
        -c build/src.linux-x86_64-2.5/fortranobject.c -o build/temp.linux-x86_64-2.5/build/src.linux-x86_64-2.5/fortranobject.o" failed with exit status 1

vous avez sûrement chargé un environnement (par exemple la variable :envvar:`LD_LIBRARY_PATH`) pointant vers
des librairies Intel, alors que vous compilez avec :program:`gfortran`.

Il est préférable de nettoyer votre environnement avant la compilation.


Quelles sont les options pour la compilation ou l'installation ?
----------------------------------------------------------------

Pour la compilation :

.. code-block:: bash

    shell> python setup.py --help build
    
Par exemple :
    
.. code-block:: bash

    shell> python setup.py install --compiler=gfortran 

Pour l'installation :

.. code-block:: bash

    shell> python setup.py --help install

Par exemple :
    
.. code-block:: bash

    shell> python setup.py install --install-lib=$HOME/lib/python
    

Je suis développeur et j'ai une erreur de type ``ImportError: libgfortran.so.2``
--------------------------------------------------------------------------------

Si vous utilisez VACUMM en tant que développeur (importation direct des sources),
une erreur typique peut se produire lors du chargement ::
    
    
    >>> import vacumm.misc
    Init? False
    libgfortran.so.2: cannot open shared object file: No such file or directory
    Trying to build it...
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/__init__.py", line 9, in <module>
        import io
    File "/home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/io.py", line 3825, in <module>
        from plot import map2, _colorbar_, savefigs as Savefigs, markers as Markers, gobjs
    File "/home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/plot.py", line 82, in <module>
        from atime import mpl,time,axis_add,compress,SpecialDateFormatter
    File "/home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/atime.py", line 2535, in <module>
        import grid as G
    File "/home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/grid/__init__.py", line 11, in <module>
        import regridding
    File "/home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/grid/regridding.py", line 87, in <module>
        exec import_interp
    File "<string>", line 1, in <module>
    ImportError: libgfortran.so.2: cannot open shared object file: No such file or directory
    >>> 

Une telle erreur se produit quand il y a incompatibilité entre votre python et 
le module :mod:`vacumm.misc.grid._interp_` compilé à partir du fichier :file:`interp.f90`.
Généralement, une version différente du compilateur à été utilisée.

Dans ce cas, allez dans le répertoire du module :mod:`vacumm.misc.grid`, et procédez ainsi ::
    
    
.. code-block:: bash

    shell> cd /home1/caparmor/sraynaud/soft/src/vacumm/lib/python/vacumm/misc/grid
    shell> rm _interp_.so
    shell> make

Le module est ainsi recompilé.

