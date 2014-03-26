.. _user.faq.install:

Installations
=============


.. _user.faq.install.intel:

"This Intel <math.h> is for use with only the Intel compilers!"
---------------------------------------------------------------


During the installation (compilation) of the library, if you encounter an error message like:


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

you surely loaded environment (eg variable :envvar:`LD_LIBRARY_PATH`) 
pointing to Intel libraries while you compile with :program:`gfortran`.

It is best to clean your environment before compilation.

What are the options for compilation or installation? 
-----------------------------------------------------

For compilation:

.. code-block:: bash

    shell> python setup.py --help build
    
For example:
    
.. code-block:: bash

    shell> python setup.py install --compiler=gfortran 

For installation:

.. code-block:: bash

    shell> python setup.py --help install

For example:
    
.. code-block:: bash

    shell> python setup.py install --install-lib=$HOME/lib/python
    

I am a developer and I have an error like ``ImportError: libgfortran.so.2``
---------------------------------------------------------------------------

If you use vacumm as a developer (direct import sources), a typical error can occur when loading::
    
    
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

This error occurs when there is an incompatibility between your python and module
:mod:`vacumm.misc.grid._interp_` compiled from :file:`interp.f90` file.
Generally, a different version of the compiler has been used.
 
In this case, go to the module directory:mod:`vacumm.misc.grid`, and proceed as:    
    
.. code-block:: bash

    $ cd lib/python/vacumm/misc/grid
    $ rm _interp_.so
    $ make

The should be recompiled.

