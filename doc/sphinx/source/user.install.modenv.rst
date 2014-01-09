.. _user.install.modenv:

Using environment modules
=========================

Introduction
------------

As it has been metionned in both previous sections, the access to the Python included in UV-CDAT needs to modify the environement variable :envvar:`PATH`, and the access to the library needs to modify the environment variable :envvar:`PYTHONPATH` if it has been installed in a different directory.

Knowing that it is possible manage and use multiple versions of CDAT and the VACUMM library, 
as needed (in an operational framework, for a given project,
for the standard users, or those users who want the latest features),
it may become difficult to administrator to ensure that users can easily access the correct versions;
as it would be unpleasant for users to have to change themselves all environment variables (ie to know which variables and know all the paths).
An easy way to handle this problem is the use of `environment modules <http://modules.sourceforge.net>`_.
This approach allows the administrator to manage multiple versions
and their dependencies transparently for the user,
who just need to insert two rows in its :file:`.cshrc` (ou :file:`.bashrc`) file 
and know :command:`module` utility to load the environment he wants.

The documentation of the utility is located here: 

    - http://modules.sourceforge.net/
    - http://modules.sourceforge.net/man/module.html
    - http://modules.sourceforge.net/man/modulefile.html
    - http://sourceforge.net/p/modules/wiki/FAQ/

Default module
--------------

The module utility makes possible to use several module version for the same software.

If we have for example:


.. code-block:: bash

    shell> ls /my/modules/vacumm
    0.1
    1.0
    trunk

To set a version as the default, you must be create in this directory a
:file:`.version` file containing:

.. code-block:: tcl

    #%Module
    set ModulesVersion "1.0"

This file indicates which module file will be loaded by default
when only the directory is specified.
The command ``module load vacumm`` 
then load the file module :file:`/my/modules/vacumm/1.0`.

The administrator side
----------------------

It is first necessary that the administrator configures the use of environment modules.

The first step is to create a directory dedicated to the module files ,
and subdirectories dedicated to versions of CDAT, vacumm and their dependencies:
    
.. code-block:: bash

    shell> mkdir /my/modules
    shell> mkdir /my/modules/cdat
    shell> mkdir /my/modules/vacumm
   
Then, the module files for CDAT and its dependencies (voir :ref:`user.install.prereq`)
should be created. 
For example, if SCRIP is a dependency and is installed in the directory :file:`/my/soft/scrip`, 
you must create the module file :file:`/my/modules/scrip` with the following content:
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "SCRIP regridder"
    }
    
    prepend-path PATH /my/soft/scrip
    
Python module :mod:`~mpl_toolkits.basemap` uses the GEOS library. 
You must create the module file :file:`/my/modules/geos/3.2` :
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "GEOS library"
    }
    
    set version 3.2
    setenv GEOS_LIB /my/soft/geos/$version/lib
    prepend-path LD_LIBRARY_PATH /my/soft/geos/$version/lib
    

**Note**: geos is now included with cdat> = 6.0, it is not necessary to have a system and a module for geos apart.

The module file of CDAT :file:`/my/modules/cdat/5.2 can then be created taking into account dependencies and conflicts with other versions:

.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "SCRIP regridder"
    }
    
    conlict cdat
    module load scrip
    module load geos/3.2
    
    set version 5.2
    set base_path /my/soft/cdat
    prepend-path PATH $base_path/$version/bin
    prepend-path LD_LIBRARY_PATH  $base_path/$version/lib
    prepend-path LD_LIBRARY_PATH  $base_path/$version/Externals/lib
    prepend-path C_INCLUDE_PATH $base_path/$version/include
 
If the default version of vacumm is installed in the python provided with CDAT,
the corresponding module file :file:`/my/modules/vacumm/1.0` then contains:
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "VACUMM python library"
    }
    
    conlict vacumm
    module load cdat

Here we only load CDAT as vacumm is installed in the default directory for installing python CDAT.

If another version of VACUMM exists 
(installed with a --prefix equal to :file:`/my/soft/vacumm-recent`), 
and using a different version of CDAT, we could have
 the file  :file:`/my/modules/vacumm/recent` :
    
.. code-block:: tcl

    #%Module
    proc ModulesHelp { } {
        puts stderr "VACUMM python library"
    }
    
    conlict vacumm
    module load cdat/6.0
    prepend-path PATH /my/soft/vacumm-recent/bin
    prepend-path PYTHONPATH /my/soft/vacumm-recent/lib/python2.7/site-packages
    
with cdat/6.0.0 module configured using this other CDAT.

The user side
-------------
    
The file :file:`$HOME/.cshrc` (ou :file:`$HOME/.bashrc`) 
will typically contain the following lines:
    
.. code-block:: bash

    # Modules
    # - initialisation
    source /usr/share/modules/init/csh # ou bash
    # - ajout du repertoire des modules
    module use /my/modules

.. warning::
    
    One should however be careful if the module is loaded cdat in a profile file as
    :file:`$HOME/.bashrc` or :file:`$HOME/.cshrc`. 
    Since these files are loaded at logon (graphic or not), 
    side effects may occur because the version of python, 
    and packages that are installed are no longer those of 
    the system and some software or even the complete graphical session 
    may be unusable.
    If this is the case, it can be solved by loading the CDAT module manually in your terminal.

The user then has an access to all environment modules defined in the previous section.
The see available modules:
    
.. code-block:: bash

    shell> module avail
     
To load the VACUMM library:
    
.. code-block:: bash

    shell> module load vacumm
    
.. note::

    It is possible to automatically load the library by inserting  ``module load vacumm``
    directly in the file  :file:`$HOME/.cshrc` after the line ``module use ...``
    
Then we check:
    
.. code-block:: bash

    shell> module list
    Currently Loaded Modulefiles:
      1) scrip      2) geos/3.2    3) cdat/5.2     4) vacumm/recent
    shell> echo $PATH
    /my/soft/cdat/5.2/bin:/my/soft/scrip:... #etc
    shell> which python
    /my/soft/cdat/5.2/bin/python
    shell> python -c "import vacumm"
    
To change the version of VACUMM in order to use, for instance,
an installation of the trunk branch that is regularly updated:
    
.. code-block:: bash

    shell> module switch vacumm vacumm/trunk
    
.. _user.install.modenv.dev:

The developer side
------------------

The checkout (svn) of VACUMM contains a module file  :file:`etc/modulefiles/vacumm`
that can be used to directly exploit the sources:

.. code-block:: bash

    shell> module use /path/to/my/vacumm/etc/modulefiles
    shell> module load vacumm
    
Not enter into conflict with (hide) another set of modules, 
you can create your personal space for modules 
and link this module file to this space:

.. code-block:: bash

    shell> mkdir -p ~/etc/modulefiles/vacumm
    shell> ln -s /path/to/my/vacumm/etc/modulefiles ~/etc/modulefiles/vacumm/dev

This module requires a module named cdat, and an error will be displayed if it is not found.

.. warning::

    If you need to specify a particular cdat module for your development version, do not to modify 
    :file:`etc/modulefiles/vacumm` (to not commit it), 
    rather create a cdat directory in your own personal module space by 
    creating a module that loads the right cdat,
    and optionally using the default module method described above.

.. code-block:: bash

    shell> mkdir -p ~/etc/modulefiles/cdat
    shell> vi ~/etc/modulefiles/cdat/dev

Learn this new module :file:`cdat/dev` with the following content:

.. code-block:: tcl

    #%Module
    module load cdat/6.0


You can now load your development version:

.. code-block:: bash

    shell> module load vacumm/dev




Et voilà, that's all folks !
    
