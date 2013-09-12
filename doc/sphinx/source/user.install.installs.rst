.. _user.install.install:

Installations
=============


For all users
-------------

This installation puts the library within the python tree.
The :program:`python` excutable must be the one provided by UV-CDAT
(see :ref:`user.install.prereq`).

One compile and install VACUMM from the main directory of the package with:

.. code-block:: bash

    shell> python setup.py install

For more options:
    
.. code-block:: bash

    shell> python setup.py --help


For our own usage
-----------------

In this case, the library is installed in the user's directory.

First, compile:

.. code-block:: bash

    shell> python setup.py build
    
Then, install in the default user's directory 
(typically  :file:`/home/<user>/.local/lib/python2.{x}/site-packages`) :
    
.. code-block:: bash

    shell> python setup.py install --user

If the library must be installed in a special directory, use:

.. code-block:: bash

    shell> python setup.py install --prefix=/my/special/dir

Then, set your environment variables to have acces both to the python library
and to the executables:

.. code-block:: bash

    shell> export PATH=/my/special/dir/bin
    shell> export PYTHONPATH=/my/special/dir/lib/python{x.x}/site-packages  # {x.x} refers to the python version


.. _user.install.install.dev:
    
Local installation for a developper
-----------------------------------

VACUMM must be used directly after an untar or a checkout by
setting the variables:

.. code-block:: bash

    shell> export PATH=~/path/to/vacumm-trunk/bin
    shell> export PYTHONPATH=~/path/to/vacumm-trunk/lib/python
    
The fortran extensions will be compilated automatically in place
during the importation of the associated module
(:mod:`~vacumm.misc.grid.regridding`).

You can also use the environment module provided
with the package as presented at section :ref:`user.install.modenv.dev`.



.. _user.install.install.config:
    
User configuration of the modules
---------------------------------

Some of the modules can be configured to change their default behaviour.
Configurations store for instance default paths.

The library is configured by default for a use on the supercomputer
from IFREMER (CAPARMOR).
If you are on your own system or you want to change your configuration,
please check the documentation section :ref:`user.install.config`.

During the installation process using the :program:`setup.py`,
you can provide a general configuration file and secondary configuration
files with :option:`--cfgfiles` option (comma separated).
All these files will be installed in the :file:`vacumm-config` directory.

For instance, the configuration of the :mod:`vacumm.bathy.bathy`
module makes a reference to a secondary configuration file
refered in the config section ``[vacumm.bathy.bathy]`` with the
key ``cfgfile_gridded``.
To alter this configuration, proceed in this way:
    
    #. Specify the name of the secondary config file in the main config file by prefixing it with ``%(conf_dir)s``, which is the directory where config files will be installed (see :ref:`user.install.config`):
        
       .. code-block:: ini
       
           [vacumm.bathy.bathy]
           cfgfile_gridded=%(conf_dir)s/bathy.gridded.cfg
           
    #. Then specify the list of config files as comandline option:
        
        
       .. code-block:: bash
       
            shell> python setup.py install --cfgfiles=myconfig.cfg,path/to/bathy.gridded.cfg
            
 
Check your installation
-----------------------

Run:
    
.. code-block:: bash

    shell> vacumm_print_config.py

    
