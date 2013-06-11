.. _user.scripts.batch:

:program:`python.batch` -- Script for batch submission
======================================================

This script helps you running a script in batch mode.

.. highlight:: bash

::

    shell> python.batch -h
    shell> python.batch myscript.py --options arguments

All the options and arguments are passed to the python script.

To properly use it, the script must load the right
python environment.
Therefore, it is suggested than you make a copy of it and edit it.

.. note::

    If you use the environment modules,
    you can set the :envvar:`VACUMM_MODULES` environment variable
    to the directory of module files that contain the 
    ::file:`vacumm` module file.    
    For instance, here we make the supposition that the
    :file:`/home/me/modulefiles/vacumm` exists:

        bash> export VACUMM_MODULES=/home/me/modulefiles



Here is a view the script for batch submission:

.. literalinclude:: ../../../bin/python.batch
    :language: bash
