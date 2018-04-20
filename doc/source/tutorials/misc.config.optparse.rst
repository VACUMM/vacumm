.. _user.tut.misc.config.optparse:

Example with :mod:`optparse`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example we'll show you how to use the :class:`~vacumm.misc.config.ConfigManager`
to handle defaults, command line arguments and a configuration file.
The command line arguments may be used through :mod:`optparse` or :mod:`argparse` module.

Specification file:

.. literalinclude:: ../../../../scripts/tutorials/misc.config.ini
    :language: ini

Configuration file:

.. literalinclude:: ../../../../scripts/tutorials/misc.config.cfg
    :language: cfg

Sample code using these specification and configuration:

.. literalinclude:: ../../../../scripts/tutorials/misc.config.py
    :language: python

The help generated is:

.. command-output:: ../../../scripts/tutorials/misc.config.py --help

The outputs of this code are:

.. command-output:: ../../../scripts/tutorials/misc.config.py

The outputs of this code with configuration setup with command line options:

.. command-output:: ../../../scripts/tutorials/misc.config.py --scalars-string="foo, 'hello world'" --lists-floats="42,3.14"

The outputs of this code with configuration setup with command line options and a configuration file :

.. command-output:: ../../../scripts/tutorials/misc.config.py --scalars-string="foo, 'hello world'" --lists-floats="42,3.14" --cfgfile=tutorials/python/misc.config.cfg
