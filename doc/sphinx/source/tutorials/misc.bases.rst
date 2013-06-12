.. _user.tut.misc.bases:

Base classes and utilities
==========================

See also the api documentations:
    - :mod:`~vacumm.misc.bases`
    - :mod:`~vacumm.misc.config`
    - :mod:`~vacumm.misc.log`

Introduction
------------

As the api documentation states, the :class:`vacumm.misc.bases.Object` provides
convenient features for logging and configuration management.

This tutorial aims to show you these features and how to customize though most
users will only need to use them as such.

In this tutorial we will define a 'MyObject' class which inherit from 
:class:`vacumm.misc.bases.Object` as you will do in your code to take advantage
of this class.

One of the mecanism of this class is to allow developpers to execute things when
the class is defined (when its module is imported). The goal of these things may
be various but the default goal is to allow a class to be parametrized using a 
configuration file.

If you also need to do things at class initialization, below is the way you may
redefine the :meth:`vacumm.misc.bases.Object.init_class` method:

.. literalinclude:: ../../../../scripts/tutorials/misc.bases.py
    :language: python

.. program-output:: ../../../scripts/tutorials/misc.bases.py


Logging features
----------------

The :class:`vacumm.misc.bases.Object` and its instances have logging capabilities
through two :class:`vacumm.misc.log.Logger` objects, one used when logging methods
are called from the class and another one dedicated to each instance of this class.

The following example introduce the basic usage of these features.

.. literalinclude:: ../../../../scripts/tutorials/misc.bases.logging.py
    :language: python

.. program-output:: ../../../scripts/tutorials/misc.bases.logging.py


Configuration features
----------------------

:class:`vacumm.misc.bases.Object` may use configurations features provided by the
:mod:`~vacumm.misc.config` module.

When your module containing the :class:`vacumm.misc.bases.Object` subclass is 
loaded, a specification file with the module name with the '.py' extension replaced
with the '.ini' extension is looked up and the default configuration is loaded
from the specification defaults.

The way you may change your class's specification file is to redefine the
:meth:`vacumm.misc.bases.Object.get_config_spec_file` method, the section name
defaults to the class name, you may also change it by redefining the
:meth:`vacumm.misc.bases.Object.get_config_section_name` method.

The :meth:`vacumm.misc.bases.Object.get_default_config` method will return the
default configuration mentionned above.

The :meth:`vacumm.misc.bases.Object.apply_config` method will be called each time a
configuration is loaded using the :meth:`vacumm.misc.bases.Object.load_config` method.

Logging configuration related to your class may be customized usgin a nested
section named 'Logger'.

Specification file:

.. literalinclude:: ../../../../scripts/tutorials/misc.bases.config.ini
    :language: ini

Configuration file:

.. literalinclude:: ../../../../scripts/tutorials/misc.bases.config.cfg
    :language: cfg

Code using these specification and configuration:

.. literalinclude:: ../../../../scripts/tutorials/misc.bases.config.py
    :language: python

The outputs of this code are:

.. program-output:: ../../../scripts/tutorials/misc.bases.config.py

Note that MyObject attributes and logger have been changed by the loaded configuration.


Debugging features
------------------

The following example shows some method calls you may use when writing your code
to trace or debug it.

.. literalinclude:: ../../../../scripts/tutorials/misc.bases.debug.py
    :language: python

.. program-output:: ../../../scripts/tutorials/misc.bases.debug.py


