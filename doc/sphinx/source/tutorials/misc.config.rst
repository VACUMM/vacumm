.. _user.tut.misc.config:

Configuration management
========================

See also the api documentation: :mod:`~vacumm.misc.config`

Introduction
~~~~~~~~~~~~

The purpose of the :mod:`vacumm.misc.config` is to allow a script to use some 
parameters stored in configuration files. These configurations are similar to 
the 'ini' files (https://en.wikipedia.org/wiki/INI_file) but in a much more 
flexible way than what the standard module :mod:`ConfigParser` provides.

:mod:`vacumm.misc.config` rely on the :mod:`~configobj`
(http://www.voidspace.org.uk/python/configobj.html) and 
:mod:`optparse` or :mod:`argparse` modules and mixes these modules to handle configurations with the 
following concepts:
    
    #) The entry point is to define a specification file which contains sections
       and options with special values defining what kind of data will be stored 
       and retrieved (data type, default value, allowed values and documentation).
    #) Then users can provide a configuration file with all or some of the sections 
       and options (missing options will use the default values defined in the 
       specifications).
    #) Finally the configuration can also be overriden by command line options 
       using the standard modules :mod:`optparse` or :mod:`argparse`.

The configuration management is made with :class:`~vacumm.misc.config.ConfigManager` .

A hierarchical configuration file can have nested sections:

.. code-block:: ini

    option1 = value1
    
    [section1]
    
        [[section11]]
        
        option11 = value11
    
    options2 = value2

The configuration object have a similar behavior as :class:`dict`, you can 
access to the ``option2`` like this: ::
    
    >>> print cfg['section1']['section11']['option11']
    value11

.. warning::
    - Although encapsulating string values with quotes in configuration files 
      and command line options is not always required, we strongly encourage to do it.
    - Note that option2 will be bound to secion11, despite indentation, any options
      for a section that appear below a nested section will be loaded as this 
      subsection option, that's why you have to always write subsections after all
      of a section options.

Example
~~~~~~~

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

