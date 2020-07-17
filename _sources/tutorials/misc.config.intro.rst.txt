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
