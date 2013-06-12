.. _user.tut.misc.config:

Configuration management
========================

See also the api documentation: :mod:`~vacumm.misc.config`

Introduction
------------

The purpose of the :mod:`vacumm.misc.config` is to allow a script to use some 
parameters stored in configuration files. These configurations are similar to 
the 'ini' files (https://en.wikipedia.org/wiki/INI_file) but in a much more 
flexible way than what the standard module :mod:`ConfigParser` provides.

:mod:`vacumm.misc.config` rely on the :mod:`configobj`
(http://www.voidspace.org.uk/python/configobj.html) and 
:mod:`optparse` modules and mixes these module to handle configurations with the 
following concepts:
    
    #) The entry point is to define a specification file which contains sections
       and options with special values defining what kind of data will be stored 
       and retrieved (data type, default value, allowed values and documentation).
    #) Then users can provide a configuration file with all or some of the sections 
       and options (missing options will use the default values defined in the 
       specifications).
    #) Finally the configuration can also be overriden by command line options.
       

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
    - Although encapsulating string values with quotes is not always required, we
      strongly encourage to do it.
    - Note that option2 will be bound to secion11, despite indentation, any options
      for a section that appear below a nested section will be loaded as this 
      subsection option, that's why you have to always write subsections after all
      of a section options.

Example
-------

Specification file:

.. literalinclude:: ../../../../scripts/tutorials/misc.config.ini
    :language: ini

Configuration file:

.. literalinclude:: ../../../../scripts/tutorials/misc.config.cfg
    :language: cfg

Sample code using these specification and configuration:

.. literalinclude:: ../../../../scripts/tutorials/misc.config.py
    :language: python

The outputs of this code are:

.. program-output:: ../../../scripts/tutorials/misc.config.py

The help generated is:

.. program-output:: ../../../scripts/tutorials/misc.config.py --help

The outputs of this code with configuration setup with command line options:

.. program-output:: ../../../scripts/tutorials/misc.config.py --scalars-string="foo, 'hello world'" --lists-floats="42,3.14"

The outputs of this code with configuration setup with command line options and a configuration file :

.. program-output:: ../../../scripts/tutorials/misc.config.py --scalars-string="foo, 'hello world'" --lists-floats="42,3.14" --cfgfile=tutorials/python/misc.config.cfg




L'exemple
---------

Pour ce tutoriel, nous créons d'abord le **fichier de spécifications** :file:`config.ini`,
définissant la nature de toutes les options, leur description
et leur valeur par défaut :

.. literalinclude:: ../../../../scripts/tutorials/config.ini
    :language: cfg

.. note:: Ce fichier doit être tout le temps accessible par le script.

Nous créons aussi le **fichier de configuration personnel** :file:`config.cfg` 
qui sera lu
à chaque exécution du script afin de définir nos propres valeurs
par défaut :

.. literalinclude:: ../../../../scripts/tutorials/config.cfg
    :language: ini

Et voici le script executable (:file:`config-plot.py`) qui fera le travail :
    
.. literalinclude:: ../../../../scripts/tutorials/misc.config.plot.py

Voici l'**aide courte** obtenue avec l'option :option:`--help` : ::
    
    Usage: config.plot.py [options] ncfile

    Script to plot a 2D netcdf variable

    Options:
      -h, --help            show a reduced help
      --long-help           show an extended help

Et voici l'**aide longue** avec l'option :option:`--long-help` : ::
    
    Usage: config.plot.py [options] ncfile

    Script to plot a 2D netcdf variable

    Options:
      -h, --help            show a reduced help
      --long-help           show an extended help

      Global options:
        General configuration options

        --cfgfile=CFGFILE   Configuration file [default: "config.cfg"]
        --var=VAR           Name of the netcdf variable [ex: "TEMP"]
        --title=TITLE       Title of the plot [ex: "%(long_name)s [TEMP]"]

      Zoom:
         zoom extensions in space

        --zoom-lon-min=ZOOM_LON_MIN
                            Min logitude [ex: "-5.5"]
        --zoom-lon-max=ZOOM_LON_MAX
                            Max longitude [ex: "-4.0"]
        --zoom-lat-min=ZOOM_LAT_MIN
                            Min latitude [ex: "47.0"]
        --zoom-lat-max=ZOOM_LAT_MAX
                            Max latitude [ex: "49.0"]


Nous exécutons alors le script avec un changement d'une des options (la latitude maximale) : 

.. code-block:: bash
    
    prompt> config.plot.py --zoom-lat-min=48 ../../../../../data/mars3d.xy.nc
    Current configuration: {'var': 'temp', 'title': '%(long_name)s [temp]', 'zoom': {'lon': {'min': -5.5, 'max': -4.0}, 'lat': {'min': 48.0, 'max': 49.0}}}
    Wrote to ./config-plot.png

Au final, le fichier :file:`config.cfg` a modifié le nom de la variable netcdf, et la ligne de commande a modifié la latitude maximale, le reste des options ayant pris les valeurs par défaut du fichier de spécifications.

La figure ainsi créée :
    
.. figure:: ../../../../scripts/tutorials/misc-config-plot.png

    Figure créé par le script :file:`misc.config.plot.py`.
