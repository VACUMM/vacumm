.. _user.install.config:

Library configuration
=====================

.. sidebar:: What is my configuration?

    If your library is not too old, type:
        
        >>> from vacumm.config import print_config
        >>> print_config()
    
    Or directly in the shell:
        
    .. code-block:: bash
    
        vacumm_print_config.py
        
    See :ref:`user.scripts.vacumm_print_config`.

Purpose
-------

Some modules need for instance to access external files (such as a bathymetric file) 
that are not included directly in the library because they are too bulky.
To make the installation generic and portable on different networks, 
it is necessary not to put these paths hardcoded in the code:
so it uses a system configuration. 
There is thus a default configuration (that of IFREMER network), 
and it is possible to alter the configuration 
for example to adapt to his or her needs or network. 
For example, Actimar has its own configuration.
And each user can also have its configuration 
which affects partly these parameters.

The concept of configuration is at two levels: 

- At the user's level who can change its personal configuration file. 
- At the developer's leve who will appeal to the API (:mod:`vacumm.config`)
  to define and use the configuration in its code.


Operating
---------

A number of modules have a default configuration in the form of a 
:file:`config.cfg` file in the same directory as the module 
or in the root directory of the library. 
The module configuration is then accessible at the section of the same name.
For example, if the :mod:`vacumm.toto.tutu` module is configurable,
its configuration file may be :file:`vacumm/toto/config.cfg`, with the content:

.. code-block:: ini

    [vacumm.toto.tutu]
    masection = mavaleur
    
    
In addition, the configuration files of support the expansion of variables in two ways:
    
    #. An options inserted in a value is automatically replaced 
       by its value (it is "extended") when it is present in the classical
       form ``%(option)s``. For example:
           
       .. code-block:: ini
       
           [vacumm.toto.tutu]
           option1 = foo
           option2 = %(option1)bar
           
       ``option2`` takes the value ``foobar``.
       There are several directory names that are currently systematically extended::
           
           - ``lib_dir`` : Root directory of the library  
             (see :func:`~vacumm.config.get_lib_dir`).
           - ``data_dir`` : Directory that stores the data used by the library and tutorials
             (see :func:`~vacumm.config.get_data_dir`).
           - ``tut_dir`` : Directory containing tutorial scripts
             (see :func:`~vacumm.config.get_tut_dir`).
           - ``dist_dir`` : Directory of packet distribution for developers 
             (see :func:`~vacumm.config.get_dist_dir`).
           - ``conf_dir`` : General directory containing the configuration files
             (see :func:`~vacumm.config.get_conf_dir`).
           - ``user_conf_dir`` : Directory containing user configuration files
             (see :func:`~vacumm.config.get_user_conf_dir`).
           - ``mod_dir`` : Directory of the module 
             (see :func:`~vacumm.config.get_mod_dir`).
           
    #. The path ``~`` is automatically replaced by the user home, 
       and environment variables that have the form ``$HOME`` ou ``${HOME}``
       are replaced by their value.
       
       .. note:: Developers that need to manage configuration can refer 
        to the documentation of the module :mod:`vacumm.config`.

Configuration values are accessible with the function 
:func:`vacumm.config.get_config_value`.


Default configuration
---------------------

The general default configuration (loaded with :func:`~vacumm.config.get_default_config`) 
is following: 

.. literalinclude:: vacumm.cfg
    :language: ini


User configuration 
------------------

It is possible to alter for a user the default configuration of one or more modules. 
The best way is to edit the file :file:`$HOME/.config/vacumm/vacumm.cfg`.
To quickly **edit** this file, use this script: :ref:`user.scripts.vacumm_edit_config`.

.. obsol
    Il existe en outre d'autres possibilités, mais qui doivent être considérées comme obsolètes.
    Les voici par ordre de priorité (see :func:`~vacumm.config.get_config_files`) :
        
        #. Créer un fichier :file:`vacumm.cfg` dans votre **répertoire courant** d'utilisation.
        #. Installer un fichier :file:`site.cfg` ou :file:`vacumm.cfg` dans le **répertoire
        principal de la librairie téléchargée** si vous l'utilisez directement sans l'installer
        (vous êtes donc un développeur !).
        #. Installer un fichier :file:`site.cfg` dans le **répertoire du module** que vous ciblez.
        #. Installer un fichier :file:`site.cfg` dans le **répertoire racine de la librairie**.
    
If you personally install the library,
you can configure it for the users using option :option:`--cfgfiles` of the :file:`setup.py` 
installation script
(see :ref:`user.install.install.config`).

Finally, you can alter this configuration online through function 
:func:`vacumm.config.set_config_value`.

