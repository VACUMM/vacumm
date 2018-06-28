.. _user.faq.config:

Configuration of the librairy
=============================


What is my configuration?
-------------------------

Type from the shell::

    $ vacum_config.py print

Where is my personnal configuration?
------------------------------------

The general config file should be there:
:file:`$HOME/.config/vacumm/vacumm.cfg`.

This file may make reference to secondary files.

You can edit it online with :func:`vacumm.config.edit_user_conf_file`:

    >>> from vacumm.config import edit_config
    >>> edit_user_conf_file()

In offline mode, use :program:`vacumm_edit_config.py`
:ref:`script <user.scripts.vacumm_edit_config>`.

I am a developer and I want to create a configuration for my module
-------------------------------------------------------------------

Suppose your module is called: :mod:`vacumm.my.module` (fichier :file:`vacumm/my/module.py`).

    #. Create the configuration file :file:`config.cfg` located in the same directory.
    #. Set your options in the ``vacumm.my.module`` section. Example:

       .. code-block:: ini

         [vacumm.my.module]
         sst_max=34

    #. Check your configuration with the :func:`~vacumm.config.get_config_value` function::

       >>> from vacumm.config import get_config_value
       >>> get_config_value('vacumm.my.module', 'sst_max')
       '34'


I am not at IFREMER and I want for instance Etopo2
--------------------------------------------------

Since its svn-583 version, the library support the automatic download with your agreement,
when it tries to access it the file.

If this is not the case, edit the your personal configuration file :file:`$HOME/.config/vacumm/vacumm.cfg`
to speficy the path to the secondary configuration file dedicated to gridded bathymetries.
This file is obviously located at the same place of the main configuration file ::
Si ce n'est pas le cas, éditer le fichier :file:`$HOME/.config/vacumm/vacumm.cfg`
For example:

.. code-block:: ini

    [vacumm.bathy.bathy]
    cfgfile_gridded=%(user_conf_dir)s/bathy.gridded.cfg


``%(user_conf_dir)s`` refers to the :file:`$HOME/.config/vacumm` directory.
You can also put an explicit path.

Then change the path to the :file:`etopo2.nc` file in :file:`bathy.gridded.cfg`:

.. code-block:: ini

    [etopo2]
    file=/path/to/etopo2.nc

