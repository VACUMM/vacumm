.. _user.faq.config:

Configurationof the librairy
=============================


What is my configuration?
-------------------------

Type :

    >>> from vacumm.config import print_config
    >>> print_config

Where is my personnal configuration?
------------------------------------

The general config file should be there: 
:file:`$HOME/.config/vacumm/vacumm.cfg`.

This file may make reference to secondary files.

You can edit it with :func:`vacumm.config.edit_config`:
    
    >>> from vacumm.config import edit_config
    >>> edit_config()

Je suis développeur et je veux créer une configuration pour mon module
----------------------------------------------------------------------

Supposons que votre module s'appelle : :mod:`vacumm.my.module` (fichier :file:`vacumm/my/module.py`).

    #. Créez le fichier de configuration :file:`config.cfg` situé dans le même répertoire.
    #. Définissez vos options dans la section ``vacumm.my.module``. Exemple :
        
       .. code-block:: ini

         [vacumm.my.module]
         sst_max=34
        
    #. Vérifier votre configuration avec :func:`~vacumm.config.get_config_value` :
        
       >>> from vacumm.config import get_config_value
       >>> get_config_value('vacumm.my.module', 'sst_max')
       '34'


Je ne suis pas à l'IFREMER et je veux Etopo2
--------------------------------------------

À partir de la version svn-583, la librairie doit prendre
en charge le téléchargement du fichier avec votre accord,
au moment où elle cherche à y accéder.

Si ce n'est pas le cas, éditer le fichier :file:`$HOME/.config/vacumm/vacumm.cfg`
pour mentionner le fichier de configuration secondaire dédié au bathymétries grillées, situé au même endroit.
Par exemple :
    
.. code-block:: ini

    [vacumm.bathy.bathy]
    cfgfile_gridded=%(user_conf_dir)s/bathy.gridded.cfg
    

``%(user_conf_dir)s`` fait référence au répertoire :file:`$HOME/.config/vacumm`.
Vous pouvez aussi mettre un chemin d'accès explicite.

Changez ensuite le chemin d'accès au fichier :file:`etopo2.nc` dans :file:`bathy.gridded.cfg`:
    

.. code-block:: ini

    [etopo2]
    file=/path/to/etopo2.nc

