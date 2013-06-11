.. _user.scripts.common_options:
    
Options communes
================

Voici une description des options que l'on peut retrouver dans divers scripts utilisateurs:

.. option::  --cfgfile

    Affiche l'usage et la liste des options possibles.

.. option::  --cfgfile

    Permet de spécifier un fichier de configuration.

.. option:: -v
.. option:: --variables

    Specifie les variables à traiter.
        - Lorsque le nommage des variables n'est pas homogène entre les fichiers, l'utilisateur peut
          spécifier des alias qui seront recherchés dans l'ordre d'appartition.
          Ainsi l'option **"-v temp,temperature"** permet d'indiquer que la variable à rechercher est **"temp"**
          et que si celle ci n'est pas trouvée, on cherchera la variable nommée **"temperature"**.
        - Cette option peut être répétée plusieurs fois afin de traiter plusieurs variables
          (ex: script.py -v temp,temperature -v sal,salinity ...)

.. option:: -o
.. option:: --output

    Motif de fichier de sortie. Les tracés sont enregistrés dans des fichiers
    dont le nom est construit à partir de ce motif. La valeur par défaut indique les substitutions
    possible. Par exemple, le script layer.py a pour motif:
      
      **"layer-%(var)s-%(depth)s-%(tmin)s-%(tmax)s.png"**
      
      Les substitutions seront:
        - **var**: la variable traitée
        - **depth**: la profondeur de la couche
        - **tmin**: la valeur minimale de temps couverte par le tracé
        - **tmax**: la valeur maximale de temps couverte par le tracé

.. option:: --show:

    Le script affiche les tracés au fur et à mesure


