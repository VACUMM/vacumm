.. _user.dev:

For developpers
***************

.. highlight:: bash

Use your own version
====================

To use the library without installing it, see doc sections :ref:`user.install.install.dev` and 
:ref:`user.install.modenv.dev`.



Coding rules
============

Please follow the :ref:`appendix.conventions`.



Architecture
============

The architecture of the library must follow the spirit
of what is explained at section :ref:`user.desc.basis.arch`.
For instance, generic stuff must go in the :mod:`vacumm.misc`
module.
When a new module is specialized, it should find a place elsewhere.

.. _user.dev.doc:

Generation of the sphinx documentation
======================================

Generation of the figures of the scripts
----------------------------------------

Some of the scripts located in the :file:`scripts` directory 
create figure that are included in the documentation.
You must run them before compiling the documentation to make
sure that all needed figures are available.
For that, just run a :program:`make` in the :file:`scripts` directory:
    
    
.. code-block:: bash

    $ cd scripts
    $ make


.. note:: It is not necessary to re-generate these figures each time you compile!


Force the regeneration of the figures of the colormaps
------------------------------------------------------

Some samples of colormap of the :mod:`vacumm.misc.color` module are displayed
in the documentation. They are usually generated automatically at compile time
by an integrated sphinx extension (the file :file:`doc/pshinx/source/sphinxext/gen_cmaps`)
each time the module is modified or when no figure is present.

You can force the regeneration of theses figure for the next compilation with::
    
    $ touch lib/python/vacumm/misc/color.py


Generate missing test files
---------------------------

The test scripts (:file:`scripts/test/test_*.py`) are presented in the appendix 
(:ref:`appendix.tests`)
thanks to corresponding rst files situated in directory :file:`doc/pshinx/source/tests`.
You can update these rst files with the :program:`make.py` present in this directory::
    
    $ cd doc/pshinx/source/tests
    $ make.py # -h for help

.. note:: The title of each of these rst files is copied from the first line of the test script.

If new files are created, add them to subversion.

    $ svn add test_newtest.py

Then make a commit:
    
    $ svn ci -m 'adding new test rst files'

Compilation
-----------

The documentation is written in `rst <http://docutils.sourceforge.net/rst.html>`_ language,
and compilated with `Sphinx <http://sphinx.pocoo.org>`_ .
Source files are located in the  :file:`doc/sphinx/source` directory.
To compile it: 

.. code-block:: bash

    $ cd doc/sphinx
    $ make      # html + pdf
    $ make html # html only
    $ make pdf  # pdf only using pdflatex

The documentation is generated in directories 
:file:`doc/sphinx/build/html` and :file:`doc/sphinx/build/latex`.


Regeneration of TikZ figures 
----------------------------

This documentation contains several figures drawn with 
:program:`pdflatex` and 
`PGF/TikZ <http://pgf.sourceforge.net>`_ (logo, architecture 
of the librairie).
The advantage is being able to put figure sources on the svn server, 
so that everyone can re-generate them.

It is necessary to have a recent version of PGF/TikZ,
you can for example get here: http://www.texample.net/tikz/builds/  
(here is `a version <http://media.texample.net/pgf/builds/pgfCVS2010-09-28_TDS.zip>`_).
For installation, proceed as follows:
    
.. code-block:: bash

    $ mkdir -p ~/texmf
    $ cd ~/texmf
    $ wget http://media.texample.net/pgf/builds/pgfCVS2010-09-28_TDS.zip
    $ unzip pgfCVS2010-09-28_TDS.zip
    $ rm pgfCVS2010-09-28_TDS.zip

.. sidebar:: What is PGF/TikZ ?

    This is a library for creating high quality figures from an TeX source code.
    The best overview is provided by the site that lists examples (tutorials): 
    http://www.texample.net/tikz/examples .
    Most of them are based on a more recent than the one installed by default on a system release.
    Its on CVS *build* that the figures in this documentation are based.

Figure TikZ can now be generated with:  

.. code-block:: bash

    $ cd doc/sphinx/sources
    $ make

The latex code is then compiled, generating a pdf which is then converted 
to ppm format and png formats.



Writing tutorials
=================

He strongly suggested to developers to create tutorials on important features they develop. 
These tutorials have two interests:
    
    - They complete the documentation.
    - They can be used to perform tests,
      thanks to the :program:`check.py` script.
      
The :program:`check.py` script is located in the tutorials directory
(:file:`scripts/tutorials`).
Its use is as follows:

.. code-block:: bash

    $ check.py [options] [pattern1 [pattern2] ...]
    
It takes as argument one (or more) global patterns to list the test scripts.
The default is ``"*.py"``.
It is then possible to exclude scripts that list with the :option:`-e` option.

This script displays to the console information about the tests, 
and store all information in the :file:`check.log` file.
The logging level at the console can be modified with the :option:`-l` option.

Examples of use:
    

.. code-block:: bash

    $ check.py -e "misc.color.py" -e misc.grid.masking.* misc.*.py
    $ check.py --loglevel=debug

    
.. rubric:: Options de :program:`check.py`
    
.. program:: check.py :

.. option:: -h, --help

    Affiche l'aide.
    
.. option:: -e, --exclude

    Adds a global pattern listing scripts to exclude tests.
    
.. option:: -l, --loglevel

    Sets the level of logging to the console. 
    This can have the following values:
        
        - ``"debug"``:  Displays standard output and standard error.
        - ``"info"``: Displays the name of the successful scripts (default).
        - ``"error"``: Displays the name of the scripts that failed.


Distributing the library as a package
=====================================

It is possible to create packages typically corresponding to specific versions (for example stability). 
The procedure is as follows:   
    
.. code-block:: bash

    $ python setup.py bdist
    
This command will then create a distributable file, whose name is close to
:file:`vacumm-0.9-svn128.linux-x86_64.tar.gz`.
This file can then be placed in the files section of the gforge site of the project 
(`à cette adresse <https://forge.ifremer.fr/frs/admin/qrs.php?package=&group_id=93>`_).


