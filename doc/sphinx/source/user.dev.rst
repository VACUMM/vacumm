.. _user.dev:

For developpers
***************

.. highlight:: bash

Use your own version
====================

To use the library without installing it, see doc sections :ref:`user.install.install` and 
:ref:`user.install.modenv`.



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

Quickly: see :ref:`user.dev.doc.comp`.


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


.. _user.dev.doc.comp:

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



.. _user.dev.tut:
    
How to add a new tutorial
=========================

It is strongly suggested to developers to create tutorials on important features they develop. 
These tutorials have two interests:
    
    - They complete the documentation: it is an example and it may add an entry
      to the gallery.
    - They can be used to perform tests,
      thanks to the :program:`check.py` script.


To add a new tutorial script:
    
1. Create your script in the :file:`scripts/tutorials`.
   
   - Its name is typically compound of subnames separated with dots.
   - If you create a figure, please use the :func:`vacumm.misc.plot.savefigs`
     function to save your figure, with the ``file__`` variable as first argument:
     this function convert dots to "-" (better with Latex) and allow to save 
     also in pdf format, which good with 1D plots or pure contours::
         
         savefigs(__file__, pdf=True)
         
2. Add a reference into the documentation.
   For that, create a new rst file in the :file:`doc/sphinx/source/tutorials`
   directory which the same root name as the script:
   
   .. code-block:: rst
       
       .. _user.tut.my.tut:

            My test script
            ==============

            Explain what it does.

            .. _fig.tut.my.tut:
            .. figure:: ../../../../scripts/tutorials/my.tut.png

                Legend of the figure.


            .. literalinclude:: ../../../../scripts/tutorials/my.tut.py

   Insert a reference to this file into the doc tree.
   Then commit the new and modified files.
       
       

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

    $ check.py -e "misc.color.py" -e misc.grid.masking.* "misc.*.py"
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

.. _user.dev.test:

How to add a new test
=====================

The goal of tests is to test a small part of the library:
we are speaking of the so-called "unit tests".
Theses tests can also be viewed as good examples of coding.
A test script can be executed directly or withing the
testing framework of the library, which is located in the :file:`test` directory.
This directory contains several modules whose name is starting with ``test_``:
they declare :class:`unittest.TestCase` classes where one classes embeds several
tests as methods starting also with ``test_`` (automatically generated).
Each of these methods execute one of the test script, and then check for the
existence of the ``result`` variable to perform more tests on this script.


To add a new one:
    
1. Create such script in the :file:`scripts/test` directory,
   with the following rules:
       
   - Makes sure the name of your script is not too generic of specific.
     Try for instance to group some these scripts under the same root name.
     If you test a function of method of a well known module or class,
     use it as a root name. for instance, if test the 
     :func:`vacumm.misc.phys.units.deg2m` function, create a script called
     :file:`test_units_deg2m.py`.
   - Insert a docstring of a single line telling what the script tests.
     See examples here: :ref:`appendix.test`.
   - Make explicit imports (without \*)::
       
       from vcmq import DS
   
   - If you create one or several figure, use the follwing syntax
     to name you file::
       
       from vcmq import code_base_name
       
       # single file
       figfile = code_base_name()
       
       # multiple files
       figfile1 = code_base_name(ext='_1.png')
       figfile2 = code_base_name(ext='_2.png')
       
   - At the end of the script, you can set the ``result`` variable.
     This variable will be latter used by the testing framework to 
     perform more checks.
     It must be a list of tuples, and the first element of each is
     of one of the follwing types:
         
     - The string ``"files"``: the second element of the tuple must
       then be a single or a list of file names. 
       The testing framework will check for the existence of these files.
       Be caferul to remove them before trying to create them in the script!
     - A string starting with ``"assert"``: It is assumed to be the name of
       a method of the class :class:`unittest.TestCase`.
       The second element of the tuple is converted as a list and passed as arguments
       to the method.
     - A callable object: The second element is passed to this object.
     
     Example::
         
         result = [
            ('files', 'myfile.nc'), ('files', ['myfile1.cfg', 'myfile2.cfg']),
            ('assertEqual', (var1, 1.5)), ('assertTrue', N.ma.allclose(var1,var2)),
            (mytestfunc, (arg1, arg2))]
            
2. Create a test module in the :file:`test` directory, if it does not exist,
   with the follwing rules:
       
   - Its name follows the same naming rules has for the test script:
     it is the root part of the names of the test scripts it will test.
     Using the example above, it must be named :file:`test_units.py`,
     and contain all the test of the :mod:`vacumm.misc.phys.units` module.
   - Its contents must be close to the following example::
       
        from utils import *


        class TestSequenceFunctions(VCTestCase):

            for test_name in [
                'test_units_deg2m',
                'test_units_m2deg',
                ]:
                exec(method_template.format(test_name))
                
        if __name__ == '__main__':
            unittest.main()
                
     The header of the :class:`TestSequenceFunctions` contains a loop
     on all the test scripts names (without suffix) that will be tested.
     If the module already exists, just add your new test script to this loop. 
      
3. If the test module does not exist, add an entry the ``TEST`` variable
  of the :file:`test/Makefile` file.
4. Update the documentation:
   
   .. code-block:: bash
   
        $ cd doc/sphinx/source/tests
        $ make.py
        
   Run the :program:`svn` commands that are displayed on the output of the 
   :program:`make.py` script. 
   If the test module was not existing, insert a reference to it in the :file:`index.rst`.
   For example:
   
   .. code-block:: rst
   
        Units
        -----

        Testing module :mod:`~vacumm.misc.phys.units`.

        .. toctree::
            :glob:

            test_units_*
    
   Then commit this file too.

Tagging of versions
===================

See also: http://en.wikipedia.org/wiki/Software_versioning

The version number must be defined in a unique location.

Vacumm set the version number into its python root package, :attr:`vacumm.__version__` module attribute (:file:`vacumm/__init__.py`)

This __version__ string must use the following conventions:
    - at least three numerical parts, separated with dots ('.'): ``__version__ = 'X.Y.Z'``
    - these three numerical parts have the following meaning:
        - X is the major version number, incremented when there is a significant change of the software (when reaching a certain maturity or when reorganization or when there are API breakage).
        - Y is the minor version number, incremented when new features have been added to the software.
        - Z is the revision number, incremented when there are some fixes to the corresponding X.Y version
    - an extra part can be present to indicate a release status: ``__version__ = '1.0.0.alpha'`` (or ``'1.0.0.beta'``, ``'1.0.0.rc1'``, ...)

When a version is ready to be tagged, proceed with the following steps:
    - check or set :attr:`vacumm.__version__`,
    - make a svn copy (which in fact only put a reference to a revision number)::
        
        $ svn copy https://gforge.ifremer.fr/svn/vacumm/trunk \
            https://gforge.ifremer.fr/svn/vacumm/tags/vacumm-X.Y.Z

Using a tag named with ``'vacumm'`` and the version string ``'X.Y.Z'`` together is a preferred way since a checkout will then result in an explicit :file:`vacumm-X.Y.Z` directory (instead of a lonely :file:`X.Y.Z`).

After a tag have been created, you can then change :attr:`vacumm.__version__` from, for example ``1.0.0`` to ``1.1.0.alpha`` to mark a difference for the people who uses the trunk.
However these users must be aware of what they're doing and then know if the ``1.0.0`` version they're using is the trunk or the real 1.0.0 version.
The most important thing to do is to correctly setup the version number before creating a new tag.

If you created the tag too quickly and need to correct something, then
    - fix it
    - remove the created faulty tag::
        
        $ svn rm https://gforge.ifremer.fr/svn/vacumm/tags/vacumm-X.Y.Z
        
    - re-create the expected tag (``svn copy``, see above).

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


