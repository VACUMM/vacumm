.. highlight:: bash

.. _user.install.prereq:

Prerequisites
=============


.. note:: It is highly suggested to install UV-CDAT and other vacumm
    dependencies using `conda <http://conda.pydata.org/docs/index.html>`_.

.. _user.prereq.list:

List of packages
-----------------

UV-CDAT
^^^^^^^

`UV-CDAT <http://uv-cdat.llnl.gov/>`_ is a meta-package containing:

    - a version of python;
    - standard and specialized python packages
      like `numpy <http://docs.scipy.org/doc/numpy/reference>`_,
      `scipy <http://docs.scipy.org/doc/scipy/reference>`_ and
      `matplotlib <http://matplotlib.org>`_;
    - python packages made by UV-CDAT providers,
      and specialized for manipulation and visualization of ocean/atmosphere data;
    - numerous other libraries and tools.


Other packages
^^^^^^^^^^^^^^

.. rubric:: Python

When you have installed UV-CDAT, here is a list of necessary or
useful  packages that can be install using :program:`pip` or
:program:`easy_install`:

.. _user.prereq.list.others.table:
.. list-table:: *Necessary or useful packages*
   :widths: 17 9 30
   :header-rows: 1

   * - Package
     - Necessary?
     - Description
   * - :pypi:`configobj (4.7.2)`
     - Yes
     - Manipulate advanced configuration files
   * - :pypi:`pytz`
     - Only for time zone support
     - Work with time zones
   * - :pypi:`paramiko`
     - Only for :mod:`~vacumm.misc.remote`
     - Use SSH2 protocol.
   * - :pypi:`seawater`
     - Yes for advanced diags.
     - Sea water properties
   * - :pypi:`xlwt`
     - Only for :func:`~vacumm.misc.misc.xls_style`.
     - Write excel files (better: :pypi:`xlutils`)
   * - :pypi:`cmocean`
     - Only if you want cmocean colormaps
     - Make these colormaps available like ``cmocean_thermal``
   * - :pypi:`PIL (1.1.7)`
     - Only sometimes for :func:`~vacumm.misc.plot.add_logo`.
     - Manipulate images
   * - :pypi:`sphinx_fortran (1.0)`
     - Doc
     - Extension to :pypi:`sphinx`: Documenting fortran codes
   * - :pypi:`sphinxcontrib-cheeseshop (0.2)`
     - Doc
     - Extension to :pypi:`sphinx`: Linking to Cheese Shop (Python Package Index) packages
   * - :pypi:`sphinxcontrib-ansi (0.6)`
     - Doc
     - Extension to :pypi:`sphinx`: Parse ANSI control sequences
   * - :pypi:`sphinxcontrib-programoutput (0.8)`
     - Doc
     - Extension to :pypi:`sphinx`: Include program output


.. rubric:: Documentation generation

These utilities are needed to fully compile the documentation.

`Graphviz <http://www.graphviz.org>`_
    Used to create hierarchical diagrams of class inheritance
    during the generation of the documentation by the
    sphinx extension :mod:`sphinx.ext.graphviz`.
    See for instance module :mod:`~vacumm.misc.core_plot`.
    The program :program:`dot` may be also needed.
    You can install both with the :pypi:`pydot` package.

`dvipng <http://savannah.nongnu.org/projects/dvipng>`_
    Used to compile latex formula of the documentation.

`enumitem <https://www.ctan.org/pkg/enumitem>`_
    This package is needed to compile the PDF version
    of the documentation.

.. _user.install.prereq.howto:

Install UV-CDAT
---------------

From conda
^^^^^^^^^^

Using the `conda installer <http://conda.pydata.org/docs/index.html>`_::

    $ conda install -c uvcdat uvcdat

From sources
^^^^^^^^^^^^

To install UV-CDAT, follow the official `instructions <http://uv-cdat.llnl.gov/install>`_.
It can be installed on both linux and mac, by compiling sources or using
availables binaries.

.. note::

    If you compile it from sources, you need `cmake <http://www.cmake.org>`_ and
    `git <http://git-scm.com>`_, and it is highly suggested to have your own
    version of `Qt4 <http://qt-project.org>`_.
    Packets are generally availables on all plateforms.

If you use UV-CDAT in operational jobs, you should install a separate version.
One way to manage several versions is to use environment modules
(see  :ref:`user.install.modenv`).


Setup the environment
---------------------

Once the installation is done, set environment variables (bash)::

    $ export PATH=/path/to/uvcdat/bin:$PATH

    # For install from sources
    $ export LD_LIBRARY_PATH=/path/to/uvcdat/Externals/lib:/path/to/uvcdat/lib:$LD_LIBRARY_PATH
    $ export C_INCLUDE_PATH=/path/to/uvcdat/Externals/include:$C_INCLUDE_PATH


Check the installation
----------------------

Then check::

    $ python -c "import cdms2"

Install other packages
----------------------

Using :program:`conda`::

    $ conda install -c conda-forge configobj PIL paramiko xlutils seawater pytz cmocean

Using pip::

    $ pip install configobj PIL paramiko xlutils seawater pytz
    $ pip install sphinxcontrib-cheeseshop sphinxcontrib-ansi sphinxcontrib-programoutput sphinx-fortran



