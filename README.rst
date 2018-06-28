VACUMM
======

.. image:: https://zenodo.org/badge/22859/VACUMM/vacumm.svg
   :target: https://zenodo.org/badge/latestdoi/22859/VACUMM/vacumm
.. image:: https://travis-ci.org/VACUMM/vacumm.svg?branch=master
    :target: https://travis-ci.org/VACUMM/vacumm

VACUMM provides generic and specialized tools for the validation of ocean models,
and more especially the MARS model from `IFREMER <http://www.ifremer.fr>`_.
The heart of VACUMM is a
`library <http://www.ifremer.fr/vacumm/library/index.html>`_  written mainly
in the `Python <http://www.python.org>`_ language,
whose `core <http://www.ifremer.fr/vacumm/library/misc.html>`_
can be used for the **preprocessing** and the
**postprocessing** of oceanic and atmospheric data coming from models or observations.
The library for instance also has specialized modules for managing outputs from
`models <http://www.ifremer.fr/vacumm/library/data/model.html>`_ and making advanced
`diagnostics <http://www.ifremer.fr/vacumm/library/diag.html>`_.

.. code-block:: python

    >>> from vcmq import *
    >>> sst = DS(data_sample('mars3d.xy.nc'), 'mars').get_sst()
    >>> map2(sst)


Features
--------

- A huge documentation with a gallery, a lot of examples and the complete API:
  http://www.ifremer.fr/vacumm
- Full UV-CDAT support and extensions.
- Matplotlib/basemap graphics with advanced plotting objects like geographical mapping tools.
- Numerous utilities for manipulating and converting time data.
- Regridding and interpolation of random or gridded data, in 1D or 2D, with curvilinear grid support.
- Helper routines for inspecting and reading NetCDF objects in single or multiple file datasets.
- Generic and specialized 1D and 2D filters working on masked variables.
- Specialized physical and numerical diagnostics, like dynamics, thermodynamics, spectral analyses, tides, etc.
- Support and extension of CF conventions for searching or formatting variables.
- Miscellaneous location utilities such as readers of sigma coordinates for ocean models, or Arakawa grid converters.
- High level generic interface for reading and post-processing NetCDF data from standard or known dataset, such as model outputs or satellite products.
- Statistical accumulator for large datasets.
- Interfaces for working with random and gridded bathymetries, and with shorelines.
- Utilities for working with masks and selection of spatial data.
- Utilities for working with input and output remote files.
- Advanced logging classes.
- Extensions to sphinx for Fortran and Python.
- A collection of scripts for some diagnostics.


Dependencies
------------

Mandatory:
`CDAT <http://uvcdat.llnl.gov>`_ (or more specifically
`cdms2 <http://uvcdat.llnl.gov>`_,
`cdutil <http://uvcdat.llnl.gov>`_,
`genutil <http://uvcdat.llnl.gov>`_ from CDAT, and
`matplotlib <https://matplotlib.org>`_,
`basemap <https://matplotlib.org/basemap>`_),
`configobj <http://www.voidspace.org.uk/python/configobj.html>`_.

Optional:
`seawater <https://pypi.python.org/pypi/seawater>`_,
`PIL <https://pypi.python.org/pypi/PIL>`_,
`pytz <http://pytz.sourceforge.net>`_,
`paramiko <http://www.paramiko.org>`_,
`xlwt <https://pypi.python.org/pypi/xlwt>`_,
`sphinx-fortran <https://pypi.python.org/pypi/sphinx-fortran>`_,
`cmocean <https://pypi.python.org/pypi/cmocean>`_.


Download
--------

To download VACUMM sources, please go to this page:
http://www.ifremer.fr/vacumm/user.install.download.html


Installation
------------

From sources::

    $ python setup.py install

Using `conda <http://conda.pydata.org/docs/index.html>`_::

    $ conda install -c vacumm -c conda-forge -c cdat  vacumm

For more information, please go to this:
http://www.ifremer.fr/vacumm/user.install.installations.html

Release notes
-------------

Release notes for each version are available here:
http://www.ifremer.fr/vacumm/appendix.release.html


Documentation
-------------

The documentation is available here:
http://www.ifremer.fr/vacumm


License
-------

VACUMM is under the :ref:`CeCiLL <appendix.license>` license,
which is compatible with well knwon GPL license.


Support
-------

You can submit `issues <https://github.com/VACUMM/vacumm/issues>`_
and `pull requests <https://github.com/VACUMM/vacumm/issues>`_
from the GitHub site.

Stephane Raynaud (raynaud (at) gmail.com),
Guillaume Charria (Guillaume.Charria (at) ifremer.fr).

See the contact page:
http://www.ifremer.fr/vacumm/contact.html


