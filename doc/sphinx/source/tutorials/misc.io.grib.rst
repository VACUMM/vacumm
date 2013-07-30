.. _tut.misc.io.grib:

Les fichiers grib
=================

Vacumm provides helper functions for loading grib1 and grib2 files data into cdms2 variables.
These functions uses the pygrib package which itself relies on the grib_api.

See:
    - pygrib: https://code.google.com/p/pygrib
    - grib_api: http://www.ecmwf.int/publications/manuals/grib_api/index.html
    - vacumm's io api: :mod:`~vacumm.misc.io`

One or more grib file(s) and parameter(s) can be loaded using the function :meth:`~vacumm.misc.io.grib_read_files` 

The script below:

.. literalinclude:: ../../../../scripts/tutorials/misc.io.grib.py
    :language: python

generate the following output:

.. program-output:: ../../../scripts/tutorials/misc.io.grib.py

