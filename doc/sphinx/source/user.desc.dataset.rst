.. _user.desc.dataset:

Using the generic data classes 
******************************

These classes are all inherited from :class:`~vacumm.data.misc.dataset.Dataset`.

Goal
====

Theses classes are generic interfaces to gridded data stored in
CF compliant netcdf files.
They provide standard method for accessing simple or complex quantities.

For instance, to access sea surface temperature, whatever the data
are satellite observations or results from models, you use::

    >>> mydata = MyData(myfile)
    >>> ss = mydata.get_sst()
    
:class:`MyData` will access the SST directly in the netcdf file,
or extract a surface slice from 3D temperature.

Specialized classes inherited from :class:`~vacumm.data.misc.dataset.Dataset`
are implemented to manage the interface to the different sources of data.

- There are some rather generic interfaces : 
  :class:`~vacumm.data.misc.dataset.OceanDataset`, 
  :class:`~vacumm.data.misc.dataset.OceanSurfaceDataset`, 
  :class:`~vacumm.data.misc.dataset.AtmosDataset`, 
  :class:`~vacumm.data.misc.dataset.AtmosSurfaceDataset`, 
  :class:`~vacumm.data.misc.dataset.Generic`.
- There are some very specialized ones :  
  :class:`~vacumm.data.model.mars3d.Mars3D`,
  :class:`~vacumm.data.model.nemo.NEMO`, 
  :class:`~vacumm.data.model.hycom.HYCOM`, 
  :class:`~vacumm.data.model.swan.SWAN`.
  This list will increase with time.

Usage
=====

You must provide at least a file name::

    >>> from vacumm.data.model.mars3d import Mars3D
    >>> mydata = Mars3D(ncfile)
   
You can also use the :func:`vacumm.data.DS` function:
  
    >>> mydata = DS(ncfile, 'mars')
    
In this latter case,  if you don't specify the dataset type (here "mars"), 
the generic class :class:`~vacumm.data.misc.dataset.Generic` is used.

Some options are supperted at initialization. For example::
    
    >>> mydata = Mars3D(ncfile, lon=(-10, 4),  logger_level='debug')
    
Then,  numerous methods are availables::
    
    >>> lon = mydata.get_lon_u()
    >>> mld = mydata.get_mld(mode='deltadens')
    

Sample of capabilities
======================

Reading multiple files
----------------------

These classes fully take advantage of the  (see documentation).
For example,  you can provide sereval files or a file name with date patterns.

    >>> Mars3D([ncfile1, ncfile2, ...])
    >>> Mars3D('myfile.%Y.nc', time=('2000', '2010'))
    
    
Searching for variables and axes
--------------------------------

The search algorithm typically exploit lists of standard names and names.
These lists are defined in the :mod:`~vacumm.data.cf` module,  
and can be completed within the specialized classes.
The functions that search for objects are 
:func:`~vacumm.misc.io.ncfind_var`and 
:func:`~vacumm.misc.io.ncfind_axis`.


Advanced variables
------------------

Generally,  a method will search for a variable which is
already in the dataset.
Other methods can also be implemented to retreive the variable.

From slices
~~~~~~~~~~~

A variable can be accesed as a slice of another variable.
For instance, SST can be defined as the surface slice of temperature
in 3D models.
This is setup within the specification of the class.

From algorithms
~~~~~~~~~~~~~~~

Sometimes,  variable can also be retreived as computation
using existing variables.
For instance, if the mixed layer depth is not store in
the netcdf file, it can be computed using temperature,  density or Kz
with the :func:`vacumm.diag.thermdyn.mixed_layer_depth`.
The :meth:`~vacumm.data.misc.dataset.Dataset.get_mld` method
will read all needed variables to compute it.


Spatial and temporal selections
-------------------------------

Spatial and temporal selections are availablesat initialization
or running time,  using the ``lon``,  ``lat``,  ``level`` and ``time``
generic keywords::
    
    >>> mars = Mars3d(ncfile,  time=('2000', '2010'))
    >>> sst = mars.get_sst(lon=(-10, 1),  lat=(43, 50))

Selection argument are those of CDAT,  and works well even
on curvilinear grids.


Vertical coordinates
--------------------

Depth can be either a 1D axis,  a volumique variable,  or
estimated from sigma coordinates (using the
:mod:`~vacumm.data.misc.sigma` module) or 
by integrated layer thicknesses.



Logging and verbosity
---------------------

Data classes integrate the logging system of module :mod:`~vacumm.misc.log`.
You can for instance define a logging level and print messages
of different level of importance:
    
    >>> mars = Mars3D(ncfile, logger_level='warning', logger_file='mars.log') 
    >>> mars.debug('debug message') # nothing is printed
    >>> mars.warning('warning message')
    



Development
===========

Adding a new class for a new dataset type
-----------------------------------------

Adding a new method
-------------------


More information
================
