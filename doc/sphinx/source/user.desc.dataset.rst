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

- There are some rather generic interfaces: 
  :class:`~vacumm.data.misc.dataset.OceanDataset`, 
  :class:`~vacumm.data.misc.dataset.OceanSurfaceDataset`, 
  :class:`~vacumm.data.misc.dataset.AtmosDataset`, 
  :class:`~vacumm.data.misc.dataset.AtmosSurfaceDataset`, 
  :class:`~vacumm.data.misc.dataset.GenericDataset`.
- There are some very specialized ones:  
  :class:`~vacumm.data.model.mars3d.Mars3D`,
  :class:`~vacumm.data.model.nemo.Nemo`, 
  :class:`~vacumm.data.model.hycom.HYCOM`, 
  :class:`~vacumm.data.model.swan.SWAN`,
  :class:`~vacumm.data.model.swan.CFSR`,
  and more.
  This list will increase with time.

Usage
=====

You must provide at least a file name::

    >>> from vacumm.data.model.mars3d import Mars3D
    >>> mydata = Mars3D(ncfile)
   
You can also use the :func:`vacumm.data.DS` function:
  
    >>> mydata = DS(ncfile, 'mars')
    
In this latter case,  if you don't specify the dataset type (here "mars"), 
the generic class :class:`~vacumm.data.misc.dataset.GenericDataset` is used.

Some options are supported at initialization. For example::
    
    >>> mydata = Mars3D(ncfile, lon=(-10, 4),  logger_level='debug')
    
Then,  numerous methods are availables::
    
    >>> lon = mydata.get_lon_u()
    >>> mld = mydata.get_mld(mode='deltadens')
    

Sample of the capabilities
==========================

Reading multiple files
----------------------

These classes fully take advantage of the  (see documentation).
For example,  you can provide sereval files or a file name with date patterns.

    >>> Mars3D([ncfile1, ncfile2, ...])
    >>> Mars3D('myfile.%Y.nc', time=('2000', '2010'))

It also accepts url, global unix patterns (``*``, ``?``, etc),
and date patterns (``%Y``, ``%m``, etc, see :func:`~vacumm.misc.atime.strftime`).
For more information on how to read a set of files, please read the 
:func:`~vacumm.misc.io.list_forecast_files` function for options.

.. note:: Files are sorted by alphabetically, unless you you pass
    the ``sort`` keyword to prevent or tune sorting.
    
    
Searching for variables and axes
--------------------------------

The search algorithm typically exploit lists of standard names and names.
These lists are defined in the :mod:`~vacumm.data.cf` module,  
and can be completed within the specialized classes.
The functions that search for objects are 
:func:`~vacumm.misc.io.ncfind_var`and 
:func:`~vacumm.misc.io.ncfind_axis`.


.. _user.desc.dataset.cap.sim:

Simple variables
----------------

You can access these variables using the generic way like the method 
:meth:`~vacumm.data.misc.dataset.OceanDataset.get_temp`.
If the variable cannot be found directly, the methods searches for it
using a name at another similar location name (no location or its physical
location according to :attr:`vacumm.data.cf.var_specs`).
Finally, the method try to get it by reading it at different 
location and interpolate it to the right one.

The method provide the **mode** keyword to choose the retreiving mode.

- ``"var"``: Read it directly from the file.
- ``"stag"``: Get it from staggered location thanks to interpolation.
- ``None``: Try them all!


.. _user.desc.dataset.cap.adv:

Advanced variables
------------------

Generally, a method will search for a variable which is
already in the dataset.
Other methods can also be implemented manually to retreive the variable,
typically from other physical variables.
You can choose which retreiving method to use thanks to the **mode** keyword:

    >>> mld = ds.get_mld(mode='var') # read it
    >>> mld = ds.get_mld(mode='deltadens') # compute it

In the second case, it is computed from the temperature, the salinity and the depths.

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
The :meth:`~vacumm.data.misc.dataset.OceanDataset.get_mld` method
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


Reshaping variables
-------------------

It can be useful sometimes to reshape a variable so that it has the same
shape as another variable, for instance to mathematical combinations
or statistics.
You can use the **asvar** keyword for such task.

    >>> bathy3d = ds.get_bathy(asvar='temp')
    
    
Interpolating a variable to another grid location
--------------------------------------------------

When a dataset is on an staggered Arakawa grid (such as a C grid),
you can use the **at** keyword to interpolate a variable
at another location.
For instance, to put sst at a U location, just use:

    >>> sst_u = ds.get_sst(at='u')
    
.. note:: Note that sometimes a method already exists to do it:

    >>> sst_u = ds.get_sst_u()


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

Adding a new get method
-----------------------

This is useful to ...

Simple variable to read
~~~~~~~~~~~~~~~~~~~~~~~

In this case, you just define an interface to read the variable in a netcdf file.

First, define a new entry in the :attr:`~vacumm.data.cf.var_specs` dictionary
(see section :ref:`user.desc.cf.specs`).
Let's name it ``myvar``.

Then you can choose between two methods:

- Add the entry name ``myvar`` to the list defined by the :attr:`auto_generic_var_names` class
  attribute (create it if no present). This method also makes the **mode** keyword available
  to choose the way you will retreive your variable :ref:`user.desc.dataset.cap.sim`.
- **Or** define the method manually with the decorator :func:`~vacumm.data.misc.dataset.getvar_fmtdoc`::

    @getvar_fmtdoc
    def get_myvar(self, **kwargs):
        return self.get_variable('myvar', **kwargs)
        


Advanced variable
~~~~~~~~~~~~~~~~~

In this case, you want to be able to read the variable or compute it
(see paragraph :ref:`user.desc.dataset.cap.adv`).

You must define the method manually and introduce the **mode** keyword,
and use the :func:`vacumm.data.misc.dataset.check_mode` function to
verify which mode the user wants.
The :func:`vacumm.data.misc.dataset.getvar_fmtdoc` format the docstring
of the method.

In the follwing example, we want to be able to read the zonal wind stress,
or compute it from the wind speed components using function :func:`~vacumm.diag.atmos.wind_stress`::

    # Define the method
    def get_taux(self, mode=None, **kwargs):
        
        # Params
        kwvar = kwfilter(kwargs, ['lon','lat','time','level','torect'], warn=False)
        kwfinal = kwfilter(kwargs, ['squeeze','order','asvar', 'at'])
        kwspeed = kwfilter(kwargs, 'speed_') # 

        # First, try to find a taux variable
        if check_mode('var', mode):
            taux = self.get_variable('taux', **kwvar)
            if taux is not None or check_mode('var', mode, strict=True): 
                return self.finalize_object(taux, **kwfinal)
            
        # Estimate it from wind speed components
        if check_mode('speed', mode):
            u10m = self.get_u10m(**kwvar)
            v10m = self.get_v10m(**kwvar)
            if u10m is not None and v10m is not None:
                taux,_ = wind_stress(u10m, v10m, format_axes=True, **kwspeed)
            if taux is not None or check_mode('speed', mode, strict=True): 
                return self.finalize_object(taux, depthup=False, **kwfinal)
      
    # Generate the docstring
    getvar_fmtdoc(get_taux, 
        mode="""Retreiving mode
    
          - ``None``: Try all modes, in the following order.
          - ``"var"``: Read it from a variable.
          - ``"speed"``: Compute it from wind speed.
          """,
        speed_rhoa="Air density (in kg.m-3)",
        speed_cd="Drag coefficient",
        )



Adding a new class for a new dataset type
-----------------------------------------

This can be useful for example to create an interface to a gridded
product such as model outputs or satellite data.

1. Create a module somewhere in package :mod:`vacumm.data` that contains a class
   that is subclass of one or several classes of :mod:`vacumm.data.misc.dataset`::
   
     class MyClass(OceanDataset):
     ...
     
2. Define important class attributes:

    - :attr:`~vacumm.data.misc.dataset.arakawa_grid_type` if the dataset is on a special grid.
    - :attr:`~vacumm.data.misc.dataset.positive` if you know that your dataset is positive up or down for levels.
    - :attr:`~vacumm.data.misc.dataset.ncobj_specs` if you want to redefine the specification 
      of some variables in the same form as in :mod:`~vacumm.data.cf.var_specs`, but with 
      two other possible keys:
      
        
        - **select**: To perform a systematic selection into the variable.
        - **squeeze**: To systematically squeeze out one or more dimensions.
      
      For instance, the SST in :class:`~vacumm.data.model.mars3d.Mars3D` is redefined
      so that it is read from the last level of 3D temperature, squeezing out the vertical
      dimension::
      
        ncobj_specs = {

            # sea surface temperature
            'sst':{
                'inherit':'temp', 
                'select':{'level':slice(-1, None)}, 
                'squeeze':'z', 
            }
            
        }
        
2. Define secondary class attributes and other special methods if needed.
3. Add an entry in :attr:`vacumm.data.dataset_specs` to create a shortcut to
   this dataset whe using :func:`vacumm.data.DS`.

