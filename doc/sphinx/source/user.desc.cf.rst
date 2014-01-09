.. _user.desc.cf:

CF conventions usage and extensions
***********************************

The module :mod:`vacumm.data.cf` defines and exploits
specifications inspired from CF conventions.
These specifications are separated in two dictionaries: 
- one for variables: :attr:`~vacumm.data.cf.var_specs`,
- one for axes: :attr:`~vacumm.data.cf.axis_specs`.

These specifications are used in the library to search
and format variables or axes.

.. _user.desc.cf.specs:

Structure of the specification dictionaries
-------------------------------------------

Here is an example::

    u3d = dict(
        names=['uz', 'u3d'],
        standard_names=['sea_water_x_velocity', 'eastward_sea_water_velocity'],
        long_names = "Sea water velocity along X", 
        units = "m s-1", 
        axes = dict(x=['lon_u'], y=['lat_u']), 
        physloc = 'u', 
        atlocs = ['t', 'u', 'v'], 
    ),


Each entry has a key that is a short name that will be used
as a generic netcdf name (here ``u3d``).
The value of an entry is a dictionary that can define de following specifications:

    - **names**: A list of lower case short name that can be found in netcdf files.
      The short name the entry (``u3d``) is insert automatically at the beginning of this list.
    - **standard_names**: A list of possible standard_names to search for the variable in
      a netcdf file. The first one is used to format variables.
    - **long_names**: One or several long_names. The first one is used to format variables.
    - **units**: One or several long_names. The first one is used to format variables.
    - **axes**: A dictionary of axis names used to format axes when they are recongnized
      as geographic axes (xyzt). Each axis name must be an entry of :attr:`vacumm.data.cf.axis_specs`.
    - **physloc**: The physical grid location of this variable on an Arakawa grid.
    - **atlocs**: A list of grid locations (see :attr:`~vacumm.data.misc.arakawa.locations`)
      that will be used to deplicate this entry to match another location 
      (with :func:`~vacumm.data.cf.specs_dup_loc`).
      For instance ``u3d`` entry will create ``u3d_t``, ``u3d_u`` and ``u3d_v`` entries.
      Then, the ``u3d`` entry will be modified depending on **physloc**:
      
        - If **physloc** is not set, ``u3d`` entry  will be merged with  
          auto-generated entries (``u3d_t``, ``u3d_u`` and ``u3d_v``):
          if you search for ``u3d``, it will also search for entries at other locations.
        - if **physloc** is set, ``u3d`` will be an alias for ``u3d_<physloc>``, here ``u3d_u``:
          if you search for ``u3d``, it will search for ``u3d_u`` if not found, and only it.
    - **inherit**: It can be used to inherit the specifications from another entry.
    - **i/jaxis**: Name of elemental axis for 2D axes.
    - The other key/val are used as attributes when formatting a variable.
    

Searching for variables and axes
--------------------------------

You can search for variables and axes in netcdf files with 
the :func:`~vacumm.misc.io.ncfind_obj` function.
The main search argument can be a simple tuple of netcdf names,
or a dictionary of names, standard_names, or even long_names and units.
This dictionary can be automatically formed by the :func:`~vacumm.data.cf.cf2search`
function for a given generic name (such as "sst").
This way of searching is used by the :class:`~vacumm.data.misc.dataset.Dataset` classes.


Formatting variables and axes
-----------------------------

You can format variables and axes to make sure they meet CF standards,
with the functions :func:`~vacumm.data.cf.format_var` and :func:`~vacumm.data.cf.format_axes`.
The attributes are defined by the :func:`~vacumm.data.cf.cf2atts` function.
  
    