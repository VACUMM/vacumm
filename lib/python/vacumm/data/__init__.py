"""
Data tools

"""

import satellite
import model
import os as _os, locale as _locale
_os.environ['LC_NUMERIC'] = 'en_US.UTF-8'
_locale.setlocale(_locale.LC_NUMERIC, 'en_US.UTF-8')
from vacumm import VACUMMError

#: Specifications of available dataset types
dataset_specs = dict(
    generic = dict(
        cls = 'vacumm.data.misc.dataset.GenericDataset', 
        desc = 'Generic dataset'
    ), 
    ocean = dict(
        cls = 'vacumm.data.misc.dataset.OceanDataset', 
        desc = 'Generic ocean dataset'
    ), 
    atmos = dict(
        cls = 'vacumm.data.misc.dataset.Atmos.Dataset', 
        desc = 'Generic atmospheric dataset'
    ), 
    mars = dict(
        cls = 'vacumm.data.model.mars3d.Mars3D', 
        desc = 'MARS3D ocean model'
    ), 
    nemo = dict(
        cls = 'vacumm.data.model.nemo.Nemo', 
        desc = 'NEMO ocean model'
    ),
    hycom = dict(
        cls = 'vacumm.data.model.hycom.HYCOM',
        desc = 'HYCOM ocean model'
   ),
    swan = dict(
        cls = 'vacumm.data.model.swan.SWAN',
        desc = 'SWAN wave model'
   ),
    cfsr = dict(
        cls = 'vacumm.data.model.cfsr.CFSR',
        desc = 'CFSR atmospheric model'
   ),
)

#: List of available dataset types
dataset_names = dataset_specs.keys()

def DS(ncfile, clsname='generic', *args, **kwargs):
    """Load a specialized :class:`vacumm.data.misc.dataset.Dataset`
    
    Available methods for retreiving data (or derived data) are the same
    for all dataset types.
    And whatever the dataset type, the outputs will all have the same
    format (names, attributes...). Formating is defined and performed
    by the :mod:`vacumm.data.cf` module.
    
    :Params:
    
        - **ncfile**: Netcdf file name(s) or pattern.
        - **clsname**: Generic class name. Please choose one of:
          %s
          
    :Example:
    
        >>> from vacumm.data import DS
        >>> mars = DS('results.nc', 'mars')
        >>> dx,dy = mars.get_resol()
        >>> sst = mars.get_sst()
    """
    if clsname is None:clsname = 'generic'
    clsname = clsname.lower()
    if clsname in dataset_names:
        modname, Clsname = _os.path.splitext(dataset_specs[clsname]['cls'])
        Clsname = Clsname[1:]
        mod = __import__(modname, fromlist=[Clsname])
        cls = getattr(mod, Clsname)
    else:
        raise VACUMMError('Wrong name of dataset type: %s. '
            ' Please choose one of the following: %s'%(clsname, ', '.join(dataset_names)))
    
    return cls(ncfile, *args, **kwargs)
   
DS.__doc__ = DS.__doc__%', '.join(dataset_names)

def setup_dataset(clsname, ncfile=None, *args, **kwargs):
    """Alias for ``DS(ncfile,  clsname, *args, **kwargs)``
    
    See :func:`DS` for arguments.
    
    This function is here mainly for backward compatibility.
    """
    return DS(ncfile,  clsname, *args, **kwargs)
