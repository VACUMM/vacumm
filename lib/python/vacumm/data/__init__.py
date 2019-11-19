# -*- coding: utf8 -*-
"""
Data tools

"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#

import os as _os, locale as _locale
try:
    _os.environ['LC_NUMERIC'] = 'en_US.UTF-8'
except:
    pass
try:
    _locale.setlocale(_locale.LC_NUMERIC, 'en_US.UTF-8')
except:
    pass

from vacumm import VACUMMError, vcwarn
import satellite

#: Specifications of available dataset types
DATASET_SPECS = {}
dataset_specs = DATASET_SPECS # compat

##: List of available dataset types
#DATASET_NAMES = DATASET_SPECS.keys()
#dataset_names = DATASET_NAMES # compat

def DS(ncfile, clsname='generic', *args, **kwargs):
    """Load a specialized :class:`vacumm.data.misc.dataset.Dataset` instance

    Available methods for retreiving data (or derived data) are the same
    for all dataset types.
    And whatever the dataset type, the outputs will all have the same
    format (names, attributes...). Formating is defined and performed
    by the :mod:`vacumm.data.cf` module.

    :Params:

        - **ncfile**: Netcdf file name(s) or pattern compatibile with
          :func:`~vacumm.misc.io.list_forecast_files`, so
          please read the documentation of this function carefully.
        - **clsname**: Generic class name. Please choose one of:
          %s
        - All other keyword are passed to the class initialization.
          See :class:`~vacumm.data.misc.dataset.Dataset` a list of these options.

    :Return: A child of :class:`vacumm.data.misc.dataset.Dataset`.

    :Example:

        >>> from vacumm.data import DS
        >>> mars = DS('results.nc', 'mars')
        >>> dx,dy = mars.get_resol()
        >>> sst = mars.get_sst()

    :See also: :ref:`user.desc.dataset`.

    """
    # Register builtin datasets
    import vacumm.data.misc.dataset
    import vacumm.data.model

    # Class name
    if clsname is None:
        clsname = 'generic'
    clsname = clsname.lower()

    # Class
    if clsname in DATASET_SPECS:
        cls = DATASET_SPECS[clsname]
    else:
        raise VACUMMError('Wrong name of dataset type: %s. '
            ' Please choose one of the following: %s'%(clsname, ', '.join(DATASET_SPECS.keys())))

    # Instantiate
    return cls(ncfile, *args, **kwargs)

DS.__doc__ = DS.__doc__%', '.join(DATASET_SPECS.keys())

def setup_dataset(clsname, ncfile=None, *args, **kwargs):
    """Alias for ``DS(ncfile,  clsname, *args, **kwargs)``

    See :func:`DS` for arguments.

    This function is here mainly for backward compatibility.
    """
    return DS(ncfile,  clsname, *args, **kwargs)

def register_dataset(cls, clsname=None, warn=True, force=True):
    """Register a new :class:`~vacumm.data.misc.dataset.Dataset` class
    that can be loaded with :func:`DS`
    """
    # Get the class when a path is given
    if isinstance(cls, str):
        modname, Clsname = _os.path.splitext(cls)
        Clsname = Clsname[1:]
        mod = __import__(modname, fromlist=[Clsname])
        cls = getattr(mod, Clsname)

    # Check inheritance
    from vacumm.data.misc.dataset import Dataset
    if not issubclass(cls, Dataset):
        raise VACUMMError('cls must be a Dataset subclass')

    # Get the name
    if clsname is None:
        clsname = cls.name
        if clsname is None:
            clsname = cls.__name__.lower()
    clsname = clsname.lower()
    if cls.name is None:
        cls.name = clsname

    # Already registered?
    if clsname in DATASET_SPECS:
        if warn:
            if force:
                ms = 'Overwriting it...'
            else:
                ms = 'Skipping...'
            vcwarn('Dataset class "{}" is already registered. '+ms)
        if not force:
            return

    # Register
    DATASET_SPECS[clsname] = cls
    return cls

def print_registered_datasets():
    """Print all registered datasets"""
    for clsname, cls in DATASET_SPECS.items():
        print clsname + ':'
        print '  description: {}'.format(cls.description or '')
        print '  class: ' + cls.__name__
