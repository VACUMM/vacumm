#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar (2010)
# 
# wilkins@actimar.fr
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


# ==============================================================================
# TODO:
# ==============================================================================
#
# - Move variable mapping feature from data.misc.profile to this module
# - Determine reference/run date from Catalog datasets, sort the loaded dataset,
#   move dataset aggregation over time into the Catalog class)
# - Add an interface (class) defining the root base feature of a Dataset
#   (aggregation using catalog, get_[time,depth,latitude,longitude,variable], get_grid, ...)
#   to provide a data format abstraction base.
# - Ensure compatibilty of possibly different axes coordinates (lat, lon, sigma, depths, ...)
#   over multiple dataset (e.g. about sigma: user.algo.readers.html)
# - Force variables order in dedicated data loaders (get_section, get_hovmoller, ...)
#
# - Check why 
#     vacumm.misc.io.NcIterBestEstimate(files = [...],time = (..., 'co'))
#   produces empty slices (e.g slice(14,2)), but doesnt if 'ccb' is used
# ==============================================================================


__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__date__ = '2011-01-17'
__doc__ = 'Common dataset utilities'


# ==============================================================================


import os, sys
from traceback import format_exc

import cdms2, MV2, numpy, pylab, seawater
from matplotlib.pyplot import colorbar
N = numpy

from vacumm.misc import auto_scale, MV2_concatenate, MV2_axisConcatenate, create_selector, \
    split_selector, dict_merge, dict_check_defaults, set_atts, broadcast
from vacumm.misc.atime import add as add_time, comptime, datetime as adatetime, Intervals
from vacumm.misc.axes import create_time, create_dep, create_lat, create_lon, guess_timeid, get_axis_type, isaxis
from vacumm.misc.bases import Object
import vacumm.misc.color as C
from vacumm.misc.io import ncget_var, ncfind_var, ncfind_obj, ncfind_axis, list_forecast_files, ncread_files, \
    ncread_best_estimate, NcIterBestEstimate, ncget_time, ncget_lon, ncget_lat, ncget_level,  \
    ncget_grid, ncget_axis
from vacumm.misc.grid.misc import meshweights, resol, create_grid, set_grid, \
    dz2depth, depth2dz, isdepthup, gridsel, makedepthup, curv2rect, isgrid,  \
    coord2slice
from vacumm.misc.grid.regridding import resol, interp1d, regrid1d, grid2xy, transect,  \
    shift1d, shift2d, extend1d, extend2d
from vacumm.misc.misc import is_iterable, kwfilter, squeeze_variable
from vacumm.misc.phys.constants import g
from vacumm.misc.plot import map2, curve2, section2, hov2, add_map_lines
from vacumm.misc.axes import islon, islat, isdep
from vacumm.data.misc.sigma import NcSigma
from vacumm.data.cf import var_specs, axis_specs, cf2search, cf2atts, generic_names, \
    generic_axis_names, generic_var_names, format_axis, format_var
from vacumm import VACUMMError
from vacumm.diag.dynamics import barotropic_geostrophic_velocity, eddy_kinetic_energy
from vacumm.diag.thermdyn import mixed_layer_depth, density

from vacumm.misc.misc import grow_variables

# Kwargs which are pop'ed from intermediate plots
# (e.g. plot section, then trajectory map, then post plot using these kwargs)
post_plot_args = ('title', 'show', 'close', 'savefig', 'savefigs')

_geonames =  {'x':'lon', 'y':'lat', 'z':'level', 't':'time'}

# Docstring formatting

_formatdoc_var_tpl = """%(func_name)s(time=None, level=None, lat=None, lon=None, squeeze=None, order=None, asvar=None, verbose=None, warn=None, **kwargs)

    Read the %(long_name)s
    
    When multiple files are accessed, data are concatenated in time in
    a best estimate fashion.
    
    :Params:
    
        - **time/level/lat/lon**, optional: For selection (tuples or slices).
        - **squeeze**, optional: Squeeze singleton dimensions 
          (see :func:`~vacumm.misc.misc.squeeze_variable`).
        - **order**, optional: Change order of axes (like ``'xty'``).
        %(other_params)s%(kwargs_params)s
        - **asvar**, optional: Grow variable to match the ``asvar`` variable,
          using :func:`~vacumm.misc.misc.grow_variables`.
    """

def formatdoc_var(func, **kwargs):
    """Format the docstring of basic variables"""
    kwargs_params = kwargs.get('_kwargs', "Other keywords are passed to :func:`~vacumm.misc.io.ncread_files`.")
    func_name = func.__name__
    if not func_name.startswith('get_'): return func
    var_name = func_name[4:]
    if not var_name in generic_var_names: return func
    long_name = var_specs[var_name]['long_names'][0]
    long_name = long_name[0].lower()+long_name[1:]
    if kwargs:
        other_params = []
        indent = '    '
        keys =  sorted(kwargs.keys())
        for key in keys:
            other_params.append('- **%s**, optional: %s\n'%(key, kwargs[key]))
        other_params = (indent*2).join(other_params)+'\n'
        other_params += indent*2
    else:
        other_params = ''
    kwargs_params = ('- '+kwargs_params) if kwargs else ''
    func.__doc__ = _formatdoc_var_tpl%locals()
    return func




class CatalogError(Exception):
    pass

class Catalog(Object):
    '''Build dataset files list based on self.get_config
    
    Variables in config (all optionnal):
        - **files**: explicit list of datasets files
        - **filepattern** , **time**: see :func:`~vacumm.misc.io.list_forecast_files`
    
    These configuration can be passed when creating the Catalog object
    like any other Object, examples:
        >>> c = Catalog(config='configfile.cfg')
        >>> c = Catalog(config=dict(files=('file1.nc', 'file2.nc')))
        >>> c = Catalog(config=dict(filepattern='model_%Y-%m-%d', time=('2010-01-01', '2011-01-01')))
        >>> c = Catalog(config=dict(filepattern='model_%Y-%m-%d', time=('2010-01-01', '2011-01-01'), files=('file1.nc', 'file2.nc')))
    '''
    
    def __init__(self, files=None, filepattern=None, time=None, **kwargs):
        self.files, self.filepattern, self.time = files, filepattern, time
        Object.__init__(self, **kwargs)
    
    @classmethod
    def find_datasets(cls, filepattern, time, **kwargs):
        '''
        Perform dataset search using :func:`~vacumm.misc.io.list_forecast_files`
        '''
        cls.info('Searching datasets using pattern: %s, time: %s', filepattern, time)
        findings = []
        if not isinstance(filepattern, (list, tuple)):
            filepattern = [filepattern]
        for fp in filepattern:
            findings.extend(list_forecast_files(fp, time, **kwargs))
        cls.info('Found %s datasets', len(findings))
        return findings
    
    def get_datasets(self):
        '''Get the logical view of datasets.
        
        Datasets are specified by this object configuration (:meth:`~vacumm.misc.bases.Object.get_config`) using:
            - files: use files specified directly
            - filepattern and time: perform dataset search using :func:`~find_datasets`
        
        :Return:
            - List of strings, the dataset files found
        
        .. todo::
            - Ensure returned dataset order (by (run)date)
            - Add a tutorial about :class:`Catalog`
        
        '''
        datasets = []
        if self.files: datasets.extend(files)
        if self.filepattern and self.time:
            datasets.extend(self.find_datasets(self.filepattern, self.time))
        return datasets

class DatasetError(Exception):
    pass

class Dataset(Object):
    '''Generic dataset representation.
    
    Handle multiple datasets (files) in a best time serie view.
    
    Variables specifications are stored in :attr:`ncobj_specs`.
    This global attribute must be refined for all new dataset.
    
    .. attribute:: ncobj_specs
    
        Dictionary describing variable specifications: final id
        of the variable as the keys, and aliases ad slice as the content.
        Here is an example::
        
            ncobj_specs = {
                'u10m':dict(
                    search=dict(names=['u', 'u10', 'uwind', '.*u.*']), 
                    select=(Ellipsis, 0)), 
                'sst':dict(
                    search=dict(standard_name="sea_surface_temperature"), 
                    select=(Ellipsis, -1)), 
                    atts=dict(units="DegC"),
    }
    
    :Params:
        - **dataset**: A dataset or list of dataset (see :meth:`load_dataset`).
        - **select**: A selector to be applied by default.
        - **time**: Time for searching files (see :func:`~vacumm.misc.io.list_forecast_files`) 
          which overwrites configuration.
        - **filepattern**, optional: File pattern for searching files 
          (see :func:`~vacumm.misc.io.list_forecast_files`), which overwrites configuration.
        - **kwargs**: passed to Object.__init__, you may then pass a config filepath/dict/ConfigObj.
    
    :Examples:
        
        >>> d = Dataset()
        >>> d.load_dataset('file1.nc')
        >>> d.load_dataset(('file2.nc', 'file3.nc'), append=True)
        
        Is equivalent to:
        
        >>> d = Dataset(('file1.nc', 'file2.nc', 'file3.nc'))
        
        
        Another example using configuration:
        
        >>> d = Dataset(config='dataset.cfg')
        
        With configuration file 'dataset.cfg':
            [Dataset]
                [[Catalog]]
                    files = 'file1.nc', 'file2.nc', 'file3.nc'
        
        .. note::
            For a subclass of Dataset, say for example MyModel, the configuration would be:
                [MyModel]
                    [[Catalog]]
                        files = 'file1.nc', 'file2.nc', 'file3.nc'
        
        .. todo::

            - Move out non generic features (mld, ped, plots) to dedicated classes (ex: Stratification, PlotBuilder classes ?)
            - Add a tutorial about :class:`Dataset`
            
            
    :Attributes:
    
    
        .. attribute:: selector
        
            A :class:`cdms2.selectors.Selector` created at initialization time.
            
        .. attribute:: selects
        
            Same as :attr:`selector` but as dictionary.
    '''
    
    ncobj_specs = {}
    
    def __init__(self, dataset=None, time=None, lon=None, lat=None, level=None, 
        ncobj_specs=None, **kwargs):
        
        # Load dataset
        if dataset is None: dataset = []
        self.dataset = []
        self.selects = dict(time=time, lon=lon, lat=lat, level=level)
        self.selector = create_selector(**self.selects)
        Object.__init__(self, **kwargs)
        self.load_dataset(dataset, time=time)
        self._ncids = {}
        
        # Load specs of variables
        self._load_ncobj_specs_(ncobj_specs)
        
    def _format_single_ncobj_specs_(self, specs, name=None):
        """Load and reformat ncobj_specs dictionary
        
        :Params:
        
            - **specs**: Dictionary of the folowing possible form::
            
                {'search': {
                    'names':[name1,...], 
                    'standard_names': [standard_name1,...],
                    'long_names': [long_name1,...]
                    },
                 'atts': {'units':units,...},
                 'select': select,
                 'squeeze': squeeze_specs,
                }
              
              If ``search`` is not found, it is taken from :mod:`vacumm.data.cf`.
              If ``'units'`` and ``'long_name'`` are not found in ``atts``,
              
              ``select`` is selector as a dictionary, cdms2 selector, tuple or list,
              passed to :func:`~vacumm.misc.io.ncread_var`.
              
            - **name**, optional: Name (id) of the variable.
        """
#        # From another generic variable
#        if 'fromobj' in specs:
#            from_specs = self._get_ncobj_specs_(specs['fromobj'])
#            specs.update(dict_merge(specs, from_specs))
#            del specs['fromobj']
        
        
        # Search specs
        search = specs.get("search", {})
        specs["search"] = search
        # - names
        names = search.get("names", [])
        if isinstance(names, tuple):
            names = list(names)
        elif not isinstance(names, list):
            names = [name]
        if name is not None:
            if name in names: names.remove(name)
        names.insert(0, name)
        search["names"] = names
        # - standard_names
        if 'standard_names' in search and not isinstance(search['standard_names'], list):
            if isinstance(search['standard_names'], tuple):
                search['standard_names'] = list(search['standard_names'])
            else:
                search['standard_names'] = [search['standard_names']]
        # - removes everything else
        for key in search:
            if key not in ['names', 'standard_names']:
                del search[key]
       
        # Attributes
        atts = specs.get("atts", {})
        specs["atts"] = atts
        # - put id
        if name is not None:
            atts['id'] = name
        # - put standard_name
        if 'standard_names' in search:
            atts.setdefault('standard_name', search['standard_names'][0])
            
        # Selection
        if 'slices' in specs:
            if "select" not in specs:
                specs['select'] = specs['slices']
            del specs['slices']
        
        # Squeeze
        if 'squeeze' in specs and not isinstance(specs['squeeze'], list):
            specs['squeeze'] = [specs['squeeze']]

        return specs

    def _load_ncobj_specs_(self, ncobj_specs=None):
        """Read the :attr:`ncobj_specs` attribute and reformat it"""
        # Get specs
        if ncobj_specs is None:
            ncobj_specs = self.ncobj_specs # get default
        else:
            self.ncobj_specs = ncobj_specs # set default
            
        # Loop on variables to format their specs
        for name, specs in ncobj_specs.items():
            self._format_single_ncobj_specs_(specs, name)
            
    def _get_ncobj_merged_specs_(self, varname, searchmode=None):
        """Get object specs merged from :mod:`vacumm.data.cf` 
        and class level specification (attribute :attr:`ncobj_specs`)
        
        Merging is done in the following order:
        
            1. Specs from :mod:`vacumm.data.cf`.
            2. Specs declared at the class definition level (attribute :attr:`ncobj_specs`).
            3. Specs from another variable specified through the 'fromobj' 
               specification key.
        
        :Params:
        
            - **varname**: Valid generic name.
        
        :Return:
        
            A dictionary of specifications.
        """
        
        # Specs from vacumm.data.cf
        fromobj = None
        if varname in generic_names:
            cf_specs = dict(search=cf2search(varname, mode=searchmode), atts=cf2atts(varname))
            genname = varname
            fromobj = cf_specs.pop('fromobj', None)
        else:
            cf_specs = None

        # Specs from current class
        if varname in self.ncobj_specs:
            local_specs = self.ncobj_specs[varname].copy()
            if genname is None: genname = varname
            fromobj = local_specs.pop('inherit', fromobj)
        else:
            local_specs = None
        
        # From other var
        if fromobj:
            from_specs = self._get_ncobj_merged_specs_(fromobj)
        else:
            from_specs = None

        # Merge
        if from_specs is not None and cf_specs is None and local_specs is None:
            raise DatasetError('Generic var/axis name "%s" has no specification defined in current class or module vacumm.data.cf'%varname)
        specs = dict_merge(cf_specs, local_specs, from_specs, mergelists=True, mergedicts=True)
        
        return specs
    
    def _get_ncobj_specs_(self, varname, attsingle=True, searchmode=None):
        """Get specs for searching, selecting and setting defaults of a variable
        
        Calls :meth:`_get_ncobj_merged_specs_` for merging specifications.
        
        :Return: genname, search, select, squeeze, atts
        """
        specs = None
        select = create_selector()
        atts = None
        genname = None
        squeeze = []
        
        # Generic name -> cf specs and/or class level specs
        if isinstance(varname , basestring) and not varname.startswith('+'): 
            
            specs = self._get_ncobj_merged_specs_(varname, searchmode=searchmode)
            genname = varname
                        
                        
        # Dict of specifications
        elif isinstance(varname, dict): 
            
            specs = self._format_single_ncobj_specs_(varname)
            if 'fromobj' in specs:
                from_specs = self._get_ncobj_merged_specs_(specs['fromobj'], searchmode=searchmode)
                specs = dict_merge(specs, from_specs)
            
            
        # Single (+varname) or list of netcdf names
        else:
            
            if isinstance(varname , basestring) and varname.startswith('+'): # Netcdf name
                varname = varname[1:]
            search = varname
            if isinstance(search, list): 
                search = tuple(search)

        # Get search, selector and attributes from specs
        if specs is not None:
            
            # Search
            search = specs['search']
            if isinstance(search['names'], list):
                search['names'] = tuple(search['names'])
            
            # Selector
            vselect = specs.get('select', None)
            if vselect is not None:
                select.refine(create_selector(vselect))
                
            # Attributes
            if atts is None:
                atts = specs.get('atts', None)
                if atts and attsingle:
                    for att, val in atts.iteritems():
                        if isinstance(val, list):
                            atts[att] = val[0]
                
            # Squeeze
            squeeze = specs.get('squeeze', [])
            if not isinstance(squeeze, list):
                squeeze = [squeeze]
            
                
        return genname, search, select, squeeze, atts
        
    def apply_config(self, config, **kwargs):
        '''Apply passed configuration (usually internal call from load_config)
        
        Call :meth:`~vacumm.misc.bases.Object.apply_config` and 
        load datasets (if config has a Catalog section)
        
        '''
        self.verbose('Applying configuration:\n  %s', '\n  '.join(config.write()))
        Object.apply_config(self, config)
        catalog = Catalog.from_config(config, nested=True)
        self.load_dataset(catalog, **kwargs)
    
    def load_dataset(self, dataset, time=None, **kwargs):
        '''Load dataset files.
        
        :Params:
            - **dataset**: can be either:
                - an instance or list of strings (filepath/filepattern/url) or 
                - an instance of :class:`Catalog`
            - **time**: used with :meth:`Catalog.find_datasets` when dataset is a string (or list of strings)
        
        :Keyword parameters:
            - **append**: keep previously loaded datasets if True
        
        :Return: The list of (newly) loaded datasets
        
        :Example:
        
            >>> d.load_dataset('../../data/model/mars/champs_%Y%m_BOBI.nc', ('2004-01-01', '2004-04-01'), patfreq=(2,'month'), verbose=True)
            [2013-01-14 09:42:43 CET Catalog INFO] Searching datasets using pattern: ['../../data/model/mars/champs_%Y%m_BOBI.nc'], time: ('2004-01-01', '2004-04-01')
            Guessing file list using:
               filepattern: ../../data/model/mars/champs_%Y%m_BOBI.nc
               time selector: ('2004-01-01', '2004-03-31')
            Found 2 files
            [2013-01-14 09:42:43 CET Catalog INFO] Found 2 datasets
            [2013-01-14 09:42:43 CET Dataset INFO] Loading datasets:
            [2013-01-14 09:42:43 CET Dataset INFO] - ../../data/model/mars/champs_200401_BOBI.nc
            [2013-01-14 09:42:43 CET Dataset INFO] - ../../data/model/mars/champs_200403_BOBI.nc
        
        '''
        append = kwargs.pop('append', None)
        # Catalog case: get its dataset files
        if isinstance(dataset, Catalog):
            datasets = dataset.get_datasets()
        # Other cases: let Catalog find files using filepaths / filepatterns & time
        else:
            if not isinstance(dataset, (list, tuple)):
                dataset = [dataset]
            datasets = Catalog.find_datasets(dataset, time, **kwargs)
        # Open new datasets
        if datasets:
            self.info('Loading datasets:')
            for i,d in enumerate(datasets):
                datasets[i] = cdms2.open(d)
                self.info('- %s', datasets[i].id)
        # If append mode, extend internal list
        if append:
            self.dataset.extend(datasets)
        # Otherwise close internal list and overwrite it
        else:
            self.close()
            self.dataset = list(datasets)
        # Return the new datasets
        return datasets
        
    def close(self):
        for ds in self.dataset:
            if ds._status_ != 'closed':
                try: ds.close()
                except: self.exception('Error closing dataset %s', ds.id)
    
    def __del__(self):
        self.close()
        
    def __len__(self):
        return len(self.dataset)
    
    def get_selector(self, time=None, lon=None, lat=None, level=None, merge=True, only=None, split=False, **kwargs):
        """Get a cdms2.selectors.Selector from specified time/lat/lon/level selection
        
        :Params:
        
            - **time/lon/lat/level**, optional: Refine or set the selector with these components. 
            - **only**, optional: Work only on one component, like "time" or "t".
            - **merge**, optional: Merge with selector created at initialization time (stored
              in attribute :attr:`selector`.)
            - **split**, optional: return a splitted selector (see :func:`~vacumm.misc.misc.split_selector`)
        """
        # Only on one component
        if only is not None:
            
            # Check "only"
            only = str(only).lower()
            _geonames = {'x':'lon', 'y':'lat', 'z':'level', 't':'time'}
            allowed = _geonames.keys()+_geonames.values()
            if not only in allowed:
                raise TypeError('Wrong selector type for "only". Please choose on of: '%', '.join(allowed))
            if only in  _geonames.keys():
                only = _geonames[only]
            
            # Get and merge
            myselect = create_selector()
            if merge:
                kw = {only: self.selects[only]}
                for key, val in kw.items():
                    if val is None: del kw[key]
                if kw:
                    myselect.refine(**kw)
            kw = {only: time}
            for key, val in kw.items():
                if val is None: del kw[key]
            if kw:
                myselect.refine(**kw)
        
            return split_selector(myselect) if split else myselect
        
        # Full selector
        myselect = create_selector()
        if merge:
            myselect.refine(self.selector)
        kw = dict(lon=lon, lat=lat, level=level, time=time)
        for key, val in kw.items():
            if val is None: del kw[key]
        if kw:
            myselect.refine(**kw)
        
        return split_selector(myselect) if split else myselect
        
    def get_axis(self, name, select=None, select2=None, dataset=None, warn=True, 
        getid=False, searchmode=None):
        """Retreive a 1D or 2D axis.
        
        :Params:
          - **name**: Generic axis name. 
          - **select** optional: Selection along this axis.
            Only slices are accepted for 2D axes.
          - **select2** optional: Selection on the other dimension
            for 2D axes.
          - **dataset** optional: find axis based on this dataset.
          - **getid**, optional: If True, only return the id (netcdf name) or None.
          - **warn**, optional: Display a warning in case of problem.
          - **searchmode**, optional: Search order (see :func:`~vacumm.misc.io.ncfind_obj`).
        
        :Return:
            - cdms2 axis or None if not found, or id
        """
        
        # Inits
        if not getid: self.debug('Searching axis %s', name)
        if dataset is None and not len(self.dataset):
            if warn: self.warning('No %s found, dataset is empty', name)
            return None
        axis = None
        if dataset is None: dataset = self.dataset[0]
        if not getid: self.debug('Using reference dataset: %s', dataset.id)
            
        # Get specs
        genname, search, base_select, squeeze, atts = self._get_ncobj_specs_(name, 
            attsingle=False, searchmode=searchmode)
        if 'long_names' in atts:
            search['long_names'] = atts['long_names']
        if 'units' in atts and 'axis' in atts and atts['axis']!='Z':
            search['units'] = atts['units']
            
        # Check cache first
        ncname = self._get_cached_ncid_(genname)
        if getid and ncname: return ncname        

        # Search id if not available
        if ncname is None:
            
            # Find
            try:
                ncname = ncfind_obj(dataset, search, regexp=True, searchmode=searchmode)
                if getid: return ncname
                if ncname is None:
                    if warn: self.warning('Generic axis not found with search specs: %s'%search)
                    return
            except:
                if not getid and warn: self.warning('Error searching for generic axis  '
                    'with specs: %s. Message: %s'%(search, format_exc()))
                return
                
        # Cache result
        self._set_cached_ncid_(genname, ncname)
        
        # Read
        if ncname in dataset.listvariables():
            axis = dataset(ncname)
        else:
            axis = dataset.getAxis(ncname)
        axis = axis.clone()
            
        # Format
        if genname is not None:
            axis = format_axis(axis, genname)
        atts = _notuplist_(atts)
        if atts:
            for key, val in atts.items():
                setattr(axis, key, val)
        
        # Select
        if select is not False and select2 is not False:
            
            # Curv case
            if len(axis.shape)==2:
                
                axtype = get_axis_type(axis, genname=True, ro=True, checkatts=False)
                if select is not None or self.selects[axtype] is not None:
                
                    if (genname and genname.startswith('lon')) or islon(axis, ro=True):
                        aname = 'lat'
                        if genname: aname += genname[3:]
                        axis = axis, getattr(self, 'get_'+aname)(False)
                        iaxis = 0
                    elif (genname and genname.startswith('lat')) or islat(axis, ro=True):
                        aname = 'lon'
                        if genname: aname += genname[3:]
                        axis = getattr(self, 'get_'+aname)(False), axis
                        iaxis = 1
                        select, select2 = select2, select
                    axis = self._select_on_axis_(axis, select, select2, 
                        global_select=True, genname=genname)[iaxis]
                    
            else:
                
                axis = self._select_on_axis_(axis, select, select2, global_select=True, genname=genname)
        
        return axis
            
    def _select_on_axis_(self, axis, select1=None, select2=None, global_select=True, genname=None):
        """Global and local selection applied to a 1D or 2D axis
        
        :Params:
        
            - **axis**: cdms2 axis on which to operate, or a tuple of them (2D case).
            - **select1**: local lon selection if axis is lon, else lat selection.
            - **select2**: local lat selection if axis is lon, else lon selection.
            - **global_select**: also apply a global selection (using self.selects).
        
        """
        if select1 is False: return axis
        if not isinstance(axis, tuple) and len(axis.shape)==1: # 1D axis
            
            selects = [select1 or select2]
            axtype = get_axis_type(axis, genname=True, ro=True, checkatts=False)
            if axtype and global_select: selects.insert(0, self.selects[axtype])
        
            for sel in selects:
                if sel is None: continue
                if isinstance(sel, slice):
                    i, j, k = sel.indices(len(axis))
                elif isinstance(sel, tuple):
                    ijk = axis.mapIntervalExt(sel)
                    if ijk is None:
                        raise DatasetError("Can't appy selection %s to axis %s"%(sel, genname or axis))
                    i, j, k = ijk
                else:
                    raise DatasetError("Selection must be a tuple or a slice: %s"%sel)
                axis = axis.subaxis(i, j, k)
            
        else:
            
            # Global selection
            if global_select:
                if self.selects['lon'] is None and self.selects['lat'] is None:
                    global_select = None
                else:
                    global_select = dict(lon = self.selects['lon'] or None, 
                        lat = self.selects['lat'] or None)
            else:
                global_select = None
            
            # Local selection
            if select1 is None and select2 is None: # No selection
                local_select = None
            else:
                lon = lat = None
                axtarget = axis[0] if isinstance(axis, tuple) else axis
                axtype = get_axis_type(axtarget, genname=True, ro=True, checkatts=False)
                if select1 is not None and select2 is not None: # Along X and Y
                    lon, lat = select1, select2
                elif select1 is not None:
                    if axtype=='lon':
                        lon = select1
                    else:
                        lat = select1
                else:
                    if axtype=='lat':
                        lon = select2
                    else:
                        lat = select2
                local_select = dict(lon=lon, lat=lat)
        
            for sel in [global_select, local_select]:
                if sel is None: continue
                if sel['lon'] is None and sel['lat'] is None: continue
                axis = gridsel(axis, **sel)
                
        return axis
        
    def _get_cached_ncid_(self, genname):
        if genname is None: return
        return self._ncids.get(genname, None)
    def _set_cached_ncid_(self, genname, ncid):
        if genname is None: return
        self._ncids[genname] = ncid
    
        
    def get_variable(self, varname, time=None, lon=None, lat=None, level=None, 
        atts=None, squeeze=False, order=None, asvar=None, torect=True, depthup=None,
        verbose=None, warn=True, searchmode=None, **kwargs):
        '''Load a variable in a best time serie fashion.
        
        :Params:
          - **varname**: Generic name, a netcdf name with a '+' a prefix, a tuple of netcdf variable names
            or dictionnary of specification names with a search argument of :func:`~vacumm.misc.io.ncfind_var` (tuple of netcdf names
            or dictionary).
          - **time/lon/lat/level**: Selector components.
          - **squeeze**: If true, call :func:`squeeze_variable` on the returned variable.
          - **order**: If not None, specify the output variable axes order.
          - **depthup**: Make depths up.
          - **torect**: Make grid rectangular if possible.
          - Other kwargs are passed to :func:`~vacumm.misc.io.ncread_files`.
        
        :Return:
            cdms2 variable or None if not found
            
        :Example:
        
            >>> get_variable('ssh', lon=(-10,0))
            >>> get_variable('+xe')
            >>> get_variable(dict(search={'standard_names':'sea_surface_height_above_sea_level'}))
        
        '''
        selstring = 'time: %(time)s, level: %(level)s, lat: %(lat)s, lon: %(lon)s'%locals()
        self.debug('Getting variable: %(varname)s,  %(selstring)s'%locals())
        if not len(self.dataset):
            if warn: self.warning('No %s variable found, dataset is empty', varname)
            return None
            
        # Specifications for searching, selecting and modifying the variable
        genname, search, select, base_squeeze, atts = self._get_ncobj_specs_(varname)
        select.refine(self.get_selector(lon=lon, lat=lat, level=level, time=time, merge=True))
        
        # Search for the variable now
        ncvarid = self._get_cached_ncid_(genname)
        if ncvarid is None:
            ncvarid = ncfind_var(self.dataset[0], search, searchmode=searchmode)
        self._set_cached_ncid_(genname, ncvarid)
        if ncvarid is None:
            if warn: 
                self.warning('Variable not found with search specs: %s (in file %s)', search, self.dataset[0])
            return
            
        # TODO: if grid not found and specs say that there must have a grid, create it with get_lon/get_lat
            
        # Curved grid case
        curvsel = CurvedSelector(self.dataset[0][ncvarid].getGrid(), select)
            
        # Read
        kwargs['torect'] = False
        try:
            var = ncread_files([d.id for d in self.dataset], ncvarid, 
                timeid=self.get_timeid(), 
    #            time=time, 
                select=select, verbose=verbose if verbose is not None else self.is_debug(), 
                squeeze=base_squeeze, **kwargs)
            if var is None:
                if warn: self.warning('No data found for variable: %s, local select: %s', varname, selstring)
                return
        except VACUMMError, e:
            if warn: self.warning('Error reading variable: %s. Message: %s', 
                varname,  e.message)
            return 
        if isinstance(var, list): var = var[0]
            
        # Post process
        var = self.finalize_variable(var, genname=genname, atts=atts, order=order, 
            squeeze=squeeze, asvar=asvar, torect=torect, depthup=depthup, curvsel=curvsel)
        self.verbose('Loaded variable: %s', self.describe(var))
            
        return var
    
    def get(self, varname, **kwargs):
        """Generic way to get a variable
        
        It first tries to find a ``get_*`` method, then call 
        :meth:`get_variable` if no method is found.
        """
        var = None
        if hasattr(self, 'get_'+varname):
            var = getattr(self, 'get_'+varname)(**kwargs)
        if var is None:
            var = self.get_variable(varname, **kwargs)
        return var
    
#    def get_time(self, *args, **kwargs):
#        return self.get_axis('t', *args, **kwargs)
    
    def get_timeid(self, warn=False):
        """Get the time id"""
        # From cache
        timeid = self._get_cached_ncid_('time')
        if timeid is not None: return  timeid
        
        if not len(self.dataset):
            if warn: self.warning('No time found, dataset is empty')
            return None
            
        # Find time id
        genname, search, _, _, _ = self._get_ncobj_specs_('time')
        timeid = ncfind_axis(self.dataset[0], search)
        self._set_cached_ncid_('time', timeid)
        return timeid
        
    
    def get_time(self, time=None, var=None, ids=None, warn=True, **kwargs):
        '''Load time axis in a best time serie fashion.
        
        :Params:
            - **time**: time selector 
        
        '''
        self.debug('Getting time, select: %s', time)
        
        # Find time id
        timeid = self.get_timeid(warn=warn)
        if timeid is None: return
                             
        # Loop on files with first level selection
        allaxes = []
        if timeid:
            
            # Pre and post time selectors
            iter_time = None
            post_times = []
            for this_time in [self.selects.get('time', None), time]:
                if isinstance(this_time, tuple) and iter_time is None:
                    iter_time = this_time
                elif this_time not in [None, False]:
                    post_times.append(this_time)
            
            # Build iterator
            iter = NcIterBestEstimate([d.id for d in self.dataset], iter_time, 
                timeid=timeid, autoclose=False)
                
            # Iterate
            for f,t in iter:
                if t is None: 
                    continue
                if f._status_=='closed':
                    f = cdms2.open(f.id)
                axis = f.getAxis(iter.timeid).clone()
                i, j, k = t.indices(len(axis))
                if i!=0 or j!=len(axis):
                    axis = axis.subaxis(i, j, k)
                allaxes.append(axis)
            iter.close()
            
            # Finalize
            if allaxes :
                
                # Concatenate
                axis =  MV2_axisConcatenate(allaxes, copy=0)
                
                # Final selections
                for this_time in post_times:
                    axis = self._select_on_axis_(axis, this_time, genname='time', 
                        global_select=False)
                
                # Format
                axis = format_axis(axis.clone(), 'time')
                
                return axis
                
            
        self.warning('No time found for select: %s'%time)
    
    def get_ctime(self, *args, **kwargs):
        """Get time axis as a list of :class:`cdtime.comptime`
        
        It is a simple shortcut to:
        
            >>> ds.get_time().asComponentTime()
            
        :Params: All arguments are passed to :meth:`get_time`
        """
        time = self.get_time(*args, **kwargs)
        if time is not None:
            return time.asComponentTime()
    
    def get_time_res(self):
        '''Get the estimated  time resolution, based on the two first time coordinates
        
        :Return:
            - resolution as datetime.timedelta
        
        '''
        time = self.get_time()
        if time is None: return
        ctime = time.asComponentTime()
        return adatetime(ctime[1]) - adatetime(ctime[0])
    
    def get_lon(self, lon=None, lat=None, **kwargs):
        '''Get longitude axis'''
        return self.get_axis('lon', lon, lat, **kwargs)
    get_longitude = get_lon
    
    def get_lon_u(self, lon=None, lat=None, **kwargs):
        '''Get longitude axis at U location'''
        return self.get_axis('lon_u', lon, lat, **kwargs)
        
    def get_lon_v(self, lon=None, lat=None, **kwargs):
        '''Get longitude axis at V location'''
        return self.get_axis('lon_v', lon, lat, **kwargs)
        
    def get_lon_f(self, lon=None, lat=None, **kwargs):
        '''Get longitude axis at F location'''
        return self.get_axis('lon_f', lon, lat, **kwargs)
        
    def get_lat(self, lat=None, lon=None, **kwargs):
        '''Get latitude axis'''
        return self.get_axis('lat', lat, lon, **kwargs)
    get_latitude = get_lat
    
    def get_lat_u(self, lat=None, lon=None, **kwargs):
        '''Get latitude axis at U location'''
        return self.get_axis('lat_u', lat, lon, **kwargs)
    
    def get_lat_v(self, lat=None, lon=None, **kwargs):
        '''Get latitude axis at V location'''
        return self.get_axis('lat_v', lat, lon, **kwargs)
        
    def get_lat_f(self, lat=None, lon=None, **kwargs):
        '''Get latitude axis at F location'''
        return self.get_axis('lat_f', lat, lon, **kwargs)
    
    
    def get_grid(self, lon=None, lat=None):
        '''Get grid'''
        axlon = self.get_lon(lon, lat)
        axlat = self.get_lat(lat, lon)
        if axlon is None or axlat is None: return None
        grid = create_grid(axlon, axlat)
        return grid
        
    def get_grid_u(self, lon=None, lat=None):
        '''Get grid at U location'''
        axlon = self.get_lon_u(lon, lat)
        axlat = self.get_lat_u(lat, lon)
        if axlon is None or axlat is None: return None
        grid = create_grid(axlon, axlat)
        return grid
        
    def get_grid_v(self, lon=None, lat=None):
        '''Get grid at V location'''
        axlon = self.get_lon_v(lon, lat)
        axlat = self.get_lat_v(lat, lon)
        if axlon is None or axlat is None: return None
        grid = create_grid(axlon, axlat)
        return grid
    
    def get_grid_f(self, lon=None, lat=None):
        '''Get grid at F location'''
        axlon = self.get_lon_f(lon, lat)
        axlat = self.get_lat_f(lat, lon)
        if axlon is None or axlat is None: return None
        grid = create_grid(axlon, axlat)
        return grid
        
    def _get_dxy_(self, xy, at='t', lon=None, lat=None, degrees=False, mode=None, 
        local=True, frompt=None, **kwargs):
        dxy = None
        atp = _at_(at, squeezet=True, focus='hor', prefix='_')
        dict_check_defaults(kwargs, warn=False, verbose=False)
        kwgeo = dict(lon=lon, lat=lat)
        kwg = kwgeo if not local else {}
        kwvar = dict(kwargs, **kwgeo)
        rmode = 'raw' if local else "median"
        if frompt is None: frompt = []
        frompt.append((at, 'centered'))
        
        # X or Y
        axis = -1 if xy=='x' else 0
        genname_deg = ['dlat', 'dlon'][axis]
        
        # In meters
        if not degrees:
        
            # Using sizes
            if check_mode('var', mode):
                dxy = self.get_variable('d'+xy+atp, **kwvar)
                if not local and dxy is not None:
                    dxy = N.median(dxy.compressed())
                if dxy is not None or check_mode('var', mode, strict=True): return dxy
        
            # Using coordinates
            if check_mode('coord', mode):
                
                # Estimate 
                dxy = self._coord2res_(axis, self._coord2dxy_, rmode, frompt, kwg, kwgeo)
                
                # Finalize
                return self._finalize_dxy_(axis, atp, dxy, 'd'+xy, kwgeo, **kwargs)
                
            return
        
        # In degrees
        # - estimate
        dxy = self._coord2res_(axis, self._coord2dlonlat_, rmode, frompt, kwg, kwgeo)
        # - finalize
        return self._finalize_dxy_(axis, atp, dxy, genname_deg, kwgeo, **kwargs)
    
    def _coord2res_(self, axis, func, rmode, frompt, kwg, kwgeo):
        """Estimate resolution in meters or degrees from coordinates at several locations"""
        dxy = None
        if frompt: # [('u', -1)]
            for pt, transforms in frompt:
                
                # Compute it
                dxy = func(axis, pt, rmode, **kwg)
                if dxy is None: continue
                if rmode!='raw': return dxy
                
                # Loop on transforms
                for transform in transforms:
                
                    # Centered estimation
                    if transform=='centered': 
                    
                        sh = list(dxy.shape)
                        sh[axis]+=1
                        dxy_ = N.zeros(sh)
                        sl = get_axis_slices(sh, axis)
                        dxy_[sl['firsts']] += shift1d(dxy, shift=-1, axis=axis)
                        dxy_[sl['last']] += shift1d(dxy, shift=1, axis=axis)[sl['last']]
                        dxy = dxy_
                    
                    else: # Using specified transforms
                        
                        transfunc, kw = transform
                        dxy = transfunc(dxy, **kw)
                        
                return dxy
                        
    
    def _coord2dxy_(self, axis, pt, rmode, **kwgeo):
        """Estimate resolution in meters from coordinates"""
        atp = _at_(pt, squeezet=True, focus='hor', prefix='_')
        grid = getattr(self, 'get_grid'+atp)(**kwgeo) # lat is also needed for projection
        if grid is None: return
        return  resol(grid, proj=True, mode=rmode, axis=axis)[1+axis]   
        
    def _coord2dlonlat_(self, axis, pt, rmode, **kwgeo):
        """Estimate resolution in degrees from coordinates"""
        atp = _at_(pt, squeezet=True, focus='hor', prefix='_')
        name = ['lat', 'lon'][axis]
        coord = getattr(self, 'get_'+name+atp)(**kwgeo)
        if coord is None: return
        return resol(coord, proj=False, mode=rmode, axis=-1)
        
    def _finalize_dxy_(self, axis, atp, dxy, genname, kwgeo, warn=True, **kwargs):
        """Reshape, put on grid and finalize as other variables"""
        # Nothing
        if dxy is None:
            if warn:
                self.warning("Can't estimate resolution from coordinates")
            return
        if N.isscalar(dxy): return dxy
        
        # Put on grid
        dxy = MV2.asarray(dxy)
        grid = getattr(self, 'get_grid'+atp)(**kwgeo)
        if dxy.ndim==1:
            if axis==0:
                dxy = MV2.transpose(N.ma.resize(dxy.reshape(1, -1), grid.shape[::-1]))
            else:
                dxy = MV2.resize(dxy, grid.shape)
        set_grid(dxy, grid)
        
        # Finalize
        return self.finalize_variable(dxy, genname=genname+atp, **kwargs)
        
    
    def get_dx(self, **kwargs):
        """Get the grid resolution along X
        
        It can be stored as an variable or computed from coordinates.
        """
        frompt = [('u', [(extend1d, {'ext':(1, 0), 'axis':-1})])]
        return self._get_dxy_('x', 't', frompt=frompt, **kwargs)
        
    def get_dx_u(self, **kwargs):
        """Get the grid resolution along X at U location"""
        frompt = [('t', [(extend1d, {'ext':(0, 1), 'axis':-1})])]
        return self._get_dxy_('x', 'u', frompt=frompt, **kwargs)
        
    def get_dx_v(self, **kwargs):
        """Get the grid resolution along X at V location"""
        frompt = [
            ('u', 
                [(extend1d, {'ext':(0, 1), 'axis':-1}), 
                (shift1d, {'shift':1, 'axis':0})]),  
            ('t', 
                ['centered', 
                (shift1d, {'shift':1, 'axis':0})]),
            ]
        return  self._get_dxy_('x', 'v', frompt=frompt, **kwargs)
        
    def get_dy(self, **kwargs):
        """Get the grid resolution along Y
        
        It can be stored as an variable or computed from coordinates.
        """
        frompt = [('u', [(extend1d, {'ext':(1, 0), 'axis':0})])]
        return self._get_dxy_('y', 't', frompt=frompt, **kwargs)
        
    def get_dy_v(self, **kwargs):
        """Get the grid resolution along Y at V location"""
        frompt = [('t', [(extend1d, {'ext':(0, 1), 'axis':0})])]
        return self._get_dxy_('y', 'v', frompt=frompt, **kwargs)
        
    def get_dy_u(self, **kwargs):
        """Get the grid resolution along Y at U location"""
        frompt = [
            ('u', 
                [(extend1d, {'ext':(0, 1), 'axis':0}), 
                (shift1d, {'shift':1, 'axis':-1})]),  
            ('t', 
                ['centered', 
                (shift1d, {'shift':1, 'axis':-1})]),
            ]
        return  self._get_dxy_('y', 'u', frompt=frompt, **kwargs)
        
    def get_resol(self, degrees=False, at='t', mode=None, **kwargs):
        """Get the horizontal grid resolutions
        
        :Params:
        
            - **degrees**, optional: In degrees or meters?
            - **local**, optional: Get resolution at each point? 
            
        :Return: ``dx,dy``
        
        :See also: :meth:`get_dx` :meth:`get_dy` :meth:`get_grid` :func:`~vacumm.misc.grid.misc.resol`
        """
        dx = dy = None
        
        # In meters
        if not degrees:
            
            # Using sizes
            kw = dict(local=False, verbose=False, mode='var', degrees=False)
            kw.update(kwargs)
            try: dx = self._get_dx_(at, **kw)
            except: pass
            try: dy = self._get_dy_(at, **kw)
            except: pass
            if dx is not None and dy is not None: return dx, dy
            
        # Using coordinates (degrees or meters)
        grid = self.get_grid(**kwargs)
        if grid is None:
            dx2 = dy2 = None
        else:
            dx2, dy2 =  resol(grid, proj=not degrees, mode="median")
        res = dx or dx2, dy or dy2
        if res[0] is None or res[1] is None:
            self.warning("Can't properly estimate resolution along both X and Y")
        return res
        
                        
        
    def get_level(self, level=None, **kwargs):
        '''Get level axis, based on :func:`get_axis`'''
        return self.get_axis('depth', level=level, **kwargs)
    
    @formatdoc_var
    def get_corio(self, **kwargs):
        '''Get Coriolis parameter'''
        return self.get_variable('corio', **kwargs)
    
    @formatdoc_var
    def get_corio_u(self, **kwargs):
        '''Get Coriolis parameter'''
        return self.get_variable('corio_u', **kwargs)
    
    @formatdoc_var
    def get_corio_v(self, **kwargs):
        '''Get Coriolis parameter'''
        return self.get_variable('corio_v', **kwargs)
    
    
    @classmethod
    def squeeze_variable(cls, var, spec=True):
        '''Squeeze a variable, preserving remaining axis
        
        
        :See also: :func:`vacumm.misc.misc.squeeze_variable`
        
        '''
        cls.debug('Squeezing variable: %s', cls.describe(var))
        
        return squeeze_variable(var, spec)
        #axes = [a for a in var.getAxisList() if a.shape[0] > 1]
        #var = var.squeeze()
        #if numpy.size(var):
            #var.setAxisList(axes)
        #return var
    
    def finalize_variable(self, var, genname=None, squeeze=False, order=None, 
        asvar=None, torect=True, atts=None, lon=None, lat=None, curvsel=None, **kwargs):
        """Finalize a variable
        
        :Params:
        
            - **genname**, optional: Generic name to format the variable using
              :func:`~vacumm.data.cf.format_var`.
            - **atts**, optional: Set these attributes to the variable.
            - **squeeze**, optional: If not False, squeeze singletons axes using
              :func:`~vacumm.misc.misc.squeeze_variable`.
            - **order**, optional: If not None, change the axes order of the variable.
              It must contains letters like 'txyz-'.
            - **asvar**, optional: Grow variable to match the ``asvar`` variable,
              using :func:`~vacumm.misc.misc.grow_variables`.
            - **asvar_<param>**: Param passed to :func:`~vacumm.misc.misc.grow_variables`.
            - **torect**, optinal: Try to convert curvilinear grid to rectangular
              grid using :func:`~vacumm.misc.grid.misc.curv2rect`.
            - **lon/lat**, optional: Additional spatial selection.
            - **curvsel**, optional: :class:`CurvedSelector` instance.
        """
        if var is None: return
        kwasvar = kwfilter(kwargs, 'asvar_')
        ax = isaxis(var)
        if atts is None: atts = {}
        
        if not isaxis(var):
            
            # Curved selector
            if curvsel is not None:
                var = curvsel.finalize(var)
            
            # Format
            if genname is not None:
                format_var(var, genname, **atts)
            elif atts:
                set_atts(var, **atts)
        
            # Squeeze singletons.
            if squeeze:
                self.debug('Squeezing variable: %s', self.describe(var))
                var = squeeze_variable(var, squeeze)
            
            # Reorder
            if order:
                self.debug('Reordering variable: %s', self.describe(var))
                var = var(order=order)
                
            # Rectangular if possible
            if torect: self._torect_(var)
            
        elif genname is not None:
            format_axis(var, genname, **atts)
        elif atts:
            set_atts(var, **atts)
            
        
        # Grow variables or axes
        if asvar is not None:
            if not cdms2.isVariable(asvar):
                asvar = self.get(asvar, warn=False)
            return grow_variables(asvar, var, **kwasvar)[1]
           
        # Post spatial selection
        if lon is not None or lat is not None:
            var = var(create_selector(lon=lon, lat=lat))
            
        return var
    
    def _torect_(self, var):
        """Convert from curv 2 rect
        
        :Params: var =  variable or grid
        """
        # Get the grid
        grid = var if isgrid(var) else var.getGrid()
        if grid is None: return var
        
        # Ids
        lon = grid.getLongitude()
        lat = grid.getLatitude()
        lonid = getattr(lon, '_oldid', lon.id)
        latid = getattr(lat, '_oldid', lat.id)
        ids = (lonid, latid)
        
        # Check cache
        if not hasattr(self, '_rgrids'):
            self._rgrids = {}
        if ids in self._rgrids:
#            print 'datset: got rect from cache'
            if isgrid(var):
                return self._rgrids[ids]
            set_grid(var, self._rgrids[ids])
            return var
                    
        # Estimate it
        var = curv2rect(var, f=self.dataset[0], mode='none')
        
        # Cache it if rectangular
        grid = var if isgrid(var) else var.getGrid()
        if isgrid(grid, curv=True):
#            print 'dataset: we cache rect'
            self._rgrids[ids] = grid
            
        return var
            
            
        
        
    ############################################################################
    ############################################################################
    ############################################################################

class AtmosSurfaceDataset(Dataset):

    ncobj_specs = {}

    @formatdoc_var
    def get_senhf(self, *args, **kwargs):
        '''Get Sensible Heat Flux'''
        return self.get_variable('senhf', *args, **kwargs)

    @formatdoc_var
    def get_lathf(self, *args, **kwargs):
        '''Get Latent Heat Flux'''
        return self.get_variable('lathf', *args, **kwargs)

    @formatdoc_var
    def get_swhf(self, *args, **kwargs):
        '''Get Short Wave Heat Flux'''
        return self.get_variable('swhf', *args, **kwargs)

    @formatdoc_var
    def get_lwhf(self, *args, **kwargs):
        '''Get Long Wave Heat Flux'''
        return self.get_variable('lwhf', *args, **kwargs)

    @formatdoc_var
    def get_evap(self, *args, **kwargs):
        '''Get Evaporation'''
        return self.get_variable('evap', *args, **kwargs)

    @formatdoc_var
    def get_rain(self, *args, **kwargs):
        '''Get Evaporation'''
        return self.get_variable('rain', *args, **kwargs)

    @formatdoc_var
    def get_lhf(self, *args, **kwargs):
        '''Get Latent Heat Flux'''
        return self.get_variable('lhf', *args, **kwargs)
    
    @formatdoc_var
    def get_taux(self, *args, **kwargs):
        '''Get zonal wind stress (westerly)'''
        return self.get_variable('taux', *args, **kwargs)
    
    @formatdoc_var
    def get_tauy(self, *args, **kwargs):
        '''Get meridional wind stress (westerly)'''
        return self.get_variable('tauy', *args, **kwargs)
    
    @formatdoc_var
    def get_u10m(self, *args, **kwargs):
        '''Get 10-m zonal wind speed (westerly)'''
        return self.get_variable('u10m', *args, **kwargs)
    
    @formatdoc_var
    def get_v10m(self, *args, **kwargs):
        '''Get 10-m meridional wind speed (westerly)'''
        return self.get_variable('v10m', *args, **kwargs)
    
class OceanSurfaceDataset(Dataset):
    
    ncobj_specs = {}
#    dict(
#        sst = dict(
#            search = dict(standard_names='sea_surface_temperature'), 
#        ), 
#    )
    @formatdoc_var
    def get_sst(self, **kwargs):
        '''Get SST'''
        return self.get_variable('sst', **kwargs)
    
    @formatdoc_var
    def get_sss(self, **kwargs):
        '''Get SSS'''
        return self.get_variable('sss', **kwargs)
    
    @formatdoc_var
    def get_ssh(self, **kwargs):
        '''Get sea surface height'''
        return self.get_variable('ssh', **kwargs)
        
    @formatdoc_var
    def get_usurf(self, **kwargs):
        '''Get U at the sea surface'''
        return self.get_variable('usurf', **kwargs)
    
    @formatdoc_var
    def get_vsurf(self, **kwargs):
        '''Get V at the sea surface'''
        return self.get_variable('vsurf', **kwargs)
    
    @formatdoc_var
    def get_hs(self, **kwargs):
        '''Get significant height of waves'''
        return self.get_variable('hs',  **kwargs)

class OceanDataset(OceanSurfaceDataset):    
    ncobj_specs = OceanSurfaceDataset.ncobj_specs.copy()
    ncobj_specs.update(dict(
    ))
    
    def finalize_variable(self, var, squeeze=False, order=None, asvar=None, torect=True, depthup=None, **kwargs):
        """Finalize a variable
        
        :Params:
        
            - **squeeze**, optional: If not False, squeeze singletons axes using
              :func:`~vacumm.misc.misc.squeeze_variable`.
            - **order**, optional: If not None, change the axes order of the variable.
              It must contains letters like 'txyz-'.
            - **asvar**, optional: Grow variable to match the ``asvar`` variable,
              using :func:`~vacumm.misc.misc.grow_variables`.
            - **asvar_<param>**: Param passed to :func:`~vacumm.misc.misc.grow_variables`.
            - **torect**, optinal: Try to convert curvilinear grid to rectangular
              grid using :func:`~vacumm.misc.grid.misc.curv2rect`.
            - **depthup**, optional: If not False, try the make depth positive up
              using :func:`~vacumm.misc.grid.misc.makedepthup`.
        """
        if var is None: return
        # Make depth positive up
        if depthup is not False:
            self._makedepthup_(var, depth=depthup)
        
        # Generic stuff
        var = Dataset.finalize_variable(self, var, squeeze=squeeze, order=order,
            torect=torect, asvar=asvar, **kwargs)
        
        return var
    

    def _get_dz_(self, at='t', warn=True, mode=None, **kwargs):
        """Get layer thickness"""
        atp = _at_(at, squeezet=True, prefix=True)
        
        # First, find a variable
        if check_mode('var', mode):
            dz = self.get_variable('dz'+atp, warn=False, **kwargs)
            if dz is not None or check_mode('var', mode, strict=True): return dz
        
        # Second, guess from vertical coordinates
        if check_mode('depth', mode):
            # - read depths
            dz2dmode = None
            if at in 'tr':
                depth = self._get_depth_('w', warn=False, mode='-dz', **kwargs) # from W points
                if depth is None:
                    depth = self._get_depth_('t', warn=False, mode='-dz', **kwargs)  # from T points
                    dz2dmode = 'center'
            elif at=='w':
                depth = self._get_depth_('t', warn=False, mode='-dz', **kwargs) # from T points
                if depth is None:
                    depth = self._get_depth_('w', warn=False, mode='-dz', **kwargs) # from W points
                    dz2dmode = 'center'
            else:
                if warn: self.warning("Computing dz from depths at points other than T and W is yet not implemented")
                return
            if depth is None:
                if warn: self.warning("Can't get depth to compute dz")
                return
            # - compute thicknesses 
            dz = format_var(depth.clone(), 'dz'+atp)
            isup = self._isdepthup_(depth)
            if at=='': # At T from W depth
                dz2dmode= 'top' if isup else 'bottom'
            else: # At W from T depth
                dz2dmode = 'bottom' if isup else 'top'
            dz[:] = depth2dz(depth, mode=dz2dmode, masked=False)
            return dz
        
    _mode_doc = """Computing mode
    
          - ``None``: Try all modes, in the following order.
          - ``"var"``: Read it from a variable.
          - ``"depth"``: Estimate from depth (see :meth:`get_depth`)
          
          You can also negate the search with
          a '-' sigme before: ``"-depth"``."""
    def get_dz(self, *args, **kwargs):
        """Get layer thickness at"""
        return self._get_dz_('t', *args, **kwargs)
    get_dz = formatdoc_var(get_dz, mode=_mode_doc)
        
    def get_dz_u(self, *args, **kwargs):
        """Get layer thickness at U location"""
        return self._get_dz_('u', *args, **kwargs)
    formatdoc_var(get_dz_u, mode=_mode_doc)
        
    def get_dz_v(self, *args, **kwargs):
        """Get layer thickness at V location"""
        return self._get_dz_('v', *args, **kwargs)
    formatdoc_var(get_dz_v, mode=_mode_doc)
        
    def get_dz_w(self, *args, **kwargs):
        """Get layer thickness at W location"""
        return self._get_dz_('w', *args, **kwargs)
    formatdoc_var(get_dz_w, mode=_mode_doc)
    
    del _mode_doc
        
    def _isdepthup_(self, depth=None):
        """Guess if depths are positive up"""
        # Cache
        if getattr(self, 'positive', None) is not None:
            return self.positive=='up'
        
        # Get depth
        if depth is None:
            depth = self.get_depth(time=slice(0, 1), warn=False)
        if depth is None: # no depth = no problem
            self.positive = 'up'
            return True
            
        # Guess
        axis = 0 if len(depth.shape)==1 else 1
        isup = isdepthup(depth, ro=False, axis=axis)
        self.positive = 'up' if isup else 'down'
        return isup
    
    def _makedepthup_(self, var, depth=None):
        """Make depths positive up"""
        isup = self._isdepthup_(depth)
        if isup: return var
        if isdep(var):
            return makedepthup(var, depth=False, strict=True)
        axis = var.getOrder().find('z')
        if axis<0: return var
        return makedepthup(var, depth=False, axis=axis, strict=True)
    
        
    def _get_depth_(self, at='t', level=None, time=None, lat=None, lon=None, 
        order=None, squeeze=None, asvar=None, torect=True, warn=True, mode=None, **kwargs):

        # Where?
        atp = _at_(at, squeezet=True, prefix=True)
        ath = _at_(at, squeezet=True, prefix=True, focus='hor')
        
        # First, try to find a depth variable
        kwfinal = dict(order=order, squeeze=squeeze, asvar=asvar, torect=torect)
        kwvar = dict(level=level, time=time, lat=lat, lon=lon, warn=False)
        kwvar.update(kwfinal)
        if check_mode('var', mode):
            depth = self.get_variable('depth'+atp, depth=False, **kwvar)
            if depth is not None or check_mode('var', mode, strict=True): 
                return self._makedepthup_(depth, depth)

        # Get selector for other tries
        selector = self.get_selector(lon=lon, lat=lat, level=level, time=time, merge=True) 
        gridmet = 'get_grid'+ath
        grid = getattr(self, gridmet)(False)
        curvsel = CurvedSelector(grid, selector)
        kwfinal['curvsel'] = curvsel
        
        # Second, find sigma coordinates
        sigma_converter = NcSigma.factory(self.dataset[0])
        if check_mode('sigma', mode):
            if sigma_converter is not None and sigma_converter.stype is None: 
                sigma_converter.close()
                sigma_converter = None
            if sigma_converter is not None:
                self.debug('Found depth refer to a sigma level, processing sigma to depth conversion')
                allvars = []
                for f,t in NcIterBestEstimate(self.dataset, self.selects['time']):
                    if t is False: continue # and when no time??? None-> ok we continue
                    #sigma = NcSigma.factory(f)
                    sel = selector(time=t)
                    self.debug('- dataset: %s: sigma: %s, select: %s', os.path.basename(f.id), sigma_converter.__class__.__name__, sel)
                    try:
                        d = sigma_converter(sel, at=at, copyaxes=True, mode='sigma')
                    except Exception, e:
                        if warn: 
                            self.warning("Can't get depth from sigma. Error: "+e.message)
                        break
                    self.debug('Sigma to depth result: %s', self.describe(d))
                    allvars.append(d)
                if allvars:
                    # Concatenate loaded depth
                    var = format_var(MV2_concatenate(allvars), 'depth'+atp)
                    return self.finalize_variable(var, depthup=var, **kwfinal)
                        
                if check_mode('sigma', mode, strict=True): return
#        if sigma_converter is not None:
#            sigma_converter.close()
            
        # Third, estimate from layer thickness
        if check_mode('dz', mode):
            
            if atp=='': # at T-location
            
                bathy = self.get_bathy(**kwvar)
                dzw = self.get_dz_w(mode='var', **kwvar)
                if bathy is not None and dzw is not None:
                    depth = dz2depth(dzw, bathy, refloc="bottom")
                    
            elif atp=='w': # at W-location
            
                ssh = self.get_ssh(**kwvar)
                dzt = self.get_dz(mode='var', **kwvar)
                if ssh is not None and dzt is not None:
                    depth = dz2depth(dzt, ssh, refloc="top") 
                    
            if depth is not None or check_mode('dz', mode, strict=True):
                var =  format_var(depth, 'depth'+atp)
                return self.finalize_variable(var, depthup=False, **kwfinal)
                        
        # Finally, find a depth axis
        if sigma_converter is None and check_mode('axis', mode): # no Z axis for sigma coordinates
            axis = self.get_axis('depth'+atp, level)
            if axis is not None:
                format_axis(axis, 'depth'+atp), 
                return self.finalize_variable(axis, depthup=axis, **kwfinal)
        if warn:
            self.warning('Found no way to estimate depths at %s location'%at.upper())
        
        
        
    _mode_doc = """Computing mode
    
          - ``None``: Try all modes, in the following order.
          - ``"var"``: Read it from a variable.
          - ``"sigma"``: Estimate from sigma coordinates.
          - ``"dz"``: Estimate from layer thinknesses (see :meth:`get_dz`)
          - ``"axis"``: Read it from an axis (if not sigma coordinates)
          
          You can specifiy a list of them: ``['dz', 'sigma']`` 
          You can also negate the search with
          a '-' sigme before: ``"-dz"``."""
    def get_depth(self, *args, **kwargs):
        return self._get_depth_('t', *args, **kwargs)
    formatdoc_var(get_depth, mode=_mode_doc)
        
    def get_depth_w(self, *args, **kwargs):
        return self._get_depth_('w', *args, **kwargs)
    formatdoc_var(get_depth_w, mode=_mode_doc)
    del _mode_doc

    @formatdoc_var
    def get_temp(self, **kwargs):
        '''Get 4D temperature'''
        return self.get_variable('temp',  **kwargs)
    
    @formatdoc_var
    def get_kz(self, **kwargs):
        '''Get 4D vertical diffusion of tracers (kz)'''
        return self.get_variable('kz',  **kwargs)
    
    @formatdoc_var
    def get_sal(self, **kwargs):
        '''Get 4D salinity'''
        return self.get_variable('sal', **kwargs)
        
    _mode_doc = """Computing mode
    
          - ``None``: Try all modes, in the following order.
          - ``"var"``: Read it from a variable.
          - ``"tempsal"``: Estimate from temperature and salinity."""
    def get_dens(self, mode=None, **kwargs):
        '''Get 4D density'''
        
        kwvar = kwfilter(kwargs, ['lon','lat','time','level','torect'], warn=False)
        kwfinal = kwfilter(kwargs, ['squeeze','order','asvar'])
        kwdens = kwfilter(kwargs, 'dens_')
        kwdepth = kwfilter(kwargs, 'depth_')
        kwdepth.update(kwvar)
        
        # First, try to find a dens variable
        if check_mode('var', mode):
            dens = self.get_variable('dens', **kwvar)
            if dens is not None or check_mode('var', mode, strict=True): 
                return self.finalize_variable(dens, **kwfinal)
            
        # Estimate it from temperature and salinity
        if check_mode('tempsal', mode):
            temp = self.get_temp(**kwvar)
            sal = self.get_sal(**kwvar)
            depth = self.get_depth(**kwdepth)
            if temp is not None and sal is not None:
                dens = density(temp, sal, depth=depth, format_axes=True, **kwdens)
            if dens is not None or check_mode('tempsal', mode, strict=True): 
                return self.finalize_variable(dens, depthup=False, **kwfinal)
    formatdoc_var(get_dens, mode=_mode_doc)
        
    @formatdoc_var
    def get_u3d(self, **kwargs):
        '''Get 4D zonal velocity'''
        return self.get_variable('u3d', **kwargs)
    
    @formatdoc_var
    def get_v3d(self, **kwargs):
        '''Get 4D meridional velocity'''
        return self.get_variable('v3d', **kwargs)
     
    @formatdoc_var
    def get_ubt(self, **kwargs):
        '''Get zonal barotropic velocity'''
        return self.get_variable('ubt', **kwargs)
    
    @formatdoc_var
    def get_vbt(self, **kwargs):
        '''Get meridional barotropic velocity'''
        return self.get_variable('vbt', **kwargs)
        
    @formatdoc_var
    def get_bathy(self, **kwargs):
        '''Get bathymetry'''
        return self.get_variable('bathy', **kwargs)

    @formatdoc_var
    def get_bathy_u(self, **kwargs):
        '''Get bathymetry'''
        return self.get_variable('bathy_u', **kwargs)
        
    @formatdoc_var
    def get_bathy_v(self, **kwargs):
        '''Get bathymetry'''
        return self.get_variable('bathy_v', **kwargs)
        

    def get_uvgbt(self, **kwargs):
        """Get zonal and meridional geostrophic velocity from SSH"""
        ssh = self.get_ssh(**kwargs)
        dxu = self.get_dx_u(degrees=False, local=True, **kwargs)
        dyv = self.get_dy_v(degrees=False, local=True, **kwargs)
        return barotropic_geostrophic_velocity(ssh, dxy=(dxu, dyv))
    
    @formatdoc_var
    def get_ugbt(self, **kwargs):
        '''Get zonal barotropic geostrophic velocity from SSH'''
        kwargs['getv'] = False
        return get_uvbt(**kwargs)
    
    @formatdoc_var
    def get_vgbt(self, **kwargs):
        '''Get meridional barotropic geostrophic velocity from SSH'''
        kwargs['getu'] = False
        return get_uvbt(**kwargs)
        
    def get_eke(self, **kwargs):
        """Get eddy kinetic energy"
        
        :See also: :func:`~vacumm.diag.dynamics.eddy_kinetic_energy`
            :meth:`get_uvgbt`
        """
        eke = self.get_variable('eke', verbose=False, warn=False, **kwargs)
        if eke is None:
            ugbt, vgbt = self.get_uvgbt(**kwargs)
            return eddy_kinetic_energy((ugbt, vgbt))
        
    
    ############################################################################
    ############################################################################
    ############################################################################
    
    def get_hsection(self, varname, depth, time=None, lat=None, lon=None, 
        timeavg=False, warn=True, **kwargs):
        """Get a horizontal section of a variable for a specified depth
        
        :Params:
        
            - **varname**: Generic var name.
            - **depth**: Target depth.
            - **timeavg**, optional: Time average of results.
            - **depth_<param>**, optional: Param is passed to :meth:`get_depth`.
            - **interp_<param>**, optional: Param is passed to 
              :meth:`~vacumm.misc.grid.regridding.interp1d`.
        """
        
        # Params
        kwargs.pop('level', None)
        kwvar = dict(time=time, lat=lat, lon=lon, torect=kwargs.pop('torect', True))
        kwdepth = kwfilter(kwargs, 'depth_')
        kwdepth.update(kwvar)
        kwinterp = kwfilter(kwargs, 'interp_', axis=1)
        
        # Get data
        var = self.get(varname, **kwvar)
        vardepth = self.get_depth(**kwdepth)
        if var is None or vardepth is None: 
            if warn: self.warning('Cannot get var or depths for hsection')
            return
        
        # Interpolate
        depth = -N.abs(depth)
        zo = create_dep([depth] if N.isscalar(depth) else depth[:])
        if len(vardepth.shape)>1:
            var, vardepth = grow_variables(var, vardepth)
            kwinterp.update(xmap=[0,2,3], xmapper=vardepth(order='...z').filled(0.))
        lvar = interp1d(var, zo, **kwinterp)
        
        # Time average
        if timeavg and lvar.getOrder().startswith('t') and lvar.shape[0]>1:
            lvar = MV2.average(lvar, axis=0)
        
        return self.finalize_variable(lvar, depthup=False, **kwargs)
        
        
    def plot_hsection(self, varname, depth, time=None, lat=None, lon=None, 
        timeavg=False, title='%(long_name)s at %(depth)im', **kwargs):
        """Plot a horizontal section of a variable for a specified depth
        
        Section is computed with :meth:`get_hsection`.
        
        :Params:
        
            - **varname**: Generic var name.
            - **depth**: Target depth.
            - **timeavg**, optional: Time average of results.
            - **depth_<param>**, optional: Param is passed to :meth:`get_depth`.
            - **interp_<param>**, optional: Param is passed to 
              :meth:`~vacumm.misc.grid.regridding.interp1d`.
            - Other arguments are passed to :func:`~vacumm.misc.plot.map2`.
        """
        # Get section
        depth = -N.abs(depth)
        var = self.get_hsection(varname, depth, time=time, lat=lat, lon=lon,
            timeavg=timeavg, squeeze=True)
            
        # Plot it
        long_name = getattr(var, 'long_name', '')
        units = getattr(var, 'units', '')
        dict_check_defaults(kwargs, bgcolor='0.5')
        return map2(var, title=title%locals(), **kwargs)
            
    
    def get_layer(self, varname, depth, select=None, timeavg=True):
        '''Get a layer of variable for a specified depth.
        
        :Params:
            - **varname**: variable to process
            - **depth**: output depth(s)
            - **select**: selector
            - **timeavg**: if true, average date along time if needed
          
        :Return:
            - var([time],latitude,longitude): layer variable
        
        '''
        self.verbose('Getting layer data'
            '\n  variable:  %s'
            '\n  depth:      %s'
            '\n  select:      %s'
            '\n  timeavg:      %s', varname, depth, select, timeavg)
        var = self.get_variable(varname, select=select)#, order='tzyx')
        if var is None: raise Exception('Variable %s not found'%varname)
        var = var.reorder('tzyx')
        #if var.getOrder() not in ('tzyx','zyx'):
        #    raise ValueError('Invalid var order: %s, [t]zyx expected'%(var.getOrder()))
        dep = self.get_depth(varname, select=select)#, order='tzyx')
        if dep is None: raise Exception('Depth not found')
        dep = dep.reorder('tzyx')
        #if dep.getOrder() not in ('tzyx','zyx'):
        #    raise ValueError('Invalid depth order: %s, [t]zyx expected'%(dep.getOrder()))
        # Time step interpolate function
        def layer(var, dep):
            self.logdesc(var, dep, title='load: ')
            # Depth interpolation axis
            odep = create_dep(is_iterable(depth) and depth or [depth])
            method='auto'
            dep = dep.reorder('yxz')
            var = var.reorder('yxz')
            var = interp1d(var, axo=odep, xmap=[0,1], xmapper=dep, method=method)
            var = var.reorder('yx')
            var = var.squeeze()
            self.logdesc(var, title='interp: ')
            return var
        tax = None
        axes = var.getAxisList()[-2:]
        if 't' in var.getOrder():
            self.verbose('Computing layers along time')
            tmp = []
            time = var.getTime()
            ctime = time.asComponentTime()
            for it in xrange(len(time)):
                self.info('- Interpolating time: %s', ctime[it])
                if 't' in dep.getOrder(): di = dep[it]
                else: di = dep
                tmp.append(layer(var[it], di))
            if timeavg:
                self.verbose('Averaging along time')
                tmp = MV2.average(tmp, axis=0)
            else: tax = time.clone()
        else:
            tmp = layer(var, dep)
        if tax is not None: axes = [tax] + axes
        var = cdms2.createVariable(tmp, id=var.id, axes=axes, attributes=var.attributes.copy())
        self.logdesc(var, title='out: ')
        return var
    
    def get_section(self, varname, xmin, ymin, xmax, ymax, select=None, timeavg=True):
        '''Get a vertical section data along a straight trajectory, not necessary zonal or meridional.
        
        :Params:
            - **varname**: variable to process
            - **xmin**: westernmost longitude coordinate
            - **ymin**: southernmost latitude coordinate
            - **xmax**: eastermost longitude coordinate
            - **ymax**: northernmost latitude coordinate
            - **select**: selector
            - **timeavg**: if true, average date along time if needed
          
        :Return: A list containing in order:
            - var(level,position): section variable
            - depth(level): depth corresponding to var's level
            - latitude(position): latitude corresponding to var's position
            - longitude(position): longitude corresponding to var's position
        
        '''
        self.verbose('Getting section data'
            '\n  variable:  %s'
            '\n  xmin:      %s'
            '\n  ymin:      %s'
            '\n  xmax:      %s'
            '\n  ymax:      %s'
            '\n  select:    %s'
            '\n  timeavg:    %s',varname, xmin, ymin, xmax, ymax, select, timeavg)
        grid = self.get_grid(varname)
        if grid is None:
            raise Exception('No grid found')
        xres, yres = resol(grid)
        select = self.get_selector(select)
        select.update(dict(longitude=(xmin-xres, xmax+xres, 'ccb'), latitude=(ymin-yres, ymax+yres, 'ccb')))
        var = self.get_variable(varname, select=select, squeeze=True, order='tzyx')
        if var is None:
            raise Exception('No %s found'%(varname))
        if var.getOrder() not in ('tzyx','zyx'):
            raise ValueError('Invalid var order: %s, [t]zyx expected'%(var.getOrder()))
        depth = self.get_depth(varname, select=select, squeeze=True, order='tzyx')
        if depth is None:
            raise Exception('No depth found')
        if depth.getOrder() not in ('tzyx','zyx'):
            raise ValueError('Invalid depth order: %s, [t]zyx expected'%(depth.getOrder()))
        #depth.setAxisList(var.getAxisList())
        # calc xy section coords
        nx = (xmax-xmin)/xres
        ny = (xmax-xmin)/xres
        nxy = int(numpy.ceil(numpy.sqrt(nx**2+ny**2)))*1.5
        xo = numpy.linspace(xmin, xmax, nxy)
        yo = numpy.linspace(ymin, ymax, nxy)
        def section(var, depth):
            # interp var
            vo = grid2xy(var, xo, yo, method='nat')
            # interp depth
            do = grid2xy(depth, xo, yo, method='nat')
            ret = vo,do
            self.verbose('\n  %s', '\n  '.join(self.describe(o) for o in ret))
            return ret
        if 't' in var.getOrder():
            self.verbose('Computing section along time')
            tmp, data = [], []
            time = var.getTime()
            ctime = time.asComponentTime()
            for it in xrange(len(time)):
                self.verbose('Interpolating time: %s', ctime[it])
                if 't' in depth.getOrder(): di = depth[it]
                else: di = depth
                tmp.append(section(var[it], di))
            for i in xrange(len(tmp[0])):
                d = [d[i] for d in tmp]
                data.append(MV2.concatenate([d]))
            vo, do = data
            tax = time.clone()
        else:
            vo, do = section(var, depth)
            tax = None
        xo, yo = create_lon(xo), create_lat(yo)
        pos = cdms2.createAxis(range(vo.shape[-1]), id='position')
        lev = cdms2.createAxis(range(vo.shape[-2]), id='level')
        lev.designateLevel()
        axes = [lev, pos]
        if tax is not None: axes = [tax] + axes
        ret = [vo, do]
        vo = cdms2.createVariable(vo, id=var.id, axes=axes, attributes=var.attributes.copy())
        do = cdms2.createVariable(do, id=depth.id, axes=axes, attributes=var.attributes.copy())
        if timeavg and 't' in var.getOrder():
            self.verbose('Averaging along time')
            for i,d in enumerate(ret):
                self.verbose('Averaging %s', d.id)
                ret[i] = MV2.average(d, axis=0)
            ret[0] = cdms2.createVariable(ret[0], id=var.id, axes=[lev,pos], attributes=var.attributes.copy())
            ret[1] = cdms2.createVariable(ret[1], id=depth.id, axes=[lev,pos], attributes=depth.attributes.copy())
        ret += [yo, xo]
        self.verbose('Final Section data:\n  %s', '\n  '.join(self.describe(o) for o in ret))
        return ret
        
        
    def get_transect(self, varname, lons, lats, 
        method='bilinear', subsamp=3, getcoords=False, timeavg=False, outaxis=None, 
        time=None, lon=None, lat=None, level=None, **kwargs):
        """Get a transect between two points 
        
        It uses :func:`~vacumm.misc.grid.regridding.transect`.
        
        :Params:
        
            - **varname**: Generic var name.
            - **lons/lats**: Specification of transect, either
            
                - Coordinates of first and last point in degrees as 
                tuples in the form ``(lon0,lon1)`` and ``(lat0,lat1)``.
                The array of coordinates is generated using :func:`transect_specs`.
                - Or explicit array of coordinates (as scalars, lists or arrays).
              
            - **subsamp**, optional: Subsampling with respect to grid cell.
            - **method**, optional: Interpolation method 
              (see :func:`~vacumm.misc.grid.regridding.grid2xy`).
            - **getcoords**, optional: Also get computed coordiantes.
            - **outaxis**, optional: Output axis.
        
                - A cdms2 axis.
                - ``None`` or ``'auto'``: Longitudes or latitudes depending 
                  on the range.
                - ``'lon'`` or ``'x'``: Longitudes.
                - ``'lat'`` or ``'y'``: Latitudes.
                - ``'dist'`` or ``'d'``: Distance in km.
            
            - **timeavg**, optional: Time average of results.
        
        :Return:
    
            ``tvar`` or ``tvar,tons,tlats``
        
        """
        # Params
        kwvar = kwfilter(kwargs, ['torect'], time=time, lat=lat, lon=lon, 
            level=level, depthup=False)
        
        # Get data
        var = self.get(varname, **kwvar)
        if var is None: return
        if len(var.shape)<2: return var
        
        # Make transect
        res = transect(var, lons, lats, subsamp=subsamp, 
            getcoords=getcoords, method=method, outaxis=outaxis)
        tvar = res[0] if getcoords else res
        
        # Time average
        if timeavg and tvar.getOrder().startswith('t') and tvar.shape[0]>1:
            tvar = MV2.average(tvar, axis=0)
        
        tvar =  self.finalize_variable(tvar, **kwargs)
        if getcoords: return (tvar, )+res[1:]
        return tvar
    
    def plot_transect(self, varname, lons, lats, 
        method='bilinear', timeavg=False, subsamp=3, outaxis=None, 
        time=None, lon=None, lat=None, level=None, 
        title='%(long_name)s along transect', minimap=None, **kwargs):
        """Plot a transect between two points
        
        :Params:
        
            - **varname**: Generic var name.
            - **lons/lats**: Specification of transect (see :meth:`get_transect`).
            - **title**, optional: Title of the figure.
            - **minimap**, optional: If True, add a minimap showing the transect
              on a map; if False, display nothing; if None, display if no minimap
              already displayed.
            - **minimap_<param>**, optional: Passed to :func:`~vacumm.misc.plot.add_map_lines`.
            - Some params are passed to :meth:`get_transect`.
            - Other params are passed to the plot function :func:`~vacumm.misc.plot.curve2`
              for 1D plots and, :func:`~vacumm.misc.plot.hov2` or 
              :func:`~vacumm.misc.plot.section2` for 2D plots.
            
        """
        
        # Get data
        kwts = dict(method=method, subsamp=subsamp, timeavg=timeavg, 
            outaxis=outaxis, time=time, lon=lon, lat=lat, level=level, 
            squeeze=True, torect=True)
        var, lons, lats = self.get_transect(varname, lons, lats, 
            getcoords=True, **kwts)
        if var is None:
            self.error("Can't get transect on variable")
            return
        if var.getLevel() is not None: 
            depth = self.get_transect('depth', lons, lats, **kwts)
            if isaxis(depth): depth = None
        else:
            depth = None
            
        # Check dimension
        if var.ndim==3:
            self.warning('Forcing time average to have a 2D transect')
            var = MV2.average(var, axis=0)
            if depth is not None and depth.getTime() is not None:
                depth = MV2.average(depth, axis=0)
            
        # Main plot
        kwmap = kwfilter(kwargs, 'minimap_')
        kwargs.update(post_plot=False, title=title)
        kwargs.setdefault('bgcolor', '0.5')
        kwargs.setdefault('ymax', 0)
        # - 1D
        if var.ndim==1:
            
            p = curve2(var, **kwargs)
                        
        # - T-
        elif 't' in var.getOrder():
            
            p = hov2(var, **kwargs)
           
        # - Z-
        else:
            kwargs.setdefault('fill', 'contourf')
            p = section2(var, yaxis=depth, **kwargs)
            
        # Minimap
        if minimap=='auto': minimap = None
        if minimap is True or (minimap is None and not getattr(p.axes, '_minimap', False)):
            kwmap.setdefault('map_zoom', 0.5)
            m = add_map_lines((lons, lats), lons, lats, **kwmap)
            p.axes._minimap = True
        
        del kwargs['title']
        p.post_plot(**kwargs)
        return p
            
    
    def get_hovmoller(self, varname, xorymin, xorymax, xory, meridional=False, select=None):
        '''Get a hovmoller(time,position) section data along a straight trajectory, zonal or meridional.
        
        Coordinates xorymin and xorymax are longitudes and xory a latitude if the section is zonal,
        latitudes and a longitude if the section is meridonal
        
        Level must be defined using the select parameter.
        
        :Params:
            - **varname**: variable to process
            - **xorymin**: westernmost longitude or southernmost latitude coordinate
            - **xorymax**: eastermost longitude or northernmost latitude coordinate
            - **xory**:  longitude or latitude coordinate
            - **meridional**:  if true, hovmoller is meridional, at a longitude and along given latitude range (default is zonal)
            - **select**: selector (should at least restrict to one level)
          
        :Return: A list containing in order:
            - var(time,position): hovmoller variable
            - latitude(position): latitude corresponding to var's position
            - longitude(position): longitude corresponding to var's position
        
        :Example:
            >>>  get_hovmoller(self, 'temp', -10, -6, 47, select=dict(level=slice(-1,None))):
        
        '''
        self.verbose('Getting hovmoller data'
            '\n  variable:    %s'
            '\n  xorymin:     %s'
            '\n  xorymax:     %s'
            '\n  xory:        %s'
            '\n  meridional:  %s'
            '\n  select:      %s', varname, xorymin, xorymax, xory, meridional, select)
        grid = self.get_grid(varname)
        if grid is None:
            raise Exception('No grid found')
        xres, yres = resol(grid)
        if select is None: select = dict()
        if meridional:
            subgrid = grid.subGridRegion((xorymin, xorymax, 'oo'), (xory, xory, 'ccb'))
            xoryax = subgrid.getLatitude()[:]
            select.update(longitude=(xory-xres, xory+xres, 'ccb'), latitude=(xorymin-yres, xorymax+yres, 'ccb'))
            yo = numpy.concatenate(([xorymin], xoryax, [xorymax]))
            xo = numpy.ones(len(yo)) * xory
        else:
            subgrid = grid.subGridRegion((xory, xory, 'ccb'), (xorymin, xorymax, 'oo'))
            xoryax = subgrid.getLongitude()[:]
            select.update(longitude=(xorymin-xres, xorymax+xres, 'ccb'), latitude=(xory-yres, xory+yres, 'ccb'))
            xo = numpy.concatenate(([xorymin], xoryax, [xorymax]))
            yo = numpy.ones(len(xo)) * xory
        
        select = self.get_selector(select)
        var = self.get_variable(varname, select=select, squeeze=True)
        if var is None:
            raise Exception('No %s found'%(varname))
        
        # interp var
        vo = grid2xy(var, xo, yo, method='nat')
        
        axes = vo.getAxisList()
        if meridional:
            axes[-1] = create_lat(yo)
        else:
            axes[-1] = create_lon(xo)
        vo.setAxisList(axes)
        
        self.verbose('Resulting variable: %s', self.describe(vo))
        self.verbose('Resulting longitude: %s', self.describe(xo))
        self.verbose('Resulting latitude: %s', self.describe(yo))
        return vo, yo, xo
    
    def get_extrema_location(self, varname, xorymin, xorymax, xory, meridional=False, extrema='min', select=None):
        '''Get positions of min/max values of variable extremas along a straight trajectory, zonal or meridional.
        
        Coordinates xorymin and xorymax are longitudes and xory a latitude if the section is zonal,
        latitudes and a longitude if the section is meridonal
        
        Level must be defined using the select parameter.
        
        :Params:

            - **varname**: variable to process
            - **xorymin**: westernmost longitude or southernmost latitude coordinate
            - **xorymax**: eastermost longitude or northernmost latitude coordinate
            - **xory**:  longitude or latitude coordinate
            - **meridional**:  if true, hovmoller is meridional, at a longitude and along given latitude range (default is zonal)
            - **extrema**: type of extrema, one of:

                - **min**: retrieve minimum values positions
                - **max**: retrieve maximum values positions

            - **select**: selector (should at least restrict to one level)
          
        :Return: A list containing in order:

            - var(time,position): loaded variable  data
            - latitude(position): latitude corresponding to var's position
            - longitude(position): longitude corresponding to var's position
        
        :Example:
            >>>  get_extrema_location(self, 'temp', -10, -6, 47, select=dict(level=slice(-1,None))):
        
        '''
        self.verbose('Getting extrema location data'
            '\n  variable:    %s'
            '\n  xorymin:     %s'
            '\n  xorymax:     %s'
            '\n  xory:        %s'
            '\n  meridional:  %s'
            '\n  extrema:     %s'
            '\n  select:      %s',varname, xorymin, xorymax, xory, meridional, extrema, select)
        var, lat, lon = self.get_hovmoller(varname, xorymin, xorymax, xory, meridional, select)
        
        ex = extrema.strip().lower()
        exfunc = getattr(numpy.ma, 'arg%s'%(ex), None)
        if ex not in ['min','max'] or exfunc is None:
            raise ValueError('Invalid extrema: %s'%(extrema))
        
        # Position of extema along time with masked values
        #iex = exfunc(var, axis=1, fill_value=-1)
        iex = exfunc(var, axis=1)
        #bad = iex == -1
        
        # Initialize localization variable
        vo = var[:,0].clone()
        if meridional:
            vo.long_name = 'Latitude of %s'%(ex)
            vo.units = 'degrees_north'
            lonorlat = var.getLatitude()
        else:
            vo.long_name = 'Longitude of %s'%(ex)
            vo.units = 'degrees_east'
            lonorlat = var.getLongitude()
        
        # Fill with lon/lat coordinates
        for j, i in numpy.ndenumerate(iex):
            vo[j] = lonorlat[i]
        
        # Remasking
        #vo[:] = MV2.masked_where(bad, vo, copy=0)
        
        self.verbose('Resulting variable: %s', self.describe(vo))
        self.verbose('Resulting longitude: %s', self.describe(lon))
        self.verbose('Resulting latitude: %s', self.describe(lat))
        return vo, lat, lon
    
    
    def get_localized_computed_values(self, varname, xorymin, xorymax, xory, meridional=False, operation='min', select=None):
        '''Get min/max/mean values of a variable along a straight trajectory, zonal or meridional.
        
        Coordinates xorymin and xorymax are longitudes and xory a latitude if the section is zonal,
        latitudes and a longitude if the section is meridonal
        
        Level must be defined using the select parameter.
        
        :Params:
            - **varname**: variable to process
            - **xorymin**: westernmost longitude or southernmost latitude coordinate
            - **xorymax**: eastermost longitude or northernmost latitude coordinate
            - **xory**:  longitude or latitude coordinate
            - **meridional**:  if true, hovmoller is meridional, at a longitude and along given latitude range (default is zonal)
            - **operation**: type of operation, one of:

                - **min**: retrieve minimum values
                - **mean**: retrieve mean values
                - **max**: retrieve maximum values

            - **select**: selector (should at least restrict to one level)
          
        :Return: A list containing in order:
            - var(time,position): loaded variable  data
            - latitude(position): latitude corresponding to var's position
            - longitude(position): longitude corresponding to var's position
        
        :Example:
            >>>  get_localized_computed_values(self, 'temp', -10, -6, 47, select=dict(level=slice(-1,None))):
        
        '''
        self.verbose('Getting localized computed data'
            '\n  variable:    %s'
            '\n  xorymin:     %s'
            '\n  xorymax:     %s'
            '\n  xory:        %s'
            '\n  meridional:  %s'
            '\n  operation:   %s'
            '\n  select:      %s',varname, xorymin, xorymax, xory, meridional, operation, select)
        var, lat, lon = self.get_hovmoller(varname, xorymin, xorymax, xory, meridional, select)
        op = operation.strip().lower()
        if op in ('mean','avg'):
            op = 'average'
        elif op not in ['min','max']:
            raise ValueError('Invalid extrema: %s'%(extrema))
        op = getattr(MV2, op, None)
        if op is None:
            raise ValueError('Invalid operation: %s'%(operation))
        vo = op(var, -1)
        vo.id = var.id
        vo.attributes.update(var.attributes)
        self.verbose('Resulting variable: %s', self.describe(vo))
        return vo, lat, lon
    
    
    def get_stratification_data(self, select, timeavg=True):
        '''Get stratification data
        
        :Params:
            - **select**: selector
            - **timeavg**: if true, average date along time if needed
          
        :Return:
            - temp, sal, dens, depth, deltadepth with shape ([time],depth,latitude,longitude)
            - densmin, densmax, densmean with shape ([time],latitude,longitude)
        
        '''
        select = self.get_selector(select=select)
        temp = self.get_temp(select=select, squeeze=True)
        sal = self.get_sal(select=select, squeeze=True)
        lat = temp.getLatitude()
        depth = self.get_depth(temp.id, select=select, squeeze=True)
        # TODO: depth may not have a time axis !
        depth.setAxisList(temp.getAxisList())
        self.verbose('Initial stratification data for select=%s:\n  %s', select, '\n  '.join(self.describe(o) for o in (temp, sal, lat, depth)))
        if temp.getOrder() not in ('tzyx','zyx'):
            raise ValueError('Invalid order: %s, [t]zyx expected'%(temp.getOrder()))
        def compute_strat(temp, sal, lat, depth):
            '''Compute strat data for one time step'''
            # Epaisseurs
            shp = depth.shape
            ddepth = meshweights(depth.reshape((depth.shape[0], numpy.product(depth.shape[1:]))), axis=0)
            ddepth = ddepth.reshape(shp)
            # Pression
            pres = seawater.csiro.pres(depth, numpy.resize(lat[:], depth.shape))
            # Densité
            dens = seawater.csiro.dens(sal, temp, pres)
            # Densité min/max
            densmin = dens.min(axis=0)
            densmax = dens.max(axis=0)
            # Densité moyenne
            densmean = MV2.average(dens, axis=0, weights=ddepth)
            ret = dens, densmin, densmax, densmean, ddepth
            self.verbose('\n  %s', '\n  '.join(self.describe(o) for o in ret))
            return ret
        if 't' in temp.getOrder():
            self.verbose('Computing stratification along time')
            tmp, data = [], []
            time = depth.getTime()
            ctime = time.asComponentTime()
            for it in xrange(len(time)):
                self.verbose('Processing time: %s', ctime[it])
                tmp.append(compute_strat(temp[it], sal[it], lat, depth[it]))
            for i in xrange(len(tmp[0])):
                d = [d[i] for d in tmp]
                data.append(MV2.concatenate([d]))
            dens, densmin, densmax, densmean, ddepth = data
        else:
            dens, densmin, densmax, densmean, ddepth = compute_strat(temp, sal, lat, depth)
        axes = temp.getAxisList()
        ddepth = cdms2.createVariable(ddepth, id='weight', axes=axes)
        dens = cdms2.createVariable(dens, id='dens', axes=axes)
        axes.pop(depth.getOrder().index('z'))
        densmin = cdms2.createVariable(densmin, id='densmin', axes=axes)
        densmax = cdms2.createVariable(densmax, id='densmax', axes=axes)
        densmean = cdms2.createVariable(densmean, id='densmean', axes=axes)
        ret = [temp, sal, dens, densmin, densmax, densmean, depth, ddepth]
        if timeavg and 't' in temp.getOrder():
            self.verbose('Averaging along time')
            for i,d in enumerate(ret):
                self.verbose('Averaging %s', d.id)
                ret[i] = MV2.average(d, axis=0)
        self.verbose('Final stratification data:\n  %s', '\n  '.join(self.describe(o) for o in ret))
        return ret
    
    _mode_doc = """Computing mode
    
          - ``None``: Try all modes, in the following order.
          - ``"var"``: Read it from a variable.
          - ``"deltatemp"``: Estimate from a difference in temperature.
          - ``"deltadens"``: Estimate from a difference in density.
          - ``"twolayers"``: Shalow water mode with two density layers.
          - ``"kz"``: Depth where ks becomes low.
          
          You can specifiy a list of them: ``['deltadens', 'deltatemp']`` 
          You can also negate the search with
          a '-' sigme before: ``"-kz"``."""
    def get_mld(self, mode=None, lat=None, **kwargs):
        
        # Params
        kwargs.pop('level', None)
        kwvar = kwfilter(kwargs, ['lon','time','torect'])
        kwvar['lat'] = lat
        kwdens = kwfilter(kwargs, 'dens_')
        kwdens.update(kwvar)
        kwdepth = kwfilter(kwargs, 'depth_')
        kwdepth.update(kwvar)
        kwfinal = kwargs.copy()
        kwfinal.pop('depthup',  None)
        kwmld = kwfilter(kwargs, ['mld_',  'deltatemp',  'deltadens',  'kzmax'])
        
        # First, try to find a MLD variable
        if check_mode('var', mode):
            mld = self.get_variable('mld', **kwvar)
            if mld is not None or check_mode('var', mode, strict=True): 
                return self.finalize_variable(mld, **kwfinal)
            
        # Estimate it
        depth = self.get_depth(**kwdepth)
        dens = None
        
        # - deltadens
        if check_mode('deltatemp', mode):           
            temp = self.get_temp(**kwvar)
            if temp is not None:
                mld = mixed_layer_depth(temp, depth=depth, mode='deltatemp', format_axes=True, **kwmld)
                if mld is not None or check_mode('deltatemp', mode, strict=True): 
                    return self.finalize_variable(mld, depthup=False, **kwfinal)
                    
        # - deltadens and twolayers (density)
        for testmode in 'deltadens', 'twolayers':
            if check_mode(testmode, mode):               
                if dens is None:
                    dens = self.get_dens(**kwdens)
                if dens is not None:
                    mld = mixed_layer_depth(dens, depth=depth, lat=lat, mode=testmode, format_axes=True, **kwmld)
                    if mld is not None or check_mode(testmode, mode, strict=True): 
                        return self.finalize_variable(mld, depthup=False, **kwfinal)
                        
        # - Kz
        if check_mode('kz', mode):
             
            kz = self.get_kz(**kwdens)
            if kz is not None:
                mld = mixed_layer_depth(kz, depth=depth, mode='kz', format_axes=True, **kwmld)
                if mld is not None or check_mode('twolayers', mode, strict=True): 
                    kwfinal['depthup'] = False
                    return self.finalize_variable(mld, **kwfinal)
             
    get_mld = formatdoc_var(get_mld, 
        mode=_mode_doc, 
        deltatemp='Temprature difference with surface.', 
        deltadens='Density difference with surface', 
        kzmax='Kz max for search for low values', 
        )
            
    
    def get_mixed_layer_depth(self, select):
        '''Get mixed layer depth
        
        MLD is computed for each time step and then averaged
        
        :Params:
            - **select**: selector with at least a time component
          
        :Return:
            - mld with shape (latitude,longitude)
        
        '''
        select = select.copy()
        time = select.pop('time')
        mlds = []
        td = self.get_time_res()
        ts = (td.days*86400+td.seconds, 'seconds')
        self.info('Detected model time step: %s', ts)
        for itv in Intervals(time, ts):
            s = select.copy()
            s['time'] = list(itv[:2])+['co']
            temp, sal, dens, densmin, densmax, densmean, depth, ddepth = self.get_stratification_data(s)
            # Profondeur max (max+demi épaisseur)
            H = depth[-1]+ddepth[-1]*.5
            mld = H*(densmax-densmean)/(densmax-densmin)
            self.info('MLD[time=%s]: %s', itv, self.describe(mld))
            mlds.append([mld])
        mlds = MV2.concatenate(mlds, axis=0)
        mld = MV2.average(mlds, axis=0)
        mld.id = 'MLD'
        mld.units = 'm'
        mld.long_name = u'Mixed layer depth'
        mld.setAxisList(temp.getAxisList()[-2:])
        self.info('AVG MLD: %s', self.describe(mld))
        return mld
    
    
    def get_potential_energy_deficit(self, select):
        '''Get potential energy deficit
        
        PED is computed for each time step and then averaged
        
        :Params:
            - **select**: selector with at least a time component
          
        :Return:
            - ped with shape (latitude,longitude)
        
        '''
        select = select.copy()
        time = select.pop('time')
        peds = []
        td = self.get_time_res()
        ts = (td.days*86400+td.seconds, 'seconds')
        self.info('Detected model time step: %s', ts)
        for itv in Intervals(time, ts):
            s = select.copy()
            s['time'] = itv
            temp, sal, dens, densmin, densmax, densmean, depth, ddepth = self.get_stratification_data(s)
            # Anomalie de densité
            danom = dens-densmean
            # Énergie potentielle disponible
            ape = danom * g
            ape *= ddepth
            # Deficit
            ped = MV2.average(ape, axis=0, weights=ddepth)
            self.info('PED[time=%s]: %s', itv, self.describe(ped))
            peds.append([ped])
        peds = MV2.concatenate(peds, axis=0)
        ped = MV2.average(peds, axis=0)
        ped.id = 'PED'
        ped.units = 'J.m^{-2}'
        ped.long_name = u"Potential energy deficit"
        ped.setAxisList(temp.getAxisList()[-2:])
        self.info('AVG PED: %s', self.describe(ped))
        return ped
    
    
    ############################################################################
    ############################################################################
    ############################################################################
    
    
    def plot_trajectory_map(self, lon, lat, **kwargs):
        ''' Plot the "legend" map of a trajectory
        
        .. todo::
            - replace this method usage by vacumm.misc.plot.traj
        
        '''
        varname=kwargs.pop('varname', None)
        grid = self.get_grid(varname)
        glob_lon, glob_lat = grid.getLongitude(), grid.getLatitude()
        mapkw = kwfilter(kwargs, 'map',
            defaults=dict(
                label=self.__class__.__name__,
                lon=(MV2.min(glob_lon), MV2.max(glob_lon)), lat=(MV2.min(glob_lat), MV2.max(glob_lat)),
                axes_rect=[0.75, 0.1, 0.2, 0.3],
                drawmeridians_size=6, drawparallels_size=6,
                show = False))
        for k in post_plot_args: mapkw.pop(k, None) # FIXME: make post_plot available for standalone trajectory plot ?
        mapkw.update(dict(vertical=True, show=False))
        # The map itself
        mp = map2(**mapkw)
        # The section line
        xx, yy = mp.map(lon, lat)
        mp.map.plot(xx, yy, 'r-')
        return mp
    
    
    def plot_layer(self, *args, **kwargs):
        '''
        Plot a layer
        
        Params:
            - **map_<keyword>**: passed to :func:`~vacumm.misc.plot.map2`
            - **plot_<keyword>**: passed to created map :func:`post_plot`
          
        Other params are passed to :func:`get_layer`
        
        '''
        mp = kwargs.pop('map', None)
        mapkw = kwfilter(kwargs, 'map', default=dict(label=self.__class__.__name__))
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        var = self.get_layer(*args, **kwargs)
        mapkw.update(show=False)
        if mp: mp(var, **mapkw)
        else: mp = map2(var, **mapkw)
        # Post plotting
        plotkw.setdefault('title', var.id)
        mp.post_plot(**plotkw)
        return mp
    
    
    def plot_section(self, varname, xmin, ymin, xmax, ymax, select=None, pmap=True, **kwargs):
        '''
        Produce a section plot.
        
        :Params: see :func:`get_section`
        
        :Plot params:
            - **sec_<keyword>**: are passed to the section plot function :func:`~vacumm.misc.plot.section2` *excepting* those about post plotting described below
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations
        
        .. todo::
            - add lat/lon position lines indicator
        
        '''
        datakw  = kwfilter(kwargs, 'data')
        # Get section data
        var,dep,lat,lon = self.get_section(varname, xmin, ymin, xmax, ymax, select, **datakw)
        pos = var.getAxis(1)
        # Compute scaling/informationnal data
#        vmin,vmax = numpy.min(var), numpy.max(var)
#        vstep = abs(vmax-vmin)/10.
        dmin,dmax = MV2.min(dep), MV2.max(dep)
#        dstep = abs(dmax-dmin)/10.
#        if dmax > 0: dmax = 0
#        var, dep = var[::-1], dep[::-1]
#        print dep[::4,0]
#        print var[::4,0]
        seckw = kwfilter(kwargs, 'sec', defaults=
            dict(
                label=self.__class__.__name__,
                figsize=(12,6),
                fill='contourf',
                xlong_name='Position', xlabel='position', xunits='pos', #xmin=, xmax=,
                #xlim=(pos[0],pos[-1]), xfmt='%s', xtickfmt='%s', 
                ylong_name='Depth', ylabel='depth', yunits=dep.units, #ymin=dmin, ymax=dmax,
                #ylim=(dmin,dmax), ytickfmt='%s', yfmt='%s',
                #vmin=vmin, vmax=vmax, units=var.units,
                #colorbar_ticks=numpy.arange(vmin, vmax+vstep, vstep),
                axes_rect=[0.1, 0.1, 0.6, 0.8]))
        mapkw = kwfilter(kwargs, 'map', keep=True, defaults=dict(varname=var.id))
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        # Plot the section
        seckw.update(yaxis=dep, levels_mode='normal', show=False)
        sc = section2(var, **seckw)
#        from vacumm.misc.plot import section
#        section(var, **seckw)
        from vacumm.misc.plot import add_grid
        add_grid((pos, dep))
        mp = None
        if pmap:
            mp = self.plot_trajectory_map(lon, lat, **mapkw)
        # Post plotting
        plotkw.setdefault('title', 'Vertical section of %s\ntime: %s\nlon: %s to %s, lat: %s to %s'%(var.id, select['time'], xmin, xmax, ymin, ymax))
        sc.post_plot(**plotkw)
        return sc
    
    
    def plot_hovmoller(self, varname, xorymin, xorymax, xory, meridional=False, select=None, pmap=True, **kwargs):
        '''
        Produce a hovmoller plot.
        
        :Params: see :func:`get_hovmoller`
        
        :Plot params:
            - **hov_<keyword>**: are passed to the section plot function :func:`~vacumm.misc.plot.hov2` *excepting* those about post plotting described below
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations
        
        '''
        datakw  = kwfilter(kwargs, 'data')
        # Get section data
        var, lat, lon = self.get_hovmoller(varname, xorymin, xorymax, xory, meridional, select, **datakw)
        # Compute scaling/informationnal data
        vmin,vmax = numpy.min(var), numpy.max(var)
        vstep = abs(vmax-vmin)/10.
        xname = meridional and 'Latitude' or 'Longitude'
        hovkw = kwfilter(kwargs, 'hov',
            defaults=dict(
                label=self.__class__.__name__,
                figsize=(12,6),
                fill='contourf',
                xlong_name=xname, xlabel=xname, xunits='', #xmin=, xmax=,
                ylong_name='Time', ylabel='time', yunits='', #ymin=dmin, ymax=dmax,
                #vmin=vmin, vmax=vmax, units=var.units,
                #colorbar_ticks=numpy.arange(vmin, vmax+vstep, vstep),
                axes_rect=[0.1, 0.1, 0.6, 0.8]))
        for k in post_plot_args: hovkw.pop(k, None)
        mapkw = kwfilter(kwargs, 'map', keep=True, defaults=dict(varname=var.id))
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        # Plot the section
        hovkw.update(time_vertical=True, show=False)
        hv = hov2(var, **hovkw)
        mp = None
        if pmap:
            mp = self.plot_trajectory_map(lon, lat, **mapkw)
        # Post plotting
        if meridional:
            hovtype = 'Meridional'
            lonlat = 'lon: %s, lat: %s to %s'%(xory, xorymin, xorymax)
        else:
            hovtype = 'Zonal'
            lonlat = 'lon: %s to %s, lat: %s'%(xorymin, xorymax, xory)
        plotkw.setdefault('title', '%s Hovmoller of %s\ntime: %s\n%s'%(hovtype, var.id, select['time'], lonlat))
        hv.post_plot(**plotkw)
    
    
    def plot_extrema_location(self, varname, xorymin, xorymax, xory, meridional=False, extrema='min', pmap=True, select=None, **kwargs):
        '''
        Produce 1D plot of min/max positions.
        
        :Params: see :func:`get_extrema_location`
        
        :Plot params:
            - **cur_<keyword>**: are passed to the section plot function :func:`~vacumm.misc.plot.curve2` *excepting* those about post plotting described below
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations
        
        '''
        datakw  = kwfilter(kwargs, 'data')
        # Get data
        var, lat, lon = self.get_extrema_location(varname, xorymin, xorymax, xory, meridional, extrema, select, **datakw)
        # Compute scaling/informationnal data
        vmin,vmax = numpy.min(var), numpy.max(var)
        vstep = abs(vmax-vmin)/10.
        xname = meridional and 'Latitude' or 'Longitude'
        curkw = kwfilter(kwargs, 'cur',
            defaults=dict(
                label=self.__class__.__name__,
                figsize=(12,6),
                xlong_name=xname, xlabel=xname, xunits=var.units, xmin=xorymin, xmax=xorymax, #xmin=vmin+vstep, xmax=vmax+vstep,
                ylong_name='Time', ylabel='time', #yunits=, #ymin=dmin, ymax=dmax,
                #vmin=vmin, vmax=vmax, units=var.units,
                axes_rect=[0.1, 0.1, 0.6, 0.8]))
        for k in post_plot_args: curkw.pop(k, None)
        mapkw = kwfilter(kwargs, 'map', keep=True, defaults=dict(varname=var.id))
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        # Plot the location curve
        curkw.update(order='td', show=False)
        cur = curve2(var, **curkw)
        if pmap:
            mp = self.plot_trajectory_map(lon, lat, **mapkw)
        # Post plotting
        if meridional:
            curtype = 'Meridional'
            lonlat = 'lon: %s, lat: %s to %s'%(xory, xorymin, xorymax)
        else:
            curtype = 'Zonal'
            lonlat = 'lon: %s to %s, lat: %s'%(xorymin, xorymax, xory)
        plotkw.setdefault('title', '%s location of %s %s\ntime: %s\n%s'%(curtype, extrema, var.id, select['time'], lonlat))
        cur.post_plot(**plotkw)
        return cur, var, lat, lon
    
    
    def plot_localized_computed_values(self, varname, xorymin, xorymax, xory, meridional=False, operation='min', pmap=True, select=None, **kwargs):
        '''
        Produce 1D plot of min/max values.
        
        :Params: see :func:`get_localized_computed_values`
        
        :Plot params:
            - **cur_<keyword>**: are passed to the section plot function :func:`~vacumm.misc.plot.curve2` *excepting* those about post plotting described below
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations
        
        '''
        datakw  = kwfilter(kwargs, 'data')
        # Get data
        var, lat, lon = self.get_localized_computed_values(varname, xorymin, xorymax, xory, meridional, operation, select, **datakw)
        # Compute scaling/informationnal data
        vmin,vmax = numpy.min(var), numpy.max(var)
        vstep = abs(vmax-vmin)/10.
        curkw = kwfilter(kwargs, 'cur',
            defaults=dict(
                label=self.__class__.__name__,
                figsize=(12,6),
                xlong_name='Time', xlabel='time', #xunits=, #xmin=, xmax=,
                ylong_name=var.id, ylabel=var.id, #yunits=var.units, ymin=vmin, ymax=vmax,
                #vmin=vmin, vmax=vmax, units=var.units,
                axes_rect=[0.1, 0.1, 0.6, 0.8]))
        for k in post_plot_args: curkw.pop(k, None)
        mapkw = kwfilter(kwargs, 'map', keep=True, defaults=dict(varname=var.id))
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        # Plot the computed curve
        curkw.update(show=False)#, order='td')
        cur = curve2(var, **curkw)
        if pmap:
            mp = self.plot_trajectory_map(lon, lat, **mapkw)
        # Post plotting
        if meridional:
            curtype = 'Meridional'
            lonlat = 'lon: %s, lat: %s to %s'%(xory, xorymin, xorymax)
        else:
            curtype = 'Zonal'
            lonlat = 'lon: %s to %s, lat: %s'%(xorymin, xorymax, xory)
        plotkw.setdefault('title', '%s %s of %s\ntime: %s\n%s'%(curtype, operation, var.id, select['time'], lonlat))
        cur.post_plot(**plotkw)
        return cur, var, lat, lon
    
    
    def plot_mld(self, select, **kwargs):
        '''
        Produce mixed layer depth map.
        
        :Params: see :func:`get_mixed_layer_depth`
        
        :Plot params:
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations
        
        '''
        mapkw = kwfilter(kwargs, 'map')
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        mld = self.get_mixed_layer_depth(select)
        l = auto_scale((mld.min(), mld.max()))
        c = C.cmap_magic(l)
        #c = C.cmap_rainbow()
        mapkw.update(dict(show=False, cmap=c))
        m = map2(mld, **mapkw)
        m.post_plot(**plotkw)
        return m,c
    
    
    def plot_ped(self, select, **kwargs):
        '''
        Produce potential energy deficit map.
        
        :Params: see :func:`get_potential_energy_deficit`
        
        :Plot params:
            - **map_<keyword>**: are passed to the map plot function :func:`~vacumm.misc.plot.map2` *excepting* those about post plotting described below
            - **plot_[show|close|savefig|savefigs]**: are passed to the post plotting function :func:`~vacumm.misc.core_plot.Plot.post_plot` at end of plotting operations
        
        '''
        mapkw = kwfilter(kwargs, 'map')
        for k in post_plot_args: mapkw.pop(k, None)
        plotkw = kwfilter(kwargs, 'plot')
        ped = self.get_potential_energy_deficit(select)
        l = auto_scale((ped.min(), ped.max()))
        c = C.cmap_magic(l)
        #c = C.cmap_rainbow()
        mapkw.update(dict(show=False, cmap=c))
        m = map2(ped, **mapkw)
        m.post_plot(**plotkw)
        return m,c
    
class AtmosDataset(AtmosSurfaceDataset):    
    ncobj_specs = AtmosSurfaceDataset.ncobj_specs.copy()
    ncobj_specs.update(dict(
    ))


def _at_(at, squeezet=False, focus=None, prefix=False):
    """Convert location letters"""
    
    if isinstance(at, basestring): at = at.lower()
    if at=='' or at is True: at = 't'
    
    # Aliases
    if at in 'rts': 
        at = 't'
        
    # Focus on hor. or ver.
    if (focus=="hor" or focus=='h' or focus=='xy') and at in 'w':
        at = 't'
    elif focus=="ver" or focus=='z':
        if at in 'uv':
            at = 't'
        elif at in 'd':
            at = "w"
    
    # Squeeze t
    if squeezet and at=='t': at = ''
            
    # Prefix
    if at != '':
        if prefix is True: 
            prefix = '_'
        if isinstance(prefix, basestring):
            at = prefix+at
    
    return at

def check_mode(validmodes, modes, strict=False):
    """Check at least one of the checked modes is valid
    
    :Params:
    
        - **validmodes**: Valid modes as a single or list of strings.
        - **modes**: Modes to check. ``None`` is accepted.
        - **strict**, optional: If a checked mode is ``None``,
          it returns ``True`` if strict is ``True`` else ``False``.
    
    :Example:
    
        >>> check_mode('var', 'var')
        True
        >>> validmodes = ['var', 'depth']
        >>> check_mode(validmodes, 'toto')
        False
        >>> check_mode(validmodes, None)
        True
        >>> check_mode(validmodes, 'var')
        True
        >>> check_mode(validmodes, 'var', 'toto')
        True
        >>> check_mode(validmodes, ['var', 'toto'])
        True
        >>> check_mode(validmodes, '-axis')
        True
    """    
    if not isinstance(validmodes, (list, tuple)):
        validmodes = [validmodes]
    if not isinstance(modes, (list, tuple)):
        modes = [modes]
    for mode in modes:
        if not isinstance(mode, (list, tuple)):
            mode = [mode]
        for m in modes:
            if m is None and not strict: return True
            rev = isinstance(m, basestring) and m.startswith('-')
            if strict and rev:
                return False
            if rev and m[1:] not in validmodes:
                return True
            elif not rev and m in validmodes:
                return True
    return False

def _notuplist_(dd):
    """Convert lists and tuples to single elements in dict (after copy)"""
    if dd is None: return
    dd = dd.copy()
    for key, val in dd.items():
        if isinstance(val, (list, tuple)):
            dd[key] = val[0]
    return dd

class GenericDataset(AtmosDataset, OceanDataset):
    """Generic :class:`Dataset` class to load everything"""
    pass
    
class CurvedSelector(object):
    """Curved grid multiple selector"""

    def __init__(self, grid, select):
    
        self.geosels = self.extract_geosels(select)
        self._post_sel = False
        if self.geosels and grid is not None and len(grid.getLongitude().shape)==2:
            islice, jslice, mask =  coord2slice(grid, *self.geosels[0])
            if islice is None: islice = ':'
            if jslice is None: jslice = ':'
            self.remove_geosels(select)
            self.xid = grid.getAxis(1).id
            self.yid = grid.getAxis(0).id
            select.refine(**{self.xid:islice, self.yid:jslice})
            self._post_sel = True
            self.mask = mask

        
    def finalize(self, var):
        if var is None or not self._post_sel: return var
        if self.mask.any():
            var[:] = MV2.masked_where(N.resize(self.mask, var.shape), var, copy=0)
        if len(self.geosels)==2:
            islice, jslice, mask =  coord2slice(var, *self.geosels[1])
            var = var(**{self.xid:islice, self.yid:jslice})
            if self.mask.any():
                var[:] = MV2.masked_where(N.resize(mask, var.shape), var, copy=0)
            return var
        return var

    @staticmethod
    def extract_geosels(select):
        """Convert select to a list of (lon,lat) selection specs"""
        lons = []
        lats = []
        for c in select.components():
            id = getattr(c, 'id', None)
            if id is None: continue
            if id=='lon':
                lons.append(c.spec)
            if id=='lat':
                lats.append(c.spec)
            continue
        n = max(len(lons), len(lats))
        if n==0: return []
        lons = broadcast(lons, n, fillvalue=None)
        lats = broadcast(lats, n, fillvalue=None)
        return zip(lons, lats)
        
    @staticmethod
    def remove_geosels(select):
        """Remove lon and lat component selections"""
        for c in select().components():
            id = getattr(c, 'id', None)
            if id is None: continue
            if id in ['lon', 'lat']:
                select._Selector__components.remove(c)
        return select
            
