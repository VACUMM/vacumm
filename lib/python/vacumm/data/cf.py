# -*- coding: utf8 -*-
"""Conventions about data formats and names"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
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

import os
from warnings import warn
from collections import OrderedDict
import string

import cdms2, MV2, re
from vacumm import VACUMMError, vcwarn
from vacumm.misc import kwfilter, dict_merge
from vacumm.misc.axes import create as create_axis, isaxis
from vacumm.misc.color import get_cmap
from vacumm.misc.grid import create_axes2d
from vacumm.misc.io import ncmatch_obj
from vacumm.misc.config import ConfigManager
from vacumm.data.misc.arakawa import ARAKAWA_LOCATIONS

__all__ = ['VAR_SPECS', 'AXIS_SPECS',
    'format_var', 'format_axis', 'format_grid', 'match_obj',
    'cf2atts', 'cf2search', 'cp_suffix', 'get_loc',
    'change_loc', 'change_loc_single', 'dupl_loc_specs', 'no_loc_single',
    'change_loc_specs', 'squeeze_loc_single', 'squeeze_loc', 'get_physloc',
    'HIDDEN_CF_ATTS', 'set_loc', 'match_known_var', 'match_known_axis',
    'CF_AXIS_SPECS', 'CF_VAR_SPECS',
    'register_cf_variable', 'register_cf_variables_from_cfg',
    'register_cf_axis', 'register_cf_axes_from_cfg',
    'CF_DICT_MERGE_KWARGS', 'get_cf_cmap', 'is_cf_known',
    'VarSpecs', 'AxisSpecs', 'CF_VAR_NAMES', 'CF_AXIS_NAMES',
]

ARAKAWA_SUFFIXES = [('_'+p) for p in ARAKAWA_LOCATIONS]
arakawa_suffixes = ARAKAWA_SUFFIXES # compat
arakawa_locations = ARAKAWA_LOCATIONS

THISDIR = os.path.dirname(__file__)

#: Joint variables and axes config specification file
CF_GCFG_INIFILE = os.path.join(THISDIR, 'cf.ini')

# Config manager for specs
CF_CFGM = ConfigManager(CF_GCFG_INIFILE)

#: Base config file for CF specifications
CF_CFG_CFGFILE = os.path.join(THISDIR, 'cf.cfg')

#: Base config object for CF specifications
CF_CFG = CF_CFGM.load(CF_CFG_CFGFILE)

CF_DICT_MERGE_KWARGS = dict(mergesubdicts=True, mergelists=True,
                         skipnones=False, skipempty=False,
                         overwriteempty=True, mergetuples=True)

class BaseSpecs(object):
    """Base class for loading variables and axes CF specifications"""
    aliases = {'id': ['ids', 'name', 'names'],
               'standard_name': ["standard_names"],
               'long_name': ['long_names']
              }
    category = None

    def __init__(self, inherit=None, names=None):
        self._dict = OrderedDict(CF_CFG[self.category])
        for name, spec in self._dict.items():
            self._dict[name] = self._dict[name].dict()
        self._names = names
        self._inherit = inherit
        self._cfgspec = CF_CFGM._configspec[self.category]['__many__'].dict()
        self._post_process_()

    @classmethod
    def get_alias_list(cls):
        al = []
        for aa in cls.aliases.values():
            al.extend(aa)
        return al

    def __getitem__(self, key):
        return self._dict[key]
    def __iter__(self):
        return iter(self._dict)
    def __len__(self):
        return len(self._dict)
    def __contains__(self, key):
        return key in self._dict
    def __str__(self):
        return str(self._dict)
    @property
    def names(self):
        return self._dict.keys()
    def items(self):
        return self._dict.items()
    def keys(self):
        return self._dict.keys()

    def register(self, name, **specs):
        """Register a new elements from its name and specs"""
        data = {self.category: {name: specs}}
        self.register_from_cfg(data)

    def register_from_cfg(self, cfg):
        """"Register new elements from a :class:`ConfigObj` instance
        or a config file"""
        if isinstance(cfg, dict) and self.category not in cfg.keys():
            cfg = {self.category:cfg}
        local_cfg = CF_CFGM.load(cfg)[self.category]
        self._dict = dict_merge(self._dict, local_cfg,
            cls=OrderedDict, **CF_DICT_MERGE_KWARGS)
        self._post_process_()

    def _post_process_(self):
        self._from_atlocs = []
        for name in self.names:
            self._check_entry_(name)
        self._add_aliases_()
        self._update_names_()

    def _check_entry_(self, name):
        """Validate an entry

        - Makes sure to have lists, except for 'axis' and 'inherit'
        - Check geo axes
        - Check inheritance
        - Makes sure that axes specs have no 'axes' key
        - Makes sure that specs key is the first entry of 'names'
        - add standard_name to list of names
        - Check duplication to other locations ('toto' -> 'toto_u')

        """
        # Wrong entry!
        if name not in self:
            return

        # Entry already generated with the atlocs key
        if name in self._from_atlocs:
            return

        # Get the specs
        if hasattr(self._dict[name], 'dict'):
            self._dict[name] = self._dict[name].dict()
        specs = self._dict[name]

        # Ids
        if name in specs['id']:
            specs['id'].remove(name)
        specs['id'].insert(0, name)
        if name=='mld':
            pass

        # Long name
        if not specs['long_name']:
            specs['long_name'].append(name.title().replace('_', ' '))

        # Physloc must be the first atlocs
        if 'physloc' in specs and 'atlocs' in specs:
            p = specs['physloc']
            if p:
                if p in specs['atlocs']:
                    specs['atlocs'].remove(p)
                specs['atlocs'].insert(0, p)

        # Geo axes (variables only)
        if 'axes' in specs:
            specs['axes'].setdefault('t', ['time'])
            suffixes = [('_'+s) for s in 'rftuv']
            specs['axes'].setdefault('y', [cp_suffix(name, 'lat', suffixes=suffixes)])
            specs['axes'].setdefault('x', [cp_suffix(name, 'lon', suffixes=suffixes)])
            for l, n in ('t', 'time'), ('y', 'lat'), ('x', 'lon'):
                if isinstance(specs['axes'][l], list) and n not in specs['axes'][l]:
                    specs['axes'][l].append(n)

        # Inherits from other specs (merge specs with dict_merge)
        if name=='mld':
            pass
        if specs['inherit']:
            from_name = specs['inherit']
            from_specs = None
            to_scan = []
            if self._inherit:
                to_scan.append(self._inherit)
            to_scan.append(self)
            for from_specs in to_scan:
                if from_name in from_specs:

                    # Merge
                    self._dict[name] = specs = dict_merge(specs, from_specs[from_name],
                        cls=dict, **CF_DICT_MERGE_KWARGS)

                    # Validate
#                    if from_specs is not self:
                    for key in specs.keys():
                        if key not in self._cfgspec:
                            del specs[key]

#                    break

        # Standard_names in ids
        if name=='mld':
            pass
        if specs['standard_name']:
            for standard_name in specs['standard_name']:
                if standard_name not in specs['id']:
                    specs['id'].append(standard_name)

        # Duplicate at other locations
        if specs['atlocs']:
            tonames =  dupl_loc_specs(self._dict, name, specs['atlocs'])
            self._from_atlocs.extend(tonames)
            specs = self._dict[name]


    def _add_aliases_(self):
        for name, specs in self._dict.items():
            for key, aliases in self.aliases.items():
                if key in specs:
                    for alias in aliases:
                        specs[alias] = specs[key]

    def _update_names_(self):
        if self._names is not None:
            while self._names:
                del self._names[0]
            self._names.extend(self.names)


class VarSpecs(BaseSpecs):
    category = 'variables'

class AxisSpecs(BaseSpecs):
    category = 'axes'


_reloc =  OrderedDict(
    [('id',  re.compile('(_([a-z]))?$', re.I).search),
    ('standard_name', re.compile('(_at_([a-z])_location)?$', re.I).search),
    ('long_name',  re.compile('( at ([a-z]) location)?$', re.I).search)],
)

#: Default location of gridded variables
DEFAULT_LOCATION = 't'
default_location = DEFAULT_LOCATION # compat

def get_loc(var, stype=None, mode='loc', default=None):
    """Get the location testing id, standard_name or long_name and other attributes

    :Params:

        - **var**: Generic id of the variable or axis OR an object CDAT object
          (variable or axis). If a CDAT object, it first searches for the
          :attr:`_vacumm_cf_location` attributes.
        - **stype**, optional: One of 'id', 'standard_name' or 'long_name'.
          If ``None``, they are all tested.
          If ``var`` is a string, location is searched interpreting it as of type ``'stype``.
          Else if ``var`` is an object ``stype`` specifies what attribute to test.
        - **mode**, optional: What you get if something found.

            - ``"loc"``: Get the location lowercase letter (like "u").
            - ``"ext"``: Same as ``"loc"`` but searches for the :attr:`position` first
              if case of a CDAT variable,
              then the :attr:`_vacumm_cf_physloc` attributes or the "physloc" specification
              if ``var`` has an entry in :attr:`CF_VAR_SPECS`, else defaults to
              :attr:`DEFAULT_LOCATION`.
            - Else, the matched string.

    :Example:

        >>> loc = get_loc('toto_u')
        >>> suffix = get_loc('toto_at_u_location', stype='standard_name', mode='full')
        >>> loc = get_loc(myvar, stype=['id','long_name'])

    :Return:

        - ``None`` if no location spec found,
        - if ``mode=="loc"``, one of possible physical
          :attr:`~vacumm.data.misc.arakawa.ARAKAWA_LOCATIONS`,
        - else, the complete matching string, like "_at_u_location".
    """
    # What to test
    if stype is None: stype = _reloc.keys() # all by default
    stypes = stype if isinstance(stype, list) else [stype]

    # CDAT object (not a string)
    if isinstance(var, list):
        if len(var)==0:
            return default
        var = var[0]
    if not isinstance(var, basestring):

        # First: _vacumm_cf_location
        if getattr(var, '_vacumm_cf_location', None):
            return var._vacumm_cf_location.lower()

        # Location from ids
        for stype in stypes:
            if stype=='id':
                loc = get_loc(var.id, stype, mode=mode)
            elif hasattr(var, stype):
                loc = get_loc(getattr(var, stype), stype, mode)
            else:
                continue
            break
        if loc is not None:
            return loc

        # Ext mode
        if mode=='ext':

            # Position attribute
            if hasattr(var, 'position'): return var.position.lower()

            # Physloc
            gid = var.id if not hasattr(var, '_vacumm_cf_name') else var._vacumm_cf_name
            gid = no_loc_single(gid, 'id')
            if gid in CF_VAR_SPECS:

                # From _vacumm_cf_physloc attribute
                if hasattr(var, '_vacumm_cf_physloc'):
                    return var._vacumm_cf_physloc

                # From specs
                if CF_VAR_SPECS[gid].get('physloc'):
                    return CF_VAR_SPECS[gid]['physloc']

        return default

    # String
    outloc = mode=='loc' or mode=='ext'
    for stype in stypes:

        # From string
        group = 2 if outloc else 1
        if stype not in _reloc: stype = 'id'
        loc = _reloc[stype](var).group(group)
        if loc and outloc: loc = loc.lower()
        if loc: return loc

        # Extended mode
        if mode=='ext':

            # From physloc
            if stype=='id' and var in CF_VAR_SPECS and CF_VAR_SPECS[var].get('physloc'):
                return CF_VAR_SPECS[var].get('physloc')

#            # Default location
#            return DEFAULT_LOCATION
    return default


def no_loc_single(name, stype):
    """Remove location specification

    :Params:

        - **name**: String to change
        - **stype**: One of 'id', 'standard_name' or 'long_name'.
    """
    loc = get_loc(name, stype, mode='full')
    if loc is not None:
        return name[:-len(loc)]
    return name

def _loc_(loc=None):
    """Get location as a single char or empty lowercase string"""
    if loc is None: loc = ''
    if loc.startswith('_'): loc = loc[1:]
    loc = loc[:1]
    if loc not in string.ascii_letters:
        raise TypeError('Wrong location: '+loc)
    return loc.lower()

def change_loc_single(name, stype, loc, squeeze=None):
    """Change location specification

    :Params:

        - **name**: String to change
        - **stype**: String type: name, standard_name or long_name.
        - **loc**: Location as a letter or None.
        - **squeeze**, optional: If specified and equal to ``loc``, location is removed
          instead of being changed.

    :Example:

        >>> change_loc_single('usurf_u', 'name', 't')
        'usurf_t'
        >>> change_loc_single('sst_at_u_location', 'standard_name', 'v')
        'sst_at_v_location'
        >>> change_loc_single('usurf_t', 'name', 'u', squeeze='u')
        'usurf'
    """
    basename = no_loc_single(name, stype)
    loc = _loc_(loc)
    if loc and (not squeeze or squeeze!=loc):
        if stype=='standard_name':
            return basename+'_at_%s_location'%loc
        if stype=='long_name':
            return basename+' at %s location'%loc.upper()
        return basename+'_'+loc
    return basename

def change_loc_specs(loc, id=None, standard_name=None, long_name=None, axes=None,
    iaxis=None, jaxis=None, squeeze=None, **kwargs):
    """Change location specification in id, standard names long name and axes id

    :Example:

        >>> specs = change_loc_specs('v', id = ['usurf_u', 'u_u'])

    :Return: A dictionary
    """
    specs = kwargs.copy()
    if 'atlocs' in specs: del specs['atlocs']

    # Attributes
    for stype in 'id', 'standard_name', 'long_name', 'iaxis', 'jaxis':
        values = locals()[stype]
        if values is None: continue
        if isinstance(values, list):
            values = list(values)
            tmp = [change_loc_single(value, stype, loc, squeeze=squeeze) for value in values]
            values = []
            for value in tmp: # unique
                if value not in values:
                    values.append(value)

        else:
            values = change_loc_single(values, stype, loc, squeeze=squeeze)
        specs[stype] = values

    # Axes
    if axes is not None:
        axes = axes.copy()
        for l in axes.keys(): # loop on x, y
            if l=='t': continue # skip time
            laxes = axes[l]
            single = not isinstance(laxes, list)
            if single: laxes = [laxes] # work on lists
            laxes = [change_loc_single(axis, 'id', loc, squeeze=DEFAULT_LOCATION) for axis in laxes]+laxes # duplicate
            lnewaxes = []
            for axis in laxes: # unique
                if axis not in lnewaxes:
                    lnewaxes.append(axis)
            if single and len(lnewaxes)==1:
                lnewaxes = lnewaxes[0]
            axes[l] = lnewaxes
        specs['axes'] = axes

    return specs

def change_loc(var, toloc, axes=True, squeeze=True):
    """Change location specifications of a variable and its axes

    It affects the id, the standard_name and the long_name when they are defined.

    :Params:

        -

    :Example:

        >>> change_loc(sst, 'u').id
        'sst_u'
        >>> change_loc(sst, 't', squeeze=True).id
        'sst'

    """

    # Squeeze physloc
    if not isinstance(squeeze, basestring):
        if squeeze:
            squeeze = get_physloc(var)
        else:
            squeeze = None

    # Change attributes
    # - get
    specs = change_loc_specs(toloc, squeeze=squeeze,
        id=var.id,
        standard_name = getattr(var, 'standard_name', None),
        long_name = getattr(var, 'long_name', None),
        )
    # - set
    var.id = specs['id']
    for att in 'standard_name', 'long_name':
        if att in specs:
            setattr(var, att, specs[att])
    # - cf name
    if hasattr(var, '_vacumm_cf_name'):
        var._vacumm_cf_name = change_loc_single(var._vacumm_cf_name, 'id', toloc,
            squeeze=squeeze)

    # Change axes and grid attributes
    if axes and cdms2.isVariable(var) or cdms2.isGrid(var):

        # Axes
        # - usual
        for meth in 'Longitude', 'Latitude', 'Time', 'Level':
            if hasattr(var, 'get'+meth):
                axismet = getattr(var, 'get'+meth) # var.getLongitude
                if axismet() is not None:
                    change_loc(axismet(), toloc, squeeze=squeeze)
        # - 2d axes
        if cdms2.isVariable(var) and isaxis(var) and var.ndim>2:
            for i in -1, -2:
                change_loc(var.getAxis(i), toloc, squeeze=squeeze)

        # Grid
        if cdms2.isVariable(var) and var.getGrid() is not None:
            change_loc(var.getGrid(), toloc, squeeze=squeeze)

    # Reference attribute
    set_loc(var, toloc)

    return var

def set_loc(var, loc, addpos=False):
    """Define (or remove) the location of a variable with :attr:`_vacumm_cf_location` attribute

    It also makes sure that the :attr:`position` attribute is the same if present.

    :Params:

        - **var**: CDAT variable.
        - **loc**: Location. If empty, remove :attr:`_vacumm_cf_location` and
          :attr:`position` attributes.
        - **addpos**, optional: Set also the :attr:`position` attribute.

    """
    # Remove
    if not loc:
        for att in '_vacumm_cf_location', 'position':
            if hasattr(var, att):
                delattr(var, att)

    # Set
    loc = loc.lower()
    var._vacumm_cf_location = loc
    if addpos:
        var.position = loc
    elif hasattr(var, 'position') and var.position.lower()!=loc:
        var.position = loc
    return var


def get_physloc(var):
    """Get the physical location of a variable

    .. note:: The physical location may be different from the current location.

    It defaults to :attr:`DEFAULT_LOCATION`
    """
    # Direct
    for att in ['_vacumm_cf_physloc']:
        if hasattr(var, att):
            return getattr(var, att)

    # Entry name in CF_VAR_SPECS (string)
    if isinstance(var, basestring):
        if var in CF_VAR_SPECS:
            var = CF_VAR_SPECS[var]
        else:
            return DEFAULT_LOCATION

    # Dict of specs
    if isinstance(var, dict):
        return var.get('physloc', DEFAULT_LOCATION)

    # Using CF_VAR_SPECS
    if hasattr(var, '_vacumm_cf_name'):
        name = var._vacumm_cf_name
    elif var.id in CF_VAR_SPECS:
        name = var.id
    else:
        return DEFAULT_LOCATION
    return CF_VAR_SPECS[name].get('physloc', DEFAULT_LOCATION)

def squeeze_loc(var):
    """Squeeze out location specification of a variable if the location is
    physloc or :attr:`DEFAULT_LOCATION`

    :Params:

        - **var**: CDAT variable.
        - **physloc**, optional: ``physloc`` specs given by the entry in :attr:`CF_VAR_SPECS`
          or :attr:`DEFAULT_LOCATION`.

    :Return: The same variable
    """
    # Physloc
    physloc = get_physloc(var)

    # Remove loc only for physloc
    change_loc(var, None, axes=True, squeeze=physloc)

def squeeze_loc_single(name, stype, physloc=None):
    """Squeeze location specs if it matches physical location"""
    if physloc is None and stype=='id' and name in CF_VAR_SPECS:
        physloc = CF_VAR_SPECS[name].get('physloc')
    if physloc is None:
        physloc = DEFAULT_LOCATION
    loc = get_loc(name, stype, mode='loc')
    if loc==physloc:
        name = change_loc_single(name, stype, None)
    return name

def dupl_loc_specs(all_specs, fromname, toloc):
    """Duplicate the specification for a variable or an axis to another or several locations

    The following rules apply:

        - If the original specifications are from a name without a specific location
          (generic), new specs (located) are appended (merged) with original specs.
        - If the specifications at target location already exist,
          it merges new specs with old specs.
        - Generic specification are systematically created or updated by merging of
          specialized ones.

    :Example:

        >>> dupl_loc_specs(CF_VAR_SPECS, 'corio', 'u') # Create the 'corio_u' entry in CF_VAR_SPECS and update 'corio'
        >>> dupl_loc_specs(CF_VAR_SPECS, 'corio_t', 'u') # Create 'corio_u' and 'corio_t'

    """
    if not fromname in all_specs:
        raise KeyError('No such entry in specifications: '+fromname)
    single = not isinstance(toloc, (list, tuple))
    if single:
        toloc = [toloc]
    tonames = []
    tomerge = []
    fromnoloc = no_loc_single(fromname, 'id')==fromname
    for loc in toloc:

        # New name (id)
        toname = change_loc_single(fromname, 'id', loc)
        tonames.append(toname)

        # New specs
        tospecs = change_loc_specs(loc, **all_specs[fromname])

        # Add a version of standard_name and long_name for T loc without location spec
        if loc=='t': # 'toto_at_t_location' -> ['toto_at_t_location','toto']
            for spname in 'standard_name', 'long_name':
                if spname not in tospecs: continue
                addspecs = change_loc_specs(None, **{spname:tospecs[spname]})
                for dd in tospecs, addspecs: # to lists
                    if not isinstance(dd[spname], list):
                        dd[spname] = [dd[spname]]
                tospecs[spname].extend(addspecs[spname])

        # For merging specs back with old specs if old has no location
#        if fromnoloc: # generic
        tomerge.append(tospecs)

        # Merge with old spec if existing
        if toname in all_specs:
            tospecs = dict_merge(tospecs, all_specs[toname], mergelists=True,
                                 cls=dict, mergedicts=True, )

#        # Default location -> add its generic version ('_t' -> + '')
#        if loc==DEFAULT_LOCATION:
#            genspecs = change_loc_specs(None, **tospecs)
#            tospecs = dict_merge(tospecs, genspecs, mergelists=True)

        # Remove atlocs attribute
        if 'atlocs' in tospecs:
            tospecs.pop('atlocs')

        # Store it
        all_specs[toname] = tospecs

    # Make a generic entry that merges all specialized ones
    if fromnoloc:
        genspecs = all_specs[fromname] # existing generic specs
    else:
        genspecs = change_loc_specs(None, **all_specs[fromname]) # remove location info
    tomerge.insert(0, genspecs)
    genname = genspecs['id'][0]
    all_specs[genname] = dict_merge(*tomerge, **CF_DICT_MERGE_KWARGS)

    if single: return tonames[0]
    return tonames

#def specs_def_loc(all_specs, names, suffixes=ARAKAWA_SUFFIXES):
#    """Make specs of a variable at a special location the specs for default location
#
#    :Example:
#
#        >>> specs_def_loc(all_specs, ['usurf_u', 'ssh_f']) # specs of 'usurf_u' is copied to 'usurf'
#    """
#    if not isinstance(names, (list, tuple)): names = [names]
#    if isinstance(suffixes, basestring): suffixes = [suffixes]
#    for name in names:
#        m = re.match('.+(%s)$'%'|'.join(suffixes), name)
#        if m is None: continue
#        dupl_loc_specs(all_specs, name, '', suffixes=ARAKAWA_SUFFIXES)


def cp_suffix(idref, id, suffixes=ARAKAWA_SUFFIXES):
    """Copy a suffix if found in an id (name) to another id"""
    if isinstance(suffixes, basestring): suffixes = [suffixes]
    m = re.match('.+(%s)$'%'|'.join(suffixes), idref)
    if m is None: return id
    return id+m.group(1)


#: List of generic variable names
CF_VAR_NAMES = []
generic_var_names = GENERIC_VAR_NAMES = GENERIC_CF_VAR_NAMES = CF_VAR_NAMES # compat

#: List of generic axis names
CF_AXIS_NAMES = []
generic_axis_names = GENERIC_AXIS_NAMES = GENERIC_CF_AXIS_NAMES = CF_AXIS_NAMES # compat

#: Specifications for variables
CF_VAR_SPECS = VAR_SPECS = var_specs = VarSpecs(names=CF_VAR_NAMES)

#: Specifications for axes
CF_AXIS_SPECS = AXIS_SPECS = axis_specs = AxisSpecs(inherit=CF_VAR_SPECS,
                                                    names=CF_AXIS_NAMES)


#: Specifications for grid formating
GRID_SPECS = {
    '': dict(lon='lon', lat='lat', level='depth'),
    't': dict(lon='lon', lat='lat', level='depth_t'),
    'u': dict(lon='lon_u', lat='lat_u', level='depth_t'),
    'v': dict(lon='lon_v', lat='lat_v', level='depth_t'),
    'f': dict(lon='lon_f', lat='lat_f', level='depth_t'),
    'w': dict(lon='lon', lat='lat', level='depth_w'),
}
#TODO: 't' must use lon_t and lat_t
GRID_SPECS['r'] = GRID_SPECS['t']
GRID_SPECS[None] = GRID_SPECS['']
grid_specs = GRID_SPECS # compat


#: List of all generic names (axes and variables): DEPRECATED!!
CF_NAMES = CF_VAR_NAMES + CF_AXIS_NAMES
generic_names = GENERIC_NAMES  = CF_NAMES# compat

def register_cf_variable(name, **specs):
    """Register a new CF variable given its generic name and its specs

    It simply calls :meth:`~VarSpecs.register` on :attr:`CF_VAR_SPECS`

    """
    CF_VAR_SPECS.register(name, **specs)
    return CF_VAR_SPECS[name]

def register_cf_variables_from_cfg(cfg):
    """Register a new CF variables from a config file or a dict

    It must contains a "variables" section.

    It simply calls :meth:`~VarSpecs.register_from_cfg` on :attr:`CF_VAR_SPECS`

    """
    CF_VAR_SPECS.register_from_cfg(cfg)

def register_cf_axis(name, **specs):
    """Register a new Cf axis given its generic name and its specs

    It simply calls :meth:`~VarSpecs.register` on :attr:`CF_AXIS_SPECS`

    """
    CF_AXIS_SPECS.register(name, **specs)
    return CF_AXIS_SPECS[name]

def register_cf_axes_from_cfg(cfg):
    """Register a new CF variables from a config file or a dict

    It must contains a "axes" section.

    It simply calls :meth:`~VarSpecs.register_from_cfg` on :attr:`CF_AXIS_SPECS`

    """
    CF_AXIS_SPECS.register_from_cfg(cfg)


def cf2search(name, mode='isa', raiseerr=True, **kwargs):
    """Extract specs from :attr:`CF_AXIS_SPECS` or :attr:`CFVAR_SPECS`
    to form a search dictionary

    :Params:

        - **name**: Generic name of an axis or a variable.
        - **mode**, optional: Search mode [default: None->``"ns"``].
          A string containg one or more of the following letters:

            - ``"n"``: Search using names (ids).
            - ``"s"``: Search using standard_name attribute.
            - ``"a"``: Search using axis attribute.
            - ``"l"``: Search using long_name attribute.
            - ``"u"``: Search using units attribute.

          The order is important.

    :Return: An :class:`colections.OrderedDict`

    :Example:

        >>> cf2search('sst', mode='isu')
        {'id':['sst'], 'id':['sst'],
        'standard_names':['sea_surface_temperature'],
        'units':['degrees_celsius']}
    """
    # Get specs
    if name in CF_VAR_SPECS:
        specs = CF_VAR_SPECS[name]
    elif name in CF_AXIS_SPECS:
        specs = CF_AXIS_SPECS[name]
    else:
        if raiseerr:
            raise VACUMMError("Wrong generic name. It should be one of: "+' '.join(CF_AXIS_SPECS.keys()+CF_VAR_SPECS.keys()))
        else:
            return

    # Form search dict
    if not isinstance(mode, basestring):
        mode = 'isa'
    mode = mode.replace('n', 'i')
    keys = []
    for m in mode:
        for key in ['id', 'standard_name', 'axis', 'long_name', 'units']:
            if key.startswith(m):
                keys.append(key)
                break
    return OrderedDict([(k, specs[k]) for k in keys if k in specs])


_attnames_exclude = ['atlocs', 'inherit', 'axes', 'physloc', 'iaxis', 'jaxis', 'select']
_attnames_firsts = ['standard_name', 'long_name', 'units', 'axis', 'valid_min', 'valid_max']
def cf2atts(name, select=None, exclude=None, ordered=True, **extra):
    """Extract specs from :attr:`CF_AXIS_SPECS` or :attr:`CF_VAR_SPECS` to form
    a dictionary of attributes (units and long_name)"""
    # Get specs
    if isinstance(name, dict):
        specs = name.copy()
    elif name in CF_VAR_SPECS:
        specs = CF_VAR_SPECS[name]
    elif name in CF_AXIS_SPECS:
        specs = CF_AXIS_SPECS[name]
    else:
        raise VACUMMError("Wrong generic name: %s. It should be one of: "%name+' '.join(CF_AXIS_SPECS.keys()+CF_VAR_SPECS.keys()))

    # Which attributes
    atts = OrderedDict() if ordered else {}
    if exclude is None: exclude = []
    elif isinstance(exclude, basestring): exclude = [exclude]
    exclude.extend(_attnames_exclude)
    for key in _attnames_firsts+specs.keys():

        # Skip aliases
        if key in BaseSpecs.get_alias_list():
            continue

        # Skip some attributes
        if (key not in specs or key in exclude or key in atts or
            (select is not None and key not in select)):
            continue

        # No lists or tuples
        value = specs[key]
        if isinstance(value, (list, tuple)):
            if len(value)==0:
                continue
            value = value[0]

        # Store it
        atts[key] = value

    # Extra
    for att, val in extra.items():
        atts[att] = val

    return atts



# Format a variable
def format_var(var, name=None, force=True, format_axes=True, order=None, nodef=True,
        mode='warn', **kwargs):
    """Format a MV2 variable according to its generic name


    :Params:

        - **var**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Generic name of variable. It should be one of
          those listed by :attr:`CF_VAR_SPECS`. If None, it is guessed
          with :func:`match_known_var`.
        - **force**, optional: Overwrite attributes in all cases.
        - **format_axes**, optional: Also format axes.
        - **nodef**, optional: Remove location specification when it refers to the
          default location (:attr:`DEFAULT_LOCATION`).
        - **mode**: "silent", "warn" or "raise".
        - Other parameters are passed as attributes, except those:

            - present in specifications to allow overriding defaults,
            - starting with 'x', 'y', or 't' which are passed
              to :func:`format_axis`.

    :Examples:

        >>> var = format_var(myarray, 'sst', valid_min=-2, valid_max=100)

    """
    # Filter keywords for axis formating
    axismeths = {'t':'getTime', 'y':'getLatitude', 'x':'getLongitude'}
    kwaxes = {}
    for k in axismeths.keys():
        kwaxes[k] = kwfilter(kwargs, k+'_')

    # Always a MV2 array
    if not cdms2.isVariable(var):
        var = MV2.asarray(var)

    # Check specs
    if name is None: # guess it
        name = match_known_var(var)
        if not name:
            if mode=='warn':
                warn("Can't guess cf name")
                return var
            elif mode=='silent':
                return var
            else:
                raise KeyError("Variable does not match any CF var")
    elif name not in CF_VAR_SPECS and name not in CF_AXIS_SPECS:
        if var.id in CF_VAR_SPECS or var.id in CF_AXIS_SPECS:
            name = var.id
        elif mode=='warn':
            warn("Generic var name not found '%s'."%name)
            return var
        elif mode=='silent':
            return var
        else:
            raise KeyError("Generic var name not found '%s'. Please choose one of: %s"%(
                name, ', '.join(CF_VAR_SPECS.keys()+CF_AXIS_SPECS.keys())))
    isaxis = name in CF_AXIS_SPECS
    if isaxis:
        specs = CF_AXIS_SPECS[name].copy()
        if 'axis' in specs:
            del specs['axis']
    else:
        specs = CF_VAR_SPECS[name].copy()
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - remove default location
    if nodef:
        refloc = specs.get('physloc', None) or DEFAULT_LOCATION
        for att in 'id', 'long_name', 'standard_name':
            if get_loc(specs[att], att)==refloc:
                specs[att] = [no_loc_single(specs[att][0], att)]
        name = specs['id'][0]
    # - id
    if ((force is True or force in [2, 'id', 'all'])
            or var.id.startswith('variable_') or
            (isaxis and var.id.startswith('axis_'))): # FIXME: use regexp
        var.id = name
    # - attributes
    forceatts = (force is True or force in ['atts', 'all'] or
        (isinstance(force, int) and force>0))
    for att, val in cf2atts(specs, **kwargs).items():
        if forceatts or not getattr(var, att, ''):
            setattr(var, att, val)
    # - physical location
    loc = get_loc(var, mode='ext')
    if not loc and 'physloc' in specs:
        loc = specs['physloc']
    if loc:
        if 'physloc' in specs and loc == specs['physloc']:
            var._vacumm_cf_physloc = loc.lower()
        set_loc(var, loc)
    # - store cf name
    var._vacumm_cf_name = name

    # Axes
    if format_axes:

        # Order
        order = var.getOrder() if not isinstance(order, basestring) else order
        if order is not None:
            if not re.match('^[xyzt-]+$', order):
                raise VACUMMError("Wrong cdms order type: "+order)
            if len(order)!=var.ndim:
                raise VACUMMError("Cdms order should be of length %s instead of %s"%(var.ndim, len(order)))

        # First check
        if 'axes' in specs:
            axspecs = specs['axes']
            formatted = []
            for key, meth in axismeths.items():
                axis = getattr(var, meth)()
                if order is not None: order.replace(key, '-')
                if axis is not None:
                    format_axis(axis, axspecs[key], **kwaxes[key])
                    formatted.append(key)

            # Check remaining simple axes (DOES NOT WORK FOR 2D AXES)
            if order is not None and order!='-'*len(order):
                for key in axismeths.keys():
                    if key in order and key not in formatted:
                        axis = var.getAxis(order.index(key))
                        format_axis(axis, axspecs[key], **kwaxes[key])


    return var

# Format an axis
def format_axis(axis, name=None, force=True, recreate=False, format_subaxes=True,
        nodef=True, mode='warn', **kwargs):
    """Format a MV2 axis according to its generic name


    :Params:

        - **axis**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Single or list of generic axis names. It should be one of
          those listed by :attr:`CF_AXIS_SPECS`.If None, it is guessed
          with :func:`match_known_axis`.
        - **force**, optional: Overwrite attributes in all cases.
        - **nodef**, optional: Remove location specification when it refers to the
          default location (:attr:`DEFAULT_LOCATION`).
        - **axes2d_<param>**, optional: <param> is passed to
          :func:`~vacumm.misc.grid.misc.create_axes2d` for 2D axes.
        - **recreate**, optional: Force recreation of axes using either

        - Other parameters are passed as attributes.

    :Examples:

        >>> axis1d = format_axis(array1d, 'lon')
        >>> axis2d = format_axis(array2d, 'lon_u', axes2d_iid = 'xi',  axes2d_jid = 'xi', recreate=True)


    """

    # Guess best generic name from a list
    if name is None: return axis
    if isinstance(name, (list, tuple)):
        for nm in name:
            if match_obj(axis, nm):
                name = nm
                break

    # Check specs
    if name is None: # guess it
        name = match_known_axis(axis)
        if not name:
            if mode=='warn':
                warn("Can't guess CF axis name")
                return axis
            elif mode=='silent':
                return axis
            else:
                raise KeyError("Variable does not match any CF axis")
    elif name not in CF_AXIS_SPECS:
        if mode=='warn':
            warn("Generic axis name not found '%s'."%name)
            return axis
        elif mode=='silent':
            return axis
        else:
            raise KeyError("Generic axis name not found '%s'. Please choose one of: %s"%(
                name, ', '.join(CF_AXIS_SPECS.keys())))
    specs = CF_AXIS_SPECS[name]
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]

    # Always a MV2 axis (1D or 2D)
    axis._oldid = axis.id
    kwaxed2d = kwfilter(kwargs, 'axes2d_')
    if not isaxis(axis) or recreate:
        if len(axis.shape)==1:
            axis = cdms2.createAxis(axis)
        else:
            xy = specs['axis'].lower()
            kwaxed2d[xy] = axis
            axis = create_axes2d(**kwaxed2d)
            return axis
    axis2d = len(axis.shape)==2

    # Apply specs
    # - merge kwargs and specs
    for key, val in kwargs.items():
        if val is None or key not in specs: continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - remove default location
    if nodef:
        for att in 'id', 'long_name', 'standard_name':
#            sname = stype+'s'
            if att not in specs: continue
            if get_loc(specs[att], att)==DEFAULT_LOCATION:
                specs[att] = [no_loc_single(specs[att][0], att)]
    # - id
    if force or axis.id.startswith('variable_') or axis.id.startswith('axis_'):
        axis.id = specs['id'][0]
    # - attributes
    for att, val in cf2atts(specs, exclude=['axis'] if axis2d else None, **kwargs).items():
        if force or not getattr(axis, att, ''):
            setattr(axis, att, val)
    # - store cf name
    axis._vacumm_cf_name = axis.id

    # Sub-axes for 2D axes
    if axis2d and format_subaxes:
        format_axis(axis.getAxis(-1), specs['iaxis'])
        format_axis(axis.getAxis(-2), specs['jaxis'])

    return axis

def format_grid(grid, pt, **kwargs):
    """Format a grid and its axes"""
    if cdms2.isVariable(grid): grid = grid.getGrid()
    if grid is None: return
    gs = GRID_SPECS[pt]
    lon = grid.getLongitude()
    lat = grid.getLatitude()
    format_axis(lon, gs['lon'], **kwargs)
    format_axis(lat, gs['lat'], **kwargs)


#: Hidden attributes of variable useful of this module
HIDDEN_CF_ATTS = ['_vacumm_cf_name', '_vacumm_cf_physloc', '_vacumm_cf_location']
hidden_cf_atts = HIDDEN_CF_ATTS # compat

def match_obj(obj, name, searchmode=None, **kwargs):
    """Check if a variable or an axis matches generic specifications

    :Params:

        - **obj**: a numpy or cdms2 axis or variable.
        - **name**: A generic names.
        - **searchmode**, optional: Passed to :func:`~vacumm.misc.io.ncmatch_obj`.
    """
    search = cf2search(name, mode=searchmode, raiseerr=False)
    if search is None: return False
    search.update(kwargs)
    return ncmatch_obj(obj, searchmode=searchmode, **search)

match_var = match_obj

def match_known_var(obj, searchmode=None, **kwargs):
    """Check if an object matches a known variable"""
    for name in CF_VAR_SPECS:
        if match_obj(obj, name, searchmode=searchmode):
            return name
    return False

def match_known_axis(obj, searchmode=None, **kwargs):
    """Check if an object matches a known axis"""
    for name in CF_AXIS_SPECS:
        if match_obj(obj, name, searchmode=searchmode):
            return name
    return False

def is_cf_known(name, type=None):
    """Is this names registered in :attr:`CF_VAR_SPECS` or
    :attr:`CF_AXIS_SPECS`"""
    assert type is None or type in ('var', 'axis')
    if type is None:
        names = CF_AXIS_NAMES + CF_VAR_NAMES
    elif type=='var':
        names = CF_VAR_NAMES
    else:
        names = CF_AXIS_NAMES
    return name in names


def get_cf_cmap(vname):
    """Get a cmap from a standard name

    cmap may be specified for a variable with the 'cmap' key of
    its entry in :attr:`CF_VAR_SPECS`.
    """
    if hasattr(vname, 'id'):
        vname = vname.id
    if vname in CF_VAR_SPECS and 'cmap' in CF_VAR_SPECS[vname]:
        cmap = get_cmap(CF_VAR_SPECS[vname]['cmap'])
        if cmap.name != CF_VAR_SPECS[vname]['cmap']:
            vcwarn("Can't get cmap '{}' for standard variable '{}'".format(
                CF_VAR_SPECS[vname]['cmap'], vname))

    return cmap
