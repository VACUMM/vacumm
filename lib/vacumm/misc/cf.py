# -*- coding: utf8 -*-
"""Conventions about data formats and names"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2018)
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
from collections import OrderedDict, UserString
import operator
import string
import re
from copy import deepcopy
from pprint import pformat

import cdms2
import MV2

from vacumm import VACUMMError, vcwarn
from .misc import kwfilter, dict_merge, match_atts, ArgList
from .config import ConfigManager
from .arakawa import ARAKAWA_LOCATIONS

__all__ = [
           'format_var', 'format_coord', 'format_grid', 'match_cf_atts',
           'cf2atts', 'cf2search', 'cp_suffix', 'get_loc',
           'change_loc', 'change_loc_single',
           'dupl_loc_specs', 'no_loc_single',
           'change_loc_specs', 'squeeze_loc_single',
           'squeeze_loc', 'get_physloc',
           'HIDDEN_CF_ATTS', 'set_loc', 'parse_loc',
           'match_known_cf_var', 'match_known_cf_coord',
           'match_cf_from_name', 'match_known_cf_obj',
           'is_known_cf_var_name', 'is_known_cf_coord_name',
           'register_cf_variable', 'register_cf_variables_from_cfg',
           'register_cf_coord', 'register_cf_coords_from_cfg',
           'CF_DICT_MERGE_KWARGS', 'get_cf_cmap',
           'is_known_cf_name', 'match_var', 'match_coord',
           'get_cf_specs', 'get_cf_standard_name',
           'CFSpecs', 'CFCatSpecs', 'standard_name_to_location',
           ]

ARAKAWA_SUFFIXES = [('_' + p) for p in ARAKAWA_LOCATIONS]
arakawa_suffixes = ARAKAWA_SUFFIXES  # compat
arakawa_locations = ARAKAWA_LOCATIONS

THISDIR = os.path.dirname(__file__)

#: Joint variables and coords config specification file
CF_GFG_INIFILE = os.path.join(THISDIR, 'cf.ini')

#: Base config file for CF specifications
CF_CFG_CFGFILE = os.path.join(THISDIR, 'cf.cfg')

#: File used for pickling :data:`CF_SPECS`
_user_cache_dir = os.path.join(os.path.expanduser('~'), '.cache')
_cache_dir = os.path.join(os.environ.get("XDG_CACHE_HOME", _user_cache_dir),
                          'vacumm')
CF_SPECS_PICKLED = os.path.join(_cache_dir, 'cf.pyk')

#: Argument passed to dict_merge to merge CF configs
CF_DICT_MERGE_KWARGS = dict(mergesubdicts=True, mergelists=True,
                            skipnones=False, skipempty=False,
                            overwriteempty=True, mergetuples=True)



_re_match_cache = {}


def _get_re_match_(pattern):
    if pattern not in _re_match_cache:
        _re_match_cache[pattern] = re.compile(pattern, re.I).match
    return _re_match_cache[pattern]


class SGLocator(UserString):
    """Staggered grid location parsing and formatting utility"""

    formats = {'name': '{root}_{loc}',
               'standard_name': '{root}_at_{loc}_location',
               'long_name': '{root} at {loc} location'}

    def __init__(self, letter=None, attrs=None, allowed='uvwtdpfr',
                 default='t', formats=None):

        # Init
        self.allowed = allowed
        self.default = default
        self.formats = self.formats.copy()

        # Build regular expressions for searching
        self.rematch = {}
        for attr, fmt in self.formats.items():
            if formats and attr in formats:
                self.formats[attr] = fmt
            pattern = fmt.format(root='(?P<root>.+)',
                                 loc='(?P<loc>[{}])'.format(self.allowed))
            self.rematch[attr] = _get_re_match_('^' + pattern + '$')

        # Parse attributes
        if attrs is None:
            attrs = {}
        else:
            loc, attrs = self.parse_attrs(attrs)
            if letter is None:
                letter = loc
        self._loc = letter or ''
        self._attrs = attrs

        UserString.__init__(self, self._loc)

    def __repr__(self):
        return '<{0}("{1}") at {2}>'.format(self.__class__.__name__,
                                            self, hex(id(self)))

    def parse_loc(self, name, attr, mode='split'):
        """Parse location from a name, standard_name or long_name

        Parameters
        ----------
        name: str
            The string to parse
        attr: str
            The attribute type within 'name', 'standard_name' and 'long_name'
        mode: str
            Ouput mode:

            - "loc": location letter
            - "re": the regular expression object
            - "match": the match string that specify the location,
              like "at_u_location")
            - "split": a dict containing the root, the location and
              the matched location string and the match regular expression
              object.
        """
        # Match
        if attr not in self.rematch:
            attr = 'name'
        m = self.rematch[attr](name)
        if mode == 're':
            return m

        # Fields
        if m is None:
            loc = None
            root = name
            match = None
        else:
            loc = m.group('loc')
            root = m.group('root')
            match = self.formats[attr].format(loc=loc, root=root)
            if loc:
                loc = loc.lower()

        # Return
        if mode == 'loc':
            return loc
        if mode == 'match':
            return match
        if mode == 'root':
            return root
        return dict(loc=loc, root=root, match=match, re=m)

    def parse_attrs(self, attrs):
        """Parse attributes to find location and roots

        The following attributes are parsed in order of priority:
        standard_name, long_name and name.

        Parameters
        ----------
        attrs: dict
            A dict of attributes.
            An attribute may be a string or list of strings.

        Return
        ------
        str or None
            The location
        attrs: dict
            A copy of the attributes dict with the location specifications
            squeezed out.
        """

        attrs = attrs.copy()
        dlocs = {}
        for attr in 'standard_name', 'long_name', 'name':
            if attr in attrs:
                al = ArgList(attrs[attr])  # to handle lists
                values = al.get()
                dlocs[attr] = []
                roots = []
                for value in values:
                    search = self.parse_loc(value, attr, mode='split')
                    if search['re'] is not None:
                        dlocs[attr].append(search['loc'].lower())
                        root = search['root']  # squeeze loc
                    else:
                        root = value
                    roots.append(root)
                attrs[attr] = al.put(roots)
        locs = []
        for ll in dlocs.values():
            locs.extend(ll)
        if not locs:
            loc = ''
        else:
            if len(locs) > 1 and any([(loc != locs[0]) for loc in locs[1:]]):
                vcwarn("Conflicting location specification in attributes: "
                       + str(dlocs))
            loc = locs[0]
        return loc, attrs

    def parse_dataarray(self, da):
        """Parse a :class:`xarray.DataArray` name and attributes"""
        attrs = da.attrs.copy()
        if da.name is not None:
            attrs['name'] = da.name
        return self.parse_attrs(attrs)

    @classmethod
    def from_name(cls, name):
        return cls(attrs={'name': name})

    @classmethod
    def from_standard_name(cls, name):
        return cls(attrs={'standard_name': name})

    @classmethod
    def from_long_name(cls, name):
        return cls(attrs={'long_name': name})

    @classmethod
    def from_dataarray(cls, da):
        attrs = da.attrs.copy()
        if da.name:
            attrs.update(name=da.name)
        return cls(attrs=attrs)

    def get_attrs(self):
        return self._attrs

    def set_attrs(self, attrs):
        _, self._attrs = self.parse_attrs(attrs)

    attrs = property(get_attrs, set_attrs, doc='Root attributes')

    def get_loc(self):
        return self._loc

    def set_loc(self, loc):
        loc = loc.lower()
        if loc:
            assert loc in self.allowed, (
                    'Invalid location: {}. Please choose one '
                    'of: {}'.format(loc, ' '.join(self.allowed)))
        else:
            loc = ''
        self._loc = loc

    loc = property(get_loc, set_loc, doc='Location in lower case')

    @property
    def LOC(self):
        """Location in upper case"""
        return self._loc.upper()

    def to(self, letter):
        """Return a copy at another location with the same attributes"""
        return SGLocator(letter, attrs=self.attrs)

    def format_attrs(self, attrs=None):
        """Format attributes with current location

        Parameters
        ----------
        attrs: dict or None

        Return
        ------
        dict
        """
        if not attrs:
            attrs = self.attrs.copy()
        else:
            _, attrs = self.parse_attrs(attrs)
        if self.loc:
            for attr, fmt in self.formats.items():
                if attr in attrs:
                    loc = self.LOC if attr == 'long_name' else self.loc
                    al = ArgList(attrs[attr])
                    values = []
                    for root in al.get():
                        values.append(fmt.format(loc=loc, root=root))
                    attrs[attr] = al.put(values)
        return attrs

    def format_dataarray(self, da, mode='keep', coords=None, dims=None):
        """Format the data array with current location

        Parameters
        ----------
        da: xarray.DataArray
        mode: str
            Mode for handling attributes

            - keep: format existing attributes only.
            - erase: use current :class:`SGLocator` attributes.
            - merge: use :class:`SGLocator` attributes when missing in data.

        Return
        ------
        xarray.DataArray
            New array with same numeric content
        """
        assert mode in ('keep', 'merge', 'erase'), 'Invalid mode'
        da = da.copy(deep=False)

        # Get attributes
        if mode != 'erase':
            _, attrs = self.parse_dataarray(da)
            if mode == 'merge':
                for attr, value in self.attrs.items():
                    if attr not in attrs:
                        attrs[attr] = value
        else:
            attrs = self.attrs

        # Format attributes
        attrs = self.format_attrs(attrs)
        da.name = attrs.pop('name', None)
        da.attrs.update(attrs)

        # Coordinates
        if coords:
            for name in coords:
                if name in da.coords and da.coords[name] is not da:
                    da.coords[name] = self.format_dataarray(
                            da.coords[name], mode=mode)
                    new_name = self.format_attrs({'name': name})['name']
                    da = da.rename({name: new_name})

        # Dims
        if dims:
            for name in dims:
                if name in da.dims and name not in coords:
                    new_name = self.format_attrs({'name': name})['name']
                    da = da.rename({name: new_name})

        return da

    def format_cf_specs(self, specs):
        """Format CF specs for a variable or a coordinates

        Parameters
        ----------
        specs: dict

        Return
        ------
        dict
        """

        # Attributes
        new_specs = self.format_attrs(specs)

        # Coords or dims
        section = 'coords' if 'coords' in new_specs else 'dims'
        new_specs[section] = new_specs[section].copy()
        for dim in 'xy':  # FIXME: and z?
            if dim in new_specs[section]:
                new_specs[section][dim] = self.format_attrs(
                        {'name': new_specs[section][dim]})['name']

        return new_specs


class CFSpecs(object):
    aliases = {'name': ['ids', 'id', 'names'],
               'standard_name': ["standard_names"],
               'long_name': ['long_names']
               }

    def __init__(self, cfg):
        self._dict = self._load_cfg_(cfg)
        self._cfgspecs = get_cf_cfgm()._configspec.dict()

        self._post_process_()

    @staticmethod
    def _load_cfg_(cfg):
        """Load and validate a configuration without post-processing"""
        if isinstance(cfg, str) and '[' in cfg:
            cfg = cfg.split('\n')
        return OrderedDict(get_cf_cfgm().load(cfg).dict())

    @property
    def categories(self):
        return self._cfgspecs.keys()

    def __getitem__(self, category):
        assert category in self.categories
        return self._cfs[category]

    def __contains__(self, category):
        return category in self.categories

    def __getattr__(self, name):
        if name in self._cfs:
            return self._cfs[name]

    def __str__(self):
        return pformat(self._dict)

    def copy_and_update(self, cfg=None, specs=None):#, names=None):
        # Copy
        obj = self.__class__(cfg={})#names=names, inherit=inherit,
#                             cfg={self.category: {}}, parent=parent)
        obj._dict = deepcopy(self._dict)

        # Update
        if cfg:
            obj.register_from_cfg(cfg)

        return obj

    def copy(self):
        return self.copy_and_update()

    __copy__ = copy


    def register_from_cfg(self, cfg):
        """"Register new elements from a :class:`ConfigObj` instance
        or a config file"""
        local_cfg = self._load_cfg_(cfg)
        self._dict = dict_merge(self._dict, local_cfg,
                                cls=OrderedDict, **CF_DICT_MERGE_KWARGS)
        print(self._dict['variables']['sst'])
        self._post_process_()

    def _post_process_(self):
        self._cfs = {}
        self._from_atlocs = {}
        for category in self.categories:
            self._cfs[category] = CFCatSpecs(self, category)
            self._from_atlocs[category] = []
            for name in list(self._dict[category].keys()):#names:
                 if isinstance(self._dict[category][name], dict):
                     self._check_entry_(category, name)
        self._add_aliases_()
#        self._update_names_()

    def _check_entry_(self, category, name):
        """Validate an entry

        - Makes sure to have lists, except for 'axis' and 'inherit'
        - Check geo coords
        - Check inheritance
        - Makes sure that coords specs have no 'coords' key
        - Makes sure that specs key is the first entry of 'names'
        - add standard_name to list of names
        - Check duplication to other locations ('toto' -> 'toto_u')

        """
        # Dict of entries for this category
        entries = self._dict[category]

        # Wrong entry!
        if name not in entries:
            return

        # Entry already generated with the atlocs key
        if name in self._from_atlocs[category]:
            return

        # Get the specs as pure dict
        if hasattr(entries[name], 'dict'):
            entries[name] = entries[name].dict()
        specs = entries[name]

        # Ids
        if name in specs['name']:
            specs['name'].remove(name)
        specs['name'].insert(0, name)

        # Long name from name or standard_name
        if not specs['long_name']:
            if specs['standard_name']:
                long_name = specs['standard_name'][0]
            else:
                long_name = name.title().replace('_', ' ')
            specs['long_name'].append(long_name.capitalize())

        # Physloc must be the first atlocs
        if 'physloc' in specs and 'atlocs' in specs:
            p = specs['physloc']
            if p:
                if p in specs['atlocs']:
                    specs['atlocs'].remove(p)
                specs['atlocs'].insert(0, p)

        # Coordinates and dimensions
        if 'coords' in specs or 'dims' in specs:
            # copy location
            key = 'coords' if 'coords' in specs else 'dims'
            specs[key].setdefault('t', ['time'])
            suffixes = [('_' + s) for s in 'rftuvw']
            for axis in 'xyz':
                if axis in specs[key]:
                    if isinstance(specs[key][axis], dict):
                        for domain in specs[key][axis]:
                            specs[key][axis][domain] = cp_suffix(
                                    name, specs[key][axis][domain],
                                    suffixes=suffixes)
                    else:
                        specs[key][axis] = cp_suffix(name, specs[key][axis],
                             suffixes=suffixes)
#            for l, n in ('t', 'time'), ('y', 'lat'), ('x', 'lon'):
#                if isinstance(specs[key][l],
#                              list) and n not in specs[key][l]:
#                    specs[key][l].append(n)

        # Inherits from other specs (merge specs with dict_merge)
        if 'inherit' in specs and specs['inherit']:
            from_name = specs['inherit']
            if ":" in from_name:
                from_cat, from_name = from_name.split(':')[:2]
            else:
                from_cat = category
            assert (from_cat != category or
                    name != from_name), 'Cannot inherit cf specs from it self'
#            from_specs_container = self[from_cat]
            # to handle high level inheritance
            self._check_entry_(from_cat, from_name)
            from_specs = None
            to_scan = []
            if 'inherit' in entries:
                to_scan.append(self._inherit)
            to_scan.append(entries)

            for from_specs in to_scan:
                if from_name in from_specs:

                    # Merge
                    entries[name] = specs = dict_merge(
                            specs,
                            from_specs[from_name], cls=dict,
                            **CF_DICT_MERGE_KWARGS)

                    # Check compatibility of keys
                    for key in list(specs.keys()):
                        if key not in self._cfgspecs[category]['__many__']:
                            del specs[key]
            specs['inherit'] = None  # switch off inheritance now

#                    break

#        # Select
#        if specs['select']:
#            for key in specs['select']:
#                try:
#                    specs['select'] = eval(specs['select'])
#                except:
#                    pass

        # Standard_names in ids
        if specs['standard_name']:
            for standard_name in specs['standard_name']:
                if standard_name not in specs['name']:
                    specs['name'].append(standard_name)

        # Duplicate at other locations
        if specs['atlocs']:
            tonames = dupl_loc_specs(entries, name, specs['atlocs'])
            self._from_atlocs[category].extend(tonames)
            specs = entries[name]

    def _add_aliases_(self):
        for entries in self._cfs.values():
            for name, specs in list(entries.items()):
                for key, aliases in list(self.aliases.items()):
                    if key in specs:
                        for alias in aliases:
                            specs[alias] = specs[key]

    def get_alias_list(cls):
        al = []
        for aa in list(cls.aliases.values()):
            al.extend(aa)
        return al


class CFCatSpecs(object):
    """Base class for loading variables and coords CF specifications"""
    aliases = {'name': ['ids', 'id', 'names'],
               'standard_name': ["standard_names"],
               'long_name': ['long_names']
               }
    category = None

    def __init__(self, parent, category):
        assert category in parent
        self.parent = parent
        self.category = category

    @property
    def _dict(self):
        return self.parent._dict[self.category]

    @classmethod
    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    def __contains__(self, key):
        return key in self._dict

    def __str__(self):
        return pformat(self._dict)

    def _assert_known_(self, name):
        assert name in self, "Invalid entry name:" + name

    @property
    def names(self):
        return list(self._dict.keys())

    def items(self):
        return list(self._dict.items())

    def keys(self):
        return list(self._dict.keys())

    def register(self, name, **specs):
        """Register a new elements from its name and specs"""
        data = {self.category: {name: specs}}
        self.parent.register_from_cfg(data)

    def _update_names_(self):
        if self._names is not None:
            while self._names:
                del self._names[0]
            self._names.extend(self.names)

    def get_atts(self, name, select=None, exclude=None, ordered=True,
                 raiseerr=True, **extra):
        self._assert_known_(name)
        return cf2atts(name, select=select, exclude=exclude, ordered=ordered,
                       raiseerr=raiseerr, category=self.category)

    def get_search_specs(self, name, mode=None, raiseerr=True, ):
        self._assert_known_(name)
        return cf2search(name, mode=mode, raiseerr=raiseerr,
                         category=self.category)


#class VarSpecs(BaseSpecs):
#    category = 'variables'
#
#
#class AxisSpecs(BaseSpecs):
#    category = 'coords'

def get_loc(var, stype=None, mode='loc', default=None):
    """Get the location testing id, standard_name or long_name and other attributes

    :Params:

        - **var**: Generic id of the variable or axis OR an object CDAT object
          (variable or axis). If a CDAT object, it first searches for the
          :attr:`_vacumm_cf_location` attributes.
        - **stype**, optional: One of 'name', 'standard_name' or 'long_name'.
          If ``None``, they are all tested.
          If ``var`` is a string, location is searched interpreting
          it as of type ``'stype``.
          Else if ``var`` is an object ``stype`` specifies
          what attribute to test.
        - **mode**, optional: What you get if something is found.

            - ``"loc"``: Get the location lowercase letter (like "u").
            - ``"ext"``: Same as ``"loc"`` but searches for the
              :attr:`position` first
              if case of a CDAT variable,
              then the :attr:`_vacumm_cf_physloc` attributes
              or the "physloc" specification
              if ``var`` has an entry in :attr:`CF_VAR_SPECS`, else defaults to
              :attr:`DEFAULT_LOCATION`.
            - Else, the matched string.

    :Example:

        >>> loc = get_loc('toto_u')
        >>> suffix = get_loc('toto_at_u_location', stype='standard_name',
                             mode='full')
        >>> loc = get_loc(myvar, stype=['name','long_name'])

    :Return:

        - ``None`` if no location spec found,
        - if ``mode=="loc"``, one of possible physical
          :attr:`~vacumm.data.misc.arakawa.ARAKAWA_LOCATIONS`,
        - else, the complete matching string, like "_at_u_location".
    """
    # What to test
    if stype is None:
        stype = list(_reloc.keys())  # all by default
    stypes = stype if isinstance(stype, list) else [stype]

    # From a numeric object
    if isinstance(var, list):
        if len(var) == 0:
            return default
        var = var[0]
    if not isinstance(var, str):

        # First: _vacumm_cf_location
        if getattr(var, '_vacumm_cf_location', None):
            return var._vacumm_cf_location.lower()

        # Location from ids
        for stype in stypes:
            if stype == 'name':
                loc = get_loc(var.name, stype, mode=mode)
            elif hasattr(var, stype):
                loc = get_loc(getattr(var, stype), stype, mode)
            else:
                continue
            break
        if loc is not None:
            return loc

        # Ext mode
        if mode == 'ext':

            # Position attribute
            if hasattr(var, 'position'):
                return var.position.lower()

            # Physloc
            gid = var.name if not hasattr(
                var, '_vacumm_cf_name') else var._vacumm_cf_name
            gid = no_loc_single(gid, 'name')
            specs = get_cf_specs(gid)
            if specs:

                # From _vacumm_cf_physloc attribute
                if hasattr(var, '_vacumm_cf_physloc'):
                    return var._vacumm_cf_physloc

                # From specs
                if specs.get('physloc'):
                    return specs['physloc']

        return default

    # String
    parse_mode = "loc" if mode == 'loc' else "match"
    for stype in stypes:

        # From string
        loc = parse_loc(var, stype, mode=parse_mode)
        if loc:
            return loc

        # Extended mode
        if mode == 'ext':

            # From physloc
            if stype == 'name':
                specs = get_cf_specs(var)
                if specs and specs.get('physloc'):
                    return specs['physloc']

#            # Default location
#            return DEFAULT_LOCATION
    return default


def standard_name_to_location(standard_name):
    """Directly access that location form the standard name"""
    return parse_loc(standard_name, stype="standard_name", mode='loc')


def no_loc_single(name, stype):
    """Remove location specification

    :Params:

        - **name**: String to change
        - **stype**: One of 'name', 'standard_name' or 'long_name'.
    """
#    loc = get_loc(name, stype, mode='full')
    return parse_loc(name, stype, mode='split')['prefix']


def _loc_(loc=None):
    """Get location as a single char or empty lowercase string"""
    if loc is None:
        loc = ''
    if loc.startswith('_'):
        loc = loc[1:]
    loc = loc[:1]
    if loc not in string.ascii_letters:
        raise TypeError('Wrong location: ' + loc)
    return loc.lower()


def format_loc(prefix, stype, loc):
    """Format location string"""
    loc = _loc_(loc)
    if stype == 'standard_name':
        return prefix + '_at_%s_location' % loc
    if stype == 'long_name':
        return prefix + ' at %s location' % loc.upper()
    return prefix + '_' + loc


def change_loc_single(name, stype, loc, squeeze=None):
    """Change location specification

    :Params:

        - **name**: String to change
        - **stype**: String type: name, standard_name or long_name.
        - **loc**: Location as a letter or None.
        - **squeeze**, optional: If specified and equal to ``loc``,
        location is removed
          instead of being changed.

    :Example:

        >>> change_loc_single('usurf_u', 'name', 't')
        'usurf_t'
        >>> change_loc_single('sst_at_u_location', 'standard_name', 'v')
        'sst_at_v_location'
        >>> change_loc_single('usurf_t', 'name', 'u', squeeze='u')
        'usurf'
    """
    prefix = no_loc_single(name, stype)
    loc = _loc_(loc)
    if loc and (not squeeze or squeeze != loc):
        return format_loc(prefix, stype, loc)
    return prefix


def change_loc_specs(loc, name=None, standard_name=None, long_name=None, coords=None,
                     idim=None, jdim=None, squeeze=None, **kwargs):
    """Change location specification in id, standard names long name and coords id

    :Example:

        >>> specs = change_loc_specs('v', name=['usurf_u', 'u_u'])

    :Return: A dictionary
    """
    name = kwargs.pop('id', name)
    specs = kwargs.copy()
    if 'atlocs' in specs:
        del specs['atlocs']

    # Attributes
    for stype in 'name', 'standard_name', 'long_name', 'idim', 'jdim':
        values = locals()[stype]
        if values is None:
            continue
        if isinstance(values, list):
            values = list(values)
            tmp = [
                change_loc_single(
                    value,
                    stype,
                    loc,
                    squeeze=squeeze) for value in values]
            values = []
            for value in tmp:  # unique
                if value not in values:
                    values.append(value)

        else:
            values = change_loc_single(values, stype, loc, squeeze=squeeze)
        specs[stype] = values

    # Axes
    if coords is not None:
        coords = coords.copy()
        for l in list(coords.keys()):  # loop on x, y
            if l == 't':
                continue  # skip time
            lcoords = coords[l]
            single = not isinstance(lcoords, list)
            if single:
                lcoords = [lcoords]  # work on lists
            lcoords = [
                change_loc_single(
                    axis,
                    'name',
                    loc,
                    squeeze=DEFAULT_LOCATION) for axis in lcoords] + lcoords  # dup
            lnewcoords = []
            for axis in lcoords:  # unique
                if axis not in lnewcoords:
                    lnewcoords.append(axis)
            if single and len(lnewcoords) == 1:
                lnewcoords = lnewcoords[0]
            coords[l] = lnewcoords
        specs['coords'] = coords

    return specs


def change_loc(var, toloc, coords=True, squeeze=True):
    """Change location specifications of a variable and its coords

    It affects the name, the standard_name and the long_name when they
    are defined.


    :Example:
    >>> change_loc(sst, 'u').name
    'sst_u'
    >>> change_loc(sst, 't', squeeze=True).name
    'sst'

    """

    # Squeeze physloc
    if not isinstance(squeeze, str):
        if squeeze:
            squeeze = get_physloc(var)
        else:
            squeeze = None

    # Change attributes
    # - get
    specs = change_loc_specs(toloc, squeeze=squeeze,
                             name=var.name,
                             standard_name=getattr(var, 'standard_name', None),
                             long_name=getattr(var, 'long_name', None),
                             )
    # - set
    var.name = specs['name']
    for att in 'standard_name', 'long_name':
        if att in specs:
            setattr(var, att, specs[att])
    # - cf name
    if hasattr(var, '_vacumm_cf_name'):
        var._vacumm_cf_name = change_loc_single(var._vacumm_cf_name, 'name',
                                                toloc,
                                                squeeze=squeeze)

    # Change coords and grid attributes
    if coords and cdms2.isVariable(var) or cdms2.isGrid(var):

        # Axes
        # - usual
        for meth in 'Longitude', 'Latitude', 'Level':
            if hasattr(var, 'get' + meth):
                axismet = getattr(var, 'get' + meth)  # var.getLongitude
                if axismet() is not None:
                    change_loc(axismet(), toloc, squeeze=squeeze)
        # - 2d coords
        from .axes import isaxis
        if cdms2.isVariable(var) and isaxis(var) and var.ndim > 2:
            for i in -1, -2:
                change_loc(var.getAxis(i), toloc, squeeze=squeeze)

        # Grid
        if cdms2.isVariable(var) and var.getGrid() is not None:
            change_loc(var.getGrid(), toloc, squeeze=squeeze)

    # Reference attribute
    set_loc(var, toloc)

    return var


def set_loc(var, loc, addpos=False):
    """Define (or remove) the location of a variable with
    :attr:`_vacumm_cf_location` attribute

    It also makes sure that the :attr:`position` attribute is the same
    if present.

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
    elif hasattr(var, 'position') and var.position.lower() != loc:
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

    # Entry name in cf specs (string)
    if isinstance(var, str):
        cf_specs = get_cf_specs(var, category='variables')
        if cf_specs is not None:
            var = cf_specs[var]
        else:
            return DEFAULT_LOCATION

    # Dict of specs
    if isinstance(var, dict):
        return var.get('physloc', DEFAULT_LOCATION)

    # Using CF_VAR_SPECS
    if hasattr(var, '_vacumm_cf_name'):
        name = var._vacumm_cf_name
    else:
        name = var.name
    cf_specs = get_cf_specs(name, category='variables')
    if cf_specs is None:
        return DEFAULT_LOCATION
    return cf_specs.get('physloc', DEFAULT_LOCATION)


def squeeze_loc(var):
    """Squeeze out location specification of a variable if the location is
    physloc or :attr:`DEFAULT_LOCATION`

    :Params:

        - **var**: CDAT variable.
        - **physloc**, optional: ``physloc`` specs given by the entry in
          :attr:`CF_VAR_SPECS`
          or :attr:`DEFAULT_LOCATION`.

    :Return: The same variable
    """
    # Physloc
    physloc = get_physloc(var)

    # Remove loc only for physloc
    change_loc(var, None, coords=True, squeeze=physloc)


def squeeze_loc_single(name, stype, physloc=None):
    """Squeeze location specs if it matches physical location"""
    if physloc is None:
        cf_specs = get_cf_specs(name, 'variables')
        if cf_specs:
            physloc = cf_specs['physloc']
    if physloc is None:
        physloc = DEFAULT_LOCATION
    loc = parse_loc(name, stype, mode='loc')
    if loc == physloc:
        name = change_loc_single(name, stype, None)
    return name


def dupl_loc_specs(fromname, toloc, category=None):
    """Duplicate the specification for a variable or an axis to another
    or several locations

    The following rules apply:

        - If the original specifications are from a name
          without a specific location
          (generic), new specs (located) are appended (merged)
          with original specs.
        - If the specifications at target location already exist,
          it merges new specs with old specs.
        - Generic specification are systematically created or
          updated by merging of
          specialized ones.

    :Example:

        >>> dupl_loc_specs('corio', 'u') # Create the 'corio_u' entry in CF_VAR_SPECS and update 'corio'
        >>> dupl_loc_specs('corio_t', 'u') # Create 'corio_u' and 'corio_t'

    """
    cf_source = get_cf_specs(name=None, category=category)
    if fromname not in cf_source:
        raise KeyError('No such entry in specifications: ' + fromname)
    single = not isinstance(toloc, (list, tuple))
    if single:
        toloc = [toloc]
    tonames = []
    tomerge = []
    fromnoloc = no_loc_single(fromname, 'name') == fromname
    for loc in toloc:

        # New name (id)
        toname = change_loc_single(fromname, 'name', loc)
        tonames.append(toname)

        # New specs
        tospecs = change_loc_specs(loc, **cf_source[fromname])

        # Add a version of standard_name and long_name for T loc without
        # location spec
        if loc == 't':  # 'toto_at_t_location' -> ['toto_at_t_location','toto']
            for spname in 'standard_name', 'long_name':
                if spname not in tospecs:
                    continue
                addspecs = change_loc_specs(None, **{spname: tospecs[spname]})
                for dd in tospecs, addspecs:  # to lists
                    if not isinstance(dd[spname], list):
                        dd[spname] = [dd[spname]]
                tospecs[spname].extend(addspecs[spname])

        # For merging specs back with old specs if old has no location
#        if fromnoloc: # generic
        tomerge.append(tospecs)

        # Merge with old spec if existing
        if toname in cf_source:
            tospecs = dict_merge(tospecs, cf_source[toname], mergelists=True,
                                 cls=dict, mergedicts=True, )

#        # Default location -> add its generic version ('_t' -> + '')
#        if loc==DEFAULT_LOCATION:
#            genspecs = change_loc_specs(None, **tospecs)
#            tospecs = dict_merge(tospecs, genspecs, mergelists=True)

        # Remove atlocs attribute
        if 'atlocs' in tospecs:
            tospecs.pop('atlocs')

        # Store it
        cf_source[toname] = tospecs

    # Make a generic entry that merges all specialized ones
    if fromnoloc:
        genspecs = cf_source[fromname]  # existing generic specs
    else:
        genspecs = change_loc_specs(
            None, **cf_source[fromname])  # remove location info
    tomerge.insert(0, genspecs)
    genname = genspecs['name'][0]
    cf_source[genname] = dict_merge(*tomerge, **CF_DICT_MERGE_KWARGS)

    if single:
        return tonames[0]
    return tonames

# def specs_def_loc(all_specs, names, suffixes=ARAKAWA_SUFFIXES):
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
    if isinstance(suffixes, str):
        suffixes = [suffixes]
    m = re.match('.+(%s)$' % '|'.join(suffixes), idref)
    if m is None or re.match('.+(%s)$' % '|'.join(suffixes), id):
        return id
    return id + m.group(1)


def _get_cf_cfg_():
    from . import cf
    if not hasattr(cf, 'CF_CFG'):
        cfgm = get_cf_cfgm()
        cf.CF_CFG = cfgm.load(CF_CFG_CFGFILE)
    return cf.CF_CFG


def get_cf_cfgm():
    """Get a :class:`~vacumm.misc.config.ConfigManager` instance to manage
    coords and variables spécifications"""
    from . import cf
    if not hasattr(cf, 'CF_CFGM'):
        cf.CF_CFGM = ConfigManager(CF_GFG_INIFILE)
    return cf.CF_CFGM

def set_cf_source(cf_source):
    """Set the current :attr:`CF_SPECS` as a :class:`CFSpecs` instance"""
    assert isinstance(cf_source, CFSpecs)
    from . import cf
    cf.CF_SPECS = cf_source
    # compat only:
    cf.CF_VAR_SPECS = cf.var_specs = cf.CF_SPECS['variables']
    cf.CF_AXIS_SPECS = cf.axis_specs = cf.CF_SPECS['coords']
    cf.CF_VAR_NAMES = cf.CF_VAR_SPECS._names
    cf.CF_AXIS_NAMES = cf.CF_AXIS_SPECS._names

def get_cf_specs(name=None, category=None):
    """Get the CF specifications for a target in a category

    Parameters
    ----------
    name: str or None
        A target name like "sst". If not provided, return all specs.
    category: str or None
        Select a category with "coords" or "variables".
        If not provided, search first in "variables", then "coords".

    Return
    ------
    dict or None
        None is return if no specs are found
    """
    # Get the source of specs
    from . import cf
    if not hasattr(cf, 'CF_SPECS'):

        # Try from cache
        import pickle
        if os.path.exists(CF_SPECS_PICKLED) and (
                os.stat(CF_CFG_CFGFILE).st_mtime <
                os.stat(CF_SPECS_PICKLED).st_mtime):
            try:
                with open(CF_SPECS_PICKLED, 'rb') as f:
                    set_cf_source(pickle.load(f))
            except Exception as e:
                vcwarn('Error while load cached cf specs: '.format(e.args))

        # Compute it from scratch
        if not hasattr(cf, 'CF_SPECS'):

            # Setup
            set_cf_source(CFSpecs(CF_CFG_CFGFILE))

            # Cache it
            try:
                cachedir = os.path.dirname(CF_SPECS_PICKLED)
                if not os.path.exists(cachedir):
                    os.makedirs(cachedir)
                with open(CF_SPECS_PICKLED, 'wb') as f:
                    pickle.dump(cf.CF_SPECS, f)
            except Exception as e:
                print(e.args)
                vcwarn('Error while caching cf specs: '.format(e.args))

    cf_source = cf.CF_SPECS

    # Select categories
    if category is not None:
        if isinstance(cf_source, CFSpecs):
            cf_source = cf_source[category]
        if name is None:
            return cf_source
        toscan = [cf_source]
    else:
        if name is None:
            return cf_source
        toscan = [cf_source['variables'], cf_source['coords']]

    # Search
    for ss in toscan:
        if name in ss:
            return ss[name]


def get_cf_coord_specs(name=None):
    """Shortcut to ``get_cf_specs(name=name, category='coords')``"""
    return get_cf_specs(name=name, category='coords')


def get_cf_var_specs(name=None):
    """Shortcut to ``get_cf_specs(name=name, category='variables')``"""
    return get_cf_specs(name=name, category='variables')


#: Specifications for grid formating
GRID_SPECS = {
    '': dict(lon='lon', lat='lat', level='depth'),
    't': dict(lon='lon', lat='lat', level='depth_t'),
    'u': dict(lon='lon_u', lat='lat_u', level='depth_t'),
    'v': dict(lon='lon_v', lat='lat_v', level='depth_t'),
    'f': dict(lon='lon_f', lat='lat_f', level='depth_t'),
    'w': dict(lon='lon', lat='lat', level='depth_w'),
}
# TODO: 't' must use lon_t and lat_t
GRID_SPECS['r'] = GRID_SPECS['t']
GRID_SPECS[None] = GRID_SPECS['']
grid_specs = GRID_SPECS  # compat


def register_cf_variable(name, **specs):
    """Register a new CF variable given its generic name and its specs

    It simply calls :meth:`~CFSpecs.register` for variables

    """
    cf_specs = get_cf_specs(category='variables')
    cf_specs.register(name, **specs)
    return cf_specs[name]


def register_cf_variables_from_cfg(cfg):
    """Register a new CF variables from a config file or a dict

    It must contains a "variables" section.

    It simply calls :meth:`~CFSpecs.register_from_cfg` for variables

    """
    cf_specs = get_cf_specs(category='variables')
    cf_specs.register_from_cfg(cfg)


def register_cf_coord(name, **specs):
    """Register a new Cf coordinate given its generic name and its specs

    It simply calls :meth:`~CFSpecs.register` for coords

    """
    cf_specs = get_cf_specs(category='coords')
    cf_specs.register(name, **specs)
    return cf_specs[name]


def register_cf_coords_from_cfg(cfg):
    """Register a new CF coordinates from a config file or a dict

    It must contains a "coords" section.

    It simply calls :meth:`~CFSpecs.register_from_cfg` for coords

    """
    cf_specs = get_cf_specs(category='coords')
    cf_specs.register_from_cfg(cfg)


def cf2search(name, mode=None, raiseerr=True, category=None, **kwargs):
    """Extract specs from :attr:`CF_AXIS_SPECS` or :attr:`CFVAR_SPECS`
    to form a search dictionary

    Parameters
    ----------
    name:
        Generic name of an coord or a variable.
    mode: optional
        Search mode [default: None->``"nsa"``].
        A string containg one or more of the following letters:

        - ``"i"`` or ``"n"``: Search using names (ids).
        - ``"s"``: Search using standard_name attribute.
        - ``"a"``: Search using axis attribute.
        - ``"l"``: Search using long_name attribute.
        - ``"u"``: Search using units attribute.

        The order is important.

    Return
    ------
    An :class:`colections.OrderedDict`

    Example
    -------

    >>> cf2search('sst', mode='isu')
    {'name':['sst'], 'name':['sst'],
    'standard_names':['sea_surface_temperature'],
    'units':['degrees_celsius']}
    """
    # Get specs
    if isinstance(name, dict):
        specs = name.copy()
    else:
        specs = get_cf_specs(name, category=category)
        if specs is None:
            if raiseerr:
                raise VACUMMError("Specifications not found for "+name)
        else:
            return

    # Form search dict
    if specs['searchmode']:
        mode = specs['searchmode']
    else:
        mode = None
    if not isinstance(mode, str):
        mode = 'nsa'
    mode = mode.replace('i', 'n')
    keys = []
    for m in mode:
        for key in ['name', 'standard_name', 'axis', 'long_name', 'units']:
            if key.startswith(m):
                keys.append(key)
                break
    return OrderedDict([(k, specs[k]) for k in keys if k in specs])


_attnames_exclude = [
    'atlocs',
    'inherit',
    'coords',
    'physloc',
    'idim',
    'jdim',
    'select',
    'searchmode',
    'cmap',
    ]
_attnames_firsts = [
    'standard_name',
    'long_name',
    'units',
    'axis',
    'valid_min',
    'valid_max',
    ]


def cf2atts(name, select=None, exclude=None, ordered=True, raiseerr=True,
            category=None, **extra):
    """Extract specs from :attr:`CF_AXIS_SPECS` or :attr:`CF_VAR_SPECS` to form
    a dictionary of attributes (units and long_name)"""
    # Get specs
    if isinstance(name, dict):
        specs = name.copy()
    else:
        specs = get_cf_specs(name, category=category)
        if specs is None:
            if raiseerr:
                raise VACUMMError("Specifications not found for "+name)
        else:
            return

    # Which attributes
    atts = OrderedDict() if ordered else {}
    if exclude is None:
        exclude = []
    elif isinstance(exclude, str):
        exclude = [exclude]
    exclude.extend(_attnames_exclude)
    for key in _attnames_firsts + list(specs.keys()):

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
            if len(value) == 0:
                continue
            value = value[0]

        # Store it
        atts[key] = value
        if key == 'name':
            atts['id'] = value

    # Extra
    for att, val in list(extra.items()):
        atts[att] = val

    return atts


def format_var(var, name=None, force=True, format_coords=True, order=None,
               nodef=True, mode='warn', **kwargs):
    """Format a MV2 variable according to its generic name


    :Params:

        - **var**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Generic name of variable. It should be one of
          those listed by :attr:`CF_VAR_SPECS`. If None, it is guessed
          with :func:`match_known_cf_var`.
        - **force**, optional: Overwrite attributes in all cases.
        - **format_coords**, optional: Also format coords.
        - **nodef**, optional: Remove location specification when it refers to the
          default location (:attr:`DEFAULT_LOCATION`).
        - **mode**: "silent", "warn" or "raise".
        - Other parameters are passed as attributes, except those:

            - present in specifications to allow overriding defaults,
            - starting with 'x', 'y', or 't' which are passed
              to :func:`format_coord`.

    :Examples:

        >>> var = format_var(myarray, 'sst', valid_min=-2, valid_max=100)

    """
    # Filter keywords for coord formating
    axismeths = {'t': 'getTime', 'y': 'getLatitude', 'x': 'getLongitude'}
    kwcoords = {}
    for k in list(axismeths.keys()):
        kwcoords[k] = kwfilter(kwargs, k + '_')

    # Always a MV2 array
    if not cdms2.isVariable(var):
        var = MV2.asarray(var)

    # Check specs
    if name is None:
        name = match_known_cf_obj(var)
    if name is not None:
        specs = get_cf_specs(name)
    else:
        specs = None
    if specs is None:
        if mode == 'warn':
            warn("Can't guess cf name")
            return var
        elif mode == 'silent':
            return var
        else:
            raise KeyError("Can't get cf specs from: "+name)
    else:
        specs = specs.copy()
    isaxis_ = 'axis' in specs
    # - merge kwargs and specs
    for key, val in list(kwargs.items()):
        if val is None or key not in specs:
            continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - remove default location
    if nodef:
        refloc = specs.get('physloc', None) or DEFAULT_LOCATION
        for att in 'name', 'long_name', 'standard_name':
            if get_loc(specs[att], att) == refloc:
                specs[att] = [no_loc_single(specs[att][0], att)]
        name = specs['name'][0]
    # - id
    if ((force is True or force in [2, 'name', 'all'])
            or var.name is None or
            (isaxis_ and var.name is None)):
        var.name = name
    # - attributes
    forceatts = (force is True or force in ['atts', 'all'] or
                 (isinstance(force, int) and force > 0))
    for att, val in list(cf2atts(specs, **kwargs).items()):
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
    if format_coords:

        # Order
        order = (var.getOrder() if not isinstance(order, str)
                 else order)
        if order is not None:
            from .axes import is_valid_order
            if not is_valid_order(order):
                raise VACUMMError("Wrong cdms order type: " + order)
            if len(order) != var.ndim:
                raise VACUMMError(
                    "Cdms order should be of length %s instead of %s" %
                    (var.ndim, len(order)))

        # First check
        if 'coords' in specs:
            axspecs = specs['coords']
            formatted = []
            for key, meth in list(axismeths.items()):
                axis = getattr(var, meth)()
                if order is not None:
                    order.replace(key, '-')
                if axis is not None:
                    format_coord(axis, axspecs[key], **kwcoords[key])
                    formatted.append(key)

            # Check remaining simple coords (DOES NOT WORK FOR 2D AXES)
            if order is not None and order != '-' * len(order):
                for key in list(axismeths.keys()):
                    if key in order and key not in formatted:
                        axis = var.getAxis(order.index(key))
                        format_coord(axis, axspecs[key], **kwcoords[key])

    return var


def format_coord(coord, name=None, force=True, recreate=False,
                format_dims=True, nodef=True, mode='warn',
                **kwargs):
    """Format a coord according to its generic name


    :Params:

        - **coord**: A :mod:`numpy` or :mod:`MV2` variabe.
        - **name**: Single or list of generic coord names.
          If None, it is guessed
          with :func:`match_known_coord`.
        - **force**, optional: Overwrite attributes in all cases.
        - **nodef**, optional: Remove location specification
          when it refers to the
          default location (:attr:`DEFAULT_LOCATION`).
        - **coords2d_<param>**, optional: <param> is passed to
          :func:`~vacumm.misc.axes.create_coords2d` for 2D coords.
        - **recreate**, optional: Force recreation of coords using either

        - Other parameters are passed as attributes.

    :Examples:

        >>> coord1d = format_coord(array1d, 'lon')
        >>> coord2d = format_coord(array2d, 'lon_u', coords2d_idim = 'xi',
                                 coords2d_jdim = 'xi', recreate=True)


    """

    # Guess best generic name from a list
    if name is None:
        return coord
    if isinstance(name, (list, tuple)):
        for nm in name:
            if match_known_cf_obj(coord, nm):
                name = nm
                break
        else:
            name = name[0]

    # Check specs
    if name is None:
        name = match_known_cf_coord(coord)
    if name is not None:
        specs = get_cf_specs(name)
    else:
        specs = None
    if specs is None:
        if mode == 'warn':
            warn("Can't guess cf name")
            return coord
        elif mode == 'silent':
            return coord
        else:
            raise KeyError("Can't get cf specs from: "+name)
    else:
        specs = specs.copy()
    # - merge kwargs and specs
    for key, val in list(kwargs.items()):
        if val is None or key not in specs:
            continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]

    # Always a coord (1D or 2D)
    coord._oldid = coord.name
    kwcoords2d = kwfilter(kwargs, 'coords2d_')
    from .axes import isaxis
    if not isaxis(coord) or recreate:
        if len(coord.shape) == 1:
            coord = cdms2.createAxis(coord) #FIXME: nocdat
        else:
            xy = specs['axis'].lower()
            kwcoords2d[xy] = coord
            from .grid import create_coords2d
            coord = create_coords2d(**kwcoords2d)
            return coord
    coords2d = len(coord.shape) == 2

    # Apply specs
    # - merge kwargs and specs
    for key, val in list(kwargs.items()):
        if val is None or key not in specs:
            continue
        # Check type
        if not isinstance(val, list) and isinstance(specs[key], list):
            val = [val]
        # Set
        specs[key] = val
        del kwargs[key]
    # - remove default location
    if nodef:
        for att in 'name', 'long_name', 'standard_name':
            #            sname = stype+'s'
            if att not in specs:
                continue
            if get_loc(specs[att], att) == DEFAULT_LOCATION:
                specs[att] = [no_loc_single(specs[att][0], att)]
    # - id
    if force or coord.name is None:
        coord.name = specs['name'][0]
    # - attributes
    for att, val in list(cf2atts(
            specs, exclude=['axis'] if coords2d else None, **kwargs).items()):
        if force or not getattr(coord, att, ''):
            setattr(coord, att, val)
    # - store cf name
    coord._vacumm_cf_name = coord.name

    # Dims for 2D coords
    if coords2d and format_dims:
        format_coord(axis.getAxis(-1), specs['idim'])
        format_coord(axis.getAxis(-2), specs['jdim'])

    return axis


def format_grid(grid, pt, **kwargs):
    """Format a grid and its coords"""
    if cdms2.isVariable(grid):
        grid = grid.getGrid()
    if grid is None:
        return
    gs = GRID_SPECS[pt]
    lon = grid.getLongitude()
    lat = grid.getLatitude()
    format_coord(lon, gs['lon'], **kwargs)
    format_coord(lat, gs['lat'], **kwargs)


#: Hidden attributes of variable useful of this module
HIDDEN_CF_ATTS = [
    '_vacumm_cf_name',
    '_vacumm_cf_physloc',
    '_vacumm_cf_location']
hidden_cf_atts = HIDDEN_CF_ATTS  # compat


def filter_search(specs, searchmode=None):
    """Order and filter attribute-based searching

    Parameters
    ----------
    specs: dict
        CF searching specs for an axis or a variable
    searchmode: None or str
        Search mode to select attributes that must be used for searching

    Return
    ------
    collections.OrderedDict
        Selected searching specs in the right order

    """
    # - get order
    all_keys = ['id', 'name', 'standard_name', 'long_name', 'units', 'axis']
    all_keys0 = [key[0] for key in all_keys]
    if searchmode is None:
        if isinstance(specs, OrderedDict):
            searchmode = ''.join([key[0] for key in list(specs.keys())])
        else:
            searchmode = ''.join(all_keys0)
#    searchmode = searchmode.replace('n', 'i')
    keys = []
    for key0 in searchmode:
        key = all_keys[all_keys0.index(key0)]
        if key0 in all_keys0 and key in specs:
            keys.append(key)
    # - reorder specs
    return OrderedDict([(key_, specs[key_]) for key_ in keys])


def match_cf_atts(obj, name=None, standard_name=None,
                  long_name=None, units=None, axis=None, ignorecase=True,
                  searchmode=None, **kwargs):
    """Check if an MV2 object (typicaly from a netcdf file) matches names,
    standard_names, etc

    It first checks the standard_name, then the names (ids),
    the axis, and finally
    the long_names and units.

    Parameters
    ----------
    obj:
        A MV2 array.
    name: optional
        Name (id) of this array, wich defaults to the id attribute.
    standard_name: optional
        List of possible standard_names.
    axis: optional
        Axis type, as one of 'x, 'y', 'z', 't'.
    long_name: optional
        List of possible long_names or callable expression
        (such as regular expression method).
    units: optional
        Same as ``long_names`` but for units.

    Example
    -------

    >>> match_cf_obj(sst, standard_name='sea_surface_temperature', id=['sst'])
    >>> import re
    >>> match_cf_obj(sst, long_name=re.compile('sea surface temp').match)
    """
    # Format
    search = OrderedDict()
    name = kwargs.pop('id', name) or name
    id = name
    for key in ('standard_name', 'id', 'name', 'axis', 'long_name', 'units'):
        val = locals()[key]
        if val is None:
            continue
        search[key] = val
        if key == 'axis':
            search[key] = val if not isinstance(key, list) else val[0]
            continue
        search[key] = val

    # Order and filter the search
    search = filter_search(search, searchmode)

    # Check attributes
    return match_atts(obj, search, ignorecase)


def match_cf_from_name(obj, name, searchmode=None, category=None,
                       **kwargs):
    """Check if a variable or an axis matches generic specifications

    :Params:

        - **obj**: a numpy or cdms2 axis or variable.
        - **name**: A generic name.
        - **searchmode**, optional: Passed to
          :func:`~vacumm.misc.match_cf_obj`.
    """
    search = cf2search(name, mode=searchmode, raiseerr=False,
                       category=category)
    if search is None:
        return False
    search.update(kwargs)
    return match_cf_atts(obj, searchmode=searchmode, **search)


def match_var(obj, name, searchmode=None, **kwargs):
    """Check if an obj matches a generic variable"""
    return match_cf_from_name(obj, name, searchmode=searchmode,
                              category='variables')


def match_coord(obj, name, searchmode=None, **kwargs):
    """Check if an obj matches a generic variable"""
    return match_cf_from_name(obj, name, searchmode=searchmode,
                              category='coords')


def match_known_cf_obj(obj, searchmode=None, **kwargs):
    """Check if an object matches a known variable or axis"""
    return (match_known_cf_var(obj, searchmode=searchmode,
                               **kwargs) or
            match_known_cf_coord(obj, searchmode=searchmode,
                                **kwargs))


def match_known_cf_var(obj, searchmode=None, **kwargs):
    """Check if an object matches a known variable"""
    for name in get_cf_var_specs():
        if match_cf_from_name(obj, name, searchmode=searchmode,
                              category='variables'):
            return name
    return False


match_known_var = match_known_cf_var


def match_known_cf_coord(obj, searchmode=None, **kwargs):
    """Check if an object matches a known axis"""
    for name in get_cf_coord_specs():
        if match_cf_from_name(obj, name, searchmode=searchmode,
                              category='coords'):
            return name
    return False


match_known_coord = match_known_cf_coord


def is_known_cf_name(name, category=None):
    """Is this name has CF specs as an axis or a variable

    Parameters
    ----------
    name: str
    category: None or str
        If a str, either "variables" or "coords"
    """
    if ((category is None or category.startswith('v')) and
            is_known_cf_var_name(name)):
        return True
    if ((category is None or category.startswith('a')) and
            is_known_cf_coord_name(name)):
        return True
    return False


def is_known_cf_var_name(name):
    """Check that a name is known as CF var"""
    return name in get_cf_var_specs()


def is_known_cf_coord_name(name):
    """Check that a name is known as CF axis"""
    return name in get_cf_coord_specs()


def get_cf_standard_name(name, category=None, mode='first'):
    """Get the standard_name of a target

    Return
    ------
    str or None
        Return None of no specs found
    """
    specs = get_cf_specs(name=name, category=category)
    if specs and specs['standard_name']:
        standard_names = specs['standard_name']
        if mode == 'first':
            standard_names = standard_names[0]
        return standard_names


def get_cf_cmap(vname):
    """Get a cmap from a standard name

    cmap may be specified for a variable with the 'cmap' key of
    its entry in :attr:`CF_VAR_SPECS`.
    """
    if hasattr(vname, 'id'):
        vname = vname.name
    elif hasattr(vname, 'name'):
        vname = vname.name
    cf_specs = get_cf_var_specs(vname)
    if cf_specs:
        from .color import get_cmap
        cmap = get_cmap(cf_specs['cmap'])
        if cmap.name != cf_specs['cmap']:
            vcwarn("Can't get cmap '{}' for standard variable '{}'".format(
                cf_specs['cmap'], vname))

    return cmap


class CFContext(object):
    """Context for changing current :class:`CFSpecs` instance

    Example
    -------
    >>> with CFContext(my_cf_source) as cfc:
    ...    cmap = get_cf_cmap('temp')

    """

    def __init__(self, cf_source):
        """
        Parameters
        ----------
        cf_source: CFSpecs
            Dict with 'variables' and 'coords' keys.
        """
        assert 'variables' in cf_source and 'coords' in cf_source
        from . import cf
        self.cf_module = cf
        self.cf_source = cf_source

    def __enter__(self):
        cf = self.cf_module
        if not hasattr(cf, 'CF_SPECS'):
            get_cf_specs()
        self.old_cf_specs = cf.CF_SPECS
        self._set_cf_source_(self.cf_source)

    def __exit__(self):
        self._set_cf_source_(self.old_cf_specs)

    def _set_cf_source_(self, cf_source):
        cf = self.cf_module
        cf.CF_SPECS = self.cf_specs = self.cf_source = cf_source
        cf.CF_VAR_SPECS = cf.VAR_SPECS = cf_source['variables']
        cf.CF_AXIS_SPECS = cf.AXIS_SPECS = cf_source['coords']
