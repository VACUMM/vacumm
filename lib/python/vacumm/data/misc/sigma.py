#!/usr/bin/env python
# -*- coding: utf8 -*-
#
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

__author__ = 'Stéphane Raynaud'
__email__ = 'raynaud@actimar.fr'

'''This module defines a class for all sigma coordinates systems'''

import re, math, numpy as N, cdms2, MV2
from traceback import format_exc
from warnings import warn
from cdms2.selectors import Selector
from vacumm.misc.misc import (selector2str, create_selector,
                              filter_level_selector)
import vacumm.data.cf as cf
from vacumm.misc.axes import axis_type
from vacumm.misc.grid import dz2depth as dz2depths
from vacumm.misc.io import NcFileObj, ncread_axis, ncread_var

RE_SN2LOC_SEARCH = re.compile(r'_at_([uvwtdfr])_location', re.I).search

def standard_name_to_location(standard_name):
    m = RE_SN2LOC_SEARCH(standard_name)
    if m is not None:
        return m.group(1)

class SigmaError(Exception):
    pass

class SigmaWarning(UserWarning):
    pass

class NcSigma(object):
    '''Abstract class for sigma coordinates interface to Netcdf files

    See :class:`NcSigmaStandard` and  :class:`NcSigmaGeneralized` for more information

    :Attributes:

        .. attribute:: f

            cdms2 file object.

        .. attribute:: levelvar

            Variable whose standard_name is of type "ocean_XX_coordinate"

        .. attribute:: formula_terms

            Dictionary to find formula variables

        .. attribute:: standard_names

            List of useful standard names

        .. attribute:: names

            Netcdf names of some variables

        .. attribute:: dzt

            Layer thickness at T point

        .. attribute:: dzu

            Layer thickness at U point

        .. attribute:: dzw

            Layer thickness at W point

    '''
    standard_names = dict(
        dz = cf.VAR_SPECS['dz']['standard_names'],
        dzu = cf.VAR_SPECS['dz_u']['standard_names'],
        dzv = cf.VAR_SPECS['dz_v']['standard_names'],
        dzw = cf.VAR_SPECS['dz_w']['standard_names'],
        depth = cf.VAR_SPECS['bathy']['standard_names'],
        depthu = cf.VAR_SPECS['bathy_u']['standard_names'],
        depthv = cf.VAR_SPECS['bathy_v']['standard_names'],
        eta = cf.VAR_SPECS['ssh']['standard_names'],
#        dz = "ocean_layer_thickness",
#        dzu = "ocean_layer_thickness_at_u_location",
#        dzv = "ocean_layer_thickness_at_v_location",
#        dzw = "ocean_layer_thickness_at_w_location",
#        depth = "model_sea_floor_depth_below_geoid",
#        depthu = "model_sea_floor_depth_below_geoid_at_u_location",
#        depthv = "model_sea_floor_depth_below_geoid_at_v_location",
#        eta = ["sea_surface_height_above_geoid"],
    )
    gridvars = dict(
        t = ['eta', 'depth'],
        u = ['depthu'],
        v = ['depthv'],
    )
    gridvars['all'] = gridvars['t']+gridvars['u']+gridvars['v']
    names = {}
    horiz_terms = ['depth', 'eta']
    sigma_name = None

    def __init__(self, ncfile, levelvars=None, formula_terms=None):

        # Open file
        self.nfo = NcFileObj(ncfile)
        self.f = self.nfo.f

        # Level variable
        if levelvars is None or isinstance(levelvars, list):
            self.levelvars = self.get_levelvars(ncfile, getcls=False, levelids=levelvars)
        else:
            self.levelvars = levelvars

        # Load formula terms
        self.load_formula_terms(formula_terms)

        # Load names from standard_names
        self._load_names_(self.f)

        # Load sigma related variables
        self.load_sigma()

#        # Load layer thickness names
#        self.load_thickness_names()

        # Init the cache
        self._cache = {}

    def close(self):
        """Closed opened file"""
        self.nfo.close()
    def __del__(self):
        self.close()

    @classmethod
    def factory(cls, f, **kwargs):
        """Get a sigma instance by scaning the standard name of axes and variables in a netcdf file

        :Parameters:

            - **ncfile**: A netcdf file name or descriptor.
            - Other keywords are used to initialize the output instance.

        :Return:

            - **sigma**: ``None`` or one of the following instance:

                - :class:`NcSigmaStandard`: Ocean sigma-coordinates
                - :class:`NcSigmaGeneralized`: Ocean s-coordinates
        """

        # Open file
        nfo = NcFileObj(f)
        f = nfo.f
        try:
            vars, sigcls = cls.get_levelvars(f, getcls=True, atlocs=False)
            if sigcls is not None:
                return sigcls(f, vars, **kwargs)
#            raise SigmaError('No valid coordinate system found in variables')
        finally:
            nfo.close()

    @classmethod
    def get_levelvars(cls, f, getcls=False, atlocs=True, levelids=None):
        """Find de level variables that have a standard name like ocean_XX_coordinate

        :Params:

             - **f**: cdms2 file object or netcdf file name.
             - **getcls**, optional: Get the associated sigma class too.
             - **levelids**, optional: Restrict search to theses ids.

        :Return: Variables as a dictionary whose keys at point positions
            such as 't', 'w'. If position is not properly estimated
            from the standard_name (using :func:`standard_name_to_location`),
            't' is assumed.

            - if not getcls: ``{'t': levelt, 'w':levelw, ...}``
            - else: ``{'t': levelt, 'w':levelw, ...}, sigcls``
        """
        nfo = NcFileObj(f)
        f = nfo.f
        targets = [(dimname, f.getAxis(dimname)) for dimname in f.listdimension()]
        targets += [var for var in f.variables.items()]
        levelvars = {}
        valid_sigcls = None
        for name,var in targets:
            if levelids is not None and name not in levelids: continue
            sigcls = cls.get_sigma_class(var)
            if sigcls is not None:
#                if atlocs:
                at = standard_name_to_location(var.standard_name)
                if at is None: at = 't'
                levelvars[at] = var
                valid_sigcls = sigcls
#                else:
#                    vars.append(var)
        if not levelvars:
            levelvars = None
        nfo.close()
        if getcls: return levelvars, valid_sigcls
        return levelvars

    @classmethod
    def get_sigma_class(cls, var):
        """Return the sigma class as identified in var or None"""
        # standard_name available ?
        if not hasattr(var, 'standard_name'): return

        # Loop on classes
        for sigcls in NcSigmaStandard, NcSigmaGeneralized:

            # Get standard names as a list
            sigma_standard_names = sigcls.standard_names[sigcls.sigma_name]+sigcls.standard_names[sigcls.sigma_name+'w']
            if not isinstance(sigma_standard_names, list):
                sigma_standard_names = [sigma_standard_names]

            # Check
            if var.standard_name in sigma_standard_names:
                return sigcls

    @classmethod
    def is_sigma(cls, *args, **kwargs):
        sc = cls.get_sigma_class(*args, **kwargs)
        return sc is not None and issubclass(sc, NcSigma)

    @classmethod
    def _load_names_(cls, f):
        """Store netcdf names of variables according to :attr:`standard_names` into :attr:`names`"""
        targets = [(dimname, f.getAxis(dimname)) for dimname in f.listdimension()]
        targets += [var for var in f.variables.items()]
        for name, standard_names in cls.standard_names.items():
            if not isinstance(standard_names, list):
                standard_names = [standard_names]
            for ncname,ncvar in targets:
                if hasattr(ncvar, 'standard_name') and ncvar.standard_name in standard_names:
                    cls.names.setdefault(name, ncname)
                    break

        return cls.names


#    def _get_from_standard_name_(self, standard_names, selector=None, mode="full"):
#        """Get a variable from its standard name
#
#        :Params:
#
#            - **standard_names**: A single standard_name or a list.
#            - **getid**, optional: If True, return the id instead of the variable.
#
#        :Return: The variabe if found, or None.
#        """
#        if not isinstance(standard_names, list):
#            standard_names = [standard_names]
#        targets = [(dimname, f.getAxis(dimname)) for dimname in f.listdimension()]
#        targets += [var for var in f.variables.items()]
#        for name,var in targets:
#            if hasattr(var, 'standard_name') and var.standard_name in standard_names:
#                if mode=="id" or mode==0:
#                    return var.id
#                if mode=="full" or mode>1:
#                    self._get_from_name_(var.id, selector)
#                return var


    def load_formula_terms(self, formula_terms=None, at=None, mode='noerr'):
        """Read formula terms at T and W points and store them in
        the dictionary attribute :attr:`formula_terms` whose keys are 't' and 'w'"""
        # Get formula
        if formula_terms is None:
            formula_terms = {}
            for at, levelvar in self.levelvars.items():
                try:
                    formula_terms[at] = levelvar.formula_terms
                except:
                    if mode!='noerr':
                        raise SigmaError('Error when reading formula_term attribute of netcdf variable or axis: '+levelvar.id)
        else:
            for at in 'wt':
                if at in formula_terms: break
            else:
                formula_terms = {'t':formula_terms}

        # Scan formulas
        self.formula_terms = {}
        for at in 'wt':
            self.formula_terms[at] = {}
            ft = formula_terms.get(at, None)
            if not ft: continue
            if not isinstance(ft, (list,tuple)):
                ft = re.split('[\W]+', ft)
            self.formula_terms[at] = {}
            for i in xrange(len(ft)/2):
                try:
                    self.formula_terms[at][ft[i*2]] = ft[i*2+1]
                except:
                    raise SigmaError('Error while scaning formula terms')

            # Patch names for eta and depth into self.names
            for name in self.horiz_terms:
                if name in ft:
                    self.names[name] = self.formula_terms[at][name]

    def _get_from_cache_(self, name, selector=None):
        """Get a variable either from the cache or
        by reading it with :func:`_get_from_name_`"""
        # Check the cache
        name = name.lower()
        if not hasattr(self, '_cache'): self._cache = {}
        selector = as_selector(selector)
        ss = selector2str(selector)
        if not self._cache.has_key(name):
            self._cache[name] = [None,None]
        if self._cache[name][0] is ss and self._cache[name][0] is not None:
            return self._cache[name][1]
        del self._cache[name]

        # Read it
        var = self._get_from_name_(name, selector)
        self._cache[name] = [ss, var]
        return var

    def clean_cache(self):
        if '_cache' in self: del self._cache

    def _get_ncname_(self, name):
        """Get netcdf var name from name"""
        if name is None: return
        name = name.lower()
        if not name.startswith('+') and name not in self.names: return
        if name.startswith('+'):
            ncname = name[1:]
        else:
            ncname = self.names[name]
        return ncname

    def _get_from_name_(self, name, selector=None):
        """Load axis or variable from its generic name

        .. note:: Tests against generic names are case insensitives.
            They are search for in :attr:`names`.

        :Return: None is not found, else a MV2 array
        """

        # Inits
        ncname = self._get_ncname_(name)
        if ncname is None: return
        selector = as_selector(selector)

        # Axes
        var = ncread_axis(self.f, ncname, mode=None)
        if var is not None:

            # Selector
            atype = axis_type(var, genname=True)[:3]
            for sname, sel in split_selector(selector)[1].items():
                if sname==var.id or (atype is not None and atype==sname[:3]):
                    if not isinstance(sel, slice):
                        ijk = var.mapIntervalExt(sel)
                    else:
                        ijk = sel.indices(len(var))
                    var = var.subaxis(*ijk)
            return var

        # Variables
        args = [ncname]
        if selector: args.append(selector)
        kwargs = {'mode':None}
        var = ncread_var(self.f, *args, **kwargs)
        if var is None: return
        if hasattr(self.f[ncname], '_FillValue') and cdms2.isVariable(var):
            var[:] = MV2.masked_values(var, self.f[ncname]._FillValue, copy=0)
        return var


    def load_sigma(self, selector=None, at=None):
        """Scan variables to find sigma (or s), and set :attr:`sigma` (or :attr:`s`) attributes"""
        if at is None: at = ['t', 'w']
        single = not isinstance(at, list)
        if single: at = [at]
        setattr(self, self.sigma_name, {})
        setattr(self, 'sigma', getattr(self, self.sigma_name)) # alias
        ss = getattr(self, self.sigma_name)
        for a in at:
            a = self._at_(a, focus='ver')
            if self.sigma_name not in self.formula_terms[a]: continue
            varname  = '+' + self.formula_terms[a][self.sigma_name]
            ncname = self._get_ncname_(varname)
            ids = [ncname] if ncname else []
            zsel = filter_level_selector(selector, ids=ids)
            var = self._get_from_cache_(varname, zsel)
            if var is None:
                raise SigmaError('Sigma variable not found for sigma coordinates: '+varname)
#            if hasattr(var, 'getValue'):
#                var = var.getValue()
#            else:
#                var = var.filled()
            ss[a] = var
            self.nz = var.shape[0]

        return ss if not single else ss[at[0]]

    def has_w(self):
        """Check if W sigma coordinates are available"""
        return 'w' in self.sigma

    def has_t(self):
        """Check if W sigma coordinates are available"""
        return 't' in self.sigma

    def get_depth(self, selector=None, at='t', mode='error'):
        """Read bottom depth at T, U or V points"""
        at = self._at_(at, squeezet=True)
        if at=="w": at = 't'
        var = self._get_from_cache_('depth'+at, selector)
        if var is None and mode=='error':
            raise SigmaError('Depth variable at %-points not found for %s coordinates: %s'(self.sigma_name, at, 'depth'+at))
        return var


    def get_eta(self, selector=None, mode='error'):
        """Scan variables to find eta and read it"""
        var = self._get_from_cache_('eta', selector)
        if var is None and mode=='error':
            raise SigmaError('Eta variable not found for %s coordinates'%self.sigma_name)
        return var

#    def load_thickness_names(self):
#        """Load name of thickness variables at T, U and W if availables
#
#        Thickness variables are identified according to their standard_name,
#        stored in attributes :attr:`standard_names`.
#
#        Results are stored in attributes :attr:`dzt_name`, :attr:`dzu_name` and :attr:`dzw_name`.
#
#        """
#        for at in 't', 'u', 'w':
#            if at=='t': at=''
#            standard_name = sef.standard_names["dz"+at]
#            self.names["dz"+at] = self._get_from_standard_name_(standard_name, mode='id')

    def get_dz(self, selector=None, at='t'):
        """Load thickness at T or W if availables

        :Params:

            - **selector**: Global selector.
            - **at**: 't' or 'w'

        :Return: The variable if available or None
        """
        at = self._at_(at, focus='ver', squeezet=True)
        return self._get_from_cache_("dz"+at, selector)



#    def create_depths(self, eta):
#        """Create an empty output depths variable"""
#        shape =  eta.shape[-2:]
#        if eta.getTime() is not None:
#            shape = (eta.shape[0],) + shape
#        depths = MV2.zeros(shape, eta.dtype)
#        return depths

    @staticmethod
    def _at_(at, squeezet=False, focus=None):
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

        return at


    def _get_dz2d_config_(self, at, selector=None):
        """Get what is needed to generate depths from layer thicknesses

        :Params:

            - **at**: Grid point name: T or W (only).
            - **selector**, optional: Global :class:`~cdms2.selectors.Selector`.

        :Returns: name, dz, ref
        """
        # Currently only at T and W points
        at = self._at_(at)
        if at not in 'tw': return

        # Check config already existing
        ss = selector2str(as_selector(selector))
        if not hasattr(self, '_dz2d_config'): self._dz2d_config={}
        if at in self._dz2d_config and ss in self._dz2d_config[at]:
            return self._dz2d_config[at][ss]

        # Generate config
        config = None
        if at=='t' and 'depth' in self.names and 'dzw' in self.names:

            config = ('dzw/depth',
                self._get_from_cache_('dzw', selector),
                self._get_from_cache_('depth', selector),
            )

        elif at=='w' and 'eta' in self.names and 'dz' in self.names:

            config = ('dz/eta',
                self._get_from_cache_('dz', selector),
                self._get_from_cache_('eta', selector),
            )

        if at not in self._dz2d_config:
            self._dz2d_config[at] = {}
        self._dz2d_config[at][ss] = config
        return config

    def dz_to_depths(self, selector=None, at='t', copyaxes=True, zerolid=False):
        """Get depths from layer thicknesses

        :Params:
        """
        # Get config
        at = self._at_(at, focus='z')
        config = self._get_dz2d_config_(at, selector)
        msg = "Can't compute depths from layer thicknesses at %s"%at
        if config is None:
            raise SigmaError(msg+": no config available")
        cfgname, dz, ref = config
        refloc = cfgname.split('/')[1]
        try:
            return dz2depths(dz, ref, refloc=refloc, zerolid=zerolid)
        except Exception, e:
            raise SigmaError(msg+': '+format_exc())

    def update_file(self, newncfile, close=False):
        """Change the netcdf file"""
        if close: self.close()
        self.nfo = NcFileObj(newncfile)
        self.f = self.nfo.f

class NcSigmaStandard(NcSigma):
    '''Ocean standard coordinates converter for netcdf files

    :Parameters:

        - **ncfile**: A netcdf file name or descriptor.
        - **levelvar**, optional: Variable or axis of levels.


    :Attributes:

    .. attribute:: sigma

        [array] Sigma coordinates

    .. attribute:: depth

        [array] Bottom depth

    '''
    standard_names = NcSigma.standard_names.copy()
    sigma_name = 'sigma'
    standard_names[sigma_name] = ["ocean_sigma_coordinate", "ocean_sigma_coordinate_at_t_location"]
    standard_names[sigma_name+'w'] = ["ocean_sigma_coordinate_at_w_location"]

    def __init__(self, ncfile, levelvars=None, formula_terms=None):

        # Basic initializations
        NcSigma.__init__(self, ncfile, levelvars=levelvars, formula_terms=formula_terms)

        # Ocean or standard sigma?
        for at in 't', 'w':
            if at in self.sigma:
                self.standard = self.sigma[at][:].sum()<0
                break
        else:
            self.standard = None
        self.stype = 'standard' if self.standard else 'ocean'

#    def load_sigma(self, selector=None):
#        """Scan variables to find sigma , and set :attr:`sigma` attribute
#
#        :Return: :attr:`sigma`
#        """
#        return NcSigma._load_sigma_('sigma', selector)


    def sigma_to_depths(self, selector=None, at='t', mode=None, copyaxes=True,
            eta=None, zerolid=None):
        """Get depths for the current state

        :Params:

            - **selector**: cdms selector or dict of selection specs.
            - **at**, optional: Where to compute them (T or Z).
            - **mode**, optional: Computational mode: "auto"/None, "sigma" or "dz".
            - **eta**: Sea surface elevation. A scalar, False, an array or defaults
              file's eta.

        """

        # Inits
        if mode is None: mode = 'auto'
        at = self._at_(at)
        atz = self._at_(at, focus='z')

        # From sigma
        if mode=='auto' or mode=='sigma':

            # Compute it
            try:

                # Read variables
                if eta is None:
                    eta = self.get_eta(selector, mode='noerr')
                elif eta is False:
                    eta = None
                depth = self.get_depth(selector)

                # Force sigma reload to allow level selection
                self.load_sigma(selector)

                # Compute it
                return sigma2depths(self.sigma[at], depth, eta, stype=self.stype,
                    copyaxes=copyaxes, zerolid=zerolid)

            except Exception, e:

                if isinstance(e, KeyError):
                    msg = "can't compute it at %s location"%e.message
                else:
                    msg = format_exc()
                msg = "Can't compute depths from sigma coordinates. Error: %s"%msg
                if mode=='sigma':
                    raise SigmaError(msg)
                else:
                    warn(msg, SigmaWarning)

        # From dz
        try:
            return self.dz_to_depths(selector=selector, at=at, copyaxes=copyaxes,
                zerolid=zerolid)
        except SigmaError, e:
            if mode=='auto':
                raise SigmaError("Can't compute depths from layer thicknesses or sigma coordinates")
            else:
                raise SigmaError(e.message)


#        # Init depths
#        notime = eta.getTime() is None
#        nt = 1 if notime else eta.shape[0]
#        nz = self.sigma.shape[0]
#        shape = (nt, nz) + eta.shape[-2:]
#        depths = MV2.zeros(shape, eta.dtype)
#        depths.long_name = 'Depths'
#        depths.units = 'm'
#        depths.id = 'depths'
#
#        # Loops
#        for it in xrange(nt):
#            for iz in xrange(nz):
#
#                # Common base
#                depths[it,iz] = self.sigma[iz] * (eta[it] + depth)
#
#                # Sigma type
#                if self.standard:
#                    depths[it,iz] += eta[it] # standard
#                else:
#                    depths[it,iz] -= depth # ocean
#
#        if notime:
#            depths = depths[0]
#
#        if copyaxes:
#            # Reaffect axes
#            # TODO:
#            # - does axes need to be copied ?
#            # - retreive z axis instead of creating it
#            # - factorize this copyaxes behavior
#            axes = []
#            if not notime:
#                axes.append(eta.getTime())
#            axes.append(cdms2.createAxis(range(self.sigma.shape[0]), id='level'))
#            axes.extend(eta.getAxisList()[-2:])
#            #print depths.shape, [type(a) for a in axes]
#            depths.setAxisList(axes)
#
#        self.f.close()
#        return depths

    __call__ = sigma_to_depths


class  NcSigmaGeneralized(NcSigma):
    '''Ocean s coordinates converter for netcdf files

    :Parameters:

        - **ncfile**: A netcdf file name or descriptor.
        - **levelvar**: Variable name of levels.


    :Attributes:

    .. attribute:: s

        [array] S coordinates

    .. attribute:: depth

        [array] Bottom depth

    .. attribute:: depth_c

        [array] Surface limit depth

    .. attribute:: a

         Surface control parameter

    .. attribute:: b

         Bottom control parameter

    .. attribute:: csu_standard_name

         Standard name of stretching at mid-layer

    .. attribute:: csw_standard_name

         Standard name of stretching at top-layer

    .. attribute:: csu

         Stretching at mid-layer

    .. attribute:: csw

         Stretching at top-layer



    '''
    standard_names = NcSigma.standard_names.copy()
    sigma_name = 's'
    standard_names.update(
        cst = "ocean_s_coordinate_function_at_midlayer",
        csw = ["ocean_s_coordinate_function_at_toplayer", "ocean_s_coordinate_function_at_interface"]
    )
    standard_names[sigma_name] = [
        "ocean_s_coordinate",
        "ocean_s_coordinate_g1",
        "ocean_s_coordinate_g2",
        "ocean_s_coordinate_g1_at_t_location"]
    standard_names[sigma_name+'w'] = [
        "ocean_s_coordinate_at_w_location",
        "ocean_s_coordinate_g1_at_w_location",
        "ocean_s_coordinate_g2_at_w_location"]
    horiz_terms = NcSigma.horiz_terms + ['depth_c']

    def __init__(self, ncfile, levelvars=None, formula_terms=None):

        # Basic initializations
        NcSigma.__init__(self, ncfile, levelvars=levelvars, formula_terms=formula_terms)
        self.stype = 'generalized'

        # Load control parameters
        self.load_controls()

        # Load stretching functions
        self.load_stretchings()

    def load_controls(self, at='t', mode='noerr'):
        """Load control parameters :attr:`a` and :attr:`b`"""
        at = self._at_(at, focus='ver')
        for term in 'a','b':
            if term not in self.formula_terms[at]:
                var = None
            else:
                varname  = self.formula_terms[at][term]
                if varname=='tetha': varname = 'theta' # Grrr
                var = self._get_from_name_('+'+varname)
                if var is None and mode!='noerr':
                    raise SigmaError('%s variable not found for %s coordinates: %s'%(term,self.sigma_name, varname))
            setattr(self, term, var)

    def load_stretchings(self):
        """Load stretching at mid and top layer and store them in :attr:`csu`"""
        self.cs = {}
        for at in 'w', 't':
            csname = "cs"+at
            if csname in self.names:
                var = self.f(self.names[csname])
            else:
                var =  None
            self.cs[at] = var


    def get_depth_c(self, selector, mode='error'):
        """Scan file for limit depth and read it"""
        var = self._get_from_cache_('depth_c', selector)
        if var is None and mode=='error':
            raise SigmaError('Depth_c variable not found for %s coordinates: '%(self.sigma_name, varname))
        return var

#    def load_sigma(self, selector=None):
#        """Scan variables to find s
#
#        :Return: :attr:`s`
#        """
#        return self._load_sigma_('s', selector)

    def sigma_to_depths(self, selector=None, at='t', mode=None, copyaxes=True,
            eta=None, zerolid=False):
        """Get depths for the current state

        :Params:

            - **selector**: cdms selector or dict of selection specs.
            - **at**, optional: Where to compute them (T or Z).
            - **mode**, optional: Computational mode: "auto"/None, "sigma" or "dz".
            - **eta**: Sea surface elevation. A scalar, False, an array or defaults
              file's eta.

        """

        # Inits
        if mode is None: mode = 'auto'
        at = self._at_(at)
        atz = self._at_(at, focus='z')

        # From sigma
        if mode=='auto' or mode=='sigma':

            try:

                # Read variables
#cval revoir le selector ici
                if eta is None:
                    eta = self.get_eta(selector, mode='noerr')
                elif eta is False:
                    eta = None
                depth = self.get_depth(selector)
                depth_c = self.get_depth_c(selector)

                # Force sigma reload to allow level selection
                self.load_sigma(selector)

                # Compute it
                return sigma2depths(self.s[atz], depth, eta, stype=self.stype,
                    cs=self.cs[atz], depth_c=depth_c, a=self.a, b=self.b,
                    copyaxes=copyaxes, zerolid=zerolid)

            except Exception, e:

                if isinstance(e, KeyError):
                    msg = "can't compute it at %s location"%e.message
                else:
                    msg = format_exc()
                msg = "Can't compute depths from sigma coordinates. Error: %s"%msg
                if mode=='sigma':
                    raise SigmaError(msg)
                else:
                    warn(msg, SigmaWarning)

        # From dz
        try:
            return self.dz_to_depths(selector=selector, at=at, copyaxes=copyaxes,
                zerolid=zerolid)
        except SigmaError, e:
            if mode=='auto':
                raise SigmaError("Can't compute depths from layer thicknesses or sigma coordinates")
            else:
                raise SigmaError(e.message)

    __call__ = sigma_to_depths




def _check_sigma_type_(stype, nogen=False):
    if stype is None:
        raise SigmaError('Sigma type undefined')
    stypes = ['standard', 'ocean']
    if not nogen: stypes.append('generalized')
    ss1 = [s[0] for s in stypes]
    if isinstance(stype, int):
        if stype<0 or stype>=len(stypes):
            raise SigmaError('Wrong sigma type id : %i. Should be > 0 and < %i'%(stype, len(stypes)))
    else:
        stype = str(stype)
        ss1 = [s[0] for s in stypes]
        if stype[:1] not in ss1:
           raise SigmaError('Wrong sigma type : %s. Should one of: '%(stype, stypes))
        stype = ss1.index(stype[:1])
    return stype



def sigma2depths(sigma, depth, eta=None, stype='standard',
        cs=None, depth_c=None, a=None, b=None, copyaxes=True,
        zerolid=False):
    """Conversion from standard or ocean sigma coordinates to depths

    :Params:

        - **sigma**: Sigma levels (abs(sigma)<1) as an 1D array.
        - **depth**: Bottom depth.
        - **eta**, optional: Sea surface elevation (with a time axis or not).
        - **stype**, optional: Sigma coordinates type

            - ``"standard"`` or ``0``: Standard.
            - ``"ocean"`` or ``1``: Ocean standard.
            - ``"generalized"`` or ``2``: Generalized (s) coordinates.

        - **cs**, optional: Stretching function (s coords only).
          If not provided, it is computed from stretching parameters.
        - **depth_c**, optional: Surface limit depth (s coords only).
        - **a**, optional: Surface control parameter (s coords only).
        - **b**, optional: Bottom control parameter (s coords only).
        - **zerolid**, optional: The surface is put at a zero depth to simulate
          observed depths. This makes the bottom to change with time if
          the sea level is varying.
    """

    # Init depths
    if eta is None: eta = 0.
    if not isinstance(eta, N.ndarray):
        eta = etam = N.ma.array(eta, dtype='d')
        withtime = False
    else:
        withtime = eta.getTime() is not None
        etam = eta.asma()
    nt = eta.shape[0] if withtime else 1
    nz = sigma.shape[0]
    shape = (nt, nz) + depth.shape
    depths = MV2.zeros(shape, eta.dtype)
    depths.long_name = 'Depths'
    depths.units = 'm'
    depths.id = 'depths'
    sigman = sigma.filled() if N.ma.isMA(sigma) else sigma
    etam = N.ma.atleast_1d(etam)
#    if not withtime:
#        etam = N.ma.resize(etam, (1, )+eta.shape)

    # Compute it
    stype = _check_sigma_type_(stype)
    if stype==2:
        if cs is None:

            if a is None or b is None or depth_c is None:
                raise SigmaError('You must prodive depth_c, and b '
                    'parameters for sigma generalized coordinates conversions')

            cs = ((1-b)*N.sinh(a*sigma)/math.sinh(a) +
                b * (N.tanh(a*(sigma+.5))-math.tanh(.5*a)) /
                (2*math.tanh(.5*a)))

        dd = depth-depth_c

    # Time loop
    for it in xrange(nt):
        for iz in xrange(nz):

            # Sigma generalized
            if stype == 2:

                depths[it, iz] = etam[it] * (1+sigma[iz])
                depths[it, iz] += depth_c * sigma[iz]
                depths[it, iz] += dd*cs[iz]

                if zerolid:
                    depths[it, iz] -= etam[it]

            else:
                # Common base
                depths[it, iz] = sigma[iz] * (etam[it] + depth)

                # Sigma type
                if stype!=0:
                    depths[it,iz] -= depth # ocean
                elif not zerolid:
                    depths[it,iz] += etam[it] # standard
        if not withtime:
            depths = depths[0]

    # Format axes
    if copyaxes:
        axes = []
        if withtime: # Time axis
            axes.append(eta.getTime())
        if isinstance(sigma, cdms2.axis.AbstractAxis): # Vertical axis
            axes.append(sigma)
        elif cdms2.isVariable(sigma):
            axes.append(sigma.getAxis(0))
        else:
            zaxis = depths.getAxis(int(withtime))
            zaxis.id = 'z'
            zaxis.long_name = 'Vertical levels'
            axes.append(zaxis)
        axes.extend(depth.getAxisList()) # Horizontal axes
        depths.setAxisList(axes) # Set axes
        grid = depth.getGrid() # Grid
        if grid is not None:
            depths.setGrid(grid)

    return depths


class SigmaStandard(object):
    """Standard or ocean sigma coordinates converters

    .. attribute:: sigma

        [array] Sigma coordinates

    .. attribute:: depth

        [array] Bottom depth

     """


    def __init__(self, sigma, depth, stype='standard'):
        self.stype = _check_sigma_type_(stype, nogen=True)
        self.depth = depth
        self.sigma = sigma
    def sigma_to_depths(self, eta=None, selector=None, copyaxes=True, zerolid=False):
        selector = as_selector(selector)
        depth = self.depth(selector)
        if eta is not None:
            eta = eta(selector)
        return sigma2depths(self.sigma, depth, eta, stype=self.stype,
            copyaxes=copyaxes, zerolid=zerolid)
    __call__ = sigma_to_depths

class SigmaGeneralized(object):
    """Ocean generalized (s) coordinates converters

    .. attribute:: s

        [array] S coordinates

    .. attribute:: depth

        [array] Bottom depth

    .. attribute:: depth_c

        [array] Surface limit depth

    .. attribute:: a

         Surface control parameter

    .. attribute:: b
         Bottom control parameter

     """
    def __init__(self, s, depth, depth_c, a, b):
        self.stype = 'generalized'
        self.depth = depth
        self.s = s
        self.depth_c = depth_c
        self.a = a
        self.b = b
    def sigma_to_depths(self, eta=None, selector=None, copyaxes=True, zerolid=False):
        selector = as_selector(selector)
        depth = self.depth(selector)
        depth_c = self.depth_c(selector)
        if eta is not None:
             eta = eta(selector)
        return sigma2depths(self.s, depth, eta, stype=self.stype, copyaxes=copyaxes,
            depth_c=depth_c, a=self.a, b=self.b, zerolid=zerolid)
    __call__ = sigma_to_depths

def as_selector(select=None):
    """Convert select to a :class:`cdms2.selectors.Selector` object"""
    if isinstance(select, cdms2.selectors.Selector):
        return select
    return create_selector(select)
#    selector = cdms2.selectors.Selector()
#    if select is None: return selector
#    if isinstance(select, dict):
#        selector.refine(**select)
#    elif isinstance(select, list):
#        selector.refine(*select)
#    else:
#        selector.refine(select)
#    return selector

if __name__=='__main__':
    a=cdms2.selectors.Selector((7,4),time=(5,6))
    print selector2str(a)
