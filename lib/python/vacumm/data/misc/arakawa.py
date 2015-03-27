#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2013-2015)
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

import cdms2
from vacumm.misc.grid.regridding import shift2d, shift1d

__all__ = ['ARAKAWA_LOCATIONS', 'ARAKAWA_POSITIONS', 'ArakawaGrid', 'CGrid', 'AGrid',
    'ArakawaGridTransfer']

#: Strings for naming standard physical ARAKAWA_LOCATIONS on grid
ARAKAWA_LOCATIONS = [
    't', # Thermodynamical quantities
    'u', # Zonal momentum quantities
    'v', # Meridional momentum quantities
    'w', # Vertical momentum quantities
    'f', # Vorticity quantities
]
locations = ARAKAWA_LOCATIONS # compat

#: Alias for :attr:`ARAKAWA_LOCATIONS`
ARAKAWA_POSITIONS = ARAKAWA_LOCATIONS
positions = ARAKAWA_POSITIONS # compat

#: Upper-case version of :attr:`ARAKAWA_LOCATIONS`
ARAKAWA_LOCATIONS_UPPER = [l.upper() for l in ARAKAWA_LOCATIONS]
Locations = ARAKAWA_LOCATIONS_UPPER# compat

# TODO: Add special ARAKAWA_LOCATIONS such as uw, vw and fw at the same place?


#: Grid types
ARAKAWA_GRID_TYPES = ['A', 'C']
grid_types = ARAKAWA_GRID_TYPES # compat

from vacumm import VACUMMError
class ArakawaGridError(VACUMMError):
    pass

class _ArakawaInterp_(object):

    def _xy_interp_(self, var, p0, p1, copy, **kwargs):
        # No grid
        if var.getGrid() is None:
            if copy: var = var.clone()
            return var

        # Delta
        dloc = self.delta_loc(p0, p1)

        # Interpolate
        var = shift2d(var, ishift=dloc[0], jshift=dloc[1], copy=True, **kwargs)
        return var

    def _z_interp_(self, var, p0, p1, copy, **kwargs):
        # No vertical axis (use vacumm.misc.grid.get_zdim?)
        if var.getLevel() is None:
            if copy: var = var.clone()
            return var

        # Z axis index
        axis = var.getOrder().index('z')

        # Delta
        dloc = self.delta_loc(p0, p1)

        # Interpolate
        var = shift1d(var, shift=dloc[2], axis=axis, copy=True, **kwargs)
        return var

    @classmethod
    def is_valid_loc(cls, p):
        """Check if a location is valid"""
        p = str(p).lower()
        if p=='': return 't'
        if p not in ARAKAWA_LOCATIONS:
            raise ArakawaGridError('Bad location in Arakawa grid: %s. '
                'Please use one of: %s'%(p, ARAKAWA_LOCATIONS))
        return p

class ArakawaGrid(_ArakawaInterp_):
    """Base class that provides a facotry and declares grid operations"""

    @staticmethod
    def factory(arg):
        """Guess the grid type and class, and instantiate it

        :Params:

            - **arg**: An explicit grid type (:attr:`ARAKAWA_GRID_TYPES`), or
              an object with a either :attr:`grid_type`, :attr:`arakawa_grid_type` or
              :attr:`grid_type` attribute.

        :Return: A :class:`ArakawaGrid` children object or ``None`` if grid type
            has not been guessed.
        """
        if arg in [AGrid, CGrid]: return arg()
        if isinstance(arg, ArakawaGrid): return arg
        if str(arg).upper() in ARAKAWA_GRID_TYPES:
            gt = arg
        else:
            gt = get_grid_type(arg)
            if gt is None: return

        return eval(gt.upper()+'Grid')()

    def __getitem__(self, p):
        p = self.is_valid_loc(p)
        return getattr(self, p)

    def __str__(self):
        return self.grid_type


    def are_same_locs(self, p0, p1):
        """Check if to ARAKAWA_POSITIONS are the same

        :Example:

            >>> mygrid.are_same_locs('t', 'u')
        """
        p0 = self.is_valid_loc(p0)
        p1 = self.is_valid_loc(p1)

        self[p0] == self[p1]

    def delta_loc(self, p0, p1):
        """Difference of relative ARAKAWA_LOCATIONS

        :Return: ``dx,dy,dz``
        """
        p0 = self.is_valid_loc(p0)
        p1 = self.is_valid_loc(p1)

        return self[p1][0]-self[p0][0], self[p1][1]-self[p0][1], self[p1][2]-self[p0][2]

    def interp(self, var, p0, p1, copy=False, mode=None, zfirst=True, **kwargs):
        """Interpolate a variable from one location to another one using
        :func:`vacumm.misc.grid.regridding.shift2d` (horizontal) and
        :func:`vacumm.misc.grid.regridding.shift1d` (vertical)

        .. note:: It does not change the attributes.

        :Params:

            - **var**: MV2 array with a grid.
            - **p0/1**: Valid ARAKAWA_LOCATIONS: 0=source, 1=destination.
              If p0 is None, it is guessed with :func:`vacumm.data.cf.get_loc`.
            - **mode**, optional: Interpolation mode at boundaries
              (see :func:`~vacumm.misc.grid.regridding.shift2d`).
            - **zfirst**, optional: Perform the vertical interpolation first,
              then the horizontal interpolation.
            - **copy**, optional: Copy the variable if same location?

        :Example:

            >>> u3d_t = CGrid().interp(u3d_u, 'u', 't', mode='extrap')
        """
        # Valid ARAKAWA_LOCATIONS
        if p0 is None:
            from vacumm.data.cf import get_loc
            p0 = get_loc(var)
#            if p0 is None:
#                raise ArakawaGridError("Can't guess location of variable: "+var.id)
        p0 = self.is_valid_loc(p0)
        p1 = self.is_valid_loc(p1)

        # Same location or no grid
        if self[p0]==self[p1]:
            if copy: var = var.clone()
            return var

        # Interpolation order
        interpmets = self._z_interp_, self._xy_interp_
        if not zfirst:
            interpmets = interpmets[::-1]

        # Interpolations
        for interpmet in interpmets:
            var = interpmet(var, p0, p1, copy, mode=mode, **kwargs)

        # Special attributes
        from vacumm.data.cf import set_loc
        set_loc(var, p1)
        set_grid_type(var, self.grid_type)

        return var
    loc2loc = interp

class AGrid(ArakawaGrid):
    """A Arakawa grid"""
    # Positions relative to T point
    t = 0, 0, 0
    u = 0, 0, 0
    v = 0, 0, 0
    w = 0, 0, 0
    f = 0, 0, 0

    grid_type = gtype = 'A'

class CGrid(ArakawaGrid):
    """C Arakawa grid"""

    # Positions relative to T point
    t = 0, 0, 0
    u = 1, 0, 0
    v = 0, 1, 0
    w = 0, 0, 1
    f = 1, 1, 0

    grid_type = gtype = 'C'


class ArakawaGridTransfer(_ArakawaInterp_):
    """To interpolate variables from one grid type to another

    .. note:: This classes does not interpolate between two grids, it just interpolates
        between relative ARAKAWA_POSITIONS. For general interpolations, please use
        :func:`~vacumm.grid.regridding.regrid2d`.

    :Example:

        >>> sst_a_v = ArakawaGridTransfer('C','A').interp(sst_c_t, 't', 'v')

    """
    def __init__(self, grid0, grid1):
        """

        :Params:

            - **grid0/1**: Arguments to :meth:`~Arakawa.factory`
        """

        # Get ArakawaGrid objects
        self.grid0 = ArakawaGrid.factory(grid0)
        self.grid1 = ArakawaGrid.factory(grid1)

    def delta_loc(self, p0, p1=None):
        """Difference of relative ARAKAWA_LOCATIONS

        :Return: ``dx,dy,dz``
        """
        p0 = self.grid0.is_valid_loc(p0)
        if p1 is None:
            p1 = p0
        else:
            p1 = self.grid0.is_valid_loc(p1)

        return (self.grid1[p1][0]-self.grid0[p0][0], self.grid1[p1][1]-self.grid0[p0][1],
            self.grid1[p1][2]-self.grid0[p0][2])


    def interp(self, var, p0=None, p1=None, copy=False, mode=None, zfirst=True, **kwargs):
        """Interpolate a variable from one location to another one using
        :func:`vacumm.misc.grid.regridding.shift2d` (horizontal) and
        :func:`vacumm.misc.grid.regridding.shift1d` (vertical)

        .. note:: It does not change the attributes.

        :Params:

            - **var**: MV2 array with a grid.
            - **p0/1**: Valid ARAKAWA_LOCATIONS: 0=source, 1=destination.
              If p0 is None, it is guessed with :func:`vacumm.data.cf.get_loc`.
              If p1 is None, it defaults to p0.
            - **mode**, optional: Interpolation mode at boundaries
              (see :func:`~vacumm.misc.grid.regridding.shift2d`).
            - **copy**, optional: Copy the variable if same location?

        :Example:

            >>> sst_a_v = ArakawaGridTransfer('C','A').interp(sst_c_t, 't', 'v')

        """

        # Valid location
        if p0 is None:
            from vacumm.data.cf import get_loc
            p0 = get_loc(var)
#            if p0 is None:
#                raise ArakawaGridError("Can't guess location of variable: "+var.id)
        p0 = self.is_valid_loc(p0)
        if p1 is None:
            p1 = p0
        else:
            p1 = self.grid0.is_valid_loc(p1)

        # Same location or no grid
        if self.grid0[p0]==self.grid1[p1]:
            if copy: var = var.clone()
            return var

        # Interpolation order
        interpmets = self._z_interp_, self._xy_interp_
        if not zfirst:
            interpmets = interpmets[::-1]

        # Interpolations
        for interpmet in interpmets:
            var = interpmet(var, p0, p1, copy, mode=mode, **kwargs)

        # Special attributes
        from vacumm.data.cf import set_loc
        set_loc(var, p1)
        set_grid_type(var, self.grid1.grid_type)

        return var

_cdms2_atts = ['_vacumm_arakawa_grid_type', '_arakawa_grid_type']
_other_atts = ['arakawa_grid_type', 'grid_type']

def get_grid_type(var):
    """Guess the Arakawa grid type

    It search for the following attributes: :attr:`arakawa_grid_type`, :attr:`grid_type`
    and :attr:`_vacumm_arakawa_grid_type`.

    :Params:

        - **var**: A :mod:`cdms2` variable or grid, an :class:`ArakawaGrid` instance or
          a :class:`~vacumm.data.misc.dataset.Dataset` instance.
          If var is a :mod:`cdms2` variable, it also check its grid if defined.

    :Return: An Arakawa grid upper-base letter, like 'C'
    """
    vv = [var]
    if cdms2.isVariable(var):
        grid = var.getGrid()
        if grid is not None:
            vv.append(grid)
    for v in vv:
        for att in _cdms2_atts+_other_atts:
            if hasattr(v, att):
                gt = getattr(a, att)
                if gt is None: return
                return str(gt).upper()

def _set_clean_atts_(var, atts, value):
    for att in atts:
        if hasattr(var, att): delattr(var, att)
    if value is not None:
        setattr(var, atts[0], value)

def set_grid_type(var, gtype):
    """Set an attribute so that var is identified as being on the specified Arakawa
    grid type.

    If var is a :mod:`cdms2` variable or grid, it sets the
    :attr:`_vacumm_arakawa_grid_type` attribute,
    else it sets the :attr:`arakawa_grid_type` attribute.

    :Params:

        - **var**:  A :mod:`cdms2` variable or grid, a
          :class:`~vacumm.data.misc.dataset.Dataset` instance.
        - **gtype**: None or one of the :attr:`grid_type` letters.

    """
    if gtype is not None:
        gtype = str(gtype).upper()
    if cdms2.isVariable(var) or cdms2.isGrid(var):
        vv = [var]
        if cdms2.isVariable(var):
            grid = var.getGrid()
            if grid is not None: vv.append(grid)
        for v in vv:
            _set_clean_atts_(v, _cdms2_atts, gtype)
    else:
        _set_clean_atts_(var, _other_atts, gtype)
    return gtype
