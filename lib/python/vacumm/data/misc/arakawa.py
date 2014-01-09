#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar (2013-2014)
# 
# raynaud@actimar.fr
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

from vacumm.misc.grid.regridding import shift2d

__all__ = ['locations', 'positions', 'ArakawaGrid', 'CGrid']

#: Strings for naming standard physical locations on grid
locations = [
    't', # Thermodynamical quantities
    'u', # Zonal momentum quantities
    'v', # Meridional momentum quantities
    'w', # Vertical momentum quantities
    'f', # Vorticity quantities
]
positions = locations # Simple alias

# TODO: Add special locations such as uw, vw and fw at the same place?


#: Grid types
grid_types = ['A', 'C']


from vacumm import VACUMMError
class ArakawaGridError(VACUMMError):
    pass

class ArakawaGrid(object):
    """Base class that provides a facotry and declares grid operations"""
    
    @staticmethod
    def factory(arg):
        """Guess the grid type and class, and instantiate it
        
        :Params:
        
            - **arg**: An explicit grid type (:attr:`grid_types`), or
              an object with a ``grid_type`` attribute.
        
        :Return: A :class:`ArakawaGrid` children object or ``None`` if grid type 
            has not been guessed.
        """
        if str(arg).upper() in grid_types:
            gt = arg
        elif hasattr(arg, 'grid_type'):
            gt = arg.grid_type
        else:
            return
        return eval(arg.upper()+'Grid')()
    
    def __getitem__(self, p):
        p = self.is_valid_loc(p)
        return getattr(self, p)
    
    @classmethod
    def is_valid_loc(cls, p):
        """Check if a location is valid"""
        p = str(p).lower()
        if p not in locations:
            raise ArakawaGridError('Bad location in Arakawa grid: %s. '
                'Please use one of: %s'%(p, locations))
        return p
    
    def are_same_locs(self, p0, p1):
        """Check if to positions are the same
        
        :Example:
        
            >>> mygrid.issame('t', 'u')
        """
        p0 = self.is_valid_loc(p0)
        p1 = self.is_valid_loc(p1)
        
        self[p0] == self[p1]
        
    def delta_loc(self, p0, p1):
        """Difference of relative locations
        
        :Return: ``dx,dy,dz``
        """
        p0 = self.is_valid_loc(p0)
        p1 = self.is_valid_loc(p1)
        
        return self[p1][0]-self[p0][0], self[p1][1]-self[p0][1], self[p1][2]-self[p0][2]
        
    def interp(self, var, p0, p1, copy=False, mode=None, **kwargs):
        """Interpolate a variable from one location to another one using 
        :func:`vacumm.misc.grid.regridding.shift2d`
        
        .. warning: It currently works only for horizontal interpolations.
        
        .. note:: It does not change the attributes.
        
        :Params:
        
            - **var**: MV2 array with a grid.
            - **p0/1**: A valid location.
            - **mode**, optional: Interpolation mode at boundaries 
              (see :func:`~vacumm.misc.grid.regridding.shift2d`).
            - **copy**, optional: Copy the variable?
            
        :Example:
        
            >>> u3d_t = CGrid().interp(u3d_u, 'u', 't', mode='extrap')
        """
        # Valid locations
        p0 = self.is_valid_loc(p0)
        p1 = self.is_valid_loc(p1)
        
        # Same location or no grid
        if self[p0]==self[p1] or var.getGrid() is None:
            if copy: var = var.clone()
            return var
            
        # Delta
        dloc = self.delta_loc(p0, p1)
        
        # Interpolate
        return shift2d(var, ishift=dloc[0], jshift=dloc[1], mode=mode, copy=copy)
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
    
