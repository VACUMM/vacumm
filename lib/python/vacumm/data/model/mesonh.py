#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
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



__author__ = 'Valerie Garnier'
__email__ = 'vgarnier@ifremer.fr'
__date__ = '2015-07-08'
__doc__ = 'MesoNH model data manipulation'


from vacumm.data import register_dataset
from vacumm.data.misc.dataset import AtmosDataset,OceanSurfaceDataset


class MesoNH(AtmosDataset,OceanSurfaceDataset):
    """:class:`~vacumm.data.misc.dataset.Dataset` class to read the
    MesoNH ocean model outputs
    
    Read the :class:`~vacumm.data.misc.dataset.Dataset` for more information
    """

    name = 'mesonh'
    domain = 'atmos'
    description = "The MesoNH non-hydrostatic mesoscale atmospheric model of the French research community"
  

# remplacer names par id a priori ?
  
    ncobj_specs = {

        # sea surface zonal wind
        'usurf':{
            'inherit':'ua',   # relatif a la methode get_ua et non au nom de variable
            'select':{'level':slice(0, None)},
            'squeeze':'z',
        }
        ,
        # sea surface meridional wind
        'vsurf':{
            'inherit':'v3d', # heritage des specs
            'select':{'level':slice(0, None)},
            'squeeze':'z',
        }

    }
    
    positive = 'up'

    def nbcells(self):

        """
        For MesoNH, the number of cells along X and Y axes must be a multiple of 2,3,5
        """ 

        list_nb_cells = []
        for l in range(10): 
            for m in range(10):
               for n in range(10):
                   list_nb_cells.append( (2**l*3**m*5**n,l,m,n) )

        return sorted(list_nb_cells)

#: Alias for :class:`MesoNH`
mesonh = MesoNH

# Register the class
register_dataset(mesonh)


if __name__ == '__main__' :

    nb_cell = mesonh.mesonh().nb_cells()
    print nb_cell[0:100]
    from vacumm.data import DS
    from vacumm.misc.plot import map2 as map
    #map(sshmod[0, Ellipsis], res=None)
    print "test utus"
    f = DS('/temp/vgarnier/AMICO/WORKENV_METHODO/AROME/./AROME_IROISE_20110901-0000_20110910-0000.nc','mesonh')
    xwind = f.get_u10m()
    xxx
    #f.get_depth()
    xxxx

