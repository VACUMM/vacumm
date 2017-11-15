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
__doc__ = 'WW3 model data manipulation'


from vacumm.data import register_dataset
from vacumm.data.misc.dataset import WaveSurfaceDataset


class WW3(WaveSurfaceDataset):
    """:class:`~vacumm.data.misc.dataset.Dataset` class to read the
    WW3 ocean model outputs
    
    Read the :class:`~vacumm.data.misc.dataset.Dataset` for more information
    """

    name = 'ww3'
    domain='wave'
    description = "The WW3 wave model from NOAA/NCEP and IFREMER"
    
    ncobj_specs = {

        # diffusivity
    #    'bathy':{'search':{'names':['dpt']}},


# remplacer names par id a priori ?
        # depth
#        'depth':{'search':{'id':['depth_t']}}, 
       
    }
    
    
    positive = 'up'

# A revoir avec un object
#    imin=1
#    imax=-1
#   jmin=1
#   jmax=-1

#: Alias for :class:`WW3`
ww3 = WW3

# Register the class
register_dataset(ww3)

if __name__ == '__main__' :
    from vacumm.data import DS
    from vacumm.misc.plot import map2 as map
    #map(sshmod[0, Ellipsis], res=None)
    print "test utus"
    f = DS('/home1/caparmor/vgarnier/ww3.190012.nc','ww3')
    stokes = f.get_utus()
    xxx
    #f.get_depth()
    xxxx

