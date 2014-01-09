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



__author__ = 'Valerie Garnier'
__email__ = 'vgarnier@ifremer.fr'
__date__ = '2012-11-14'
__doc__ = 'NEMO model data manipulation'


from vacumm.data.misc.dataset import OceanDataset


class Nemo(OceanDataset):
    '''
    .. todo::
      - move specific code from Dataset to NEMO (variable mapping)
    '''
    
    
    ncbj_specs = {

        # depth
        'depth':{'search':{'names':['depth']}}, 
       
        # salinity
        'sal':{'search':{'names':['vosaline']}}, 
        
        # sea surface height
        'ssh':{'search':{'names':['sossheig']}}, 
    
        # temperature
        'temp':{'search':{'names':['votemper']}}, 
    
        # salinity
        'sal':{'search':{'names':['vosaline']}}, 
    
        # zonal current (3d)
        'u3d':{'search':{'names':['vozocrtx']}}, 
        
        # zonal current (barotrope)
        #'ubt':{'search':{'names': ['vosaline']}},
        
        # meridional current (3d)
        'v3d':{'search':{'names':['vomecrty']}}, 
    
        # meridional current (2d)
        #'vbt':{'search':{'names':['vomecrty']}}
        
	# -- Atmosphere --
        # net downward heat flux
	'sohefldo':{'search':{'names':['sohefldo']}},

	# cloud cover
	'soccov':{'search':{'names':['soccov']}},

	# surface heat flux: damping
	'sohefldp':{'search':{'names':['sohefldp']}},

        # specific humidity
	'sohumspe':{'search':{'names':['sohumspe']}},

        # latent downward heat flux
        'lathf':{'search':{'names':['solhflup']}},

	# Longwave Downward Hear flux
	'lwhf':{'search':{'names':['solwfldo']}},

        # Sensible Downward Heat Flux
	'senhf':{'search':{'names':['sosbhfup']}},
	
        # Shortwave Radiation
	'soshfldo':{'search':{'names':['soshfldo']}},

	# Air temperature at 2m
	'sotemair':{'search':{'names':['sotemair']}},

	# Concentration/dilution water flux
	'sowaflcd':{'search':{'names':['sowaflcd']}},

	# Surface Water Flux: Damping
	'sowafldp':{'search':{'names':['sowafldp']}},

	# Net Upward Water Flux
	'sowaflup':{'search':{'names':['sowaflup']}},

	# Total Precip
	'sowapre':{'search':{'names':['sowapre']}},

	# Wind speed module at 10m
	'sowindsp':{'search':{'names':['sowindsp']}},

	# -- Rivers --
	# River runoffs
	'sornf':{'search':{'names':['sornf']}},

        # SST
        'sst':{
            'inherit':'temp', 
            'select':{'level':slice(0, 1)}, 
            'squeeze':'z', 
        }
        ,
        # sea surface salinity
        'sss':{
            'inherit':'sal',
            'select':{'level':slice(0, 1)},
            'squeeze':'z',
        }
        ,
        # sea surface zonal velocity
        'usurf':{
            'inherit':'u3d',
            'select':{'level':slice(0, 1)},
            'squeeze':'z',
        }
        ,
        # sea surface meridional velocity
        'vsurf':{
            'inherit':'v3d',
            'select':{'level':slice(0, 1)},
            'squeeze':'z',
        }
    }
    
    
    positive = 'down'
    
if __name__ == '__main__' :
    from vacumm.data import setup_dataset
    from vacumm.misc.plot import map2 as map
    f = setup_dataset('nemo', dataset='/tmp/oo7/oo/mfstep_sys4/best_estimate/2004/INGV_MFSTEP-SYS4_20040702_ASLV.nc')
    sshmod = f.get_ssh()
    #map(sshmod[0, Ellipsis], res=None)
    f = setup_dataset('nemo', dataset='/tmp/oo7/oo/mfstep_sys4/best_estimate/2004/INGV_MFSTEP-SYS4_20040702_TEMP.nc')
    temp = f.get_temp()
    xxx
    #map(temp[0, 0, Ellipsis], res=None)  # axes verticaux inverses dans NEMO !!!!!!!
    f.get_depth()
    xxxx

