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



__author__ = 'Valerie Garnier'
__email__ = 'vgarnier@ifremer.fr'
__date__ = '2012-11-14'
__doc__ = 'NEMO model data manipulation'


from vacumm.data import register_dataset
from vacumm.data.misc.dataset import OceanDataset


class Nemo(OceanDataset):
    """:class:`~vacumm.data.misc.dataset.Dataset` class to read the
    NEMO ocean model outputs

    Read the :class:`~vacumm.data.misc.dataset.Dataset` for more information
    """
    name = 'nemo'
    domain='ocean'
    description = "The NEMO ocean model"


    ncobj_specs = {

        # depth
        'depth':{'search':{'id':['depth']}},

        # salinity
        'sal':{'search':{'id':['vosaline']}},

        # sea surface height
        'ssh':{'search':{'id':['sossheig']}},

        # temperature
        'temp':{'search':{'id':['votemper']}},

        # salinity
        'sal':{'search':{'id':['vosaline']}},

        # zonal current (3d)
        'u3d':{'search':{'id':['vozocrtx']}},

        # zonal current (barotrope)
        #'ubt':{'search':{'id': ['vosaline']}},

        # meridional current (3d)
        'v3d':{'search':{'id':['vomecrty']}},

        # meridional current (2d)
        #'vbt':{'search':{'id':['vomecrty']}}

        # -- Atmosphere --
            # net downward heat flux
        'sohefldo':{'search':{'id':['sohefldo']}},

        # cloud cover
        'soccov':{'search':{'id':['soccov']}},

        # surface heat flux: damping
        'sohefldp':{'search':{'id':['sohefldp']}},

            # specific humidity
        'sohumspe':{'search':{'id':['sohumspe']}},

        # latent downward heat flux
        'lathf':{'search':{'id':['solhflup']}},

        # Longwave Downward Hear flux
        'lwhf':{'search':{'id':['solwfldo']}},

            # Sensible Downward Heat Flux
        'senhf':{'search':{'id':['sosbhfup']}},

            # Shortwave Radiation
        'soshfldo':{'search':{'id':['soshfldo']}},

        # Air temperature at 2m
        'sotemair':{'search':{'id':['sotemair']}},

        # Concentration/dilution water flux
        'sowaflcd':{'search':{'id':['sowaflcd']}},

        # Surface Water Flux: Damping
        'sowafldp':{'search':{'id':['sowafldp']}},

        # Net Upward Water Flux
        'sowaflup':{'search':{'id':['sowaflup']}},

        # Total Precip
        'sowapre':{'search':{'id':['sowapre']}},

        # Wind speed module at 10m
        'sowindsp':{'search':{'id':['sowindsp']}},

        # -- Rivers --
        # River runoffs
        'sornf':{'search':{'id':['sornf']}},

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


#: Alias for :class:`Nemo`
NEMO = Nemo

# Register the class
register_dataset(NEMO, warn=False)




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

