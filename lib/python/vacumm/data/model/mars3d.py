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



__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__date__ = '2011-01-17'
__doc__ = 'MARS3D model data manipulation'


from vacumm.misc import dict_merge
from vacumm.data import register_dataset
from vacumm.data.misc.dataset import OceanDataset, AtmosSurfaceDataset
#from vacumm.data.cf import specs_def_loc


class Mars3D(OceanDataset,AtmosSurfaceDataset):
    """:class:`~vacumm.data.misc.dataset.Dataset` class to read the
    MARS3D ocean model (IFREMER) outputs

    Read the :class:`~vacumm.data.misc.dataset.Dataset` for more information
    """
    name = 'mars'
    domain='ocean'
    description = "The MARS3D ocean model from IFREMER"

    #: Grid type
    arakawa_grid_type = 'C'

    #: Radius of earth (m)
    earth_radius = 6.367e+6

    #: Gravity (m s-2)
    gravity = 9.81

    #: Positive up?
    positive = 'up'

    #: Name of time axis
    taxis_name = 'time'

    #: Name of vertical axis
    zaxis_name = 'levels'

    #: Name of meridional / j axis
    yaxis_name = 'nj'

    #: Name of zonal / i axis
    xaxis_name = 'ni'

    # Local specs
    ncobj_specs = {

        # diffusivity
        'kz':{'search':{'id':['kzm']}},

        # sea surface temperature
        'sst':{
            'inherit':'temp',
            'select':{'level':slice(-1, None)},
            'squeeze':'z',
        }
        ,
        # sea surface salinity
        'sss':{
            'inherit':'sal',
            'select':{'level':slice(-1, None)},
            'squeeze':'z',
        }
        ,
        # sea surface zonal velocity
        'usurf':{
            'inherit':'u3d',   # relatif a la methode get_u3d et non au nom de variable
            'select':{'level':slice(-1, None)},
            'squeeze':'z',
        }
        ,
        # sea surface meridional velocity
        'vsurf':{
            'inherit':'v3d', # heritage des specs
            'select':{'level':slice(-1, None)},
            'squeeze':'z',
        }
    }
#    ncobj_specs = dict_merge(_local_obj_specs, OceanDataset.ncobj_specs)
#    ncobj_specs = dict_merge(ncobj_specs, AtmosSurfaceDataset.ncobj_specs)
#    del _local_obj_specs


#: Alias for :class:`Mars3D`
MARS3D = Mars3D

# Register the class
register_dataset(MARS3D)

