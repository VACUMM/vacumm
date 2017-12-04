#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
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


from vacumm.data import register_dataset
from vacumm.data.misc.dataset import OceanDataset, AtmosSurfaceDataset
from vacumm.misc import dict_merge
#from vacumm.data.cf import specs_def_loc


class GlobCurrent(OceanDataset,AtmosSurfaceDataset):
    """:class:`~vacumm.data.misc.dataset.Dataset` class to read the
    GlobCurrent satellite-derived currents

    Read the :class:`~vacumm.data.misc.dataset.Dataset` for more information
    """
    name = 'globcurrent'
    domain='ocean'
    description = 'The GlobCurrent currents dataset (http://www.globcurrent.org/)'

    #: Grid type
    arakawa_grid_type = 'A'

    # Local specs
    ncobj_specs = {

        'ugbt': {'search': {'id': ['eastward_geostrophic_current_velocity']}},
        'vgbt': {'search': {'id': ['northward_geostrophic_current_velocity']}},
    }

# Register the class
register_dataset(GlobCurrent, warn=False)




