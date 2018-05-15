"""Module for quickly loading VACUMM content"""
from __future__ import absolute_import
from __future__ import print_function

# Generic
from vcmqm import *

# Data

from vacumm.data import setup_dataset, DS, register_dataset


from vacumm.data.misc.dataset import (Dataset, OceanDataset, AtmosDataset,
                                      GenericDataset)



# Diag

from vacumm.diag.thermdyn import density, mixed_layer_depth

from vacumm.diag.dynamics import (barotropic_geostrophic_velocity,
                                  coriolis_parameter, kinetic_energy)

