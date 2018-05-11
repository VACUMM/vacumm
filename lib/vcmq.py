"""Module for quickly loading VACUMM content"""
from __future__ import absolute_import
from __future__ import print_function

# Generic
from vcmqm import *

# Data

from vacumm.data import setup_dataset, DS, register_dataset

#from vacumm.data.cf import (
#    format_var, format_axis, format_grid, change_loc, get_cf_cmap,
#    set_loc, get_loc, GENERIC_NAMES, GENERIC_AXIS_NAMES, GENERIC_VAR_NAMES,
#    AXIS_SPECS, VAR_SPECS, GRID_SPECS, ARAKAWA_SUFFIXES, HIDDEN_CF_ATTS,
#    match_known_cf_var, match_known_cf_axis, match_known_cf_obj,
#    CF_AXIS_SPECS, CF_VAR_SPECS,
#    register_cf_variable, register_cf_variables_from_cfg,
#    register_cf_axis, register_cf_axes_from_cfg,
#    CF_VAR_NAMES, CF_AXIS_NAMES, is_cf_known,
#    )

#match_obj = match_known_cf_obj
#match_known_var = match_var = match_known_cf_var
#match_known_axis = match_known_cf_axis

#from vacumm.data.misc.sigma import NcSigma, sigma2depths, SigmaGeneralized, SigmaStandard

#from vacumm.data.misc.arakawa import (
#    ArakawaGrid, ARAKAWA_LOCATIONS as arakawa_locations,
#    CGrid, AGrid, ArakawaGridTransfer
#    )

from vacumm.data.misc.dataset import (
    Dataset, OceanDataset, AtmosDataset, GenericDataset)

# - diag

from vacumm.diag.thermdyn import density, mixed_layer_depth

from vacumm.diag.dynamics import (
    barotropic_geostrophic_velocity, coriolis_parameter, kinetic_energy
    )

