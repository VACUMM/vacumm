"""Module for quickly loading VACUMM content"""
from __future__ import absolute_import
from __future__ import print_function

# Generic
from vcmqm import *

# Data
from vacumm.data import *
from vacumm.data.misc.dataset import *

# Diag
from vacumm.diag.thermdyn import *
from vacumm.diag.dynamics import *
from vacumm.diag.atmos import *
from vacumm.diag.dynamics import *
from vacumm.diag.thermdyn import *
from vacumm.diag.spectrum import *

# Bathymetry and shoreline
from vacumm.bathy.bathy import *
from vacumm.bathy.shorelines import *

# Tide
from vacumm.tide.filters import *
from vacumm.tide.station_info import *
from vacumm.tide.marigraph import *
