"""Module for quickly loading vacumm.misc content"""
from __future__ import absolute_import
from __future__ import print_function



# System

import os
import sys
import shutil
import re
import glob
import gc
import warnings
from collections import OrderedDict
import argparse
from traceback import format_exc

# Numeric

import numpy
npy = N = np = numpy
testing = N.testing
assert_allclose = testing.assert_allclose
assert_raises = testing.assert_raises

# CDAT

import MV2
import cdms2
import cdtime
import cdutil
from cdutil import averager
from genutil import minmax

# Matplotlib

import matplotlib.pyplot as plt
#import pylab as P
mp = P = plt
from matplotlib import use, rc, rcdefaults, rcParams
from matplotlib.dates import (DayLocator, HourLocator, WeekdayLocator,
                              MinuteLocator, YearLocator, MonthLocator,
                              DateFormatter,  SecondLocator,
                              date2num, num2date)
from matplotlib.ticker import (AutoLocator, MaxNLocator, FixedLocator,
                               IndexLocator, LinearLocator, LogLocator,
                               MultipleLocator, AutoMinorLocator)

# VACUMM

# - root

from vacumm import (VACUMMError, VACUMMWarning, help as vchelp,
                    vacumm_warning, vacumm_warn, vcwarn, VACUMM_CFG)

from vacumm.config import *


# - misc level 0 (no dep or vacumm root)
from vacumm.misc.misc import *

from vacumm.misc.constants import *

from vacumm.misc.units import *

from vacumm.misc.file import *

from vacumm.misc.exception import *

from vacumm.misc.config import *
getspec = get_spec

from vacumm.misc.log import *

from vacumm.misc.namelist import *

# from vacumm.misc.docstrings import *

from vacumm.misc.weakrefset import *

from vacumm.misc.xml import *

from vacumm.misc.elapsed_time import *

from vacumm.misc.arakawa import *

from vacumm.misc.bases import *

from vacumm.misc.remote import *

from vacumm.misc.poly import *
create_polygons = polygons = as_polygons

# - misc level 1 (dep on misc)

from vacumm.misc.cf import *
match_obj = match_known_cf_obj
match_known_var = match_var = match_known_cf_var
match_known_axis = match_known_cf_axis


from vacumm.misc.axes import *

from vacumm.misc.geo import *

from vacumm.misc.poly import *
envelop = convex_hull

from vacumm.misc.filters import *

from vacumm.misc.color import *
discretize_cmap = discretise_cmap
pastelize = pastelise
land_color = land
ocean_color = ocean
sea_color = sea

from vacumm.misc.basemap import *
merc = get_merc

# - misc level 2  (dep on level 1)

from vacumm.misc.grid import *


# - misc level 3 (dep on levels 2 and 3)

from vacumm.misc.kriging import *
krign = krig

from vacumm.misc.atime import *
datetime_ = adatetime = datetime
datetime64_ = adatetime64 = datetime64
#del datetime, datetime64
add_time = add
is_time_interval = is_interval
tcompress = compress
round_time_interval = round_interval
round_time = round_date

from vacumm.misc.io import *

from vacumm.misc.sigma import *

from vacumm.misc.stats import *

# - misc level 4

from vacumm.misc.masking import *

from vacumm.misc.core_plot import *

from vacumm.misc.plot import *
curve = curve2
map = map2
hov = hov2
stick = stick2
section = section2
bar = bar2

from vacumm.misc.regridding import *

from vacumm.misc.sdata import *



