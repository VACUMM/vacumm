"""
Misc tools

Misc submodule can be used with vacumm.misc.*
instead of vacumm.misc.misc.*.
"""
from misc import *
import color
import io
import plot
import axes
import atime
import stats
import filters
import phys
import grid
import math

import os as _os, locale as _locale
_os.environ['LC_NUMERIC'] = 'en_US.UTF-8'
try:
    _locale.setlocale(_locale.LC_NUMERIC, 'en_US.UTF-8')
except:
    pass

