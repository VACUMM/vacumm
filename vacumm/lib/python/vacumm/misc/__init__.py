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
#import archive
import phys
#import axml
import grid
import math
#try:
#    import easyPypar
#except:
#    pass
#import elapsed_time
import os as _os, locale as _locale
_os.environ['LC_NUMERIC'] = 'en_US.UTF-8'
try:
    _locale.setlocale(_locale.LC_NUMERIC, 'en_US.UTF-8')
except:
    pass

#try:
#   if pf: print 'code 1'
#   import code
#   if pf: print 'code 0'
#except:
#   pass
#
