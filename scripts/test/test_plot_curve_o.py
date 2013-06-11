"""Test :func:`~vacumm.misc.plot.curve2` with an arbitrary axis"""

# Imports
from vcmq import *

# Init
var = MV2.arange(5.)
var.units = '$Kg$'
var.long_name = 'Flour'
axis = var.getAxis(0)
axis.units = 'm'
axis.long_name = 'Distance'


# Plot
figfile = code_base_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True)

# Unittest
result = {"files":figfile}
