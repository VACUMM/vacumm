"""Test :func:`~vacumm.misc.plot.curve2` with an arbitrary axis"""

# Imports
from vcmq import MV2, code_file_name, os, code_file_name, curve2

# Init
var = MV2.arange(5.)
var.units = '$Kg$'
var.long_name = 'Flour'
axis = var.getAxis(0)
axis.units = 'm'
axis.long_name = 'Distance'


# Plot
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True)

