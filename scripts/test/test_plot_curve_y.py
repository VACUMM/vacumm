"""Test :func:`~vacumm.misc.plot.curve2` with a latitude axis"""

# Imports
from vcmq import MV2, code_file_name, os, code_file_name, curve2, create_lat

# Init
var = MV2.arange(5.)
var.units = 'mm'
var.long_name = 'Precipitation'
create_lat(var.getAxis(0))

# Plot
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True, latex_units=True)

