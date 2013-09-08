"""Test :func:`~vacumm.misc.plot.curve2` with a longitude axis"""

# Imports
from vcmq import MV2, code_base_name, os, code_base_name, curve2, create_lon

# Init
var = MV2.arange(5.)
var.units = r'm s^{-1}'
var.long_name = 'Speed'
var.setAxis(0, create_lon(var.getAxis(0)))

# Plot
figfile = code_base_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True)

# Unittest
result = {"files":figfile}
