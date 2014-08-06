"""Test :func:`~vacumm.misc.plot.curve2` with a depth axis"""

# Imports
from vcmq import MV2, code_file_name, os, code_file_name, curve2, create_dep

# Init
var = MV2.arange(5.)
var.units = '${}^{\circ}C$'
var.long_name = 'Temperature'
depth = create_dep((-80,0.1,20.))
var.setAxis(0, depth)

# Plot
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True)

