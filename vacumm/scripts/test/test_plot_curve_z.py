"""Test :func:`~vacumm.misc.plot.curve2` with a depth axis"""

# Imports
from vcmq import *

# Init
var = MV2.arange(5.)
var.units = '${}^{\circ}C$'
var.long_name = 'Temperature'
depth = create_dep((-80,0.1,20.))
var.setAxis(0, depth)

# Plot
figfile = 'test_plot_curve_z.png'
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True)

# Unittest
result = {"files":figfile}
