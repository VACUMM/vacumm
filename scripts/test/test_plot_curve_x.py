"""Test :func:`~vacumm.misc.plot.curve2` with a longitude axis"""

# Imports
from vcmq import *

# Init
var = MV2.arange(5.)
var.units = r'm s^{-1}'
var.long_name = 'Speed'
var.setAxis(0, create_lon(var.getAxis(0)))

# Plot
figfile = 'test_plot_curve_x.png'
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True)

# Unittest
result = {"files":figfile}
