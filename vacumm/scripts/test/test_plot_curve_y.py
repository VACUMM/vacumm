"""Test :func:`~vacumm.misc.plot.curve2` with a latitude axis"""

# Imports
from vcmq import *

# Init
var = MV2.arange(5.)
var.units = 'mm'
var.long_name = 'Precipitation'
create_lat(var.getAxis(0))

# Plot
figfile = 'test_plot_curve_y.png'
if os.path.exists(figfile): os.remove(figfile)
curve2(var, savefig=figfile, show=False, close=True, latex_units=True)

# Unittest
result = {"files":figfile}
