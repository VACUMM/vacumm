"""Test :func:`~vacumm.misc.plot.curve2` with a latitude axis"""

# Imports
from vcmq import MV2, curve2, create_lat

# Init
var = MV2.arange(5.)
var.units = 'mm'
var.long_name = 'Precipitation'
create_lat(var.getAxis(0))

# Plot
curve2(var, latex_units=True)

