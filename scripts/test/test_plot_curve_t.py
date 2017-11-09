"""Test :func:`~vacumm.misc.plot.curve2` with a time axis"""

# Imports
from vcmq import MV2, curve2

# Init
var = MV2.arange(5.)
var.units = 'W'
var.long_name = 'Power'
axis = var.getAxis(0)
axis.units = 'months since 2013-01-01'
axis.long_name = 'Time'
axis.axis = 'T'


# Plot
curve2(var)

