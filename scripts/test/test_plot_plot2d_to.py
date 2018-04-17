"""Test :func:`~vacumm.misc.plot.plot2d` with a time axis"""

# Imports
from vcmq import MV2, plot2d

# Init
var = MV2.reshape(MV2.arange(10*8), (8,10))
var.units = 'Count'
var.long_name = 'Carottes'
x = var.getAxis(1)
x.units = 'months since 2000'
x.long_name = 'Time'
y = var.getAxis(0)
y.units = 'Hz'
y.long_name = 'Frequency'

# Plot
plot2d(var, fill='pcolor', order='t-')
