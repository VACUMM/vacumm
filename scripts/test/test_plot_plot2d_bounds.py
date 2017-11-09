"""Test :func:`~vacumm.misc.plot.plot2d` with specified cell bounds"""

# Imports
from vcmq import MV2, N, plot2d

# Init
var = MV2.reshape(MV2.arange(4*3), (3, 4))
var.units = 'Count'
var.long_name = 'Navets'
x = var.getAxis(1)
x.units = 'km'
x.long_name = 'Distance'
y = var.getAxis(0)
y.units = 'm'
y.long_name = 'Height'
y.designateLevel()

# Bounds
y[:] = [1.5, 3.5, 4.5]
y2db = N.array([0., 3, 4, 5])

# Plot
plot2d(var, savefig=figfile, y2db=y2db, fill='pcolor',
    xmin=-.5, xmax=3.5, ymin=0, ymax=5, cmap='jet')

