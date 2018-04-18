"""Test :func:`~vacumm.misc.plot.curve2` with a depth axis"""

# Imports
from vcmq import MV2, curve2, create_dep

# Init
var = MV2.arange(5.)
var.units = '${}^{\circ}C$'
var.long_name = 'Temperature'
depth = create_dep((-80,0.1,20.))
var.setAxis(0, depth)

# Plot
curve2(var, show=False)

