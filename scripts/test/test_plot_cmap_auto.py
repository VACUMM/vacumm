"""Test :func:`~vacumm.misc.plot.plot2d` with automatic colormaps"""

# Imports
from vcmq import MV2, code_file_name, plot2d, os

# Init
var = MV2.reshape(MV2.arange(10*8), (8,10))
var.id = 'navets'
var.getAxis(0).id = 'frequency'
var.getAxis(1).id = 'distance'

# Plot
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
# - positive
plot2d(var, subplot=326, show=False, figsize=(6, 8), cmap='auto',
    title='Positive?')
# - negative
var[:] *= -1
plot2d(var, subplot=325, show=False, cmap='auto', title='Negative?')
# - force normal
var[:] *= -1
plot2d(var, subplot=322, show=False, cmap='auto', levels_mode='normal',
    title='Force normal')
# - force cmap
plot2d(var, subplot=323, show=False, cmap='turbid', title='Force cmap: turbid')
# - normal
var[:] += 500
plot2d(var, subplot=321, show=False, cmap='auto', title='Normal?')
# - symetric
var[:] -= var.mean()
plot2d(var, subplot=324, savefig=figfile, show=True, cmap='auto',
    tight_layout=True, close=True, title='Symetric?')


