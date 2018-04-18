"""Test :func:`~vacumm.misc.plot.add_logo`"""

# Imports
from vcmq import P, add_logo, data_sample

# Inits
logofile = data_sample('logo_ifremer.png')
P.figure()
P.plot([2, 6])

# Default
add_logo(logofile, scale=1)

# Upper right / no rescale
add_logo(logofile, loc='upper right')

# Rescale
add_logo(logofile, loc='upper left', scale=2)

# Alpha
add_logo(logofile, loc='lower right', alpha=0.2)

