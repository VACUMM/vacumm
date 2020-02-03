"""Test :func:`~vacumm.misc.plot.add_logo`"""

# Imports
from vcmq import N, P, add_logo, os, data_sample
import matplotlib.image as mpimg

# Inits
logofile = data_sample('logo_ifremer.png')
P.plot([2, 6])

# Default
add_logo(logofile, scale=1)

# Upper right / no rescale
add_logo(logofile, loc='upper right')

# Rescale
add_logo(logofile, loc='upper left', scale=2)

# Alpha
add_logo(logofile, loc='lower right', alpha=0.2)

# Save
figfile = __file__[:-2]+'png'
if os.path.exists(figfile): os.remove(figfile)
P.savefig(figfile)
