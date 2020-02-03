"""Test :func:`~vacumm.misc.plot.map` with an arcgisimage as background"""

# Imports
from vcmq import map, os

# Plot
figfile = __file__[:-2]+'png'
if os.path.exists(figfile): os.remove(figfile)
map(lon=(-7, -3), lat=(46, 49), arcgisimage='ocean',
    show=True, close=True,
    epsg=4326, savefig=figfile, res=None)

