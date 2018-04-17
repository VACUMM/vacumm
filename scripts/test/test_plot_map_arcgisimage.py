"""Test :func:`~vacumm.misc.plot.map` with an arcgisimage as background"""

# Imports
from vcmq import map2

# Plot
map2(lon=(-7, -3), lat=(46, 49), arcgisimage='ocean',
    show=True, epsg=4326, res=None)

