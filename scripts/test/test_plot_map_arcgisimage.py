"""Test :func:`~vacumm.misc.plot.map` with an arcgisimage as background"""

# Imports
from vcmq import map

# Plot
map(lon=(-7, -3), lat=(46, 49), arcgisimage='natgeo',
    show=True, epsg=4326, res=None)

