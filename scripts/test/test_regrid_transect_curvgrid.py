"""Test :func:`~vacumm.misc.regridding.transect` on a curvilinear grid"""

# Inits
ncfile = "swan.four.nc"
lon0 = -5.1
lon1 = -4.9
lat0 = 48.35
lat1 = 48.55


# Imports
from vcmq import (cdms2, data_sample, transect, curve2, add_map_lines)

# Read data
f = cdms2.open(data_sample(ncfile))
hs = f('HS', squeeze=1)
f.close()

# Compute transect
hst, lons, lats = transect(hs, (lon0,lon1), (lat0,lat1), getcoords=True)

# Plot along time
s = curve2(hst, figsize=(8,3), title='Spatial transect on curved grid', 
    show=False, top=0.85)

# Add a small map to show the transect positions
o = add_map_lines(hs, lons, lats, map_bgcolor='w', map_bbox= [.6, .2, .3, .5], map_anchor='W')

