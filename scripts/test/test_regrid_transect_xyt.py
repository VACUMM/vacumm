"""Test :func:`~vacumm.misc.grid.regridding.transect` with interpolation in both time and space"""

# Inits
ncfile = "mars2d.xyt.nc"
lon0 = -5.4
lon1 = -4.7
lat0 = 48.1
lat1 = 48.6
time0 = "2008-08-15 07:00"
time1 = "2008-08-15 16:00"


# Imports
from vcmq import (cdms2, data_sample, transect, stick2, 
    transect_specs, add_map_lines)

# Read data
f = cdms2.open(data_sample(ncfile))
u = f('u')
v = f('v')
f.close()

# Transect specs
tlons, tlats = transect_specs(u.getGrid(), lon0, lat0, lon1, lat1)
ttimes = (time0, time1)

# Compute transect
tu = transect(u, tlons, tlats, times=ttimes)
tv = transect(v, tlons, tlats, times=ttimes)

# Plot along time
s = stick2(tu, tv, figsize=(8,3), title='Space-time transect of speed',
    show=False, top=0.85, quiver_width=0.002)

# Add a small map to show the transect positions
add_map_lines(u.getGrid(), tlons, tlats, map_zoom=0.5)
