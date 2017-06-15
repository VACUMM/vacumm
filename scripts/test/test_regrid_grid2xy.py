"""Test :func:`~vacumm.misc.grid.regridding.grid2xy` with interpolation in X, Y, Z and T"""

# Inits
np = 100
nx = 20
ny = 15
nz = 10
nt = 5
lon0 = -5.4
lon1 = -4.7
lat0 = 48.1
lat1 = 48.6
dep0 = -200
dep1 = 0
time0 = "2008-08-15 07:00"
time1 = "2008-08-15 16:00"
ne = 3



# Imports
from vcmq import (N, MV2, code_file_name, os, P, create_lon, create_lat, create_dep,
                  create_time, lindates, create_axis, reltime, grid2xy,
                  comptime)

# Rectangular xyzt with 1d z
# - data
lon = create_lon(N.linspace(lon0, lon1, nx))
lat = create_lat(N.linspace(lat0, lat1, ny))
dep = create_dep(N.linspace(dep0, dep1, nz))
time = create_time(lindates(time0, time1, nt))
extra = create_axis(N.arange(ne), id='member')
vi = MV2.array(N.arange(nx*ny*nz*nt*ne, dtype='d').reshape(ne, nt, nz, ny, nx),
                 axes=[extra, time, dep, lat, lon], copy=False,
                 fill_value=1e20)
N.random.seed(0)
xo = N.random.uniform(lon0, lon1, np)
yo = N.random.uniform(lat0, lat1, np)
zo = N.random.uniform(dep0, dep1, np)
to = comptime(N.random.uniform(reltime(time0, time.units).value,
                      reltime(time1, time.units).value, np),
                      time.units)
# - interpolation
vo = grid2xy(vi, xo=xo, yo=yo, zo=zo, to=to, method='linear')

P.scatter(xo, vo[0])
P.scatter(xo, vo[2])
p.show()

## Plot along time
#s = stick2(tu, tv, figsize=(8,3), title='Space-time transect of speed',
#    show=False, top=0.85, quiver_width=0.002)
#
## Add a small map to show the transect positions
#add_map_lines(u.getGrid(), tlons, tlats, map_zoom=0.5)
#
## Save
#figfile = code_file_name(ext='png')
#if os.path.exists(figfile): os.remove(figfile)
#s.savefig(figfile, pdf=True)
