"""Test :func:`~vacumm.misc.grid.regridding.grid2xy` with interpolation in X, Y, Z and T"""

# Inits
np = 100
nx = 20
ny = 15
nz = 10
nz = 5
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
from vcmq import (N, code_file_name, os, P, create_lon, create_lat, create_dep,
                  create_time)

# Rectangular xyzt with 1d z
# - data
lon = create_lon(N.linspace(lon0, lon1, nx))
lat = create_lat(N.linspace(lat0, lat1, ny))
dep = create_dep(N.linspace(dep0, dep1, nz))
time = create_time(N.linspace(time0, time1, nt))
extra = create_axis(N.arange(ne), id='member')
vi = MV2.asarray(N.arange(nx*ny*nz*nt*ne).reshape(ne, nt, nz, nt, nx),
                 axes=[extra, time, dep, lat, lon])
N.random.seed(0)
xo = N.random.uniform(lon0, lon1, np)
yo = N.random.uniform(lat0, lat1, np)
zo = N.random.uniform(dep0, dep1, np)
to = N.random.uniform(time0, time1, np)
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
