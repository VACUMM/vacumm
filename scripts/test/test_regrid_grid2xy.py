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
data = N.resize(lat[:], (ne, nt, nz, nx, ny)) # function of y
data = N.moveaxis(data, -1, -2)
#data = N.arange(nx*ny*nz*nt*ne, dtype='d').reshape(ne, nt, nz, ny, nx)
vi = MV2.array(data,
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
assert vo.shape==(ne, np)
try:
    N.testing.assert_allclose(vo[0], yo)
except:
    pass
    xxxx
# - plot
P.figure()
P.scatter(yo, vo[0])
#P.scatter(yo, vo[2])
P.xlabel('y')
P.ylabel('vo')
P.title('grid2xy: rectangular XYZT with 1d Z = f(Y)')
P.savefig(code_file_name(ext='.rxyzt.png'))

# Rectangular xyt only
vi_xyt = vi[:, :, 0]
vo = grid2xy(vi_xyt, xo=xo, yo=yo, to=to, method='linear')
assert vo.shape==(ne, np)
N.testing.assert_allclose(vo[0], yo)
# - plot
P.figure()
P.scatter(yo, vo[0])
#P.scatter(yo, vo[2])
P.xlabel('y')
P.ylabel('vo')
P.title('grid2xy: rectangular XYT = f(Y)')
P.savefig(code_file_name(ext='rxzt.png'))

# Zo present but not requested
vo = grid2xy(vi, xo=xo, yo=yo, to=to, method='linear')
assert vo.shape==(ne, nz, np)
N.testing.assert_allclose(vo[0, 0], yo)

#P.show()
P.close()

