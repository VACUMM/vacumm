"""Test :func:`~vacumm.misc.regridding.grid2xy` with interpolation in X, Y, Z and T"""

# Imports
from vcmq import (np, MV2, plt, create_lon, create_lat, create_dep,
                  create_time, lindates, create_axis, reltime, grid2xy,
                  comptime, set_grid, rotate_grid, add_grid)

# Inits
npt = 100
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
ne = 4
nez = 2


# Rectangular xyzt with 1d z data and coords
# - data
lon = create_lon(np.linspace(lon0, lon1, nx))
lat = create_lat(np.linspace(lat0, lat1, ny))
dep = create_dep(np.linspace(dep0, dep1, nz))
time = create_time(lindates(time0, time1, nt))
extra = create_axis(np.arange(ne), id='member')
data = np.resize(lat[:], (ne, nt, nz, nx, ny))  # function of y
data = np.moveaxis(data, -1, -2)
#data = np.arange(nx*ny*nz*nt*ne, dtype='d').reshape(ne, nt, nz, ny, nx)
vi = MV2.array(data, axes=[extra, time, dep, lat, lon], copy=False,
               fill_value=1e20)
np.random.seed(0)
xo = np.random.uniform(lon0, lon1, npt)
yo = np.random.uniform(lat0, lat1, npt)
zo = np.random.uniform(dep0, dep1, npt)
to = comptime(np.random.uniform(reltime(time0, time.units).value,
                                reltime(time1, time.units).value, npt),
                time.units)

# Rectangular xyzt with 1d z
vo = grid2xy(vi, xo=xo, yo=yo, zo=zo, to=to, method='linear').asma()
von = grid2xy(vi, xo=xo, yo=yo, zo=zo, to=to, method='nearest').asma()
assert vo.shape==(ne, npt)
np.testing.assert_allclose(vo[0], yo)
kwp = dict(vmin=vi.min(), vmax=vi.max())
plt.figure(figsize=(6, 3))
plt.subplot(121)
plt.scatter(xo, yo, c=vo[0],  cmap='jet', **kwp)
add_grid(vi.getGrid())
plt.title('linear4d')
plt.subplot(122)
plt.scatter(xo, yo, c=von[0], cmap='jet', **kwp)
add_grid(vi.getGrid())
plt.title('nearest4d')
plt.figtext(.5, .98, 'grid2xy in 4D', va='top', ha='center', weight='bold')
plt.tight_layout()
plt.show()

# Reversed z and y
vi_revz = vi[:, :, ::-1, ::-1, :]
vo = grid2xy(vi_revz, xo=xo, yo=yo, zo=zo, to=to, method='linear').asma()
np.testing.assert_allclose(vo[0], yo)


# Rectangular xyt only
vi_xyt = vi[:, :, 0]
vo = grid2xy(vi_xyt, xo=xo, yo=yo, to=to, method='linear').asma()
assert vo.shape == (ne, npt)
np.testing.assert_allclose(vo[0], yo)

# Rectangular xy only
vi_xy = vi[:, 0, 0]
vo = grid2xy(vi_xy, xo=xo, yo=yo, method='linear').asma()
assert vo.shape==(ne, npt)
np.testing.assert_allclose(vo[0], yo)

# Rectangular xyzt with 5d z
zi_5d = np.resize(dep[:], (nez, nt, ny, nx, nz))
zi_5d = np.moveaxis(zi_5d, -1, 2)
vo = grid2xy(vi, zi=zi_5d, xo=xo, yo=yo, zo=zo, to=to, method='linear').asma()
assert vo.shape==(ne, npt)
np.testing.assert_allclose(vo[0], yo)

# Reversed 5d z
zi_5d_rev = zi_5d[:, :, ::-1, :, :]
vo = grid2xy(vi_revz, zi=zi_5d_rev, xo=xo, yo=yo, zo=zo, to=to,
             method='linear').asma()
np.testing.assert_allclose(vo[0], yo)

# Zi present but not requested
vo = grid2xy(vi, xo=xo, yo=yo, to=to, method='linear').asma()
assert vo.shape==(ne, nz, npt)
np.testing.assert_allclose(vo[0, 0], yo)

# Zi and Ti present but not requested
vo = grid2xy(vi, xo=xo, yo=yo, method='linear').asma()
assert vo.shape==(ne, nt, nz, npt)
np.testing.assert_allclose(vo[0, 0, 0], yo)

# Curvilinear xy only
vi_xyc = vi[:, 0, 0]
gridc = rotate_grid(vi_xyc.getGrid(), 30)
set_grid(vi_xyc, gridc)
vi_xyc[:] = np.ma.resize(gridc.getLatitude()[:], vi_xyc.shape)
vo = grid2xy(vi_xyc, xo=xo, yo=yo, method='linear').asma()
assert vo.shape==(ne, npt)
np.testing.assert_allclose(vo[0], yo)

