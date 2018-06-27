"""Test :func:`~vacumm.misc.grid.misc.GriddedSelector`"""

# %% Imports
from vcmq import N, MV2, create_grid, rotate_grid, GriddedSelector

# %% Rectangular grid
rgrid = create_grid((-14., -4, .25), (43., 52, .25), 'rect')
var = MV2.reshape(MV2.arange(2*rgrid.shape[0]*rgrid.shape[1]),
                  (2,) + rgrid.shape)
gs = GriddedSelector(rgrid, lon=(-9, -6, 'co'), lat=(44.1, 48.2))
rvar = gs(var)
assert rvar.shape == (2, 16, 12)
assert rvar.getLongitude()[-1] == -6.25
assert rvar.getLatitude()[-1] == 48

# %% Curvilinear
cgrid = rotate_grid(rgrid, 30.)
gs = GriddedSelector(cgrid, lon=(-9, -6, 'co'), lat=(44.1, 48.2))
cvar = gs(var)
assert cvar.shape == (2, 20, 18)
assert cvar.getLongitude().asma()[~gs.mask].min() == -9
assert cvar.getLatitude().asma()[~gs.mask].max() <= 48.2

# %% Unstructured
N.random.seed(0)
ugrid = create_grid(N.random.uniform(-14, -4, 100),
                    N.random.uniform(43, 52, 100), 'unstruct')
var = MV2.resize(ugrid.getLatitude().asma(), (2,) + ugrid.shape)
gs = GriddedSelector(ugrid, lon=(-9, -6, 'co'), lat=(44.1, 48.2))
uvar = gs(var)
assert uvar.getGrid().shape == (34,)
assert uvar.getLongitude().min() >= -9
assert uvar.getLatitude().max() < 52
