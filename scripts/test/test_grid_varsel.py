"""Test :func:`~vacumm.misc.grid.misc.varsel`"""

# %% Imports
from vcmq import N, MV2, create_grid, varsel, set_grid

# %% With specified grid
rgrid = create_grid((-14., -4, .25), (43., 52, .25), 'rect')
var_ = MV2.reshape(MV2.arange(2*rgrid.shape[0]*rgrid.shape[1]),
                   (2,) + rgrid.shape)
var = var_.clone()
cache = {}
rvar0 = varsel(var, grid=rgrid, lon=(-9, -6, 'co'), lat=(44.1, 48.2),
               cache=cache)
assert rvar0.shape == (2, 16, 12)
assert rvar0.getLongitude()[-1] == -6.25
assert rvar0.getLatitude()[-1] == 48
set_grid(var, rgrid)

# %% Without grid specification
rvar1 = varsel(var,  lon=(-9, -6, 'co'), lat=(44.1, 48.2))
N.testing.assert_allclose(rvar0.asma(), rvar1.asma())

# %% With cached gridded selector
assert 'gridded_selector' in cache
rvar2 = varsel(var_.clone(),  lon=(-9, -6, 'co'), lat=(44.1, 48.2),
               cache=cache)
N.testing.assert_allclose(rvar0.asma(), rvar2.asma())
