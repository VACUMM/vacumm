"""Test :func:`~vacumm.misc.grid.misc.rotate_grid`"""

# %% Imports
from vcmq import N, create_grid, rotate_grid

# %% Create grid
rgrid = create_grid((0., 25, 1.), (0., 10, 1.), 'rect')
cgrid = rotate_grid(rgrid, 30.)

# %% check
assert cgrid.getLongitude().ptp(axis=0).mean() == 4.5


