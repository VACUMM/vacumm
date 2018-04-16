"""Test method :meth:`~vacumm.data.misc.dataset.Dataset.get_grid` on MENOR"""

# Inits
ncfile="menor.nc"

# Imports
from vcmq import DS, data_sample

# Open
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')

# Read grids
grid = ds.get_grid()
grid_t = ds.get_grid_t()
grid_u = ds.get_grid_u()
grid_v = ds.get_grid_v()


# Read with selection
subgrid = ds.get_grid(lon=(4,5))
