"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_hsection` with MFS"""

# Imports
from vcmq import DS, data_sample

# Inits
ncfile = "mfs.nc"
depth = -1000.

# Setup dataset
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')

# Plot hsection
ds.plot_hsection('temp', depth, fill='contourf', show=False,
                 res='i',  cmap='cmocean_thermal', colorbar_shrink=.8,
                 right=1)
