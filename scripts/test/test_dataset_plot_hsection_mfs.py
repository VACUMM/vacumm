"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_hsection` with MFS"""

# Inits
ncfile = "mfs.nc"
depth = -1000.

# Imports
from vcmq import DS, data_sample

# Setup dataset
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')

# Plot hsection
ds.plot_hsection('temp', depth, fill='contourf', show=False, close=True)

