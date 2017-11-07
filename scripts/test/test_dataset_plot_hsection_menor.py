"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_hsection` with MENOR"""

# Inits
ncfile = "menor.nc"
depth = -1000.

# Imports
from vcmq import DS, data_sample, os

# Setup dataset
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')

# Plot hsection
ds.plot_hsection('temp', depth, savefig=figfile, fill='contourf')

