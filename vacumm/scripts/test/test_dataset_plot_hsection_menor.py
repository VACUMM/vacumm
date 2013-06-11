"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_hsection` with MENOR"""

# Inits
ncfile = "menor.nc"
depth = -1000.

# Imports
from vcmq import *

# Setup dataset
ds = setup_dataset('mars', data_sample(ncfile))

# Plot hsection
figfile = 'test_dataset_plot_hsection_menor.png'
if os.path.exists(figfile): os.remove(figfile)
ds.plot_hsection('temp', depth, savefig=figfile, show=False, close=True)

# For unittest
result = dict(files=figfile)
