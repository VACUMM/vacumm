"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_hsection` with MFS"""

# Inits
ncfile = "mfs.nc"
depth = -1000.

# Imports
from vcmq import *

# Setup dataset
ds = setup_dataset('nemo', data_sample(ncfile))

# Plot hsection
figfile = 'test_dataset_plot_hsection_mfs.png'
if os.path.exists(figfile): os.remove(figfile)
ds.plot_hsection('temp', depth, savefig=figfile, show=False, close=True)

# Unittest
result = dict(files=figfile)
