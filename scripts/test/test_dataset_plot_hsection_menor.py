"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_hsection` with MENOR"""

# Inits
ncfile = "menor.nc"
depth = -1000.

# Imports
from vcmq import DS, data_sample, os, os, code_base_name

# Setup dataset
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')

# Plot hsection
figfile = code_base_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
ds.plot_hsection('temp', depth, savefig=figfile, show=False, close=True)

