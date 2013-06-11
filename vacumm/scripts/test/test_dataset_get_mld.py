"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_mld` with deltadens in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import *

# Read data
ds = setup_dataset('mars', data_sample(ncfile))
mld = ds.get_mld(mode='deltadens', squeeze=True)

# Plot surface density
figfile = 'test_dataset_get_mld.png'
if os.path.exists(figfile): os.remove(figfile)
map2(mld, savefig=figfile, close=True, show=False, vmax=600)

# For unittest
result = {'files':figfile}
