"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_dens` in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import *

# Read data
ds = setup_dataset('mars', data_sample(ncfile))
dens = ds.get_dens(squeeze=True)

# Plot surface density
figfile = 'test_dataset_get_dens.png'
if os.path.exists(figfile): os.remove(figfile)
map2(dens[-1], savefig=figfile, close=True, vmin=1025.,show=False)

# For unittest
result = {'files':figfile}
