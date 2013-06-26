"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_mld` with deltadens in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, os, map2, data_sample, code_base_name

# Read data
ds = DS(data_sample(ncfile), 'mars')
mld = ds.get_mld(mode='deltadens', squeeze=True)

# Plot surface density
figfile = code_base_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
map2(mld, savefig=figfile, close=True, show=False, vmax=600)

# For unittest
result = {'files':figfile}
