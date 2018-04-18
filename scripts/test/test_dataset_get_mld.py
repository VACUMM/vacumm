"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_mld` with deltadens in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, map2, data_sample

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
mld = ds.get_mld(mode='deltadens', deltadens=0.06, squeeze=True)

# Plot surface density
map2(mld, show=False, fill='pcolormesh', cmap='vacumm_previmer', contour=False)
