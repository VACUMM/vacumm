"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_mld` with deltadens in MENOR"""

# Imports
from vcmq import DS, map2, data_sample

# Inits
ncfile = "menor.nc"

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
mld = ds.get_mld(mode='deltadens', deltadens=0.06, squeeze=True)

# Plot surface density
map2(mld, show=True, fill='pcolormesh', cmap='cmocean_dense',
     contour=False)
