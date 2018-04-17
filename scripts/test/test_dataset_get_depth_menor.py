"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_depth` in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, map2, data_sample

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', logger_level='critical')
depth = ds.get_depth(squeeze=True)

# Plot midlevel
map2(depth[depth.shape[0]/2], vmax=0, show=False,
    colorbar_shrink=0.6, title='Testing Dataset.get_depth()')
