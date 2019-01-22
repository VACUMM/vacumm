"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_depth` in MFS"""

# Imports
from vcmq import DS, data_sample, isdep

# Inits
ncfile = "mfs.nc"

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'nemo', logger_level='critical')
depth = ds.get_depth(squeeze=True)

# Checks
assert isdep(depth)
assert depth[0] < depth[-1] < 0
