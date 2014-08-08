"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_depth` in MFS"""

# Inits
ncfile = "mfs.nc"

# Imports
from vcmq import DS, os, map2, data_sample, code_file_name, isdep

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'nemo', logger_level='critical') 
depth = ds.get_depth(squeeze=True)

# For unittest
result = [
    (isdep, depth), 
    ('assertTrue', depth[0]<depth[-1]<0), 
]
