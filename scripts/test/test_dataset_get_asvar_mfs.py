"""Test the asvar keyword of :meth:`~vacumm.data.misc.dataset.Dataset.finalize_object` on MFS"""

# Imports
from vcmq import DS, data_sample

# Inits
ncfile = "mfs.nc"

# Read data
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')
temp = ds.get_temp()
depth = ds.get_depth(asvar=temp)

# Checks
assert depth.shape == temp.shape
