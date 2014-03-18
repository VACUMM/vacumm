"""Test the asvar keyword of :meth:`~vacumm.data.misc.dataset.Dataset.finalize_object` on MFS"""

# Inits
ncfile = "mfs.nc"

# Imports
from vcmq import DS, data_sample

# Read data
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')
temp = ds.get_temp()
depth = ds.get_depth(asvar=temp)

# For unittest
result = {'assertTupleEqual':[depth.shape, temp.shape]}
