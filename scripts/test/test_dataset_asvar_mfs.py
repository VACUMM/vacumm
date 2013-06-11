"""Test the asvar keyword of :meth:`~vacumm.data.misc.dataset.Dataset.finalize_variable` on MFS"""

# Inits
ncfile = "mfs.nc"

# Imports
from vcmq import *

# Read data
ds = setup_dataset('nemo', data_sample(ncfile))
temp = ds.get_temp()
depth = ds.get_depth(asvar=temp)

# For unittest
result = {'assertTupleEqual':[depth.shape, temp.shape]}
