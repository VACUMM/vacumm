"""Test the asvar keyword of :meth:`~vacumm.data.misc.dataset.Dataset.finalize_object` on MENOR"""

# Inits
ncfile="menor.nc"

# Imports
from vcmq import DS, data_sample

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
temp = ds.get_temp()
bathy = ds.get_bathy(asvar=temp)

# For unittest
result = {'assertTupleEqual':[bathy.shape, temp.shape]}


