"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_temp` on MFS"""

# Imports
from vcmq import DS, data_sample

# Inits
ncfile = "mfs.nc"

# Read data
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')
temp = ds.get_temp()
lon = temp.getLongitude()

# Checks
assert lon is not None
