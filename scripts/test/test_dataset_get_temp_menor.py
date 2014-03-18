"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_temp` on MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, data_sample

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
temp = ds.get_temp()
lon = temp.getLongitude()

# For unittest
result = {'assertIsNotNone':lon}
