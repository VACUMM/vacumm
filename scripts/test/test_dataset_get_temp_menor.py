"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_temp` on MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import *

# Read data
ds = setup_dataset('mars', data_sample(ncfile))
temp = ds.get_temp()
lon = temp.getLongitude()

# For unittest
result = {'assertIsNotNone':lon}
