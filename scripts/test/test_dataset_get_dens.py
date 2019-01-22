"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_dens` in MENOR"""

# Imports
from vcmq import DS, os, map2, data_sample

# Inits
ncfile = "menor.nc"

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
dens = ds.get_dens(squeeze=True)

# Plot surface density
map2(dens[-1], close=True, vmin=1025., show=False, res='i')
