"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_dens` in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, os, map2, data_sample, code_file_name

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
dens = ds.get_dens(squeeze=True)

# Plot surface density
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
map2(dens[-1], savefig=figfile, close=True, vmin=1025.,show=False)

# For unittest
result = {'files':figfile}
