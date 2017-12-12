"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_depth` in MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, os, map2, data_sample, code_file_name

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', logger_level='critical')
depth = ds.get_depth(squeeze=True)

# Plot midlevel
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
map2(depth[depth.shape[0]/2], savefig=figfile, close=True, vmax=0, show=False,
    colorbar_shrink=0.6, title='Testing Dataset.get_depth()')

# For unittest
result = [
    ('files', figfile),
    ]
