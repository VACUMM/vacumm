"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_depth` in MENOR"""

# %% Imports
from vcmq import DS, map2, data_sample

# %% Inits
ncfile = "menor.nc"

# %% Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', logger_level='critical')
depth = ds.get_depth(squeeze=True)

# %% Plot midlevel
map2(depth[int(depth.shape[0]/2)], vmax=0, show=False, res='i',
     colorbar_shrink=0.6, title='Testing Dataset.get_depth()',
     cmap='cmo.deep_r')
