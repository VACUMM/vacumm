"""Test :meth:`~vacumm.data.misc.dataset.Dataset.plot_transect` on MFS"""

# Imports
from vcmq import DS, data_sample

# Inits
ncfile = "mfs.nc"
lons = (3.125, 5.5)
lats = (43.4375, 42.125)

# Setup dataset
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')

# Plot transect
ds.plot_transect('temp', lons, lats, outaxis='dist', fill='contourf',
                 figsize=(6, 4), top=0.9, right=0.89, show=False)
