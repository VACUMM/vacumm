"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_transect` with MENOR"""

# Inits
ncfile = "menor.nc"
lons = (3.76, 4.7)
lats = (43.1, 42)

# Imports
from vcmq import DS, data_sample

# Setup dataset
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')

# Plot transect
ds.plot_transect('temp', lons, lats, figsize=(6,4), fill='contourf', nmax=30, 
    top=0.9, right=0.89, show=False)
