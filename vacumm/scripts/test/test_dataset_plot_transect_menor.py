"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_transect` with MENOR"""

# Inits
ncfile = "menor.nc"
lons = (3.76, 5.2)
lats = (43.1, 42)

# Imports
from vcmq import *

# Setup dataset
ds = setup_dataset('mars', data_sample(ncfile))

# Plot transect
figfile = 'test_dataset_plot_transect_menor.png'
if os.path.exists(figfile): os.remove(figfile)
ds.plot_transect('temp', lons, lats, figsize=(6,4), 
    savefig=figfile, top=0.9, right=0.89, show=False, close=True)

# For unittest
result = dict(files=figfile)
