"""Test :meth:`~vacumm.data.misc.dataset.Dataset.plot_transect` on MFS"""

# Inits
ncfile = "mfs.nc"
lons = (3.125,5.5)
lats = (43.4375,42.125)

# Imports
from vcmq import *

# Setup dataset
ds = setup_dataset('nemo', data_sample(ncfile))

# Plot transect
figfile = 'test_dataset_plot_transect_mfs.png'
if os.path.exists(figfile): os.remove(figfile)
ds.plot_transect('temp', lons, lats, outaxis='dist', show=False,
    figsize=(6,4), savefig=figfile, top=0.9, right=0.89, close=True)

# For unittest
result = dict(files=figfile)
