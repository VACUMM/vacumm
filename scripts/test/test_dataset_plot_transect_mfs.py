"""Test :meth:`~vacumm.data.misc.dataset.Dataset.plot_transect` on MFS"""

# Inits
ncfile = "mfs.nc"
lons = (3.125,5.5)
lats = (43.4375,42.125)

# Imports
from vcmq import DS, data_sample, os, code_file_name

# Setup dataset
ds = DS(data_sample(ncfile), 'nemo', logger_level='critical')

# Plot transect
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
ds.plot_transect('temp', lons, lats, outaxis='dist', show=False, fill='contourf',
    figsize=(6,4), savefig=figfile, top=0.9, right=0.89, close=True)

