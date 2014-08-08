"""Test :meth:`~vacumm.data.misc.OceanDataset.plot_transect` with MENOR"""

# Inits
ncfile = "menor.nc"
lons = (3.76, 4.7)
lats = (43.1, 42)

# Imports
from vcmq import DS, data_sample, os, code_file_name

# Setup dataset
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')

# Plot transect
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
ds.plot_transect('temp', lons, lats, figsize=(6,4), fill='contourf', nmax=30, 
    savefig=figfile, top=0.9, right=0.89, show=False, close=True)
