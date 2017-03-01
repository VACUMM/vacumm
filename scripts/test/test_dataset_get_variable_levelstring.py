"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_variable` with level as a string"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, os, map2, data_sample, code_file_name

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', logger_level='critical')
tbot = ds.get_temp(level='bottom', squeeze=True)
tsurf = ds.get_temp(level='surf', squeeze=True)

# Plot bottom
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
kw = dict(contour=True, fill='contourf',  colorbar_shrink=0.8,
    cmap='cmocean_thermal', linewidth=.3)
m = map2(tsurf, subplot=211, close=False, show=False, figsize=(5, 7),
    title='Testing Dataset.get_temp(level="surf")',  **kw)
map2(tbot, subplot=212, close=True, show=True, m=m,
    title='Testing Dataset.get_temp(level="bottom")',
    savefig=figfile, **kw)

