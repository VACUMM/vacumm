"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_variable` with level as a string"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, os, map2, data_sample

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', logger_level='critical')
tbot = ds.get_temp(level='bottom', squeeze=True)
tsurf = ds.get_temp(level='surf', squeeze=True)

# Plot bottom
kw = dict(contour=True, fill='contourf',  colorbar_shrink=0.8,
    cmap='cmocean_thermal', linewidth=.3)
m = map2(tsurf, subplot=211, show=False, figsize=(5, 7),
    title='Testing Dataset.get_temp(level="surf")',  **kw)
map2(tbot, subplot=212, m=m,
    title='Testing Dataset.get_temp(level="bottom")',
    **kw)

