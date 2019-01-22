"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_variable` with level as a string"""

# Imports
from vcmq import DS, map2, data_sample

# Inits
ncfile = "menor.nc"

# Read data
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', logger_level='critical')
tbot = ds.get_temp(level='bottom', squeeze=True)
tsurf = ds.get_temp(level='surf', squeeze=True)

# Plot bottom
kw = dict(contour=True, fill='contourf',  colorbar_shrink=0.8,
          cmap='cmocean_thermal', linewidth=.3, show=False,
          vmin=min(tbot.min(), tsurf.mean()),
          vmax=max(tbot.max(), tsurf.max()))
m = map2(tsurf, subplot=211, figsize=(5, 7),
         title='Testing Dataset.get_temp(level="surf")',  **kw)
map2(tbot, subplot=212, m=m,
     title='Testing Dataset.get_temp(level="bottom")',
     **kw)

# Checks
assert tbot.mean() < tsurf.mean()
