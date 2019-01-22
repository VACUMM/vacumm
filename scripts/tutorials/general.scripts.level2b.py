from __future__ import print_function
from vcmq import DS, map2, data_sample

# Parameters
ncfile = "menor.nc"
lon=(3.,4.5)
lat=(42, 42.8)

# Read temperature and compute MLD directly
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', lon=lon, lat=lat)
mld = ds.get_mld(mode='deltatemp', deltatemp=0.1, squeeze=1)

# Plot
map2(mld, proj='merc', figsize=(6, 4), autoresize=0, cmap='cmocean_dense',
     fill='contourf', show=False, colorbar_shrink=0.7, right=1, close=False)
