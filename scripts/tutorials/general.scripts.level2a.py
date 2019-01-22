from __future__ import print_function
from vcmq import DS, map2, data_sample, mixed_layer_depth

# Parameters
ncfile = "menor.nc"
lon=(3.,4.5)
lat=(42, 42.8)

# Read temperature and compute depth directly
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', lon=lon, lat=lat)
temp = ds.get_temp(squeeze=1)
depth = ds.get_depth(squeeze=1)

# Compute MLD
mld = mixed_layer_depth(temp, depth, mode='deltatemp', deltatemp=0.1)

# Plot
map2(mld, proj='merc', figsize=(6, 4), autoresize=0, cmap='cmocean_dense',
     fill='contourf', show=False, colorbar_shrink=0.7, right=1, close=False)
