# Parameters
ncfile = "menor.nc"
lon=(3.,4.5)
lat=(42, 42.8)

# Imports
from vcmq import DS, map2, data_sample
from vacumm.diag.thermdyn import mixed_layer_depth

# Read temperature and depth
ncfile = data_sample(ncfile)
ds = DS(ncfile, 'mars', lon=lon, lat=lat)
temp = ds.get_temp(squeeze=1)
depth = ds.get_depth(squeeze=1)

# Compute MLD
mld = mixed_layer_depth(temp, depth, mode='deltatemp')

# Plot
map2(mld, proj='merc', figsize=(6, 6), autoresize=0,
    colorbar_shrink=0.7, right=1, savefigs=__file__, show=False,
    close=True)


print 'Done'
