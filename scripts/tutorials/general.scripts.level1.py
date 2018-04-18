from __future__ import print_function

# Parameters
ncfile = "menor.nc"
lon=(3.,4.5)
lat=(42, 42.8)

# Imports

import cdms2
from vcmq import data_sample, mixed_layer_depth, map2, NcSigma, curv2rect

# Read temperature
print('Read')
ncfile = data_sample(ncfile)
f = cdms2.open(ncfile)
temp = curv2rect(f('TEMP'))(lon=lon,  lat=lat,  squeeze=1)

# Compute depth
print('depth')
s = NcSigma.factory(f)
depth = curv2rect(s())(lon=lon,  lat=lat, squeeze=1)
f.close()

# Compute MLD
print('mld')
print(temp.shape, depth.shape)
mld = mixed_layer_depth(temp, depth, mode='deltatemp', deltatemp=0.1)

# Plot
print('plot')
map2(mld, proj='merc', figsize=(6, 4), autoresize=0,
    fill='pcolormesh',  contour=False, show=False,
    colorbar_shrink=0.7, right=1, close=False)


print('Done')
