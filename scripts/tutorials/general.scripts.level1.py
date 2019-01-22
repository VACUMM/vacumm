from __future__ import print_function
import cdms2
from vcmq import data_sample, mixed_layer_depth, map2, NcSigma, curv2rect

# %% Parameters
ncfile = "menor.nc"
lon = (3., 4.5)
lat = (42, 42.8)

# %% Read temperature
ncfile = data_sample(ncfile)
f = cdms2.open(ncfile)
temp = curv2rect(f('TEMP'))(lon=lon,  lat=lat,  squeeze=1)

# %% Compute depth
s = NcSigma.factory(f)
depth = s()
f.close()
depth.getLevel().units = ''  # workaround for cdat bug
depth = curv2rect(depth)(lon=lon, lat=lat, squeeze=1)

# %% Compute MLD
mld = mixed_layer_depth(temp, depth, mode='deltatemp', deltatemp=0.1)

# %% Plot
map2(mld, proj='merc', figsize=(6, 4), autoresize=0, cmap='cmocean_dense',
     fill='contourf', show=False, colorbar_shrink=0.7, right=1, close=False)
