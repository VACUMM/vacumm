# Parameters
ncfile = "menor.nc"
lon=(3.,4.5)
lat=(42, 42.8)

# Imports
from vcmq import setup_dataset, map2, data_sample

# Read temperature and depth
ncfile = data_sample(ncfile)
ds = setup_dataset('mars', ncfile) # lon=lon, lat=lat)
mld = ds.get_mld(mode='deltatemp', squeeze=1)

# Plot
map2(mld, proj='merc', figsize=(6, 6), autoresize=0, 
    colorbar_shrink=0.7, right=1, savefigs=__file__)


print 'Done'
