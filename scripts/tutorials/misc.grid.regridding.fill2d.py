"""Fill 2D masked data and mask on land"""
from vcmq import (cdms2, MV2, plt, data_sample, map2, fill2d, masked_polygon)

# Read SST
f = cdms2.open(data_sample('mars3d.xy.nc'))
temp = f('temp', lon=(-5.5, -4), lat=(47, 49))
f.close()
temp.long_name = 'Masked'

# Create an hole
temp[15:60, 40:80] = MV2.masked

# Fill it
tempf = fill2d(temp, method='carg')
tempf.long_name = 'Filled'

# Mask on GSHHS land
tempf[:] = masked_polygon(tempf, 'i', copy=0)

# Plot
plt.rc('font', size=9)
kw = dict(vmin=9, vmax=13, proj='merc', res='i', show=False,
          nmax_levels=10, colorbar=False)
map2(temp, subplot=121, hspace=.2,
     left=.05, right=.95, figsize=(6, 4.5), **kw)
map2(tempf, subplot=122, **kw)
plt.rcdefaults()
