from __future__ import print_function
from matplotlib import rc;rc('font', size=11)
import cdms2, MV2
from vcmq import (data_sample, curv_grid, set_grid, rotate_grid,
    map2, add_grid, regrid2d)
import pylab as P

# Read data
zone = dict(yc=slice(0, 40), xc=slice(10, 40))
f = cdms2.open(data_sample('swan.four.nc'))
lon2d = f('longitude', **zone)
lat2d = f('latitude', **zone)
hs = MV2.masked_values(f('HS', squeeze=1, **zone), f['HS']._FillValue)
dlon = float(f.longitude_resolution)
dlat = float(f.latitude_resolution)
f.close()

# Assign the curvilinear grid to the variable
cgridi = curv_grid(lon2d, lat2d)
hs = set_grid(hs, cgridi)

# Rotate the grid
cgrido = rotate_grid(hs.getGrid(), -30)


print('Regridding')
print(' - with SCRIP/conservative')
hs_scripcons = regrid2d(hs, cgrido, 'conservative')
print(' - with SCRIP/bilineare')
hs_scripbilin = regrid2d(hs, cgrido, 'bilinear')


print('Plots')
P.rcParams['font.size'] = 10
P.figure(figsize=(5, 8))
P.subplots_adjust(hspace=.28, bottom=.07, left=.08, right=.98)
kwplot = dict(show=False, colorbar=False, vmin=hs.min(), vmax=hs.max(),
    drawparallels_size=8, drawmeridians_size=8, drawmeridians_rotation=45.,
    xhide='auto',  yhide='auto')
m = map2(hs, title='Original', subplot=311, **kwplot)
add_grid(cgrido, lw=.7, alpha=.3)
map2(hs_scripcons, title='SCRIP : cellave', subplot=312, m=m, **kwplot)
add_grid(cgridi, lw=.7, alpha=.3)
map2(hs_scripbilin, title=u'SCRIP : bilinear', subplot=313, m=m, **kwplot)
add_grid(cgridi, lw=.7, alpha=.3)
P.rcdefaults()
