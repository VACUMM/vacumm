# -*- coding: utf8 -*-
from __future__ import print_function
from matplotlib import rc;rc('font', size=11)
import cdms2, MV2
from vcmq import (data_sample, curv_grid, set_grid, rotate_grid, 
    map2, add_grid, regrid2d)
import pylab as P

# Lecture des donnees
zone = dict(yc=slice(0, 40), xc=slice(10, 40))
f = cdms2.open(data_sample('swan.four.nc'))
lon2d = f('longitude', **zone)
lat2d = f('latitude', **zone)
hs = MV2.masked_values(f('HS', squeeze=1, **zone), f['HS']._FillValue)
dlon = float(f.longitude_resolution)
dlat = float(f.latitude_resolution)
f.close()

# Affectation de la grille curvilineaire à la variable
cgridi = curv_grid(lon2d, lat2d)
hs = set_grid(hs, cgridi)

# Rotation de la grid
cgrido = rotate_grid(hs.getGrid(), -30)


print('Regrillage')
print(' - par SCRIP/conservative')
hs_scripcons = regrid2d(hs, cgrido, 'conservative')
print(' - par SCRIP/bilineaire')
hs_scripbilin = regrid2d(hs, cgrido, 'bilinear')


print('Plots')
P.rcParams['font.size'] = 10
P.figure(figsize=(4, 7))
P.subplots_adjust(hspace=.28, bottom=.07, left=.08, right=.98)
kwplot = dict(show=False, colorbar=False, vmin=hs.min(), vmax=hs.max(),
    drawparallels_size=8, drawmeridians_size=8, drawmeridians_rotation=45.,
    xhide='auto',  yhide='auto')
m = map2(hs, title='Original', subplot=311, **kwplot)
add_grid(cgrido, lw=.7, alpha=.3)
map2(hs_scripcons, title='SCRIP : remapping', subplot=312, m=m, **kwplot)
add_grid(cgridi, lw=.7, alpha=.3)
map2(hs_scripbilin, title=u'SCRIP : bilinéaire', subplot=313, m=m, **kwplot)
add_grid(cgridi, lw=.7, alpha=.3)
P.rcdefaults()
