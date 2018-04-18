# -*- coding: utf8 -*-
import cdms2, MV2, pylab as P
from vcmq import data_sample, map2

# Lecture de la temperature
f = cdms2.open(data_sample('mars3d.xy.nc'))
temp = f('temp', lon=(-5.5, -4), lat=(47, 49))
f.close()
temp.long_name = 'Original'

# On crée un trou
temp[15:60, 40:80] = MV2.masked

# On rempli des trous
from vacumm.misc.grid.regridding import fill2d
tempf = fill2d(temp, method='carg')
tempf.long_name = 'Rempli'

# On masque la terre
from vacumm.misc.grid.masking import masked_polygon
tempf[:] = masked_polygon(tempf, 'h', copy=0)

# Plots
P.rc('font', size=9)
kw = dict(vmin=9, vmax=13, show=False)
map2(temp, subplot=211, hspace=.2, bottom=.05,
    left=.08, top=.97, figsize=(4.5, 8), nmax=10, **kw)
map2(tempf, subplot=212, **kw)
P.rcdefaults()
