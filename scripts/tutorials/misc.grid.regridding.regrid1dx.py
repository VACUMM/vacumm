# -*- coding: utf-8 -*-
import cdms2, numpy as N, pylab as P
from vcmq import data_sample, create_dep, regrid1d, add_grid, section2

# Lecture de la temperature
f = cdms2.open(data_sample('mars3d.xz.nc'))
t = f('temp', lon=(-4.9, -4.43))
h0 = f('h0', lon=(-4.9, -4.43)).filled(0.)
f.close()
t.long_name = 'Original'

# Creation d'un axe de profondeur
ddep = 5.
dep = create_dep((-h0.max(),0.+ddep , ddep))

# Creation de l'axe etendu (taille (ndep,nx))
depths = -N.outer(t.getAxis(0)[::-1], h0)
dep[-1] = depths[-1].max()

# Regrillage lineaire
tr = regrid1d(t, dep, 'linear', axi=depths, axis=0)
tr.long_name = 'Interpolated'
t.getAxis(0).designateLevel()

# Plot
P.figure(figsize=(5.5, 6))
kwplot = dict(show=False, bgcolor='.5', ylim=(-80, 0))
# - original
section2(t, yaxis=depths, subplot=211, hspace=.3, **kwplot)
add_grid((t.getLongitude(), depths[:]), linewidth=.3)

# - regridded
section2(tr, subplot=212, **kwplot)
add_grid((tr.getLongitude(), dep[:]), linewidth=.3)

