# -*- coding: utf8 -*-
from vacumm.misc.grid.regridding import regrid2d
import MV2
from vacumm.misc.grid import create_grid
from vacumm.misc.plot import map2
gmax=40
ggi = create_grid((0., gmax, 1.), (0., gmax, 1.))
ggo = create_grid((0.01, gmax, .5), (0.01, gmax, .5))

vari = MV2.reshape(MV2.arange(ggi.size(), dtype='f'), ggi.shape)
vari.setGrid(ggi)
varo = regrid2d(vari, ggo, 'bilinear')
map2(varo, clabel_glow=True, resolution=None, show=False, savefigs=__file__,
    close=True)
