# -*- coding: utf8 -*-
from vcmq import regrid2d, create_grid, map2, MV2

gmax=40
ggi = create_grid((0., gmax, 1.), (0., gmax, 1.))
ggo = create_grid((0.01, gmax, .5), (0.01, gmax, .5))

vari = MV2.reshape(MV2.arange(ggi.size(), dtype='f'), ggi.shape)
vari.setGrid(ggi)
varo = regrid2d(vari, ggo, 'bilinear')
map2(varo, clabel_glow=True, resolution=None, show=False,
    close=False)
