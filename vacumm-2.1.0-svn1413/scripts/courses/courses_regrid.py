#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Le regrillage"""


from vcmq import N, cdms2, curve2, map2, regrid1d, regrid2d, data_sample, \
    create_grid, regrid_method, CDATRegridder, resol, add_grid, P


# Lecture
f = cdms2.open(data_sample('mfs.nc'))
select = dict(lat=(42.5, 43.5), lon=(3, 4))
ssti = f('votemper', level=slice(0, 1), squeeze=1, **select)
f.close()
gridi = ssti.getGrid()
dxi, dyi = resol(gridi)


# Nouvelle grille
factor = 0.634                                                # -> CHANGEZ LE FACTEUR
dxo = dxi*factor
dyo = dyi*factor
grido = create_grid((3.13, 3.8, dxo), (42.62, 43.4, dyo))


# Quelle méthode ?
method = regrid_method(gridi, grido)
print method


# Regrillage
ssto = regrid2d(ssti, grido, method=method)                 # -> CHANGEZ LA METHODE
# -> UTILISEZ USECDR=TRUE


# Regrilleur de CDAT
regridder = CDATRegridder(gridi, grido, method='cellave')   # -> CHANGEZ LA METHODE
ssto2 = regridder.regrid(ssti)
ssto3, cdr = regrid2d(ssti, grido, cdr=regridder, getcdr=True)  
# -> VERIFIEZ LE TYPE


# Plots
kwp = dict(vmin=ssti.min(), vmax=ssti.max(), res=None, show=False, colorbar=False, 
    drawmeridians_linewidth=0, drawparallels_linewidth=0, grid=False, yhide='auto', **select)
map2(ssti, title='Original',figsize=(15, 5), subplot=131, **kwp)
add_grid(gridi, alpha=1)
map2(ssto, title='Avec regrid2d', subplot=132, **kwp)
add_grid(grido, alpha=1)
add_grid(gridi, edges=False, centers=True, alpha=.5, markersize=3, marker='o')
map2(ssto2, title='Avec CDATRegridder', subplot=133, **kwp)
add_grid(grido, alpha=1)
add_grid(gridi, edges=False, centers=True, alpha=.5, markersize=3, marker='o')
P.tight_layout()
P.show()


