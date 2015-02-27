#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Masquages évolués

Modules: :mod:`vacumm.misc.grid.masking`
"""

from vcmq import *
from vacumm.misc.grid.masking import polygon_mask, erode_coast, get_coastal_indices, \
    get_coastal_indices, envelop, polygons, polygon_select, get_dist_to_coast, \
    convex_hull


# Creation d'un masque via un trait de cote
grid = create_grid((-6., -4., 0.1), (47.5, 49.1, 0.1))  # grille de l'iroise
resolution = 'h'                                        # resolution gshhs : testez 'i'
mask = polygon_mask(grid, resolution, thresholds=0.5)   # testez threshold = 0.1
mvmask = MV2.array(mask)
set_grid(mvmask, grid)
map2(mvmask, fillcontinents=False, drawcoastlines_color='r', drawcoastlines_zorder=100,
    res=resolution, proj='merc', fill='pcolor', drawcoastlines_linewidth=2,
    cmap=cmap_rs(['1', '0.6']), contour=False, colorbar=False,
    savefig ='courses_masking_0.png', show=True)


# Cote via mask
indices = get_coastal_indices(mask)
dist = get_dist_to_coast(grid, mask)


# Erosion de la cote
f = cdms2.open(data_sample('mars3d.xy.nc'))
temp = f('temp')
f.close()
temp2 = erode_coast(temp, maxiter=3)                # changez maxiter et tracez


# Enveloppe de points
xy = P.randn(2, 500)*0.5+N.array([-5,48]).reshape(2,-1)
xe, ye = convex_hull(xy, method='delaunay')         # changez la methode


# Creer des plolygones
xpoly = [-5.8, -5, -5, -5.8]
ypoly = [48, 48, 48.5, 48.5]
polys = polygons([[xpoly, ypoly]])                  # meme chose avec des min/max + tracez
# -> appliquez a polygon_mask


# Selection de points
xsel, ysel = polygon_select(xy[0], xy[1], polys)    # essay l'option mask
P.plot(xy[0], xy[1], 'o')
P.plot(xpoly, ypoly, 'k-', lw=2)
P.plot(xe, ye)
P.plot(xsel, ysel, 'ro')
P.savefig('courses_masking_1.png')
P.show()


# Detection de lacs et iles
# -> testez GetLakes

