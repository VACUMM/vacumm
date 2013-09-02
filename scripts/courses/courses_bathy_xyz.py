#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Manipulation de bathymetrie en semis de points (:mod:`vacumm.bathy.bathy`)"""

from vcmq import N, os, merc, create_grid, resol, map2, P
from vacumm.bathy.bathy import XYZBathy, XYZBathyMerger

# Creation de fausses bathymetries xyz
# - fonction generatrice
def gene_bathy(xc, yc, xr, yr, n=500, amp=30.):
    noise = N.random.random(n)
    a = N.random.random(n)*N.pi*2
    r = N.random.random(n)
    x = xc + xr*r*N.cos(a)
    y = yc + yr*r*N.sin(a)
    return N.asarray([x, y, N.exp(-(x-xc)**2/xr**2-(y-yc)**2/yr**2)*amp+noise])
# - creations
xyz1 = XYZBathy(gene_bathy(-5, 48.3, .3, .15)) # top right
xyz2 = XYZBathy(gene_bathy(-5.45, 48.1, .45, .25, n=1000, amp=20), transp=False) # bot left

# Recuperation de donnees
x = xyz1.x
y = xyz1.y
z = xyz1.z

# Exclusions/selections
xyz1.select([[-5.4, 48.1], [-4.8, 48.1], [-5.1, 48.5]])
xyz1.exclude([-5.2, 48., -5, 48.25])


# Infos
print xyz1.xmin                                         # -> essayer get_xmin avec mask
print xyz1.resol()                                      # -> essayer avec deg=...

# Plot
kwp = dict(size=40, map_res=None, masked_alpha=.1)
xyz1.plot(mode='both', **kwp)                           # -> essayer autre mode

# Fichier
fbathy = __file__[:-2]+'xyz'
xyz1.save(fbathy)
xyz3 = XYZBathy(fbathy)                             
xyz3.plot(title='XYZ3', **kwp)


# Fusions
xyz = xyz1 + xyz2
merger = XYZBathyMerger()
merger += xyz1                                          # -> supprimer/ajouter
merger.append(xyz2)
for i, b in enumerate(merger):
    b.long_name = 'Niveau : %i'%i # Pour les legendes
    b.set_transp(False) # Opacite
merger.plot(mode='cluster', **kwp)                      # -> essayer autre mode
xyz = merger.get_xyz(long_name='Merged')
xyz.plot(**kwp)


# Regrillage
grid_auto = xyz.grid
grid_manual = create_grid((-5.3, -4.91, .01), (48.1, 48.41, .01))
gridded = xyz.togrid(grid_manual)
map2(gridded, res=None)

