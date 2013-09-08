#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Interpolation depuis des positions aleatoires"""

# Generation des donnees
# - sur grille
import numpy as N
xr = N.arange(20.)
yr = N.arange(10.)
xxr, yyr = N.meshgrid(xr, yr)
zzr = (N.sin(xxr*N.pi/6)*N.sin(yyr*N.pi/6) + \
    N.exp(-((xxr-7.)**2+(yyr-7.)**2)/4.**2))*100.
zzr -= zzr.min()
vminmax = dict(vmin=zzr.min(), vmax=zzr.max())
# - un jeu alatoire (input)
ij = N.unique((N.random.rand(50)*zzr.size).astype('i'))
xi, yi, zi = xxr.flat[ij], yyr.flat[ij], zzr.flat[ij]
# - un autre jeu aleatoire (output)
ij2 = N.unique((N.random.rand(100)*zzr.size).astype('i'))
xo, yo = xxr.flat[ij2], yyr.flat[ij2]


# Sans VACUMM

# - CDAT natgrid
from nat import Natgrid
I = Natgrid(xi, yi, xr, yr)
I.ext = 1                                               # Changer
I.igr = 1                                               # les parametres
zzoc = I.rgrd(zi).T

# - Matplotlib / R. Kern
from matplotlib.delaunay import Triangulation
tri = Triangulation(xi,yi)
I = tri.nn_extrapolator(zi)                             # Tester les autres methodes
zzork = I(xxr,yyr)  # vers grille
zork = I(xo,yo) # vers points


# Avec VACUMM
from vacumm.misc.grid.regridding import GridData, griddata, xy2xy
zzov = griddata(xi, yi, zi, (xr, yr), method='nat', ext=True, sub=10)     
# -> Tester la methode "carg"
# -> Testez parametre sub=...
# -> Essayer avec GridData
zov2 = xy2xy(xi, yi, zi, xo, yo)

# Krigeage
from vacumm.misc.grid.kriging import krig
zzok = krig(xi, yi, zi, xxr.ravel(), yyr.ravel(), nproc=1).reshape(zzr.shape)
# -> Tester nproc et npmax


# Plots
from vcmq import meshbounds, P
xxrb, yyrb = meshbounds(xr, yr)
P.figure(figsize=(10, 8))
axis = [xxrb.min(), xxrb.max(), yyrb.min(), yyrb.max()]
#
P.subplot(332)
P.pcolormesh(xxrb, yyrb, zzr, **vminmax)
P.scatter(xi, yi, c='k')
P.title('Original')
P.axis(axis)
#
P.subplot(334)
P.pcolormesh(xxrb, yyrb, zzoc, **vminmax)
P.title('CDAT/Natgrid')
P.axis(axis)
#
P.subplot(335)
P.pcolormesh(xxrb, yyrb, zzork, **vminmax)
P.title('MPL/R.K. grid')
P.axis(axis)
#
P.subplot(336)
P.pcolormesh(xxrb, yyrb, zzr, **vminmax)
P.scatter(xo, yo, c=zork, s=50, **vminmax)
P.title('MPL/R.K. random')
P.axis(axis)
#
P.subplot(337)
P.pcolormesh(xxrb, yyrb, zzov, **vminmax)
P.title('VACUMM/nat grid')
P.axis(axis)
#
P.subplot(338)
P.pcolormesh(xxrb, yyrb, zzr, **vminmax)
P.scatter(xo, yo, c=zov2, s=50, **vminmax)
P.title('VACUMM/nat random')
P.axis(axis)
#
P.subplot(339)
P.pcolormesh(xxrb, yyrb, zzok, **vminmax)
P.title('Kriging')
P.axis(axis)
#
P.tight_layout()
P.show()

# Transect
