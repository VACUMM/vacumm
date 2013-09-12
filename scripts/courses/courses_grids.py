#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Axes et grilles avec VACUMM"""

from vcmq import N, MV2, create_lon, create_lat, create_grid, isgrid, isrect, islon, set_grid,  get_grid, get_axis, varsel, resol, curv2rect, get_xy, meshgrid, meshcells, create_dep, isregular, P, rotate_grid, shiftgrid, extendgrid, create_axes2d, isdepthup, coord2slice, monotonic, xshift, depth2dz, get_closest

# Créer

# - axes
lon = create_lon((2., 11, 2.))  # -> SPECIFIEZ LE LONG_NAME
lat = create_lat(N.arange(43, 50.))
dep = create_dep((0., 10))
# -> AFFICHEZ LES INFOS
xx, yy = N.meshgrid(N.arange(5.), N.arange(4.))
lon2d, lat2d = create_axes2d(xx, yy)
ii = lon2d.getAxis(1)

# - grille
grid = create_grid(lon, lat)    # -> ESSAYEZ AVEC LON EXPLICITE
gridc = create_grid(lon2d, lat2d)


# Verifier
print islon(lon)
print isgrid(grid)              # -> TEST PARAM CURV=...
print isrect(gridc)             # -> CREEZ GRILLE NON RECT ET RETESTER
print isdepthup(dep)            # -> TESTEZ EN CHANGEANT ATTRIBUT POSITIVE ET VALEURS
print isregular(lon)            


# Affecter
var = MV2.ones(grid.shape)
set_grid(var, grid)
varc = MV2.ones(gridc.shape)
set_grid(varc, gridc)            # -> VERIFIEZ ID DES AXES


# Récupérer
mygrid = get_grid(gridc)        # -> TESTEZ AVEC (LON,LAT) ET PARAMS STRICT ET INTERCEPT
mylon2d = get_axis(gridc, 1)    # -> COMPAREZ AVEC .GETAXIS()


# Sélectionner
print coord2slice(lon, lon=(2, 4.5))# -> COMPAREZ AVEC .GETINTERVALEXT(...)
print coord2slice(grid, lon=(4, 8), lat=(44, 46)) # -> TESTEZ SUR GRIDC
# -> TESTEZ VARSEL


# Transformer
gridcr = curv2rect(gridc)       # -> TESTEZ AVEC VOTRE GRILLE NON RECT 
ax = create_lon([350, 0, 10.])
print monotonic(ax)
print xshift(var, 2)[0]         # -> COMPAREZ A VAR ET VERIFIEZ AXE
# -> TESTEZ XEXTEND
grido = rotate_grid(grid, 30)   # -> TRACEZ LES LONGITUDE (PCOLOR)
grids = shiftgrid(grid, 1, -1)  # -> VERIFIEZ PUIS TESTEZ SUR GRIDC + PASSER DE T À U
gride = extendgrid(grid, iext=(2, 3))   # -> VERIFIER PUIS TESTEZ LE MODE 


# Exploiter
print resol(grid)               # -> TESTEZ EN METRES ET SUR AXE
print depth2dz(dep)
print get_closest(lon2d, lat2d, 2.3, 1.2)


# Utilitaires sut les coordonnées
xx, yy = meshgrid(xx, yy[:, 0])
xxb, yyb = meshcells(xx, yy)      # -> EN 1D?
