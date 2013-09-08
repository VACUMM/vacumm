#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Filtres 1D et 2D"""

from vcmq import cdms2, curve2, map2, generic1d, generic2d, shapiro2d, bartlett1d, data_sample, curve2, shapiro1d, gaussian1d, gaussian2d

# Lectures
f = cdms2.open(data_sample('mars3d.t.nc'))
sst1d = f('temp')
f.close()
f = cdms2.open(data_sample('menor.nc'))
sst2d = f('TEMP', time=slice(0, 1), level=slice(-1, None), squeeze=1)
f.close()



# 1D

# - filtrages
sst1d_gen13 = generic1d(sst1d, 13)
sst1d_gen3 = generic1d(sst1d, [1., 1., 1.])
sst1d_sha = shapiro1d(sst1d)        # -> SHAPIRO AVEC GENERIC1D
sst1d_bar13 = bartlett1d(sst1d, 13)
# -> TESTEZ GAUSSIEN SUR 13 HEURES

# - plots
curve2(sst1d, 'k', label='Original', show=False, figsize=(12, 5))
curve2(sst1d_gen13, 'r', label='Generic 13 pts', show=False)
curve2(sst1d_gen3, 'g', label='Generic 3 pts', show=False)
curve2(sst1d_sha, 'b', label='shapiro', show=False)
curve2(sst1d_bar13, 'm', label='Bartlett', legend=True)


# -> MASQUEZ UNE PARTIE DES DONNEES ET JOUEZ AVEC LE PARAMETRE MASK=

# -> LISEZ UN BLOCK 3D ET FILTREZ LE SUIVANT LE TEMPS


# 2D

# - filtrage
sst2d_gen13 = generic2d(sst2d, 13)
sst2d_gau13 = gaussian2d(sst2d, 13)

# - plots
kw = dict(vmin=sst2d.min(), vmax=sst2d.max(), colorbar=False, nmax=18)
map2(sst2d, title='Original', figsize=(13, 3), subplot=131, show=False, **kw)
map2(sst2d_gen13, title='Generic 13', subplot=132, show=False, **kw)
map2(sst2d_gau13, title='Gauss 13', subplot=133, show=True, **kw)


# -> JOUEZ AVEC MASK=
