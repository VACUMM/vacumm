#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Diagnostics thermodynamiques (:mod:`vacumm.diag.thermdyn`)"""

from vcmq import DS, data_sample, NcSigma, map2, density, mixed_layer_depth

# Lecture des données de base
ds = DS(data_sample('menor.nc'))
temp = ds.get_temp(squeeze=1)
sal = ds.get_sal(squeeze=1)
depth = ds.get_depth(squeeze=1)


# Densité
dens_surf = density(temp[-1], sal[-1])                      # -> CALCULER DENSITE 3D
map2(dens_surf)


# Couche mélangée
mld = mixed_layer_depth(temp, sal, depth, mode='deltadens') # -> ESSAYER AUTRES MODES
map2(mld)
