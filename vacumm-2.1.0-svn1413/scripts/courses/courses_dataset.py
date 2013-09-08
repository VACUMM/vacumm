#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Exploitation generique de donnes grilles (:mod:`vacumm.data.misc.dataset`)"""

from vcmq import DS, data_sample, map2

# Initialisation
ds = DS(data_sample('menor.nc'), 'mars', lon=(4, 5), lat=(42.5, 43.5), log_level='debug')
print ds.dataset


# Logging
ds.info('Pratique')                                 # -> ESSAYER WARNING
ds.set_loglevel('debug')


# Coordonnees
lon = ds.get_lon()                                  # -> ESSAYER AU POINT U
grid = ds.get_grid()


# Variables
temp = ds.get_temp(level=slice(-1, None), squeeze=True)
sst = ds.get_sst(squeeze=True)                      # -> VERIFIER ATTRIBUTS
ds.info('Plot SST')
map2(sst)
sal = ds.get_sal(lon=(4., 4.5))
print temp.shape, sal.shape

# Generique
ds2 = DS(data_sample('mfs.nc'), 'nemo', lon=(4, 5), lat=(42.5, 43.5))
sst2 = ds2.get_sst(squeeze=True)                      # -> VERIFIER ATTRIBUTS

# Plus evolue
depth = ds.get_depth() # sigma
print depth.shape
mld = ds.get_mld()
ds.plot_transect('temp', lons=(4.1, 4.9), lats=(43.3, 42.6), figsize=(12, 6)) # -> OUTAXIS=
