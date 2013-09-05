#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Gestion des coordonnées sigma et interpolation 1D"""

from vcmq import NcSigma, SigmaGeneralized, sigma2depths, cdms2, N, section2, create_depth, regrid1d, data_sample

  
# Lecture de la température
f =cdms2.open(data_sample('mars3d.tsigyx.nc'))
data_in = f('TEMP') # T-YX (- = level)


# Détermination des profondeurs d'après les sigma
sigma_class = NcSigma.factory(f)                # détection auto de la classe
# -> VERIFIER QUE sigma_class EST BIEN SigmaGeneralized
sigma_converter = sigma_class(copyaxes=True)    # initialisation du convertisseur
depths_in = sigma_converter().filled()          # lecture eta, etc + conversion
f.close()


# Creation de l'axe des profondeurs cibles
depths = N.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
    22,24,26,28,30,32,34,36,38,40,45,50,55,60,65,70,75,80,85,90,95,100,120,140,160])
depths = -depths[::-1] # croissantes et négatives
depth_out = create_depth(depths)


# Interpolation
xmap = (0, 2, 3)                                # la profondeur varie en T/Y/X
xmapper = N.rollaxis(depths_in, 1, 4)           # profondeur = dernier axe
data_out = regrid1d(data_in, depth_out, axis=1, method='linear', 
    xmap=xmap, xmapper=xmapper, extrap=1)


# Plot
kw = dict(vmin=10, vmax=14, xhide='auto', add_grid=True, ymax=0, fill='contourf')   # FILL=PCOLOR ?
section2(data_in[0, :, 10], yaxis=depths_in[0, :, 10], subplot=211, title='Sigma', show=False, **kw)
s = section2(data_out[0, :, 10], subplot=212, title='Z', **kw)


