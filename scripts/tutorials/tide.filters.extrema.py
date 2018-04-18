# -*- coding: utf8 -*-
import cdms2
from vcmq import data_sample, curve2, extrema, zeros

# Lecture d'une serie 1D de niveau de la mer du modele
f = cdms2.open(data_sample('tide.sealevel.BREST.mars.nc'))
sea_level = f('sea_level', ('2006-10-01', '2006-10-02'))[1::4] # Toutes les heures
f.close()

# Recuperation des pleines mers, basses mers et zeros
bm, pm = extrema(sea_level, spline=True, reference='mean')
zz = zeros(sea_level, ref='mean')

# Plots
curve2(sea_level, 'ko', markersize=3, figsize=(6, 4), show=False)
curve2(zz, 'go', linewidth=0, show=False, xstrict=False)
curve2(pm, 'ro', linewidth=0, show=False, xstrict=False)
curve2(bm, 'bo', linewidth=0, xstrict=False, title="Niveau de la mer",
   show=False)
