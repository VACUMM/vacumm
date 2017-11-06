# -*- coding: utf8 -*-
# Lecture d'une serie 1D de niveau de la mer du modele
import cdms2
from vacumm.config import data_sample
f = cdms2.open(data_sample('tide.sealevel.BREST.mars.nc'))
sea_level = f('sea_level', ('2006-10-01', '2006-10-02'))[1::4] # Toutes les heures
f.close()

# Recuperation des pleines mers, basses mers et zeros
from vacumm.tide.filters import extrema, zeros
bm, pm = extrema(sea_level, spline=True, reference='mean')
zz = zeros(sea_level, ref='mean')

# Plots
from vacumm.misc.plot import curve2 as curve
curve(sea_level, 'ko', markersize=3, figsize=(6, 4), show=False)
curve(zz, 'go', linewidth=0, show=False, xstrict=False)
curve(pm, 'ro', linewidth=0, show=False, xstrict=False)
curve(bm, 'bo', linewidth=0, xstrict=False, title="Niveau de la mer",
   savefigs=__file__, savefigs_pdf=True, show=False, close=True)
