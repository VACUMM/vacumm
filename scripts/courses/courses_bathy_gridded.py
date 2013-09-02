#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Manipulation de bathymetrie sur grille (:mod:`vacumm.bathy.bathy`)"""

from vcmq import cdms2, data_sample, create_grid
from vacumm.bathy.bathy import GriddedBathy, NcGriddedBathy, bathy_list, GriddedBathyMerger, plot_bathy

# Lecture d'une variable dans un fichier
cdms2.axis.latitude_aliases.append('y')
cdms2.axis.longitude_aliases.append('x')
f = cdms2.open(data_sample('ETOPO2v2g_flt.grd'))
var = f('z', lon=(-6.1, -4), lat=(47.8, 48.6))
f.close()

# Creation de l'objet de bathy grillee
bathy = GriddedBathy(var)

# Versions masquees
bathy_masked1 = bathy.masked_bathy()                        # -> nb de points masques ?
bathy.set_shoreline('f')                                    # -> essayer 'i'
bathy_masked2 = bathy.masked_bathy()                        # -> nb de points masques ?

# Plot
bathy.plot(title='Direct')


# Bathy standard
from vacumm.bathy.bathy import NcGriddedBathy, bathy_list
print bathy_list().keys()
sbathy = NcGriddedBathy(lon=(-5.2, -3), lat=(48., 48.9), name='etopo2')
sbathy.plot(title='Bathy par defaut')


# Fusion de bathy
mgrid = create_grid((-6.1, -3, 0.01), (47.8, 48.9, .006))   # -> changer la grille
from vacumm.bathy.bathy import GriddedBathyMerger
merger = GriddedBathyMerger(mgrid)
merger += sbathy
merger += bathy
merger.plot(title='Fusion / plot direct')

# Plot via plot_bathy
mbathy = merger.merge()
plot_bathy(mbathy, title='Fusion / plot via plot_bathy')


# Regrillage
rgrid = create_grid((-6.1, -3, 0.05), (47.8, 48.9, .02))
Mbathy = GriddedBathy(mbathy)
rbathy = Mbathy.regrid(rgrid)
plot_bathy(rbathy, title='Regrillee')

