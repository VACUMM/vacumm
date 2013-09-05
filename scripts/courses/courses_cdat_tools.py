#!/usr/bin/env python
# -*- coding: utf8 -*- 
"""
Utilitaires intégrés à UVCDAT

Modules: :mod:`cdtime`, :mod:`cdutil`, :mod:`genutil`
"""


# Le module cdtime: la gestion du temps

# - temps absolu
import cdtime
ctime = cdtime.comptime(2000, 1, 1, 6, 23)              # mettre moins d'arguments
print ctime.hour                                        # verifier les autres composantes
print cdtime.comptime(2000) > cdtime.comptime(1999)
print cdtime.comptime(2000).cmp(cdtime.comptime(1999))
print ctime.add(1, cdtime.Hour).hour                    # essayez 70s

# - temps relatif
rtime = cdtime.reltime(23, 'minutes since 2000-1-1 6')
print rtime.value, rtime.units
print rtime.torel('hours since 2000').value

# - conversions
print rtime.tocomp()
print ctime.torel('months since 2000')                  # essayez les "weeks'
print cdtime.s2c('2000-10').month                       # testez cdtime.s2r, cdtime.r2r...





# Le module cdutil : utilitaires orientes climat

# - contenu
import cdutil
print dir(cdutil)

# - chargement des données (vent Pacifique central sur plusieurs années)
from vcmq import *
f = cdms2.open(data_sample('uv_pacific.nc'))
u = f('uwnd')
f.close()

# - construire une climatologie mensuelle et des anomalies
cdutil.setTimeBoundsMonthly(u)                          # importance des bounds (autres ?)
uclim = cdutil.ANNUALCYCLE.climatology(u)               # climato
uanom = cdutil.ANNUALCYCLE.departures(u, ref=uclim)     # anomalies
print uclim.std(), uanom.std()
djf = cdutil.times.Seasons('DJF')                       # creation d'une saison
udjf = djf(u)                                           # extraction
dfj = cdutil.DJF                                        # des saisons existent déjà

# - averager
ut = cdutil.averager(u, axis='yx',  weights=cdutil.area_weights(u)) # moyenne spatiale
help(cdutil.averager)
#  -> essayez la moyenne temporelle

# - regions et selecteurs
equator = cdutil.region.domain(lat=(-2, 2))
select = cdms2.selectors.Selector(lon=slice(0, 3), time=('1950', cdtime.comptime(1960)))
print u(equator)(select).shape
#  -> appliquez à la lecture du fichier






# Le module genutil : utilitaires generiques

# - contenu
import genutil
print dir(genutil)

# - statistics standards
print dir(genutil.statistics)
ustd = genutil.statistics.std(u, axis=0)                # testez les options
ustd.info()
print N.ma.allclose(ustd, u.std(axis=0))

# - salstat
print dir(genutil.salstat)
print genutil.salstat.kurtosis(u, axis='t')             # essayez d'autres fonctions

# - divers
_, uu = genutil.grower(u, u[0])                         # testez d'autres slices
print uu.shape
print genutil.minmax(u, u**2)
print u[:, 0, 0](genutil.picker(time=['1960-01-03', '1965-05-03'])) # testez option match
