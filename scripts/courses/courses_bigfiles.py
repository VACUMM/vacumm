#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Calculs sur les gros fichiers"""

from vcmq import cdms2, Intervals, data_sample, StatAccum, cdtime, curve2, N


# Initialisations

# - ouverture du fichier et infos temporelles
f = cdms2.open(data_sample('mars3d.xt.xe.nc'))
ctimes = f['xe'].getTime().asComponentTime()

# - passage par un accumulateur de stats
sa = StatAccum(tmean=True)


# Boucle sur des intervals journaliers
for itv in Intervals((ctimes[0], ctimes[-1], 'cc'), 'day'): # -> ESSAYER 2 JOURS

    # Lecture
    print itv
    tmp = f('xe', time=itv)
    
    # Accumulation
    sa += tmp
    del tmp
    
# Finalisation

# - restitution
xem = sa.get_tmean()
print N.ma.allclose(xem.asma(),f('xe').mean(axis=0))
curve2(xem)
del sa

# - fermeture
f.close()

