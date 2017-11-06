# -*- coding: utf8 -*-
import pylab as P

# Definition de la region de travail Bretagne
zone = (-6.5, 47.2, -2, 49.05)

# Recuperation du trait Europeen
from vacumm.bathy.shorelines import GSHHS
tc = GSHHS(input='h', clip=zone)

# Interpolation du niveau de la mer sur le trait
xyz = tc.bathy()

# Plot
xyz.plot(size=10, savefigs=__file__, show=False)#, savefigs_pdf=True)
P.close()
