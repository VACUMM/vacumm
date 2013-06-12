# Lecture des divers traits de cote
from vacumm.bathy.shorelines import *
zone = (-5.15, 48.42, -5.03, 48.49)
gmt = GSHHS(clip=zone)
euro = EUROSION(clip=zone)
thc = Histolitt(clip=zone)

# Trace
kwpt  = dict(fill=False, points=True, s=2., alpha=.7, points_linewidth=0, show=False)
thc.plot(color='r', zorder=12, m_left=.1, 
    label='Histolitt (SHOM/IGN)', m_figsize=(5.5, 6), **kwpt)
euro.plot(color='g', zorder=11, label='EUROSION', **kwpt)
gmt.plot(fill=False, color='k', linewidth=1.5, zorder=10, label='GSHHS', show=False)

# Fin de plot
from pylab import show, legend, title
from vacumm.misc.plot import savefigs
title("Traits d'Ouessant")
savefigs(__file__, pdf=True)

