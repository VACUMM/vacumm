from vcmq import GSHHS, EUROSION, Histolitt
from pylab import legend, title

# Lecture des divers traits de cote
zone = (-5.15, 48.42, -5.03, 48.49)
gmt = GSHHS(clip=zone)
euro = EUROSION(clip=zone)
thc = Histolitt(clip=zone)

# Trace
kwpt  = dict(fill=False, points=True, s=2., alpha=.7, points_linewidth=0, show=False)
thc.plot(color='r', zorder=12, m_left=.1,
    label='Histolitt (SHOM/IGN)', m_figsize=(5.5, 6), **kwpt)
euro.plot(color='g', zorder=11, label='EUROSION', m='auto', **kwpt)
gmt.plot(fill=False, color='k', linewidth=1.5, zorder=10, label='GSHHS',
    m='auto', show=False)

# Fin de plot
legend()
title("Ushant shorelines")

