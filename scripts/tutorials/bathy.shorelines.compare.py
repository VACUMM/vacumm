"""Compare two shorelines"""
from vcmq import GSHHS, Histolitt
from pylab import legend, title

# %% Read shapes
zone = (-5.15, 48.42, -5.03, 48.49)
gmt = GSHHS('h', clip=zone)
thc = Histolitt(clip=zone)

# %% Plot
kw  = dict(fill=False, show=False, linewidth=1.5)
thc.plot(color='tab:blue', zorder=12, m_left=.1, m='auto',
    label='Histolitt (SHOM/IGN)', m_figsize=(5.5, 6), **kw)
gmt.plot(color='k', zorder=10, label='GSHHS (h)',
    m='auto', **kw)
legend(loc="lower right")
title("Ushant shorelines")
