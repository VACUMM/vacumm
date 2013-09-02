#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Manipulation de traits de cote (:mod:`vacumm.bathy.shorelines`)"""

from vcmq import N, os, merc, create_grid, resol, map2, P
from vacumm.bathy.shorelines import GSHHS, Histolitt

# Lecture de divers traits sur zone reduite
zone = (-6.5, 47.2, -2, 49.05)
gf = GSHHS('h', clip=zone)
gl = GSHHS('l', clip=zone)                      # -> essayer Histolitt

# Trace
kwpt  = dict(fill=False, points=True, s=2., alpha=.7, points_linewidth=0)
gf.plot(color='r', zorder=12, m_left=.1, label='High resolution', show=False, m_figsize=(5.5, 6), **kwpt)
gl.plot(color='g', zorder=11, label='Low resolution', **kwpt)

# Infos
print gf.resol()
print gl.resol()
print gf.xmax
xy = gl.xy                                      # -> essayer get_point() et get_xy()
print xy.shape

# Bathy cotiere
xyz = gf.bathy()
xyz.plot(size=10)

