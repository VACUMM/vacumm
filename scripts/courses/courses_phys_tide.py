#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Diagnostics de maree"""

from vcmq import cdms2, P, curve2, savefigs, data_sample
from vacumm.tide.filters import demerliac, godin
from vacumm.tide.filters import extrema, zeros
from vacumm.tide.marigraph import Marigraph
from vacumm.tide.station_info import StationInfo

# Stations
station = StationInfo('Brest')
print station.attributes()
print station.name, station.longitude
print 'Niveau moyen a Brest:', station.nm


# Read sea level at Brest
f = cdms2.open(data_sample("tide.sealevel.BREST.mars.nc"))
sea_level = f('sea_level')
f.close()


# Surcotes/decotes
cotes, tide = demerliac(sea_level, get_tide=True)            # -> ESSAYER GODIN
kwp = dict(date_fmt='%b', date_locator='month')
curve2(sea_level, 'b', show=False, figsize=(15, 4), **kwp)
curve2(tide, 'r', **kwp)
curve2(cotes, figsize=(15, 4), **kwp)


# Extremas
slzoom1 = sea_level(('2006-10-01', '2006-10-02'))[::4] # Toutes les heures
bm, pm = extrema(slzoom1, spline=True, ref='mean')           # -> SANS SPLINES
zz = zeros(slzoom1, ref='mean')                              # -> AUTRES REFERENCE ?
curve2(slzoom1, 'ko', markersize=3, figsize=(6, 4), show=False)
curve2(zz, 'go', linewidth=0, show=False, xstrict=False)
curve2(pm, 'ro', linewidth=0, show=False, xstrict=False)
curve2(bm, 'bo', linewidth=0, xstrict=False, title="Niveau de la mer")


# Outil marégraphique
slzoom2 = sea_level(('2006-10', '2006-11'))[::4] # Toutes les heures
mg = Marigraph(slzoom2, verbose=True)
tide = mg.tide(tide_filter='demerliac')                     # -> ESSAYER COTES/HIGH...
mg.plot()                                                   # -> SELECTION
