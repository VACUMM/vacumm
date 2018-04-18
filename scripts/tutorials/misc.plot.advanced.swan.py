# -*- coding: utf8 -*-
# Lecture et masquage de Hs
from vcmq import  cdms2, data_sample

f = cdms2.open(data_sample('swan.four.nc'))
lon = f('longitude')
lat = f('latitude')
missing_value = f['HS']._FillValue
hs = f('HS', squeeze=1) 
f.close()
hs[:] = cdms2.MV.masked_object(hs, missing_value, copy=0)

# Trace sur grille irrégulière
from vacumm.misc.plot import map2
map2(hs, xaxis=lon, yaxis=lat,  figsize=(6, 4), 
    long_name = 'Significant wave height', fill='pcolormesh', contour=False, 
    left=.12, right=1, show=False)

