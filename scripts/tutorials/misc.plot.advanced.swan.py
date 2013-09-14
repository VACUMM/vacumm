# -*- coding: utf8 -*-
# Lecture et masquage de Hs
import cdms2, MV2
from vacumm.config import data_sample
from vcmq import code_base_name

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
    savefigs=code_base_name(ext='png'), left=.12, right=1, show=False)

