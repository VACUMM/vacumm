#! /usr/bin/env python
# -*- coding: utf-8 -*- 
# 
# Extraction des courants juste en-dessous de la couche d'Ekman
# La profondeur de la couche d'Ekman est calculee a partir des forcages 
# atmospheriques au pas de temps considere.
#
# G. Charria (01/2011)
# ________________________________________________________________

from vacumm.misc.axes import create_time, create_dep, create_lat, create_lon
from vacumm.misc.io import ncread_best_estimate, NcIterBestEstimate, list_forecast_files
import cdms2, os
import cdtime, MV2
from vacumm.misc.grid.regridding import interp1d, grid2xy, interp2d
from vacumm.misc.grid.masking import erode_coast
from vacumm.data.misc.sigma import Sigma
from vacumm.misc.atime import strftime, daily
import numpy as np

# Pour produire du NetCDF3
cdms2.setNetcdfShuffleFlag(0); cdms2.setNetcdfDeflateFlag(0); cdms2.setNetcdfDeflateLevelFlag(0)

#dir_model = '/home2/creizic/pgarreau/jonsmod_015'
dir_model = '/temp/gcharria'

#print os.listdir(dir_model)

ctdeb=cdtime.comptime(2007,10,10,0,0,0)
#ctfin=cdtime.comptime(2007,10,10,3,0,0)
ctfin=cdtime.comptime(2007,12,7,0,0,0)

ctest = ctdeb
# Boucle sur la periode consideree.
while ctest < ctfin:
    print ctest

    c1 = ctest
    c2 = c1.add(3,cdtime.Hours)
    U = ncread_best_estimate('UZ',os.path.join(dir_model,"champs_V8.06_%Y%m%d%H0000_Z.nc"), (c1, c2))
    V = ncread_best_estimate('VZ',os.path.join(dir_model,"champs_V8.06_%Y%m%d%H0000_Z.nc"), (c1, c2))

    
    # -- Profondeurs d'interpolation
    lo = U.getLongitude()
    la = U.getLatitude()
    


    # Longitudes toute les 1 minute
    loi = np.arange(lo[0],lo[-1],1./60.)
    # Latitudes toute les 1 minute
    lai = np.arange(la[0],la[-1],1./60.)

    xo = cdms2.createAxis(loi)
    yo = cdms2.createAxis(lai)


    LONGITUDE = create_lon(loi, id='longitude', attributes=dict(long_name='Longitude of each location',standard_name='longitude',units='degrees_east',valid_min='-180.',valid_max='180.',axis='X'))
    LATITUDE = create_lat(lai, id='latitude', attributes=dict(long_name='Latitude of each location',standard_name='latitude',units='degrees_north',valid_min='-90.',valid_max='90.',axis='Y'))

    z = U.getLevel()
    DEPTH = create_dep(z,attributes=dict(long_name='sea water depth',standard_name='depth',units='m',valid_min='0.',valid_max='12000.'))

    #axes = [TIME,DEPTH,LATITUDE,LONGITUDE]
    axes = [DEPTH,LATITUDE,LONGITUDE]

    Uis=np.arange(U.getLevel().__len__()*yo.__len__()*xo.__len__()).reshape(U.getLevel().__len__(),yo.__len__(),xo.__len__())
    Vis=np.arange(U.getLevel().__len__()*yo.__len__()*xo.__len__()).reshape(U.getLevel().__len__(),yo.__len__(),xo.__len__())


    Uis = cdms2.createVariable(Uis, typecode='f',id='UZ', axes=axes, attributes=dict(long_name='3d zonal velocity',standard_name='eastward_sea_water_velocity',units='m.s-1',valid_min='-100.',valid_max='100.'))
    Vis = cdms2.createVariable(Vis, typecode='f',id='VZ', axes=axes, attributes=dict(long_name='3d meridional velocity',standard_name='northward_sea_water_velocity',units='m.s-1',valid_min='-100.',valid_max='100.'))



    for iz, dep in enumerate(U.getLevel()):
        Ui2 = interp2d(U[0,iz,:,:], (xo,yo), method='bilinear')
        Vi2 = interp2d(V[0,iz,:,:], (xo,yo), method='bilinear')
 
        Uis[iz,:,:]=Ui2
        Vis[iz,:,:]=Vi2



    date=strftime('_%Y%m%d%H0000',c1)
    file_out = 'MARS_MENOR_'+date+'.nc'
    f = cdms2.open(os.path.join(dir_model,file_out), 'w')
    #f.write(LATITUDE)
    #f.write(LONGITUDE)
    #f.write(TIME)
    #f.write(DEPTH)


    f.write(Uis) # ecriture d'une variable
    f.write(Vis)
    f.history = 'Created with '+__file__.encode('utf8') # attribut global
    f.close() # fermeture

    del Uis, Vis

    ctest = c2

