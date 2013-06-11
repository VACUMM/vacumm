# Creation d'un axe de temps
from vacumm.misc.axes import *
import numpy as N, cdms2
time_axis = create_time(N.arange(10.),
    'days since 2006-10-01',long_name='Mon axe de temps')

# Creation des axes geographiques
# - longitude: on change l'id de 'lon' a 'longitude'
lon_axis = create_lon(N.arange(5)-5.,id='longitude')
# - latitude
lat_axis = create_lat(N.arange(10)*.5+44.)
# - profondeur
dep_axis = create_dep(N.arange(0,-400.,50),units='m')

# Un axe quelconque
bad_axis = cdms2.createAxis([1],id='pipo')

# Verification des types d'axes avec 'axis_type'
# - affichage pour tous
print ' | '.join(['%s:%s'%(axis.id,axis_type(axis)) 
    for axis in time_axis,lon_axis,lat_axis,dep_axis,bad_axis])
#  -> time:t | other_time:t | longitude:x | lat:y | depth:z | pipo:-
# - verification ponctuele
print islon(lon_axis),islon(lat_axis),is_geo_axis(lat_axis)
#  -> True False False

# Reformattage des axes pour deviner identifier
#  ceux geographiques (lon,lat,dep,time)
# Le reformattage peut changer : id, units, long_name.
# - un axe pourri mais avec 'long_name' explitcite
time_axis2 = cdms2.createAxis([6],id='other_time')
time_axis2.long_name = 'Time'
check_axis(time_axis2) 
print istime(time_axis2) # en fait, 'check_axis' appelle 'istime'
#  -> True
# - une variable tiree d'un fichier (voir le tutoriel (*@\ref{lst:misc.io.netcdf}@*))
from vacumm.config import data_sample
f = cdms2.open(data_sample('mars3d.xt.xe.nc'))
var = f('xe',time=slice(0,1), squeeze=1)
f.close()
check_axes(var)
lon_axis = var.getAxis(0)
print lon_axis.axis, islon(lon_axis)
#  -> X True
