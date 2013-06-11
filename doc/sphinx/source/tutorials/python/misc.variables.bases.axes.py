import numpy as N, cdms2, MV2

# Creation d'un axe de longitude
# - base
lon = cdms2.createAxis([-5.,-4.,-3.],id='lon')
lon.long_name = 'Longitude'
lon.units = 'degree_east'
# - ajout de lon.axis='X' et lon.modulo = '360.'
lon.designateLongitude()

# Latitude
lat = cdms2.createAxis([46.,47.,48.],id='lat')
lat.long_name = 'Latitude'
lat.units = 'degree_north'
lat.designateLatitude() # lat.axis = 'Y'

# Pronfondeur
depth = cdms2.createAxis([-200.,-100.,-0.],id='depth')
depth.long_name = 'Depth'
depth.units = 'm'
depth.designateLevel() # depth.axis = 'Z'

# Temps
# - creation
time = cdms2.createAxis([0.,1.,2.],id='time')
time.long_name = 'Time'
time.units = 'days since 2006-08-01'
time.designateTime() # time.axis = 'T'
# - verif
ctime = time.asComponentTime()
print ctime,ctime[1].day
#  -> [2006-8-1 0:0:0.0, 2006-8-2 0:0:0.0, 2006-8-3 0:0:0.0] 2
rtime = time.asRelativeTime()
print rtime,rtime[1].value
#  -> [0.00 days since 2006-08-01, 1.00 days since 2006-08-01,
#      2.00 days since 2006-08-01] 1.0


# Maintenant, on cree une variable avec ces axes
#- methode direct
temp1 = cdms2.createVariable(N.ones((3,3,3,3)),typecode='f',id='temp',
    fill_value=1.e20,axes=[time,depth,lat,lon],copyaxes=0,
    attributes=dict(long_name='Temperature',units='degC'))
# - remarque
print cdms2.createVariable is MV2.array
#  -> True         (ce sont les meme fonctions !)
# - une methode indirecte
#   . initialisation
temp2 = MV2.array(N.ones((3,3,3,3))).astype('f')
#   . attributs
temp2.id = 'temp'
temp2.long_name = 'Temperature'
temp2.units = 'degC'
temp2.set_fill_value(1.e20) # <=> temp2.setMissing(1.e20)
#   . axes
temp2.setAxisList([time,depth,lat,lon])
#   . ou par exemple individuellement
temp2.setAxis(1,depth)
