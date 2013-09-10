#!/usr/bin/env python
# -*- coding: utf8 -*- 
"""UV-CDAT - Axes and grids (creating, modifying, axes ...)"""

# Imports
import numpy as N, cdms2, MV2

# -----------------------------------------------------------------------------------------------------------
# Longitude axis creation
# - base
lon = cdms2.createAxis([-5.,-4.,-3.],id='lon')
lon.long_name = 'Longitude'
lon.units = 'degree_east'
# - add lon.axis='X' and lon.modulo = '360.'
lon.designateLongitude()

# Latitude
lat = cdms2.createAxis([46.,47.,48.],id='lat')
lat.long_name = 'Latitude'
lat.units = 'degree_north'
lat.designateLatitude() # lat.axis = 'Y'

# Depth
depth = cdms2.createAxis([-200.,-100.,-0.],id='depth')
depth.long_name = 'Depth'
depth.units = 'm'
depth.designateLevel() # depth.axis = 'Z'

# Time
# - creation
time = cdms2.createAxis([0.,1.,2.],id='time')
time.long_name = 'Time'
time.units = 'days since 2006-08-01'
time.designateTime() # time.axis = 'T'
# - check
ctime = time.asComponentTime()
print ctime,ctime[1].day
#  -> [2006-8-1 0:0:0.0, 2006-8-2 0:0:0.0, 2006-8-3 0:0:0.0] 2
rtime = time.asRelativeTime()
print rtime,rtime[1].value
#  -> [0.00 days since 2006-08-01, 1.00 days since 2006-08-01,
#      2.00 days since 2006-08-01] 1.0


# Now, we create a variable using these axes.
#- straightforward method
temp1 = cdms2.createVariable(N.ones((3,3,3,3)),typecode='f',id='temp',
    fill_value=1.e20,axes=[time,depth,lat,lon],copyaxes=0,
    attributes=dict(long_name='Temperature',units='degC'))
print cdms2.isVariable(temp1)
# - Remark
print cdms2.createVariable is MV2.array
#  -> True         (These are the same functions !)
# - other methos
#   . initialization
temp2 = MV2.array(N.ones((3,3,3,3))).astype('f')
#   . attributes
temp2.id = 'temp'
temp2.long_name = 'Temperature'
temp2.units = 'degC'
temp2.set_fill_value(1.e20) # <=> temp2.setMissing(1.e20)
#   . axes
temp2.setAxisList([time,depth,lat,lon])
#   . or for example for each axis individually
temp2.setAxis(1,depth)

# Selection as for file
print temp2(time=("2006-08-01", "2006-08-03", "co")).shape

# The grid itself
# - get it
grid = temp2.getGrid()
# - check axes
grid.getLongitude() is temp2.getLongitude() is grid.getAxis(1) is temp2.getAxis(3)
# - create it!
grid2 = cdms2.createGenericGrid(lat,lon)
# - set it
temp2.setGrid(grid2)
print cdms2.isGrid(grid2)
