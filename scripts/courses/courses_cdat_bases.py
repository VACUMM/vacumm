#!/usr/bin/env python
# -*- coding: utf8 -*- 
"""UV-CDAT - First commands (MV2 variables, grid, links with numpy arrays)"""

# Imports
import matplotlib
matplotlib.use('qt4agg')
import cdms2, MV2, sys
import numpy as N
from vacumm.misc.axes import create_lon,create_lat, create_time, create_dep, create_depth, islon, isdep, istime
from vacumm.misc.plot import map2 as map, hov2, section2
from vacumm.config import data_sample

# -----------------------------------------------------------------------------------------------------------
# ---- From a NetCDF file ...
print 10*'-'+' ... From a NetCDF File ... '+10*'-'
f = cdms2.open(data_sample('mars3d.xy.nc'))
ncarr=f('temp', lat=slice(5,-1), lon=(-6.2, -4))
f.close()
print 'Array dimension :', ncarr.shape
print 'Array type - ncarr - :',type(ncarr)
# print 'Longitude :',ncarr.getLongitude()
# print 'Latitude :',ncarr.getLatitude()

print '2D Grid: ',ncarr.getGrid()

# ---- Copy or not copy !!!!
print 'An example: ',ncarr[50,50]

# ==> Practice: test and understant the differences WITH and WITHOUT copy.
newnccarr = ncarr.clone() # WITH COPY
newnccarr = ncarr # WITHOUT COPY
print '+ nccarr: ',ncarr[50,50]
print '+ newnccarr: ',newnccarr[50,50]
newnccarr[50,50] = 25.
print '+ newnccarr modified: ',newnccarr[50,50]
print '+ nccarr after modification of newnccarr: ',ncarr[50,50]

# Back to numpy array ...
nparr=ncarr.getValue()
print 'Array type - nparr - :',type(nparr)

sys.exit() # End of the run
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

# ---- From a Numpy array ...
print 10*'-'+' ... From a numpy array ... '+10*'-'

# ---- A 2-D numpy array ...
print 10*'-'+' Numpy array '+10*'-'
arr = N.random.rand(10,15)
print 'Array dimension :', arr.shape
print 'Array type - arr - :',type(arr)
# print arr

# ---- From numpy array to MV2 (Masked array) ...
print 10*'-'+' MV2 array '+10*'-'
marr = MV2.array(arr)
print 'Array type - marr - :',type(marr)
# print marr

# Mask value lower than 0.5 ...
# More details on masked array operations: http://docs.scipy.org/doc/numpy/reference/routines.ma.html
marr=MV2.masked_where((marr < 0.5),marr)
# print marr.mask
# print marr
print 'Array dimension :', marr.shape
print 'Fill value: ',marr.fill_value

# Back to numpy array ...
nparr=marr.getValue()
# print nparr

# ---- A Cdms2 object = MV2 + Axes + Attributes ...
print 10*'-'+' CDMS2 '+10*'-'
# cdms2.createVariable()

# Geographic axis creation
# - longitude: changing 'lon' id to 'longitude'
ax1 = create_lon(N.arange(10)-7.,id='longitude')
# - latitude
ax2 = create_lat(N.arange(15)*.5+44.)

# ==> Practice: Create cdarr with depth/lat axes. - see doc Vacumm -
# # - depth
# ax1 = create_dep(N.linspace(-1000,0,10))
# # - latitude
# ax2 = create_lat(N.arange(15)*.5+44.)


# ==> Practice: Create cdarr with time/lat axes. - see doc Vacumm -
# # - time
# ax1 = create_time(N.arange(10.),
    # 'days since 2006-10-01',long_name='Mon axe de temps')
# # - latitude
# ax2 = create_lat(N.arange(15)*.5+44.)

# - cdms2 variable creation
cdarr = cdms2.createVariable(marr,axes=[ax1,ax2],id='test',attributes=dict(long_name='the main test',units='$m^3 s^{-1}$'))
print 'Array type - cdarr - :',type(cdarr)
print 'Variable :',cdarr.id
# print 'Longitude :',cdarr.getLongitude()
# print 'Latitude :',cdarr.getLatitude()
# print cdarr

# - Result snapshot
if islon(cdarr.getAxis(0)):
    map(cdarr, contour=False)
if isdep(cdarr.getAxis(0)):
    section2(cdarr, contour=False)
if istime(cdarr.getAxis(0)):
    hov2(cdarr, contour=False)
