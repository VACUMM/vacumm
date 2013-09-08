#!/usr/bin/env python
# -*- coding: utf8 -*- 
""" UV-CDAT - Time (time, datetime, cdtime, ...) """

# Imports
import matplotlib
matplotlib.use('qt4agg')

# import time, datetime, cdtime
import numpy as N
from vacumm.config import data_sample
# from StringIO import StringIO
import cdms2 as cdms
from vacumm.misc.atime import strftime, ch_units, strptime
import scipy.io as SC
import matplotlib.pyplot as mp

# Inits
drifterfile = "drifter.txt"
adcpfile = "adcp.txt"


# -----------------------------------------------------------------------------------------------------------
# ---- Read an ASCII file ----
print 10*'-'+' ASCII file '+10*'-'
# We load a cdms variable
f = cdms.open(data_sample('mars3d.t.nc'))
temp = f('temp')
u =  f('u')
v =  f('v')
f.close()

# We get the Time
time =  temp.getTime() # axis
ctime = time.asComponentTime() # time cdtime.comptime()

# Create the file with the format : YYYY/MM/DDZHH:MM TEMP U V
f = open('misc.io.ascii.1.dat', 'w')
f.write('# Ligne de commentaire\n')
for it in xrange(len(temp)):
    t = strftime('%Y/%m/%dZ%H:%M',  ctime[it])
    f.write('%s %.4f %f %f\n' % (t, temp[it], u[it], v[it]))
f.close()
# => Practice: Read just a small part of the file (10 last lines and/or 10 first lines)

# Quick check (two first lines)
f = open('misc.io.ascii.1.dat')
print ''.join(f.readlines()[:3])
f.close()

# Quick writing using Numpy
# - Creation
time_units = 'hours since %s'%ctime[0]
newtime = ch_units(time, time_units)[:]
data = N.array([newtime,temp.filled(999.), 
    u.filled(999.), v.filled(999.)],copy=0)
f = open('misc.io.ascii.2.dat', 'w')
f.write('# Ligne de commentaire\n')
N.savetxt(f, data.transpose(), fmt='%.3f', delimiter='\t')
#   Note : we can enter f or the file name in N.savetxt()
# - Check
f = open('misc.io.ascii.2.dat')
print ''.join(f.readlines()[:2])
f.close()

# Advanced reading
# - Time converter
def convtime(s):
    return strptime(s, '%Y/%m/%dZ%H:%M').torel(time_units).value

# - uncomplete load
tt, uu, vv = N.loadtxt('misc.io.ascii.1.dat',  comments='#', 
    usecols=[0, 2, 3], converters={0:convtime}, unpack=True)
# - check
print tt[0], uu[0], vv[0]

# => Practice: Explore N.genfromtxt and try to read misc.io.ascii.1.dat with this function.

# -----------------------------------------------------------------------------------------------------------
# ---- Read/Write an binary MAT file (Matlab) ----
print 10*'-'+' .MAT file '+10*'-'

# -- Example on ASPEX (Scanfish) data
# - Load a mat file
mat = SC.loadmat(data_sample('SectionG_simplified.mat'))
T=mat['T']
lat=mat['lat']
P=mat['P']

print '* Close the figure to execute following actions ...'
# - Figure
mp.figure()
mp.scatter(lat,-P,s=20,c=T,edgecolors='none')
mp.colorbar()
mp.xlabel('Latitude')
mp.ylabel('Depth')
mp.ylim((-110,0))
mp.title('Scanfish ASPEX Section G')
mp.grid()
mp.show()

# Practice: Write only the salinity in another .MAT file.

# -----------------------------------------------------------------------------------------------------------
# ---- Read/Write an binary npy file (Numpy) ----
print 10*'-'+' .NPY file '+10*'-'
arr = N.ones(10)
# - Write
N.save('test.npy',arr)

# - Load
a1=N.load('test.npy')
print a1

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

