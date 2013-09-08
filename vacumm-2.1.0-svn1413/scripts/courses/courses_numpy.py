#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Les bases numpy"""

import numpy as N, scipy, scipy.io


# Créations de base
a = N.array([1, 2, 3])
print a.dtype
a = N.ones(3)
a = N.zeros(3, dtype='i')
a = N.zeros((2,2), '?')
print N.empty(2)
a = N.array(['AZERT','AZERT','AZERT'], dtype='|S5') 
a = N.zeros(2,dtype='f8,i4,a5')
print a.dtype
# -> CHANGER LE TYPE AVEC .ASTYPE()


# Plus évoluées
a = N.identity(2)
a = N.eye(2,3)
a = N.diag([1,2])
a = N.tri(2)
a = N.arange(5)
a = N.linspace(-5,5,11)
#a = N.r_[-5:5:11j]
#a = N.mgrid[-5:5:11j]
a, b = N.meshgrid([0,1,2],[3,4])
# -> ESSAYER MGRID AVEC 2D


# Générateurs
a = N.random.random((2,2)) # -> TESTER RANDINT, NORMAL, WEIBULL


# Copie (IMPORTANT)
a = N.array([1, 2, 3])
b = N.array(a)      # -> ESSAYER AVEC .COPY()
c = a
d = a[:2]
a[0] = 100
print b[0], c[0], d[0]


# Matrices
a2 = N.random.random((2,2))
m2 = N.matrix(a2)   # -> ESSAYER AVEC .ASMATRIX()
# -> UTILISER MAT()
# -> MULTIPLIER (CARRE)
print m2.I


# Type enregistrement
rec = N.rec.array(
    [(53,'Jean'),(22,'Pierre')],
    names='age,name')
# -> ESSAYER AVEC REC.FROMARRAYS()
print rec.age
print rec[0]
# -> ESSAYER SORT() AVEC ORDER=...


# Transformations
a = N.arange(6)
a.shape = 2,3                   # -> ESSAYER .RESHAPE()
print a
print a.ravel()
print a.T                       # -> ESSAYER .TRANSPOSE()
c = N.zeros((2,3,4,5))
print N.swapaxes(c,1,2).shape   # -> ESSAYER ROLLAXIS
a = N.diag([1.,2.,3.])
print N.fliplr(a)
print N.rot90(a)
# -> ESSAYER ROLL()


# Split
print N.split(N.arange(6),2)
print N.tile(b, [2,3])
print 'b',b
x = N.array([[1,2],[3,4]])
print N.repeat(x, [3,2], axis=0)


# Join
a = N.arange(6).reshape(2,3)
b = a+10
print N.concatenate((a,b), axis=1)
print N.column_stack((a,b))
print N.hstack((a,b))           # -> ESSAYER VSTACK()


# Extractions
a = N.arange(200).reshape((10,20))
print a[2::3, ::-1]
print a[1, slice(None, -1)]
print a[:, [2, 5]]
a = N.ones(3)
print a[a>0]
print N.argmax([3,5,1,-3]) # -> TROUVER LE MAX AVEC CE RESULTAT

# Affectation
a = N.arange(6).reshape(2, 3)
N.place(a, a>2, [44, 55])
print a
N.put(a, [0,3], [100, 200])
print a
N.putmask(a, a<10, -a)
print a


# Lecture fichier 

# - ascii
x,y = N.loadtxt('courses_numpy.txt',comments='#', usecols=(1,2),unpack=True)
data = N.loadtxt('courses_numpy.txt',    usecols=(0,1,2),skiprows=2,
    dtype={'names':('date','temp','psal'), 'formats':('S10','f','f')})

print data[0]
import pylab
from datetime import datetime
def s2d(s):
    return pylab.date2num(datetime.strptime(s,'%Y-%m-%d'))
data = N.recfromtxt('courses_numpy.txt',skiprows=2,
    usecols=(0,1,2),converters={0:s2d}, dtype={'names':('date','temp','psal'),
    'formats':('f','f','f')})
print data.date
N.savetxt('courses_numpy2.txt.gz', data, fmt='%d')

# - binaire
N.save('courses_numpy.npy', data)
data2 = N.load('courses_numpy.npy')

# - matlab
vect = N.arange(10)
num = 5.
scipy.io.savemat('courses_numpy.mat', {'vect':vect})
# -> ESSAYER LOADMAT()


# Arithmetique et trigo
a = N.arange(3.)
b = N.array([10., 20, 20])
print a+b, N.add(a, b)
print b**a-10
a += 1
print a
a = N.random.randn(10,3)
a -= a.mean(axis=0) # 3D-2D=3D
print N.log(N.exp(2.))
print N.angle(1+1j, deg=True)
print N.clip([-0,5,10], 2, 8)
# -> ESSAYER MAXIMUM
print N.arctan2(1.,1.)/N.pi
print N.degrees(N.pi/2)


# Precision
a = N.array([0.1254, 158.23])
print N.round(a, decimals=2)    # -> ESSAYER DECIMALS=-1
print N.rint(a)
b = N.array([2.1, 2.9, -2.1, -2.9])
print N.ceil(b)                 # -> ESSAYER FLOOR() ET FIX()


# Calculs axiaux
a = N.arange(6.).reshape(2,3)*2
print a.sum(axis=1)
print N.diff(a, axis=1)
print N.gradient(a)
print N.trapz([1,2,3])
print N.convolve([1.,2,3],[1,2,1])
print N.interp(2.5,[1,2,3],[3,2,0])


# Statistiques
a = N.arange(12).reshape(3,4)
print a.max(axis=0)
print a.mean(axis=1)
print N.average(a,axis=1, weights = [2,1,1,0])
print a.std(axis=0,ddof=1)
a = N.random.randn(10,3)
print N.corrcoef(a,rowvar=0)
print N.cov(a,rowvar=0).shape
print N.correlate([1,2,3],[0,1,0.5],"same")
count, bins = N.histogram(a,3,normed=True)


# Tableaux masqués
a = N.arange(5)
print N.ma.masked_greater(a, 3)
b = N.ma.masked_values(a,2)
print b.mask
b.set_fill_value(-1)
print b.filled()    # -> PASSER UN ARGUMENT (NOMBRE)
c = N.ma.asarray(a)
print c.mask, c.mask is N.ma.nomask, N.ma.isMA(c)
print N.ma.getmaskarray(c)
a = N.ma.array([0,2,4,6,8], mask=[0,1,0,0,0])
b = N.ma.masked_where(a==8,a,copy=0) 
# -> CHANGER A ET VERIFIER B
print a[~a.mask]    # -> COMPARER AVEC .COMPRESSED()


# Nan et autres
a = N.array([0,3,N.nan,2])
print N.nanmax(a), N.nanargmax(a)
print N.nansum(a)
print N.isnan(a)
# -> COMPARER NAN ET NAN
a[0] = N.inf
# -> TESTER ISINF, ISNEGINF, ISFINITE
print N.nan_to_num(a)
a = N.array([0,3,N.nan])
print a.mean()
m = N.ma.masked_where(N.isnan(a),a)
print m.mean()


# Algèbre
a = N.array([[1, 0.], [0, 1]])
b = N.array([[4, 1.], [2, 2]])
print N.dot(a,b), N.inner(a,b)
x = N.linalg.solve(a, b)
U,S,VH = N.linalg.svd(a)
print N.allclose(N.linalg.pinv(a), N.matrix(a).I)
x,residues,rank,s =N.linalg.lstsq(a, b)
# VOIR SCIPY AUSSI


# FFT
y = N.sin(N.arange(100.)) + N.random.randn(100)
sp = N.fft.fft(y)
freq = N.fft.fftfreq(y.shape[0])
print N.allclose(y, N.fft.ifft(sp).real)


# Polynomes
x = N.array([0.0, 1.0, 2.0, 3.0,  4.0,  5.0])
y = N.array([0.0, 0.8, 0.9, 0.1, -0.8, -1.0])
z = N.polyfit(x, y, 3)







