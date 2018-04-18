"""Test external griddata and natgrid tools"""

from vcmq import P, N
from nat import Natgrid



# Generate data
xr = N.arange(20.)
yr = N.arange(10.)
xxr, yyr = N.meshgrid(xr, yr)
zzr = (N.sin(xxr*N.pi/6)*N.sin(yyr*N.pi/6) + \
    N.exp(-((xxr-7.)**2+(yyr-7.)**2)/4.**2))*100.
zzr -= zzr.min()
vminmax = dict(vmin=zzr.min(), vmax=zzr.max())
ij = N.unique((N.random.rand(50)*zzr.size).astype('i'))
xi, yi, zi = xxr.flat[ij], yyr.flat[ij], zzr.flat[ij]

ij2 = N.unique((N.random.rand(50)*zzr.size).astype('i'))
xo2, yo2 = xxr.flat[ij2], yyr.flat[ij2]


# CDAT natgrid
I = Natgrid(xi, yi, xr, yr)
I.ext = 1
I.igr = 1
zzoc = I.rgrd(zi).T

# Matplotlib / R. Kern
from matplotlib.delaunay import Triangulation
tri = Triangulation(xi,yi)
I = tri.nn_extrapolator(zi)
zzok = I(xxr,yyr)


# Plot
P.figure(figsize=(4, 11))
P.subplot(311)
P.pcolormesh(zzr, **vminmax)
P.title('Original')
stdref = zzr.std()
P.subplot(312)
P.pcolormesh(zzoc, **vminmax)
P.title('CDAT Natgrid: err=%i%%'%((zzr-zzoc).std()*100/stdref))
P.subplot(313)
P.pcolormesh(zzok, **vminmax)
P.title('MPL/R. Kern : err=%i%%'%((zzr-zzok).std()*100/stdref))



