# -*- coding: utf8 -*-
import numpy as N, cdms2, MV2, pylab as P
from vcmqm import meshbounds, interp2d, add_grid

# Input variables
xi = N.arange(20.)
yi = N.arange(10.)
xxi, yyi = N.meshgrid(xi, yi)
xib, yib = meshbounds(xi, yi)
vari = (N.sin(xxi*N.pi/6)*N.sin(yyi*N.pi/6) +
    N.exp(-((xxi-7.)**2+(yyi-7.)**2)/4.**2))*100.
vminmax=dict(vmin=vari.min(), vmax=vari.max())
vari = cdms2.createVariable(vari)
vari.setAxis(-2, cdms2.createAxis(yi))
vari.setAxis(-1, cdms2.createAxis(xi))
vari[3:4, 3:7] = MV2.masked

# Output grid
xo = cdms2.createAxis(N.linspace(-3., 23., 70))
yo = cdms2.createAxis(N.linspace(-3., 13., 40))

# Interpolation
# - bilinear
varob = interp2d(vari, (xo, yo), method='bilinear')
# - nearest
varon = interp2d(vari, (xo, yo), method='nearest')

# Plot
xob, yob = meshbounds(xo[:], yo[:])
lims = [xob.min(), xob.max(), yob.min(), yob.max()]
# -
P.figure(figsize=(5.5, 7))
P.subplots_adjust(bottom=.07, hspace=.35)
P.subplot(311)
P.pcolor(xib, yib, vari, **vminmax)
P.axis(lims)
add_grid((xo, yo), linewidth=.3)
P.title('Original')
# -
P.subplot(312)
P.pcolor(xob, yob, varob, **vminmax)
P.axis(lims)
add_grid((xi, yi), linewidth=.3)
P.title('Bilinear')
# -
P.subplot(313)
P.pcolor(xob, yob, varon, **vminmax)
P.axis(lims)
add_grid((xi, yi), linewidth=.3)
P.title('Nearest')

