"""Test class :func:`~vacumm.misc.grid.kriging.OrdinaryCloudKriger`"""

npi = 50

from vcmq import P, savefigs, code_file_name
from vacumm.misc.grid.kriging import (gridded_gauss3, random_gauss3,
    random_points, OrdinaryCloudKriger)

# Generate random field
xg, yg, zzg = gridded_gauss3()
xi, yi, zi = random_gauss3(np=npi)

# Setup kriger
ock = OrdinaryCloudKriger(xi, yi, zi, n=0) # nugget = 0

# Interpolate to grid and get error
xo, yo = P.meshgrid(xg[::4], yg[::4])
xo.shape = yo.shape = xo.size
zo, zoe = ock(xo, yo, geterr=True)

# Relative error
zoemax = ock.variogram_func(P.inf)
zoe /= zoemax
zoe *= 100

# Variogram evolution
rr = P.linspace(0, P.sqrt(xg.ptp()**2+yg.ptp()**2), 100)
vv = ock.variogram_func(rr)

# Plot
P.figure(figsize=(8, 8))
axis = [xg.min(), xg.max(), yg.min(), yg.max()]
kwmm = dict(vmin=zi.min(), vmax=zi.max())
kwsc = dict(lw=0.2, **kwmm)
ax = P.subplot(221, aspect=1)
P.imshow(zzg,  extent=axis, interpolation='bilinear', origin='lower', alpha=.2, **kwmm)
P.scatter(xi, yi, c=zi, s=20, lw=0.2, **kwmm)
P.axis(axis)
P.title('Input points')
P.subplot(222, aspect=1)
P.scatter(xo, yo, c=zo, s=40, lw=0.2, **kwmm)
P.axis(axis)
P.title('Interpolated points')
P.subplot(223, aspect=1)
P.scatter(xo, yo, c=zoe, s=40, lw=0.2, vmin=0, vmax=100)
P.axis(axis)
P.colorbar(shrink=.8)
P.title('Relative interpolation error')
P.subplot(224)
P.plot(rr, vv)
P.xlabel('Distance')
P.ylabel('Error')
P.title('Fitted variogram')
P.tight_layout()
savefigs(code_file_name(), verbose=False)
P.close()
