"""Test function :func:`~vacumm.misc.grid.kriging.krig`"""

npi = 500
npo = 200

from vcmq import P, savefigs, code_file_name
from vacumm.misc.grid.kriging import gridded_gauss3, random_gauss3, random_points, krig

# Generate random field
xg, yg, zzg = gridded_gauss3()
xi, yi, zi = random_gauss3(np=npi)

# Interpolate to random points
xo, yo = random_points(np=npo)
zo = krig(xi, yi, zi, xo, yo)

# Plot
# - source data
axis = [xg.min(), xg.max(), yg.min(), yg.max()]
kw = dict(vmin=zzg.min(), vmax=zzg.max())
kwim = dict(extent=axis, interpolation='bilinear', origin='lower', alpha=.2, **kw)
kwsc = dict(lw=0.2, **kw)
P.figure(figsize=(6, 3.5))
P.subplot(121)
P.title('Source field')
P.imshow(zzg, **kwim)
P.scatter(xi, yi, c=zi, **kwsc)
P.axis(axis)
# - interpolated data
P.subplot(122)
P.title('Interpolated points')
P.imshow(zzg, **kwim)
P.scatter(xo, yo, c=zo, marker='s', s=30, **kwsc)
P.axis(axis)
P.tight_layout()
savefigs(code_file_name(), verbose=False)
P.close()
