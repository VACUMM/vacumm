"""Test function :func:`~vacumm.misc.grid.kriging.OrdinaryKriger` in parallel mode"""
 
npi = 2000
npo = 200
npmax = 500
nproc = 2


from vcmq import P, savefigs, code_base_name
from vacumm.misc.grid.kriging import gridded_gauss3, random_gauss3, random_points, OrdinaryKriger
from time import time

# Random and gridded input fields
xg, yg, zzg = gridded_gauss3()
xi, yi, zi = random_gauss3(npts=npi)

# Init kriger
kriger = OrdinaryKriger(xi, yi, zi, npmax=npmax, nproc=nproc)

# Output points
xo, yo = random_points(np=npo)

# Interpolate
t0 = time()
zo = kriger(xo, yo)
#print '%.2fs'%(time()-t0)

# Plot
P.figure(figsize=(6, 8))
axis = [xg.min(), xg.max(), yg.min(), yg.max()]
kw = dict(vmin=zi.min(), vmax=zi.max())
kwim = dict(extent=axis, interpolation='bilinear', origin='lower', alpha=.2, **kw)
kwsc = dict(lw=0.2, **kw)
P.subplot(211)
P.imshow(zzg,  **kwim)
P.scatter(xi, yi, c=zi, s=20, **kwsc)
P.title('Input points')
P.subplot(212)
P.imshow(zzg, **kwim)
P.scatter(xo, yo, c=zo, s=40, **kwsc)
P.title('Interpolated points')
P.tight_layout()
savefigs(code_base_name(), verbose=False)
P.close()

