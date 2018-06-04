"""Test class :func:`~vacumm.misc.grid.kriging.SimpleCloudkriger` with observation errors"""

from vcmq import SimpleCloudKriger, N, P, variogram_model, code_file_name, OrdinaryCloudKriger

# Kriging params
sill = 12.
range = 40

# Input
xi = [20., 26.,  50., 70]
yi = [20, 26., 70., 50]
zi = [2., 6., 15., 4.]
ei = [.3, .6, 3., .3]

# Output positions
nx = ny = 101
xg = N.linspace(1, 100, 101)
yg = N.linspace(1, 100, 101)

# Some inits
xi = N.array(xi)
yi = N.array(yi)
zi = N.array(zi)
ei = N.array(ei)
zi -= zi.mean()
xxg, yyg = N.meshgrid(xg, yg)
xo = xxg.ravel()
yo = yyg.ravel()
vgm = variogram_model('linear', n=0, s=sill, r=range)

# Setup the krigers
sck = SimpleCloudKriger(xi, yi, zi, vgf=vgm)#, farvalue=0)
scke = SimpleCloudKriger(xi, yi, zi, vgf=vgm, e=ei)#, farvalue=0)
#sck = OrdinaryCloudKriger(xi, yi, zi, vgf=vgm)
#scke = OrdinaryCloudKriger(xi, yi, zi, vgf=vgm, e=ei)

# Interpolate
zo = sck(xo, yo)
zoe = scke(xo, yo)
zzg = zo.reshape(ny, nx)
zzge = zoe.reshape(ny, nx)


# Plot

vmin = min(zi.min(), zo.min(), zoe.min())
vmax = max(zi.max(), zo.max(), zoe.max())
vmax = max(abs(vmin), abs(vmax))
vmin = -vmax
cmap = 'cmocean_tempo'
cmap = 'cmocean_delta'
kw = dict(vmin=vmin, vmax=vmax)
P.figure(figsize=(6, 5.5))

ax = P.subplot(221)
P.pcolormesh(xxg, yyg, zzg, cmap=cmap, **kw)
P.colorbar()
P.scatter(xi, yi, c=zi, s=100, cmap=cmap, **kw)
P.title('Without obs error')

P.subplot(222, sharex=ax, sharey=ax)
P.pcolormesh(xxg, yyg, zzge, cmap=cmap, **kw)
P.colorbar()
P.scatter(xi, yi, c=zi, s=100, cmap=cmap, **kw)
#P.axis('image')
P.title('With obs error')

P.subplot(223, sharex=ax, sharey=ax)
vmax = abs(zzge-zzg).max()
vmin = -vmax
P.pcolormesh(xxg, yyg, zzge-zzg, cmap=cmap, vmin=vmin, vmax=vmax)
P.colorbar()
P.scatter(xi, yi, c='k', s=100, cmap=cmap)
#P.axis('image')
P.title('With minus without')

P.subplot(224, sharex=ax, sharey=ax)
P.scatter(xi, yi, c=ei, s=100, cmap='cmocean_amp', vmin=0)
P.colorbar()
P.title('Observation errors')

P.tight_layout()
P.savefig(code_file_name(ext='.png'))
P.show()
P.close()
