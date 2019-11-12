"""Test the fortran function :f:func:``"""
from utils import np, plt, meshcells
from vacumm.fortran.interp import cellave1dxx

nx = 17
nyi = 20
nyo = 12
u, v = np.mgrid[-3:3:nx*1j, -3:3:nyi*1j]-2
vari = np.asarray(u**2+v**2)
yib = np.linspace(-1000., 0., nyi+1)
yob = np.linspace(-1200, 200, nyo+1)
vari[int(nx/3):int(2*nx/3), int(nyi/3):int(2*nyi/3)] = np.nan
xb = np.arange(nx+1)
yyib = np.resize(yib, (nx, nyi+1))
yyob = np.resize(yob, (nx, nyo+1))
dyi = (yib[1]-yib[0])*0.49
dyo = (yob[1]-yob[0])*0.49
yyib += np.random.uniform(-dyi, dyi, yyib.shape)
yyob += np.random.uniform(-dyo, dyo, yyob.shape)

varoa = cellave1dxx(vari, yyib, yyob, conserv=0, extrap=0)
varoe = cellave1dxx(vari, yyib, yyob, conserv=0, extrap=2)
varoc = cellave1dxx(vari, yyib, yyob, conserv=1, extrap=0)

sumi = np.nansum(vari*np.resize(np.diff(yyib, axis=1), vari.shape), axis=1)
sumo = np.nansum(varoc*np.resize(np.diff(yyob, axis=1), varoc.shape), axis=1)
np.testing.assert_allclose(sumi, sumo)

yyib = meshcells(yyib, axis=0).T
yyob = meshcells(yyob, axis=0).T
kw = dict(vmin=np.nanmin(vari), vmax=np.nanmax(vari))
axlims = [xb[0], xb[-1], yob[0], yob[-1]]
fig, [[ax0, ax1], [ax2, ax3]] = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
ax0.pcolor(xb, yyib, vari.T, **kw)
ax0.axis(axlims)
ax0.set_title('Original')
ax1.pcolor(xb, yyob, varoa.T, **kw)
ax1.axis(axlims)
ax1.set_title('Cell average')
ax2.pcolor(xb, yyob, varoe.T, **kw)
ax2.axis(axlims)
ax2.set_title('With extrap')
ax3.pcolor(xb, yyob, varoc.T, **kw)
ax3.axis(axlims)
ax3.set_title('Conservative')
plt.tight_layout()
