"""Test the fortran function :f:func:`cellave1d`"""
from utils import np, plt
from vacumm.fortran.interp import cellave1d

nx = 17
nyi = 20
nyo = 12
u, v = np.mgrid[-3:3:nx*1j, -3:3:nyi*1j]-2
vari = np.asarray(u**2+v**2)
yib = np.linspace(-1000., 0., nyi+1)
yob = np.linspace(-1200, 200, nyo+1)
vari[int(nx/3):int(2*nx/3), int(nyi/3):int(2*nyi/3)] = np.nan
xb = np.arange(nx+1)
xxib, yyib = np.meshgrid(xb, yib)
xxob, yyob = np.meshgrid(xb, yob)

varoa = cellave1d(vari, yib, yob, conserv=0, extrap=0)
varoe = cellave1d(vari, yib, yob, conserv=0, extrap=2)
varoc = cellave1d(vari, yib, yob, conserv=1, extrap=0)

sumi = np.nansum(vari*np.resize(np.diff(yib), vari.shape), axis=1)
sumo = np.nansum(varoc*np.resize(np.diff(yob), varoc.shape), axis=1)
np.testing.assert_allclose(sumi, sumo)

kw = dict(vmin=np.nanmin(vari), vmax=np.nanmax(vari))
axlims = [xb[0], xb[-1], yob[0], yob[-1]]
fig, [[ax0, ax1], [ax2, ax3]] = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
ax0.pcolor(xxib, yyib, vari.T, **kw)
ax0.axis(axlims)
ax0.set_title('Original')
ax1.pcolor(xxob, yyob, varoa.T, **kw)
ax1.axis(axlims)
ax1.set_title('Cell average')
ax2.pcolor(xxob, yyob, varoe.T, **kw)
ax2.axis(axlims)
ax2.set_title('With extrap')
ax3.pcolor(xxob, yyob, varoc.T, **kw)
ax3.axis(axlims)
ax3.set_title('Conservative')
plt.tight_layout()
