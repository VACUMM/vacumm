"""Test the fortran function :f:func:`interp1dx`"""
from utils import np, plt, meshcells
from vacumm.fortran.interp import interp1dx

nx = 17
nyi = 10
nyo = 30
u, v = np.mgrid[-3:3:nx*1j, -3:3:nyi*1j]-2
vari = np.asarray(u**2+v**2)
yi = np.linspace(-1000., 0., nyi)
yo = np.linspace(-1200, 100, nyo)
vari[int(nx/3):int(2*nx/3), int(nyi/3):int(2*nyi/3)] = np.nan
x = np.arange(nx)
dyi = (yi[1]-yi[0])*0.49
yyi = np.resize(yi, vari.shape)+np.random.uniform(-dyi, dyi, vari.shape)
xxib, yyib = meshcells(x, yyi.T)
xxob, yyob = meshcells(x, yo)

varon = interp1dx(vari, yyi, yo, 0, extrap=0)
varol = interp1dx(vari, yyi, yo, 1, extrap=0)
varoh = interp1dx(vari, yyi, yo, 3, extrap=0)

kw = dict(vmin=np.nanmin(vari), vmax=np.nanmax(vari))
axlims = [x[0], x[-1], yo[0], yo[-1]]
fig, [[ax0, ax1], [ax2, ax3]] = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))
ax0.pcolor(xxib, yyib, vari.T, **kw)
ax0.axis(axlims)
ax0.set_title('Original')
ax1.pcolor(xxob, yyob, varon.T, **kw)
ax1.axis(axlims)
ax1.set_title('Nearest1dx')
ax2.pcolor(xxob, yyob, varol.T, **kw)
ax2.axis(axlims)
ax2.set_title('Linear1dx')
ax3.pcolor(xxob, yyob, varoh.T, **kw)
ax3.axis(axlims)
ax3.set_title('Hermit1dx')
plt.tight_layout()
