"""Test the fortran function :f:func:`bilin`"""
from utils import np, plt, meshcells
from vacumm.fortran.interp import linear2d

# %% Input
nxi = 15
nyi = 10
u, v = np.mgrid[-3:3:nyi*1j, -3:3:nxi*1j]-2
vari = np.asarray(u**2+v**2)
xi = np.arange(nxi)
yi = np.arange(nyi)
vari[int(nyi*0.4), int(nxi*0.4):int(nxi*0.4)+2] = np.nan
xxib, yyib = meshcells(xi, yi)

# %% Output grid
nxo = 40
nyo = 25
xo = np.linspace(int(nxi*0.2), int(nxi*1.2), nxo)
yo = np.linspace(int(-nyi*0.2), int(nyi*0.8), nyo)
xxob, yyob = meshcells(xo, yo)

# %% Interp
vari.shape = (1, )+vari.shape
varo = linear2d(vari, xi, yi, xo, yo, 0)

# %% Plot
kw = dict(vmin=np.nanmin(vari), vmax=np.nanmax(vari))
axlims = [min(xi.min(), xo.min()), max(xi.max(), xo.max()),
          min(yi.min(), yo.min()), max(yi.max(), yo.max())]
plt.figure(figsize=(8, 4))
plt.subplot(121)
plt.pcolor(xxib, yyib, vari[0], **kw)
plt.axis(axlims)
plt.title('Original')
plt.subplot(122)
plt.pcolor(xxob, yyob, varo[0], **kw)
plt.axis(axlims)
plt.title('Bilinear')
plt.tight_layout()
