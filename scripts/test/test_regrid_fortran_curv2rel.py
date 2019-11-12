"""Test the fortran function :f:func:`curv2rel`"""
from utils import np, plt, meshcells
from vacumm.fortran.interp import curv2rel


# Input grid
x0, y0 = 0., 2.
nxi = 3
nyi = 3
dxi = (2., 2.)
dyi = (-2., 2.)
xxi = np.zeros((nyi, nxi))
xxi[0] = x0 + np.arange(nxi) * dxi[0]
for j in range(1, nyi):
    xxi[j] = xxi[j-1] + dyi[0]
yyi = np.zeros((nyi, nxi))
yyi[:, 0] = y0 + np.arange(nyi) * dyi[1]
for j in range(1, nyi):
    yyi[:, j] = yyi[:, j-1] + dxi[1]
xxbi, yybi = meshcells(xxi, yyi)
relpos2index = lambda fi, fj: fj * nyi + fi
ii, jj = np.meshgrid(np.arange(nxi)+.5, np.arange(nyi)+.5)
iib, jjb = meshcells(ii, jj)
zzi = jj * nyi + ii  # relpos2index(ii, jj)
zzbi = jjb * nyi + iib  # ny+1? relpos2index(iib, jjb)

# Input random points
np.random.seed(0)
npt = 100
xxo = np.random.random(npt)*(xxbi.max()-xxbi.min()) + xxbi.min()
yyo = np.random.random(npt)*(yybi.max()-yybi.min()) + yybi.min()
xxo = np.concatenate((xxo, xxi.ravel()))
yyo = np.concatenate((yyo, yyi.ravel()))

# Convert to relative indices
pp, qq = curv2rel(xxbi, yybi, xxo, yyo)
bad = pp < 0
zzo = (qq-1) * nyi + pp-1  # relpos2index(pp-1, qq-1)
zzo[bad] = -1

# Plot
vmin = zzbi.min()
vmax = zzbi.max()
plt.pcolor(xxbi, yybi, zzi, vmin=vmin, vmax=vmax)
plt.scatter(xxo, yyo, c=zzo, vmin=vmin, vmax=vmax, s=80)
plt.grid()
plt.title('curv2rel')
