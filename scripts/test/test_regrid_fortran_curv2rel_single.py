"""Test the fortran function :f:func:`curv2rel_single`"""
from utils import np, plt, meshcells
from vacumm.fortran.interp import curv2rel_single

# %% Input grid
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
relpos2index = lambda fi, fj, nyi: fj * nyi + fi
ii, jj = np.meshgrid(np.arange(nxi)+.5, np.arange(nyi)+.5)
iib, jjb = meshcells(ii, jj)
zzi = relpos2index(ii, jj, nyi)
zzbi = relpos2index(iib, jjb, nyi)

# %% Input random points
np.random.seed(0)
npt = 100
xxo = np.random.random(npt)*(xxbi.max()-xxbi.min()) + xxbi.min()
yyo = np.random.random(npt)*(yybi.max()-yybi.min()) + yybi.min()
xxo = np.concatenate((xxo, xxi.ravel()))
yyo = np.concatenate((yyo, yyi.ravel()))
npt = xxo.size

# %% Convert to relative indices
zzo = np.zeros(npt)-1
for i, (xo, yo) in enumerate(zip(xxo, yyo)):
    p, q = curv2rel_single(xxbi, yybi, xo, yo)
#    print xo, yo, '|', p, q, (q-1) * nxi + p - 1
    if q >= 0:
        zzo[i] = relpos2index(p-1, q-1, nyi)

# %% Plot
vmin = np.nanmin(zzbi)
vmax = np.nanmax(zzbi)
plt.pcolor(xxbi, yybi, zzi, vmin=vmin, vmax=vmax)
plt.scatter(xxo, yyo, c=zzo, vmin=vmin, vmax=vmax, s=80)
plt.grid()
plt.title('curv2rel_single')
