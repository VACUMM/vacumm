"""Test the fortran function :f:func:`closest2d`"""
from utils import np, plt, rotate_grid, plot_grid, meshcells
from vacumm.fortran.interp import closest2d


# %% Input grid
xxi, yyi = rotate_grid(np.arange(5), np.arange(4), 0)

# %% Input random points
np.random.seed(0)
npt = 100
xxo = np.random.random(npt)*(xxi.max()-xxi.min()) + xxi.min()
yyo = np.random.random(npt)*(yyi.max()-yyi.min()) + yyi.min()

# %% Closest and plot
for xo, yo in zip(xxo, yyo):
    i, j = closest2d(xxi, yyi, xo, yo, nogeo=0)
    plt.plot([xo, xxi[j-1, i-1]], [yo, yyi[j-1, i-1]], 'k')
xxib, yyib = meshcells(xxi, yyi)
plot_grid(xxib, yyib)
plt.scatter(xxi, yyi, s=30, marker='o', c='k')
plt.scatter(xxo, yyo, c='r', s=30, zorder=11, marker='s')
plt.title('closest2d')
plt.axis('image')
