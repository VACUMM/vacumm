"""Test the fortran function :f:func:`nearest2d`"""
from utils import np, plt, rotate_grid, plot_grid, meshcells
from vacumm.fortran.interp import nearest2d

# %% Input grid
xxi, yyi = rotate_grid(np.arange(5), np.arange(4), 30)
vari = np.resize(yyi, (20, ) + yyi.shape)
nb = 10
xxbi, yybi = meshcells(xxi, yyi)

# %% Output grid
xxo, yyo = rotate_grid(np.linspace(0, 6, 50)-1, np.linspace(0, 4, 35)+1., -20)
xxbo, yybo = meshcells(xxo, yyo)

# %% Nearest
varo = nearest2d(vari, xxi, yyi, xxo, yyo, nb)

# %% Plot
vmin = varo.min()
vmax = varo.max()
plt.figure(figsize=(8, 4))
plt.subplot(121, aspect=1)
plt.pcolormesh(xxbi, yybi, vari[0], vmin=vmin, vmax=vmax)
plot_grid(xxbo, yybo)
plt.title('original')
plt.subplot(122, aspect=1)
plt.pcolormesh(xxbo, yybo, varo[0], vmin=vmin, vmax=vmax)
plot_grid(xxbi, yybi)
plt.title('nearest2d')
plt.axis('image')
plt.tight_layout()
