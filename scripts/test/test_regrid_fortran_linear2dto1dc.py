"""Test fortran function :f:func:`linear2dto1dc`"""

from utils import plt, np, rotate_grid
from vacumm.fortran.interp import linear2dto1dc


# %% Input grid and data
nxy = 15
xi = np.arange(nxy*1.)
yi = np.arange(nxy*1.)
xxi, yyi = rotate_grid(xi, yi, 30)
zzi = np.array(yyi)
zzi[int(nxy*0.3):int(nxy*0.8), int(nxy*0.3):int(nxy*0.8)] = np.nan
zzi.shape = 1, nxy, nxy

# %% Output positions
no = 1000
xo = np.random.uniform(-nxy/4., nxy+nxy/4., no)
yo = np.random.uniform(-nxy/4., nxy+nxy/4., no)

# %% Interpolate
zo = linear2dto1dc(xxi, yyi, zzi, xo, yo)

# %% Plot
kw = dict(vmin=np.nanmin(zzi), vmax=np.nanmax(zzi))
plt.figure(figsize=(6, 6))
plt.subplot(111, aspect=1)
plt.contourf(xxi, yyi, zzi[0], **kw)
plt.scatter(xxi, yyi, c='k')
plt.scatter(xo, yo, c=zo[0], s=50, **kw)
plt.title('linear2dto1dc')
