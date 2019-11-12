"""Test fortran function :f:func:`linear2dto1d`"""

from utils import plt, np
from vacumm.fortran.interp import linear2dto1d


# %% Input grid and data
nxy = 10
xi = np.arange(nxy*1.)
yi = np.arange(nxy*1.)
xxi, yyi = np.meshgrid(xi, yi)
zi = np.array(yyi)
zi[int(nxy*0.3):int(nxy*0.6), int(nxy*0.3):int(nxy*0.6)] = np.nan
zi.shape = 1, nxy, nxy

# %% Output positions
no = 500
xo = np.random.uniform(-nxy/4., nxy+nxy/4., no)
yo = np.random.uniform(-nxy/4., nxy+nxy/4., no)

# %% Interpolate
zo = linear2dto1d(xi, yi, zi, xo, yo)

# %% Plot
kw = dict(vmin=np.nanmin(zi), vmax=np.nanmax(zi))
plt.figure(figsize=(6, 6))
plt.subplot(111, aspect=1)
plt.contourf(xxi, yyi, zi[0], **kw)
plt.scatter(xxi, yyi, marker='o', color='k', s=3)
plt.scatter(xo, yo, c=zo[0], s=50, **kw)
plt.title('linear2dto1d')
