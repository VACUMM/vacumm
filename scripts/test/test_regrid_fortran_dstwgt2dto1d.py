"""Test fortran function :f:func:`dstwgt2dto1d`"""

from utils import plt, np
from vacumm.fortran.interp import dstwgt2dto1d

# %% Input grid and data
nxy = 15
xi = np.arange(nxy*1.)
yi = np.arange(nxy*1.)
xxi, yyi = np.meshgrid(xi, yi)
zi = np.array(yyi)
zi[int(nxy*0.3):int(nxy*0.8), int(nxy*0.3):int(nxy*0.8)] = np.nan
zi.shape = 1, nxy, nxy

# %% Output positions
no = 1000
xo = np.random.uniform(-nxy/4., nxy+nxy/4., no)
yo = np.random.uniform(-nxy/4., nxy+nxy/4., no)

# %% Interpolate
zo = dstwgt2dto1d(xi, yi, zi, xo, yo)

# %% Plot
kw = dict(vmin=np.nanmin(zi), vmax=np.nanmax(zi))
plt.figure(figsize=(6, 6))
plt.subplot(111, aspect=1)
plt.contourf(xxi, yyi, zi[0], **kw)
plt.scatter(xxi, yyi, color='k')
plt.scatter(xo, yo, c=zo[0], s=50, **kw)
plt.title('dstwgt2dto1d')
