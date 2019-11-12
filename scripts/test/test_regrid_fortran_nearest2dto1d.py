"""Test fortran function :f:func:`nearest2dto1d`"""

from utils import np, plt
from vacumm.fortran.interp import nearest2dto1d


# Input grid and data
nxy = 15
xi = np.arange(nxy*1.)
yi = np.arange(nxy*1.)
xxi, yyi = np.meshgrid(xi, yi)
zi = np.array(yyi)
zi[int(nxy*0.3):int(nxy*0.8), int(nxy*0.3):int(nxy*0.8)] = np.nan
zi.shape = 1, nxy, nxy

# Output positions
no = 1000
xo = np.random.uniform(-nxy/4., nxy+nxy/4., no)
yo = np.random.uniform(-nxy/4., nxy+nxy/4., no)

# Interpolate
zo = nearest2dto1d(xi, yi, zi, xo, yo)

# Plot
kw = dict(vmin=np.nanmin(zi), vmax=np.nanmax(zi))
plt.figure(figsize=(6, 6))
plt.subplot(111, aspect=1)
plt.contourf(xxi, yyi, zi[0], **kw)
#add_grid((xi, yi), edges=False, centers=True, marker='o')
plt.scatter(xo, yo, c=zo[0], s=50, **kw)
plt.title('nearest2dto1d')
plt.tight_layout()
