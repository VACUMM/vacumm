"""Test fortran function :f:func:`dstwgt2dto1d`"""

from vcmq import P, N, add_grid
from vacumm.fortran.interp import dstwgt2dto1d


# Input grid and data
nxy = 15
xi = N.arange(nxy*1.)
yi = N.arange(nxy*1.)
xxi, yyi = N.meshgrid(xi, yi)
zi = N.ma.array(yyi)
zi[int(nxy*0.3):int(nxy*0.8), int(nxy*0.3):int(nxy*0.8)] = N.ma.masked
zi.shape = 1, nxy, nxy

# Output positions
no = 1000
xo = N.random.uniform(-nxy/4., nxy+nxy/4., no)
yo = N.random.uniform(-nxy/4., nxy+nxy/4., no)

# Interpolate
mv = zi.get_fill_value()
zo = dstwgt2dto1d(xi,yi,zi.filled(mv),xo,yo,mv)
zo = N.ma.masked_values(zo, mv)

# Plot
kw = dict(vmin=zi.min(), vmax=zi.max())
P.figure(figsize=(6, 6))
P.subplot(111, aspect=1)
P.contourf(xxi, yyi, zi[0], **kw)
add_grid((xi, yi), edges=False, centers=True, marker='o')
P.scatter(xo, yo, c=zo[0], s=50, **kw)
P.title('dstwgt2dto1d')

