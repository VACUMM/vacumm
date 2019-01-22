"""Test fortran function :f:func:`bilin2dto1dc`"""

from vcmqm import P, N, add_grid, rotate_grid
from vacumm.fortran.interp import bilin2dto1dc


# Input grid and data
nxy = 15
xi = N.arange(nxy*1.)
yi = N.arange(nxy*1.)
gridi = rotate_grid((xi, yi), 30)
xxi = gridi.getLongitude().getValue()
yyi = gridi.getLatitude().getValue()
zzi = N.ma.array(yyi)
zzi[int(nxy*0.3):int(nxy*0.8), int(nxy*0.3):int(nxy*0.8)] = N.ma.masked
zzi.shape = 1, nxy, nxy

# Output positions
no = 1000
xo = N.random.uniform(-nxy/4., nxy+nxy/4., no)
yo = N.random.uniform(-nxy/4., nxy+nxy/4., no)

# Interpolate
mv = zzi.get_fill_value()
zo = bilin2dto1dc(xxi,yyi,zzi.filled(mv),xo,yo,mv)
zo = N.ma.masked_values(zo, mv)

# Plot
kw = dict(vmin=zzi.min(), vmax=zzi.max())
P.figure(figsize=(6, 6))
P.subplot(111, aspect=1)
P.contourf(xxi, yyi, zzi[0], **kw)
add_grid((xxi, yyi), edges=False, centers=True, marker='o')
P.scatter(xo, yo, c=zo[0], s=50, **kw)
P.title('bilin2dto1dc')


