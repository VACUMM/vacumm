"""Test the class :class:`~vacumm.misc.regridding.CurvedInterpolator`"""

from vcmq import (plt, np, set_grid, MV2, add_grid, create_time,
                  CurvedInterpolator, rotate_grid)


# Curved grid
nxy = 10
nt = 5
lon = np.arange(nxy*1.)
lat = np.arange(nxy*1.)
time = create_time((nt, ), 'years since 2000')
gridi = rotate_grid((lon, lat), 30)
xxi = gridi.getLongitude()[:].filled()
yyi = gridi.getLatitude()[:].filled()
vari = MV2.resize(yyi, (nt, nxy, nxy))
vari.setAxis(0, time)
set_grid(vari, gridi)
kw = dict(vmin=vari.min(), vmax=vari.max())
plt.figure(figsize=(10, 3.5))
plt.subplot(131, aspect=1)
plt.contourf(xxi, yyi, vari[0].asma(), **kw)
add_grid(gridi, edges=False, centers=-1)
xylims = (xxi.min(), xxi.max(), yyi.min(), yyi.max())
plt.axis(xylims)
plt.title('Curved grid')

# Interpolate to grid
xg, yg = np.meshgrid(np.arange(-3.5, 14.5), np.arange(-3.5, 14.5))
nxyg = xg.shape
cig = CurvedInterpolator(gridi, (xg, yg), g2g=True)
varog = cig(vari)
plt.subplot(132, aspect=1)
plt.scatter(xg, yg, c=varog[0].asma(), s=120, linewidth=0, **kw)
add_grid(gridi, edges=False, centers=-1)
xylims = (xxi.min(), xxi.max(), yyi.min(), yyi.max())
plt.axis(xylims)
plt.title('Interp to grid')

# Interpolate to random
nr = 40
np.random.seed(0)
xr = np.random.uniform(size=nr)*nxy
yr = np.random.uniform(size=nr)*nxy
cir = CurvedInterpolator(gridi, (xr, yr))
varor = cir(vari)
plt.subplot(133, aspect=1)
plt.scatter(xr, yr, c=varor[0].asma(), s=120, linewidth=0, **kw)
add_grid(gridi, edges=False, centers=-1)
xylims = (xxi.min(), xxi.max(), yyi.min(), yyi.max())
plt.axis(xylims)
plt.title('Interp to random')
plt.tight_layout()

