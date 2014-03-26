"""Test :func:`~vacumm.misc.grid.misc.coord2slice`"""

from vcmq import code_base_name, P, os, N, create_lon
from vacumm.misc.grid import coord2slice, create_grid, create_axes2d, meshbounds
from vacumm.misc.plot import add_grid

figfiles = []
figfile = code_base_name(ext=False)+'_%i.png'
def plot(xx, yy, target, label, figfiles, lon=None, lat=None, show=False):
    xs, ys, mask = coord2slice(target, lon=lon, lat=lat)
    P.figure(figsize=(6, 3.5))
    P.title('Target=%(label)s / select: lon=%(lon)s, lat=%(lat)s'%locals())
    add_grid((xx, yy))
    xx = xx.asma()
    yy = yy.asma()
    if isinstance(lon, tuple): 
        P.axvline(lon[0], color='m', ls='--', lw=2)
        P.axvline(lon[1], color='m', ls='--', lw=2)
    elif isinstance(lon, slice):
        i, j, k = lon.indices(xx.shape[1])
        P.plot(xx[:, i], yy[:, i], 'c--', lw=2)
        P.plot(xx[:, j-1], yy[:, j-1], 'c--', lw=2)
    if isinstance(lat, tuple): 
        P.axhline(lat[0], color='m', ls='--', lw=2)
        P.axhline(lat[1], color='m', ls='--', lw=2)
    elif isinstance(lat, slice):
        i, j, k = lat.indices(yy.shape[0])
        P.plot(xx[i], yy[i], 'c--', lw=2)
        P.plot(xx[j-1], yy[j-1], 'c--', lw=2)
    P.xticks(N.arange(xx.min()-1, xx.max()+1))
    P.yticks(N.arange(yy.min()-1, yy.max()+1))
    xxi, yyi = xx, yy
    xx = xx[ys, xs]
    yy = yy[ys, xs]
#    mask = mask[ys, xs]
    xxb, yyb = meshbounds(xx, yy)
    P.pcolor(xxb, yyb, mask, shading='faceted')
    P.scatter(xx.ravel(), yy.ravel(), c=(0, 1, 0))
    P.grid('on')
    P.axis('image')
    P.tight_layout()
    i = len(figfiles)
    savefig = figfile%i
    if os.path.exists(savefig): os.remove(savefig)
    P.savefig(savefig)
    figfiles.append(savefig)
    if show: P.show()
    else: P.close()
    
print '1D axis'
lon1d = create_lon((0, 10))
print coord2slice(lon1d, lon=(2.5, 4., 'cc'))
print coord2slice(lon1d, lon=(2.5, 4., 'ccb')) 
print coord2slice(lon1d, lon=slice(3, 6))
print coord2slice(lon1d, lat=(6, 8))
print coord2slice(lon1d, lon=(60, 70))

print 'Rect grid'
grid = create_grid((0, 10.), (20, 30.))
print coord2slice(grid, lon=(0., 3.5), lat=slice(3, 5))
print coord2slice(grid, lat=(21,21, 'ccb')) 

print '2D axis'
lon2d = N.empty((10, 10.))
for i in xrange(10): 
    lon2d[i] = lon1d[:]+i
lat2d = N.resize((N.arange(10)+20), (10, 10)).T
lon2d, lat2d = create_axes2d(lon2d, lat2d)
kw = dict(show=False)
plot(lon2d, lat2d, lon2d, 'lon2d', figfiles, lon=(2, 4), **kw)
plot(lon2d, lat2d, lon2d, 'lon2d', figfiles, lon=(2, 4), lat=slice(0, 2), **kw)
plot(lon2d, lat2d, lat2d,  'lat2d', figfiles, lat=(22, 26.6,'ccb'), **kw)

print 'Curv grid'
grid = create_grid(lon2d, lat2d)
plot(lon2d, lat2d, grid, 'grid', figfiles, lon=(8, 11, 'cc'), lat=(21.9, 26., 'cc'), **kw)
plot(lon2d, lat2d, grid, 'grid', figfiles, lon=slice(2, 5), lat=(23.4, 23.6, 'ccb'), **kw)
print coord2slice(grid,lon=(8,8,'ccb'),lat=(24,24,'ccb'))
