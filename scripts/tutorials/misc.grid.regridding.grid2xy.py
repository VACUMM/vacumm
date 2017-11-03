# Lecture des donnnees
import cdms2, MV2, numpy as N
from vacumm.config import data_sample
select=dict(lon=(-5.3, -4.72), lat=(47.9, 48.8), time=slice(0, 24))
f = cdms2.open(data_sample('mars2d.xyt.nc'))
v =  MV2.masked_values(f('v',**select), 0., copy=False)
f.close()

# Diagonale
lon = v.getLongitude()
lat = v.getLatitude()
nd = N.sqrt(len(lon)**2.+len(lat)**2)/2.
xo = N.linspace(lon[0], lon[-1], nd)
yo = N.linspace(lat[0], lat[-1], nd)

# Interpolation
from vacumm.misc.grid.regridding import grid2xy
vo = grid2xy(v, xo, yo, method='bilinear')

# Plot
# - variable interpolee
from vacumm.misc.plot import hov2 as hov, savefigs, map2 as map
hov(vo, cmap='jet', show=False,  top=.9, date_fmt='%H',
    colorbar_shrink=.5,  left=.13)
# - carte + diagonal
import pylab as P
m = map(v[0],  xhide=True, yhide=True, contour=False,
    title=False, autoresize=0, cmap='cmap_bwr',
    colorbar=False, axes_rect=[.78, .78, .2, .2], show=False)
m.map.plot(xo, yo, 'r-', lw=2)
savefigs(__file__)
P.close()
