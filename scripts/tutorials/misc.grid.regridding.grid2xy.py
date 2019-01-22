"""Interpolate gridded data to random positions"""
from vcmqm import cdms2, MV2, np, data_sample, hov2, map2, grid2xy

# %% Read data
select=dict(lon=(-5.3, -4.72), lat=(47.9, 48.8), time=slice(0, 24))
f = cdms2.open(data_sample('mars2d.xyt.nc'))
v = MV2.masked_values(f('v',**select), 0., copy=False)
f.close()
v.long_name = 'Meridional velocity'

# %% Transect
lon = v.getLongitude()
lat = v.getLatitude()
nd = int(round(np.sqrt(len(lon)**2.+len(lat)**2)/2.))
xo = np.linspace(lon[0], lon[-1], nd)
yo = np.linspace(lat[0], lat[-1], nd)

# %% Interpolation
vo = grid2xy(v, xo, yo, method='bilinear')

# %% Plot
# - interpolated data
hov2(vo, cmap='delta', show=False,  top=.9, date_fmt='%H',
     colorbar_shrink=.5,  left=.13)
# - map + transect
m = map2(v[0],  xhide=True, yhide=True, contour=False, proj='merc',
         title=False, autoresize=0, cmap='delta', res='l',
         colorbar=False, axes_rect=[.78, .78, .2, .2], show=False)
m.add_lines(xo, yo, color='tab:red', linewidth=2)
