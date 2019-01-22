"""Compare cella average versus linear interpolation"""
from vcmq import (cdms2, MV2, plt, np, data_sample, create_lon, curve2,
                  regrid1d, yscale)

# Read temperature
f = cdms2.open(data_sample('mars3d.x.uv.nc'))
v = f('v', squeeze=1, lon=(-5.08, -4.88))
f.close()
v = MV2.masked_values(v, 0.)
v.long_name = 'Original'

# Add some noise
np.random.seed(1)
v[:] += np.random.random(len(v))

# Create longitudes
lon = v.getLongitude().getValue()
dlon = np.diff(lon).mean()
# - low resolution
lon_lr = create_lon((lon.min()+dlon*0.7, lon.max()+dlon, dlon*5.3))
# - high resolution
lon_hr = create_lon((lon.min()+dlon, lon.max()+dlon, dlon/3.3))
# - dict
lons = dict(low=lon_lr, high=lon_hr)

# Interp and plot
plt.rcParams['font.size'] = 9
plt.figure(figsize=(5.5, 6))
kwplot = dict(show=False, vmin=v.min(), vmax=v.max(), alpha=.7, linewidth=.6)
for ilh, resdst  in enumerate(lons):

    # Interp
    vlinear = regrid1d(v, lons[resdst], 'linear')
    vremap = regrid1d(v, lons[resdst], 'cellave')

    # Plots
    plt.subplot(2, 1, ilh+1)
    curve2(v, 'o-', markersize=3, color='k', label=u'Original',
           hspace=.3, **kwplot)
    curve2(vremap, 'o-', markersize=2, label=u'Cellave',
           color='tab:blue', **kwplot)
    curve2(vlinear, 'o-', markersize=2, label=u'Linear',
           color='tab:red', **kwplot)
    yscale(1.1)
    plt.title('Toward {} resolution'.format(resdst))
    if not ilh:
        plt.legend(loc='lower left')#.legendPatch.set_alpha(.6)
