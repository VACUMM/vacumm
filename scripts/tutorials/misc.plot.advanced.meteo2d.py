# %% Imports
import cdms2, MV2
from vacumm.misc.plot import map2
from vacumm.misc.color import cmap_custom
from vacumm.misc.filters import generic2d
from vacumm.misc.phys.units import ms2kt
from vacumm.config import data_sample
from matplotlib import rc
rc('font', size=9)

# %% Read data
f = cdms2.open(data_sample('wrf_2d.nc'))
select = dict(lat=slice(4, -4), lon=slice(4, -4), squeeze=1)
slp = f('pslvl', **select)
rain = f('rain', **select)
u = f('u10m', **select)
v = f('v10m', **select)
f.close()

# %% Surface pressure
slp[:] = generic2d(slp*0.01, 5)
map2(slp, fill=False, contour=True, show=False,  figsize=(6, 5.5),
    title='WRF Bretagne',
    fillcontinents=True, zorder=5, projection='merc', right=1, bottom=.07,
    lowhighs=True, lowhighs_smooth=9, lowhighs_zorder=5, fillcontinents_color='.95',
    lowhighs_color=(0, .5, 0), contour_colors=[(0, .5, 0)], drawcoastlines_linewidth=.6)

# %% Rain
cmap_rain = cmap_custom([('0.9', 0), ('b', .8), ('r', 1.)])
rain[:] = MV2.masked_less(rain, 0.1)
map2(rain, cmap=cmap_rain, vmin=0.,
    fill='pcolor', fillcontinents=False, show=False,
    shadow_xoffset=4, shadow_yoffset=-4,shadow_width=4, colorbar_shrink=.7,
    alpha=.7, shadow=True, shadow_alpha=.3, contour=False, zorder=15)

# %% Wind
u[:] = ms2kt(u*5)
v[:] = ms2kt(v*5)
m = map2((u[::3, ::3], v[::3, ::3]), fill=False, contour=False, barbs=True, projection='merc',
    quiver_sizes=dict(height=.2, spacing=.15), quiver_linewidths=.8, zorder=10,
    shadow=True, quiver_alpha=.5, savefigs=__file__, show=False,
    savefigs_pdf=False,
    fillcontinents=False)
