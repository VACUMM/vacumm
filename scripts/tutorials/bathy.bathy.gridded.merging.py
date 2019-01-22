# -*- coding: utf8 -*-
from vcmq import (cdms2, data_sample, regrid2d, resol, create_grid,
                  GriddedBathy, GriddedBathyMerger, add_grid)

# %% Create from Smith and Sandwell
cdms2.axis.longitude_aliases.append('x')
cdms2.axis.latitude_aliases.append('y')
f = cdms2.open(data_sample('ETOPO2v2g_flt.grd'))
# - large
var_large = f('z', lon=(-7, -1), lat=(46, 49))
# - small
var_small = f('z', lon=(-6, -.5), lat=(44.5, 46.5))
f.close()
# - regrid large toward low resolution
xr, yr = resol(var_large.getGrid())
grid_large = create_grid((-7., -1, xr*4.5), (46., 49, yr*4.5))
var_large = regrid2d(var_large, grid_large)

# %% Add shoreline for masking purpose
bathy_large = GriddedBathy(var_large, shoreline='l')
bathy_small = var_small

# %% Create an intermediate resolution final grid
final_grid = create_grid((-6., 0.1, xr*2.5), (45, 47.1, yr*2.5))

# %% Setup the merger
merger = GriddedBathyMerger(final_grid)
# - add the low resolution bathy first
merger += bathy_large
# - then the high resolution one
merger += bathy_small

# %% Add shoreline and max value to the merger for masking purpose
merger.set_shoreline('l')

# %% Merge
bathy = merger.merge()

# %% Plot
merger.plot(show=False, map_res='l', xymasked=False, map_proj='cyl',
            nmax_levels=16, vmax=0, levels_mode='negative', right=1,)
kwgrid = dict(linewidth=.6, alpha=.7, samp=2, zorder=10)
add_grid(grid_large, color='r', **kwgrid)
add_grid(var_small.getGrid(), color='#00ff00', **kwgrid)
