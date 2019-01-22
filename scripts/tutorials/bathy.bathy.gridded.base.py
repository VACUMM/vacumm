# -*- coding: utf8 -*-

from vcmq import (cdms2, data_sample, GriddedBathy, create_grid, P,
                  plot_bathy, auto_cmap_topo, land)

# %% Read bathymetry
cdms2.axis.latitude_aliases.append('y')
cdms2.axis.longitude_aliases.append('x')
f = cdms2.open(data_sample('ETOPO2v2g_flt.grd'))
var = f('z', lon=(-6.1, -3), lat=(47.8, 48.8))
f.close()

# %% Create a Bathy object
bathy = GriddedBathy(var)

# %% Shoreline resolution for masking
shoreline_resolution = 'l'
bathy.set_shoreline(shoreline_resolution)

# %% Get masked and unmasked bathies
bathy_orig = bathy.bathy(mask=False)
bathy_masked = bathy.bathy(mask=True)

# %% On new grid
new_grid = create_grid((-5.5, -4.6, 0.09), (48.25, 48.6, .06))

# %% Regrid
interp_bathy_masked = bathy.regrid(new_grid)

# %% Regrid and mask
# - interpolation
interp_bathy_orig = bathy.regrid(new_grid, mask=False)
# - masking
interp_bathy_orig_masked = GriddedBathy(
        interp_bathy_orig, shoreline=shoreline_resolution).bathy()

# %% Plots
P.rc('font', size=9)
P.rc('axes', titlesize=9)
kwplot = dict(resolution=None, show=False, colorbar=False, contour=False,
              top=.97, hspace=.25, bottom=.03, right=.98, left=0.08,
              autoresize=False, vmax=bathy_orig.max(), vmin=bathy_orig.min())
# - colormap
kwplot['cmap'] = auto_cmap_topo(bathy_orig)
kwplot['bgcolor'] = land
# - direct
bathy.plot(title='Original', mask=False, figsize=(7, 4), subplot=221,
           key=1, **kwplot)
bathy.plot(title='Masked', mask=True, subplot=222, key=2,
           yhide=True, **kwplot)
# - indirect
plot_bathy(interp_bathy_masked, title='Masked->Interpolated', subplot=223,
           key=3, **kwplot)
plot_bathy(interp_bathy_orig_masked, title='Interpolated->Masked',
           yhide=True, subplot=224, key=4, **kwplot)
P.rcdefaults()
