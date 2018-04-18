# -*- coding: utf8 -*-
from vcmq import (cdms2, data_sample, regrid2d, resol, create_grid, 
    GriddedBathy, GriddedBathyMerger, add_grid)

# Création de bathymétries fictives à partir de Smith and Sandwell
cdms2.axis.longitude_aliases.append('x')
cdms2.axis.latitude_aliases.append('y')
f = cdms2.open(data_sample('ETOPO2v2g_flt.grd'))
# - large
var_large = f('z', lon=(-7, -1), lat=(46, 49))
# - petite
var_small = f('z', lon=(-4.5, -.5), lat=(44.5, 46.5))
f.close()
# - regrillage de la large vers une grille moins fine
xr, yr = resol(var_large.getGrid())
grid_large = create_grid((-7., -1, xr*4.5), (46., 49, yr*4.5))
var_large = regrid2d(var_large, grid_large)

# On ajoute un traît de côte pour le masquage
bathy_large = GriddedBathy(var_large, shoreline='i')
bathy_small = var_small

# Création de la grille finale de résolution intermédiaire
final_grid = create_grid((-6., .5, xr*2.5), (45, 47., yr*2.5))

# On crée maintenant le merger
merger = GriddedBathyMerger(final_grid)
# - ajout de la bathy basse resolution en premier (en dessous)
merger += bathy_large
# - puis ajout de celle haute résolution
merger += bathy_small

# On définit le traît de côte pour le masquage
merger.set_shoreline('h')

# Fusion vers la grille finale
bathy = merger.merge()

# Plot
merger.plot(show=False)
kwgrid = dict(linewidth=.5, alpha=.5, samp=2)
add_grid(grid_large, color='r', **kwgrid)
add_grid(var_small.getGrid(), color='#00ff00', **kwgrid)
