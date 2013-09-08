# Creation de la grille (15x15)
nx, ny = 15, 10
import numpy as N
from vacumm.misc.axes import create_lon, create_lat
lon = create_lon(N.linspace(-5.2, -4.4, nx))
lat = create_lat(N.linspace(48, 48.6, ny))

# Parametres de plot
import pylab as P
import vacumm.misc.color as C
P.figure(figsize=(6, 5))
P.subplots_adjust(top=.92,left=.03,bottom=.03,wspace=.05,hspace=.2)
kwplot = dict(colorbar=False,cmap=C.cmap_linear([(.6, .8, 1), C.land]), \
    contour=False,show=False,fillcontinents=False,fill='pcolor',
    drawmeridians=False, drawparallels=False)

# Creation du mask avec trait de cote fin ('f')
from vacumm.misc.grid.masking import polygon_mask
from vacumm.misc.plot import map2
import MV2 as MV
# - mask si point central sur terre (mode = 0 = 'inside')
lon2,lat2 = N.meshgrid(lon,lat)
mask0 = polygon_mask((lon2, lat2), 'f', mode=0, ocean=False)
mask0 = MV.array(mask0, axes=[lat, lon]) # Forme cdms
mask0.long_name = "'inside'"
map2(mask0, subplot=221, key=1, **kwplot)
# - mask si point central sur terre (mode = 1 = 'intersect')
seuils = .5
mask1 = polygon_mask((lon2, lat2), 'f', mode=1, thresholds=seuils, ocean=False)
mask1 = MV.array(mask1, axes=[lat, lon]) # Forme cdms
mask1.long_name = "'intersect' : seuil=%g" %seuils
map2(mask1, subplot=222, key=2, **kwplot)
# - mask si point central sur terre (mode = 1 = 'intersect')
seuils = .3
mask1 = polygon_mask((lon2, lat2), 'f', mode=1, thresholds=seuils, ocean=False)
mask1 = MV.array(mask1, axes=[lat, lon]) # Forme cdms
mask1.long_name = "'intersect' : seuil=%g" %seuils
map2(mask1, subplot=223, key=3, **kwplot)
# - mask si point central sur terre (mode = 1 = 'intersect')
seuils = (.3, .7)
mask1 = polygon_mask((lon2, lat2), 'f', mode=1, thresholds=seuils, ocean=False)
mask1 = MV.array(mask1, axes=[lat, lon]) # Forme cdms
mask1.long_name = "'intersect' : seuils=(%g/%g)" %seuils
map2(mask1, subplot=224, key=4, **kwplot)

# Figure
from vacumm.misc.plot import savefigs
savefigs(__file__, pdf=True)
P.show()

