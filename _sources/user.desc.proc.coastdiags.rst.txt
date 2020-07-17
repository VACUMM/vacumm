.. _user.desc.proc.coastdiags:


Diagnostics des processus côtiers
*********************************

.. _user.desc.proc.coastdiags.bour:
    
Bourrelet d’eau froide
======================


Cette tâche consiste à représenter le bourrelet d’eau froide à travers 
une section verticale pouvant n’être
ni zonale, ni méridienne.
La procédure est la suivante :
    
- La première étape consiste à calculer la résolution de la grille.
- Celle-ci est utilisé ensuite pour définir la section et le domaine
  spatial de lecteur des données.
- Le bloc spatial 3D contenant la future section est ensuite lu.
- Si le champs est en coordonnées sigma, il convient aussi de 
  récupérer les profondeurs associées à chaque points de grille.
- Les champs 3D (quantité physique, et profondeur si besoin) 
  sont alors interpolés spatialement sur la section (utilisation de :func:`~vacumm.misc.grid.regridding.grid2xy`).
- Le champs n'est finalement pas interpolé verticalement mais tracé directement 
  sur sa grille verticale potentiellement irrégulière.

On définit les limites de la section : ::
    
    xmin, xmax = -5., -3.
    ymin, ymax = 43, 45
    
On récupère la résoltion de la grille avc un lecteur :class:`MARS3D` : ::
    
    # Lecture de la grille
    mars = MARS3D(config, dates)
    grid = mars.get_grid()
    
    # Calcul de la résolution en degrés
    from vacumm.misc.grid import resol
    xres, yres = resol(grid)
    
On calcul les coordonnées de la section par exemple en
évaluant approximativement le nombre de points de grille,
et en multipliant cette valeur par un facteur pour être sûr
d'avoir assez de point de section : ::
    
    
    nx = (xmax-xmin)/xres
    ny = (xmax-xmin)/xres
    nxy = int(N.ceil(N.sqrt(nx**2+ny**2)))*1.5
    xx = N.linspace(xmin, xmax, nxy)
    yy = N.linspace(ymin, ymax, nxy)
    
    
On récupére le block 3D, avec une marge déduite de la résolution afin d'être sûr
d'englober la section avec les points de grilles : ::
    
    # Initialisation
    mars = MARS3D(config, date, 
        lon=(xmin-xres, xmax+xres, 'ccb'), 
        lat=(xmax-xres, xmax+xres, 'ccb'))
        
    # Lecture de la température
    temp = mars.get_temp()
    
    # Lecture des profondeurs en sigma
    with_sigma = hasattr(mars, get_sigma_depths)
    if with_sigma:
        depths = temp.get_depths()
 
    
On interpole sur la section : ::
    
    from vacumm.misc.grid import grid2xy
    xytemp = grid2xy(temp, xx, yy)
    if  with_sigma:
        xydepths = grid2xy(depths, xx, yy)


On trace la section avec indication de sa position : ::

    # Section
    from vacumm.misc.plot import section2, map2
    if not with_sigma:
        section2(xytemp, show=False)
    else:
        section2(xytemp, yaxis=xydepths, show=False)
        
    # Carte
    m = map2(lon=(xmin-dx, xmax+dx), lat=(ymin-dy, ymax+dy),
         axes_rect = [0.8, 0.8, 0.98, 0.98],
         drawmeridians_size=6, drawparallels_size=6,
         show = False)
    xx, yy = m.map(xx, yy)
    m.map.plot(xx, yy, 'r-')
    m.show()



Front d’Ouessant
================

Le front d’Ouessant est diagnostiqué à travers l’utilisation de la SST à une latitude et dans un inter-
valle de longitudes donnés. 
Les graphiquues suivants sont générés :
– Un diagramme de Hovmoller.
– Un évolution zonale du minimum de SST.
– La valeur du minimum de température sur cette section en fonction du temps.
 
Afin de permettre les comparaisons de plusieurs sources de données,
la SST est interpolée sur la section.

Les paramètres d'entrée sont les suivants : ::
    
    lon_min = -6.
    lon_max = -5.5
    lat = 47.5

Section à partir d'une grille rectangulaire
-------------------------------------------

On estime la résolution de la grille comme dans le cas de la section :ref:`user.desc.proc.coastdiags.bour`,
et on obtient ``xres``, ``yres`` et ``grid``.
La section zonale est alors définie à partir des points de la grille et des extrémités spécifiées de la section.

Les points de la grille sont obtenus ainsi : ::
    
    subgrid = grid.subGridRegion((lat, lat, 'ccb'), (lon_min, lon_max, 'oo'))
    xgrid = subgrid.getLongitude()[:]

On ajoute les extrémités pour définir la section : ::
    
    xx = N.concatenate(([lon_min], xgrid, [lon_max]))
    yy = N.ones(len(xx))*lat

Le modèle est lu par exemple de la manière suivante : ::
    
    mars = MARS3D(cfg, 
        lon=(lon_min-xres, lon_max+xres, 'ccb'), 
        lat=(lat-yres, lat_max+yres, 'ccb'))
    sst2d = mars.get_sst()
    
On interpole alors sur la section comme dans le cas du bourrelet : ::

    from vacumm.misc.grid.regridding import grid2xy
    sst = grid2xy(sst2d, xx, yy)
    
    
Section à partir d'une grille curvilinéaire
-------------------------------------------

TODO ?
    
Diagnostics et graphiques
-------------------------

Evolution de la température le long d’une section : ::
    
    from vacumm.misc.plot import hov2
    hov2(sst)
    
Localisation méridienne du minimum de température à un degré de latitude donné : ::
    
    # Fonction de récupération de la position des min
    def get_minx(sst):
        
        # Position min à chaque instant avec valeurs masquées
        imin = N.ma.argmin(sst, axis=1, fill_value=-1)
        bad = imin == -1 
    
        # Initialisation variable
        sst_minx = sst[:,0].clone()
        sst_minx.long_name = 'Longitude of minimum'
        sst_minx.units = 'degrees_east'
        
        # Remplissage
        for j, i in N.ndenumerate(imin):
            sst_minx[j] = sst.asma()[j, i]
        
        # Masquage
        sst_minx[:] = MV2.masked_where(bad, sst_minx, copy=0)
        
    # Estimations de la position minimum
    minx_mars = get_minx(sst_mars)
    minx_sat = get_minx(sst_sat)
    
    # Plots
    from vacumm.misc.plot import curve2
    curve2(minx_mars, 'b', show=False, label='MARS')
    curve2(minx_sat, 'ok', show=False, label='Satellite')
    import matplotlib.pyplot as P
    P.legend(loc='best')


Valeur du minimum de température à un degré de latitude donné : ::

        
    # Estimations de la position minimum
    min_mars = MV2.min(sst_mars, axis=1)
    min_sat = MV2.min(sst_sat, axis=1)
    
    # Plots
    from vacumm.misc.plot import curve2
    curve2(min_mars, 'b', show=False, label='MARS')
    curve2(min_sat, 'ok', show=False, label='Satellite')
    import matplotlib.pyplot as P
    P.legend(loc='best')


