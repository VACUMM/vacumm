.. _user.desc.proc.warmpool:


Warm pool
*********

L'objectif est d'estimer l'accumulation d'eau à la fin de l'été au fond du Golfe de Gascogne,
dans le modèles et les observations satellites.

Préalables
==========

Un préalable essentiel aux diagnostics est la colocalisation des SST du modèle et satellite,
lue toutes les deux vers minuit.
Cela implique deux considérations importantes :

    - Les scripts doivent être capables de sappliquer à des données satellites sur grille curvilinéaire.
    - La colocalisation implique l'interpolation spatiale d'un des deux champs vers la grille de plus basse résolution.


La première étape consiste donc à comparer les résolutions des champs en entrée : ::
    
    from vacumm.misc.grid import resol
    from vacumm.misc.grid.basemap import merc
    
    # Résolution du modèle
    grid_mod = MARS3D(cfgmod, dates).get_grid()
    xres_mod, yres_mod = resol(grid_mod)
    
    # Résolution du modèle
    grid_sat = Satellite(cfgsat, dates).get_grid()
    xres_sat, yres_sat = resol(grid_sat)
    
    # Conversions en mètres
    xres_mod, yres_mod = merc(xres_mod, yres_mod)
    xres_sat, yres_sat = merc(xres_sat, yres_sat)
    
    # Resolution commune X & Y
    res_mod = N.sqrt(xres_mod**2+yres_mod**2)
    res_sat = N.sqrt(xres_sat**2+yres_sat**2)


Nous devons maintenant choisir quelle sera la grille cible (modèle ou satellite ?).
On se base pour cela sur un facteur limite à définir,
que l'on utilise pour comparer la résolution des deux grilles,
en favorisant le type rectangulaire rectangulaire : ::
    
    # Plus le facteur est faible, plus la résolution de la grille
    # satellite curvilinéaire doit être fine pour que cette dernière soit choisie
    ratio_curvgrid_choice = 0.8
    
    # Maintenant, le choix
    from vacumm.misc.grid import isgrid
    if isgrid(grid_sat, curv=True):  # La grille satellite est bien curvilinéaire
    
        choose_sat = (res_mod / res_sat) < ratio_grid_choice
        
    else: # Cas de deux grilles rectangulaires
    
        choose_sat = res_sat > res_mod 
    

Nous devons maintenant colocaliser les champs : ::
    
    # Minuit
    from vacumm.misc.atime import midnight_interval
    midnight = midnight_interval(date)    
    cdms2.setAutoBounds(1)
    
    # Lecture des données vers minuit.
    mod = MARS3D(cfgmod, midnight).load_sst()
    sat = Satellite(cfgsat, midnight).load_sst()
    
    # Regrillage
    from vacumm.misc.grid.regridding import regrid2d
    if choose_sat: # Vers grille satellite
        
        mod = regrid2d(mod, grid_sat)
        
    else: # Vers grille modèle
        
        sat = regrid2d(sat, grid_mod)
        
        
.. note::
    
    Par défaut, 
    la méthode d'interpolation est choisie par la fonction
    :func:`~vacumm.misc.grid.regridding.regrid_method` en se basant
    sur le rapport des résolutions des deux grilles.
    Pour plus d'information, voir la section ":ref:`user.desc.proc.regridding`".
        
.. warning::
    
    Le regrillage vers une grille satellite curvilinéaire nécessite
    d'avoir l'exécutable :program:`scrip` accessible dans la variable d'environnement
    :envvar:`PATH`.

On suppose par la suite que les SST sont colocalisées.
  
    
.. _user.desc.proc.warmpool.extent:
    
Extension de la *Warm pool* à une date donnée
=============================================

La *warm pool* satellite est tracée avec une palette grise
sur un fond de SST du modèle, dont la *warm pool*  est elle aussi 
mise en valeur (avec un contour, par exemple rouge) : ::

    # Inits
    from vacumm.misc.plot import map2

    # Définition de la warm pool
    wp_temp = 20
    wp_color = 'r'
    wp_linewidth = 1.5

    # SST du modèle
    levels_mod = N.arange(N.floor(mod.min()), N.ceil(mod.max()))
    map2(mod, levels=levels_mod, show=False)
    if (mod.asma().min()<=wp_temp and mod.asma().max()>=wp_temp):
        map2(mod, levels=wp_temp, linewidth=wp_linewidth, 
            contour_colors=wp_color, fill=False, show=False)
    else:
        print u'Pas de warm pool dans le modèle'
        
    # SST satellite seulement sur la warm pool
    sat = MV2.masked_less(sat, wp_temp) # On masque ailleurs
    levels_sat = N.arange(N.floor(sat.min()), N.ceil(sat.max()))
    if (sat.asma().min()<=wp_temp and sat.asma().max()>=wp_temp):
        map2(mod, levels=levels_sat, cmap='cmap_grey')
    else:
        print u'Pas de warm pool dans les données satellites'


Évolution temporelle de la *warm pool*
======================================

Le but est de comparer une mesure de l'évolution temporelle de la 
*warm pool* du modèle et des données satellites.
Cette dernière est estimée en $km^2$ dans trois cas différents, à partir
des SST colocalisées :
    
    - Dans le modèle.
    - Dans le modèle, mais avec un masque commun aux observations.
    - Dans les observations, avec un masque commun au modèle.
    
À titre indicatif, une estimation du recouvrement des 
*warm pool* modélisée et observée est fournie à travers deux indices
au cours du temps, en se basant sur les champs colocalisés avec même masque :
    
    - La fraction de *warm pool* commune rapportée à la *warm pool* modélisée.
    - La fraction de *warm pool* commune rapportée à la *warm pool* observée.
    



Pour procéder, nous devons effectuer une boucle temporelle sur les dates
à minuit de l'interval de temps considéré.
Afin d'optimiser un peu les traitements, nous pouvons créer une classe dédiée
à la colocalisation des SST, dont le coeur se base ce qui a été présenté précédemment : ::
    
    class Coloc(object):
        """Colocalize model and observation SST
        
        :Example:
            
            >>> coloc = Coloc(grid_mod, grid_sst, ratio = ratio_curvgrid_choice)
            >>> sst_mod, sst_sat = coloc(sst_mod, sst_sat)
        """
        
        def __init__(self, grid_mod, grid_sst, ratio=0.8):
            
            self.grid_mod = grid_mod
            self.grid_sat = grid_sat
            
            # Cf. plus haut
            [...]
            self.choose_sat = ...
            
        def __call__(self, sst_mod, sst_sat):
            
            if self.choose_sat:
                sst_mod = regrid2(sst_mod, self.grid_sat)
            else:
                sst_sat = regrid2(sst_sat, self.grid_mod)
            return sst_mod, sst_sat

            
La boucle temporelle peut s'écrire ainsi en supposant 
que l'ont traite les données par blocs de ``nmaxdays`` : ::
    
    # Taille des blocs XYT en jours
    nmaxdays = 30
    
    # Colocalisateur
    coloc = Coloc(grid_mod, grid_sat, ratio = ratio_curvgrid_choice)
    
    
    # Initialisations
    ext = dict(MOD=[], mod=[], sat=[]) # Extensions
    rec = dict(mod=[], sat=[]) # Recouvrements
    ctimes = []
    areas = None
    
    # On boucle sur les jours
    from vacumm.misc.atime import IterDates, Intervals
    for interval in Intervals((date_min, date_max), (nmaxdays, 'day')):
        
        # Lecture
        block_mod = MARS3D(cfgmod, interval).load_sst()
        block_sat = Satellite(cfgsat, interval).load_sst()
        
        # Colocalisation
        mod, sat = coloc(mod, sat)
        
        # Surface de chaque cellule
        if areas is None:
            from vacumm.misc.grid import meshweights, get_xy
            xx, yy = get_xy(mod.getGrid(), num=True)
            areas = meshweights(xx, yy, proj=True) # proj -> m2
            areas *= 1e-6 # km2
        
        # Boucle à minuit
        for date in IterDates(midnight_interval(interval, (1, 'day'))):
            
            # Extraction
            try:
                mod = block_mod(time=(date, date, 'ccb'), squeeze=1)
            except:
                continue
            try:
                sat = block_sat(time=(date, date, 'ccb'), squeeze=1)
            except:
                continue
            
            # Calcul des masques
            wp_mod = (mod.asma()>wp_temp).astype('f').filled(0.)
            wp_sat = (sat.asma()>wp_temp).astype('f').filled(0.)          
            
            # Aire totale dans le modèle
            ext['MOD'].append((wp_mod*areas).sum())
            
            # Fusion des masques modèle/satellite
            mask = (sat.mask | mod.mask).astype('f')
            
            # Aire du modèle masqué
            ext['mod'].append((wp_mod*mask*areas).sum())
            
            # Aire des données satellites
            ext['sat'].append((wp_sat*mask*areas).sum())
            
            # Recouvrement
            intersect = (wp_mod*wp_sat).sum()
            rec['mod'].append(-1. if ext['mod']==0. else intersect/ext['mod'])
            rec['sat'].append(-1. if ext['sat']==0. else intersect/ext['sat'])
            
            # Dates
            ctimes.extend(mod.getTime().asComponentTime())
            
    # Concatenations et formattages
    from vacumm.misc.axes import create_time
    taxis = create_time(ctimes)
    for dd in ext, rec:
        for key, values in dd.items():
            dd[key] = MV2.array(values)
            dd[key].setAxis(0, taxis)
    for var in ext.values():
        var.units = 'km^2'
    for var in rec.values():
        var[:] = MV2.masked_less(var, 0.)
        var.units = '%%'
        var[:] *= 100.
    ext['MOD'].long_name = 'Full model warm pool extent'
    ext['mod'].long_name = 'Model warm pool extent restricted to satellite'
    ext['sat'].long_name = 'Satellite warm pool extent'
    rec['mod'].long_name = 'Fraction of model warm pool common to satellite'
    rec['sat'].long_name = 'Fraction of satellite warm pool common to model'
    
    # Tracé des extensions
    from vacumm.misc.plot import curve2
    curve2(ext['MOD'], 'r--', vmin=0, show=False)
    curve2(ext['mod'], 'r-', show=False)
    curve2(ext['sat'], 'b-', show=False).legend(loc='upper left')
    
    # Tracé des fractions de recouvrement
    curve2(rec['mod'], 'r+', twin='x', vmin=0, vmax=100, show=False)
    curve2(rec['sat'], 'bs').legend(loc='upper right')
    
.. note::
    
    On suppose ici que ``date_max`` désigne une borne ouverte de l'intervalle de temps :
    elle n'est donc jamais atteinte explicitement.
    Ainsi, l'intervalle ``('2000', '2001')`` ne bouclera que sur l'année 2000.