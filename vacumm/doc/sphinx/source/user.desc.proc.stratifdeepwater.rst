.. _user.desc.proc.stratifdeepwater:


Stratification : comparaisons en eaux profondes
***********************************************

Le but est ici de comparer les couches mélangées modélisée et observée en eaux profondes.
L'approche est équivalente à celle du calcul de profondeur de la section ":ref:`user.desc.proc.stratifinsitu`" en eau peu profonde,
mais avec un algorithme différent pour le calcul de profondeur 
(cf. ":ref:`user.desc.proc.stratifinsitu.mld`").

L'approche la plus pertinente et simple est celle basée sur un différentiel en température ou en densité
entre les eaux de proche surface (par exemple 10m) et celles de la base de la couche mélangée.
On peut pour cela utiliser les critères de [deBoyerMontegutetal2003]_.

Pour l'implémentation, nous allons raisonner sur la densité.
Nous supposons avoir calculé cette dernière comme en eau peu profonde,
en supposant le dernier niveau le plus proche de la surface.

La première étape consiste à calculer la densité de référence en proche surface : ::
    
    # Profondeur de référence
    dep = -10
    
    # Axe pour interpolation
    from vacumm.misc.axes import create_dep
    depaxis = create_dep([dep])

    # Interpolations
    from vacumm.misc.grid.regridding import interp1d
    dens_obs_ref = interp1d(dens_obs, depaxis, xmap=-1, deps_obs)


La deuxième étape vise à calculer la profondeur ayant une densité plus importante
que celle de rérence d'une valeur prédéfinie.
Pour cela, on calcul se sert du masque des données situées dans
la couche mélangée (critère en densité), et du masque de celles en dehors (profondes).
Une interpolation linéaire de la prodondeur est alors faite entre
le niveau de fond de la couche mélangée (densité max et profondeur négative min) 
et le niveau supérieur des eaux profondes (densité min et profondeur négative max) : ::
    
    # Valeur du différentiel de densité (cf. de Boyer Montégut et all, 2003)
    delta_dens = 0.03 

    # Densité cible
    dens_obs_target = dens_obs_ref+0.03 
    
    # Masques de la couche mélangée et des eaux profondes
    dens_obs_target3d = MV2.resize(dens_obs_target, dens_obs.shape)
    dens_obs_good = (dens_obs.asma() <= dens_obs_target3d.asma()).filled(False)
    dens_obs_bad = (dens_obs.asma() > dens_obs_target3d.asma()).filled(False)
    
    # Profondeurs juste au dessus et en dessous de la MLD
    deps_obs_above = MV2.masked_where(dens_obs_bad, deps_obs).min(axis=0)
    deps_obs_below = MV2.masked_where(dens_obs_good, deps_obs).max(axis=0)
    
    # Masques associés
    from vacumm.misc import closeto
    mask_above = closeto(deps_obs, MV2.resize(deps_obs_above, deps_obs.shape))
    mask_below = closeto(deps_obs, MV2.resize(deps_obs_below, deps_obs.shape))
    
    
    # Densités juste au dessus et en dessous de la MLD
    dens_obs_above = MV2.masked_where(dens_obs, mask_above).max(axis=0)
    dens_obs_below = MV2.masked_where(dens_obs, mask_below).min(axis=0)
    
    # Interpolation
    dens_obs_delta = dens_obs_above-dens_obs_below
    deps_obs_mld = (dens_obs_target-dens_obs_above)*deps_obs_above
    deps_obs_mld += (dens_obs_below-dens_obs_target)*deps_obs_below
    deps_obs_mld = MV2.where(dens_obs_delta.mask, MV2.masked, deps_obs_mld/dens_obs_delta)
    
    
    
