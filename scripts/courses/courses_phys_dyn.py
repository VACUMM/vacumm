#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Diagnostics thermodynamiques et dynamiques (:mod:`vacumm.diag.thermdyn` et :mod:`vacumm.diag.dynamics`)"""

from vcmq import DS, data_sample, map2, density, mixed_layer_depth,  barotropic_geostrophic_velocity, kinetic_energy, shapiro2d

# Lecture des données de base
ds = DS(data_sample('menor.nc'))
temp = ds.get_temp(squeeze=1)
sal = ds.get_sal(squeeze=1)
depth = ds.get_depth(squeeze=1)
ssh = ds.get_ssh(squeeze=1)


# Densité
dens_surf = density(temp[-1], sal[-1])                      # -> CALCULER DENSITE 3D
map2(dens_surf)


# Couche mélangée
mld = mixed_layer_depth((temp, sal), depth, mode='deltadens') # -> ESSAYER AUTRES MODES
map2(mld, vmax=600.)


# Vitesse geostrophique
# - vitesses
ug, vg = barotropic_geostrophic_velocity(ssh)
# - energie cinetique
ke = kinetic_energy((ug, vg))                                 # -> TESTER AVEC SSH
ke.long_name = 'Surface geostrophique kinetic energy'
ke = shapiro2d(ke)
# - plot
map2((ke, ug, vg), fill='pcolormesh', vmax=0.2, quiver_vmax=1.5, figsize=(10, 8), 
    quiver_samp=2, quiver_width=2, contour=False, cmap_tint=0.3, 
    quiverkey_value=1, quiver_scale=13)
    
    
    
