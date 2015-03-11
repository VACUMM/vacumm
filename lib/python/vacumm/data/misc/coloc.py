#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#


__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__date__ = '2011-01-17'
__doc__ = 'Dataset Colocalizations'


# ==============================================================================


import os, sys

import cdms2, MV2, numpy, pylab, seawater
from matplotlib.pyplot import colorbar

from vacumm.data.misc.sigma import NcSigma
from vacumm.misc import auto_scale
from vacumm.misc.atime import add as add_time, comptime, datetime as adatetime, Intervals
from vacumm.misc.axes import create_time, create_dep, create_lat, create_lon
from vacumm.misc.bases import Object
from vacumm.misc.color import cmap_magic
from vacumm.misc.io import ncget_var, ncfind_var, list_forecast_files, ncread_best_estimate, NcIterBestEstimate
from vacumm.misc.grid.misc import meshweights, resol
from vacumm.misc.grid.regridding import resol, interp1d, regrid1dold, grid2xy
from vacumm.misc.misc import is_iterable, kwfilter
from vacumm.misc.phys.constants import g
from vacumm.misc.plot import map2, curve2, section2, hov2


class Colocator(Object):

#    def __init__(self, data1, data2):
#        Object.__init__(self)
#        self.data1 = data1
#        self.data2 = data2

    def coloc_mod_on_pro(self, model, profiles, varnames, select=None, method='nearest'):
        '''Colocalize model on profile data.

        Load model data corresponding to the selected profiles positions and time.

        Returns loaded model longitudes, latitudes, depths and requested variable(s)

        :Params:
            - **model**: model data :class:`~vacumm.data.misc.dataset.Dataset`
            - **profiles**: profile data :class:`~vacumm.data.misc.profile.ProfilesDataset`
            - **varnames**: variables to load (ex: ('temp','sal') or (('temp','temperature'),('sal','salinity'))
            - **select**: selector
            - **method**: coloc method (**nearest** or **interp**)

        :Return:
            - **lons_mod**: model longitude coordinates, shape: (profile)
            - **lats_mod**: model latitude coordinates, shape: (profile)
            - **deps_mod**: model depth coordinates, shape: (level,profile)
            - **var1**: requested variables, shape: (level,profile)
            - ...
            - **varN**

        .. todo::
            - also load and return profile data here
            - exclude coords where profile data is masked (no data for specified depth)
            - return time coordinates
            - return depth and vars with shape (profile,level)

        '''

        self.verbose('Colocalizing %s on %s\nvarnames: %s\nselect: %s\n method: %s', model.__class__.__name__, profiles.__class__.__name__, varnames, select, method)
        prof_pro = profiles.get_axis('profile', select=select)
        if prof_pro is None or not len(prof_pro):
            raise Exception('No profiles found, aborting')
        lev_pro = profiles.get_axis('level', select=select)
        time_pro = profiles.get_variable('time', select=select)
        lons_pro = profiles.get_variable('longitude', select=select)
        lats_pro = profiles.get_variable('latitude', select=select)
        dates = create_time(time_pro).asComponentTime()

        self.info('Number of profiles: %s', len(dates))
        self.info('Profiles time coverage: %s to %s', dates[0], dates[-1])

        # Init model
        td = model.get_time_res()
        dtmax = (td.days*86400+td.seconds, 'seconds')
        self.info('Detected model time step: %s', td)
        grid_mod = model.get_grid()
        xres, yres = resol(grid_mod)
        time_mod = model.get_time()
        ctime_mod = time_mod.asComponentTime()
        self.info('Model time coverage: %s to %s', ctime_mod[0], ctime_mod[-1])

        level_mod = model.get_level(select=select)
        lons_mod = MV2.zeros((len(prof_pro),))+MV2.masked
        lats_mod = lons_mod.clone()
        deps_mod = MV2.zeros((len(level_mod), len(prof_pro)))+MV2.masked
        deps_mod.setAxis(1, prof_pro)
        lons_mod.id, lats_mod.id, deps_mod.id = 'longitude', 'latitude', 'depth'

        # Creation des variables demandees
        variables = []
        for n in varnames:
            v = MV2.zeros((len(level_mod), len(prof_pro)))+MV2.masked
            v.setAxis(1, prof_pro)
            v.id = is_iterable(n) and n[0] or n
            variables.append(v)

        cdms2.setAutoBounds(1) # ???

        # Boucle temporelle
        for ip, date in enumerate(dates):
            try:
                # Limites spatiales
                lon = lons_pro[ip]
                lat = lats_pro[ip]
                lon_min = lon-2*xres
                lon_max = lon+2*xres
                lat_min = lat-2*yres
                lat_max = lat+2*yres
                date_interval = (add_time(date, - dtmax[0], dtmax[1]), add_time(date, dtmax[0], dtmax[1]), 'ccb')

                self.info('Colocalizing data for date %s, lon: %s, lat: %s', date, lon, lat)

                # Methode 1 : donnees les plus proches
                if method == 'nearest':
                    sel = dict(time=(date, date, 'ccb'), longitude=(lon, lon, 'ccb'), latitude=(lat, lat, 'ccb'))
                    # Verifier la disponibilite des donnees
                    if time_mod.mapIntervalExt(sel['time']) is None:
                        self.warning('Time interval %s not found', sel['time'])
                        continue
                    if grid_mod.getLatitude().mapInterval(sel['latitude']) is None:
                        self.warning('Latitude coordinate %s not found', sel['latitude'])
                        continue
                    if grid_mod.getLongitude().mapInterval(sel['longitude']) is None:
                        self.warning('Longitude coordinate %s not found', sel['longitude'])
                        continue

                    # Load tmp depth to get lon & lat coordinates
                    #tmp = model.get_depth(select=sel, squeeze=False) # tmp squeezed !!! see sigma ?
                    tmp = model.get_variable(varnames[0], select=sel, squeeze=False)
                    lons_mod[ip] = tmp.getLongitude()[0]
                    lats_mod[ip] = tmp.getLatitude()[0]
                    deps_mod[:, ip] =  model.get_depth(select=sel, squeeze=True)
                    for iv,vn in enumerate(varnames):
                        variables[iv][:,ip] = model.get_variable(vn, select=sel, squeeze=True)

                # Methode 2 : interpolation
                elif method == 'interp':
                    sel = dict(time=date_interval, longitude=(lon_min, lon_max), latitude=(lat_min, lat_max))
                    if time_mod.mapIntervalExt(sel['time']) is None:
                        self.warning('Time interval %s not found', sel['time'])
                        continue
                    if grid_mod.getLatitude().mapInterval(sel['latitude']) is None:
                        self.warning('Latitude coordinate %s not found', sel['latitude'])
                        continue
                    if grid_mod.getLongitude().mapInterval(sel['longitude']) is None:
                        self.warning('Longitude coordinate %s not found', sel['longitude'])
                        continue

                    # Lectures
                    order = 'tzyx'
                    # lon & lat du profile car interp sur cette position
                    lons_mod[ip], lats_mod[ip] = lon, lat
                    deps_mod_tzyx = model.get_depth(select=sel, order=order, squeeze=True)
                    tmp_tzyx = []
                    for iv,vn in enumerate(varnames):
                        tmp_tzyx.append(model.get_variable(vn, select=sel, order=order, squeeze=True))

                    # Interpolations temporelles
                    mctime = tmp_tzyx[0].getTime()
                    mrtime = mctime.asRelativeTime()
                    d0 = date.torel(mctime.units).value - mrtime[0].value
                    d1 = mrtime[1].value - date.torel(mctime.units).value
                    f0 = d0 / (d0 + d1)
                    f1 = d1 / (d0 + d1)
                    deps_mod_zyx = f0 * deps_mod_tzyx[0] + f1 * deps_mod_tzyx[1]
                    tmp_zyx = []
                    for iv,vn in enumerate(varnames):
                        tmp_zyx.append(f0 * tmp_tzyx[iv][0] + f1 * tmp_tzyx[iv][1])
                    del tmp_tzyx

                    # Interpolations spatiales
                    deps_mod[:,ip] = numpy.squeeze(grid2xy(deps_mod_zyx, numpy.array([lon]), numpy.array([lat]), method='nat'))
                    for iv,vn in enumerate(varnames):
                        variables[iv][:,ip] = numpy.squeeze(grid2xy(tmp_zyx[iv], numpy.array([lon]), numpy.array([lat]), method='nat'))
                    del tmp_zyx

                else:
                    raise ValueError('Invalid colocation method: %s'%(method))
            except:
                self.exception('Failed to colocalize data for date %s', date)

        for v in [deps_mod] + variables:
            v.getAxis(0).id = 'level'
            v.getAxis(0).designateLevel()

        data = tuple([lons_mod, lats_mod, deps_mod] + variables)
        self.verbose('Colocalized data:\n  %s', '\n  '.join(self.describe(o) for o in data))
        return data


    def coloc_strat_mod_on_pro(self, model, profiles, select, **kwargs):
        '''Get colocalized stratification data of a model on profiles

        See: :func:`coloc_mod_on_pro`

        :Return:
            - **lons_mod, lats_mod, deps_mod**:
            - **temp_mod, sal_mod**: temperature and salinity
            - **pres_mod, dens_mod**: pressure and density

        '''
        # Coloc donnees modele sur profiles
        varnames = (('temp','temperature'), ('sal', 'psal','salinity'))
        lons_mod, lats_mod, deps_mod, temp_mod, sal_mod = \
            self.coloc_mod_on_pro(model, profiles, varnames, select, **kwargs)
        # Calcul pression & densite
        pres_mod = MV2.zeros(temp_mod.shape)+MV2.masked
        pres_mod.setAxisList(temp_mod.getAxisList())
        dens_mod = pres_mod.clone()
        for ip in xrange(len(lons_mod)):
            pres_mod[:, ip] = seawater.csiro.pres(deps_mod[:, ip], numpy.resize([lats_mod[ip]], deps_mod[:, ip].shape))
            dens_mod[:, ip] = seawater.csiro.dens(sal_mod[:, ip], temp_mod[:, ip], pres_mod[:, ip])
        pres_mod.id = 'pressure'
        dens_mod.id = 'density'
        return lons_mod, lats_mod, deps_mod, temp_mod, sal_mod, pres_mod, dens_mod


    def plot_layer_mod_on_pro(self, model, profiles, varname, depth, select=None, **kwargs):
        '''Get a layer of variable for a specified depth.

        :Params:

          - **varname**: variable to process
          - **depth**: output depth(s)

          Other params, see: :func:`coloc_mod_on_pro`

        '''
        self.verbose('Plotting layer of colocalized %s on %s\nvarname: %s\ndepth: %s\nselect: %s\nkwargs: %s',
            model.__class__.__name__, profiles.__class__.__name__, varname, depth, select, kwargs)

        lons_mod, lats_mod, deps_mod, var_mod = \
            self.coloc_mod_on_pro(
                model, profiles,
                # If varname is a list of possible names, wrap it into a
                # tuple of one element as we are requiring only one variable
                len(numpy.shape(varname)) and (varname,) or varname,
                select, **kwargs)

        odep = create_dep(is_iterable(depth) and depth or [depth])
        var_mod = var_mod.reorder('-z')
        var_mod = interp1d(var_mod, axo=odep, xmap=0, xmapper=deps_mod) # Required order: ...z
        var_mod = var_mod.reorder('z-')
        model.verbose('layer:  %s', var_mod)

        lons_pro, lats_pro, var_pro = profiles.get_layer(varname, depth, select=select)
        profiles.verbose('layer:  %s', var_pro)

        # =====
        # Tracés
        vmin = var_pro.min()
        vmax = var_pro.max()
        if var_mod.count():
            vmin = min(var_mod.min(), vmin)
            vmax = max(var_mod.max(), vmax)
        else:
            self.warning('No model data')
        if 'latitude' in select:
            lat_min, lat_max = select['latitude'][:2]
        else:
            la = model.get_latitude()
            lat_min, lat_max = min(la), max(la)
        if 'longitude' in select:
            lon_min, lon_max = select['longitude'][:2]
        else:
            lo = model.get_longitude()
            lon_min, lon_max = min(lo), max(lo)

        levels = auto_scale((vmin, vmax))
        vmin = levels[0]
        vmax = levels[-1]
        # Création de la palette
        cmap = cmap_magic(levels)
        # On trace de champ du modèle moyené sur la période choisie
        m = map2(lon=(lon_min, lon_max), lat=(lat_min, lat_max), show=False)
        # On trace le modèle en fond

        # =====
        msca = None
        try: msca = m.map.scatter(lons_mod, lats_mod, s=100, c=var_mod, vmin=vmin, vmax=vmax, cmap=cmap, label='model')
        except: self.exception('Failed plotting model data')
        # =====
        # Triangulation du modèle
#        import matplotlib.tri as tri
#        triang = tri.Triangulation(lons_mod, lats_mod)
#        # On trace le modèle en fond
#        mod = m.axes.tripcolor(triang, var_mod, vmin=vmin, vmax=vmax, cmap=cmap)
        # =====

        # On trace les observations avec des points plus petits
        psca = None
        try: psca = m.map.scatter(lons_pro, lats_pro, s=25, c=var_pro, vmin=m.vmin, vmax=m.vmax, cmap=m.cmap, label='profiles')
        except: self.exception('Failed plotting profile data')

        # Colorbar
        if msca is not None:
            colorbar(msca)
        elif psca is not None:
            colorbar(psca)

    def plot_mld_mod_on_pro(self, model, profiles, select=None, deep=False, **kwargs):
        '''Plot mixed layer depth of model data correspponding to profiles position

        :Params:

          - **deep**: deep water computation mode if true

          Other params, see: :func:`coloc_mod_on_pro`

        '''
        self.verbose('Plotting MLD of colocalized %s on %s\nselect: %s\ndeep: %s\nkwargs: %s',
            model.__class__.__name__, profiles.__class__.__name__, select, deep, kwargs)

        # =====
        # Lecture des donnees de stratification
        lons_mod, lats_mod, deps_mod, temp_mod, sal_mod, pres_mod, dens_mod = \
            self.coloc_strat_mod_on_pro(model, profiles, select, **kwargs)

        # =====
        # Calcul mld modele

        # MLD en eaux peu profondes
        if not deep:
            ddeps = meshweights(deps_mod, axis=0)
            # Densité min/max
            dmin = dens_mod.min(axis=0)
            dmax = dens_mod.max(axis=0)
            # Densité moyenne
            dmean = numpy.ma.average(dens_mod, axis=0, weights=ddeps)
            # Profondeur max (max+demi épaisseur)
            H = deps_mod[-1]+ddeps[-1]*.5
            # MLD
            mld_mod = H*(dmax-dmean)/(dmax-dmin)
        # MLD en eaux profondes
        else:
            # Profondeur de référence
            dep = -10
            # Axe pour interpolation
            from vacumm.misc.axes import create_dep
            depaxis = create_dep([dep])
            # Interpolations
            dens_mod_ref = dens_mod_ref.reorder('-z')
            dens_mod_ref = regrid1dold(dens_mod, axo=depaxis, xmap=0, xmapper=deps_mod) # Required order: ...z
            dens_mod_ref = dens_mod_ref.reorder('z-')
            # Valeur du différentiel de densité (cf. de Boyer Montégut et all, 2003)
            delta_dens = 0.03
            # Densité cible
            dens_mod_target = dens_mod_ref+0.03
            # Masques de la couche mélangée et des eaux profondes
            dens_mod_target3d = MV2.resize(dens_mod_target, dens_mod.shape)
            dens_mod_good = (dens_mod.asma() <= dens_mod_target3d.asma()).filled(False)
            dens_mod_bad = (dens_mod.asma() > dens_mod_target3d.asma()).filled(False)
            # Profondeurs juste au dessus et en dessous de la MLD
            deps_mod_above = MV2.masked_where(dens_mod_bad, deps_mod).min(axis=0)
            deps_mod_below = MV2.masked_where(dens_mod_good, deps_mod).max(axis=0)
            # Masques associés
            from vacumm.misc import closeto
            mask_above = closeto(deps_mod, MV2.resize(deps_mod_above, deps_mod.shape))
            mask_below = closeto(deps_mod, MV2.resize(deps_mod_below, deps_mod.shape))
            # Densités juste au dessus et en dessous de la MLD
            dens_mod_above = MV2.masked_where(dens_mod, mask_above).max(axis=0)
            dens_mod_below = MV2.masked_where(dens_mod, mask_below).min(axis=0)
            # Interpolation
            dens_mod_delta = dens_mod_above-dens_mod_below
            deps_mod_mld = (dens_mod_target-dens_mod_above)*deps_mod_above
            deps_mod_mld += (dens_mod_below-dens_mod_target)*deps_mod_below
            deps_mod_mld = MV2.where(dens_mod_delta.mask, MV2.masked, deps_mod_mld/dens_mod_delta)
            mld_mod = deps_mod_mld
        # Finalize
        mld_mod = MV2.array(mld_mod)
        mld_mod.units = 'm'
        mld_mod.long_name = u'Profondeur de la couche de melange'
        mld_mod.setAxisList(lons_mod.getAxisList())
        model.verbose('MLD:  %s', mld_mod)

        # =====
        # Calcul mld profiles
        mld_pro, lats_pro, lons_pro = profiles.get_mld(select)
        profiles.verbose('MLD:  %s', mld_pro)

        # =====
        # Tracés
        vmin = mld_pro.min()
        vmax = mld_pro.max()
        if mld_mod.count():
            vmin = min(mld_mod.min(), vmin)
            vmax = max(mld_mod.max(), vmax)
        else:
            self.warning('No model data')
        if 'latitude' in select:
            lat_min, lat_max = select['latitude'][:2]
        else:
            la = model.get_latitude()
            lat_min, lat_max = min(la), max(la)
        if 'longitude' in select:
            lon_min, lon_max = select['longitude'][:2]
        else:
            lo = model.get_longitude()
            lon_min, lon_max = min(lo), max(lo)

        levels = auto_scale((vmin, vmax))
        vmin = levels[0]
        vmax = levels[-1]
        # Création de la palette
        cmap = cmap_magic(levels)
        # On trace de champ du modèle moyené sur la période choisie
        m = map2(lon=(lon_min, lon_max), lat=(lat_min, lat_max), show=False)
        # On trace le modèle en fond

        # =====
        msca = None
        try: msca = m.map.scatter(lons_mod, lats_mod, s=100, c=mld_mod, vmin=vmin, vmax=vmax, cmap=cmap, label='model')
        except: self.exception('Failed plotting model data')
        # =====
        # Triangulation du modèle
#        import matplotlib.tri as tri
#        triang = tri.Triangulation(lons_mod, lats_mod)
#        # On trace le modèle en fond
#        mod = m.axes.tripcolor(triang, mld_mod, vmin=vmin, vmax=vmax, cmap=cmap)
        # =====

        # On trace les observations avec des points plus petits
        psca = None
        try: psca = m.map.scatter(lons_pro, lats_pro, s=25, c=mld_pro, vmin=m.vmin, vmax=m.vmax, cmap=m.cmap, label='profiles')
        except: self.exception('Failed plotting profiles data')

        # Colorbar
        if msca is not None:
            colorbar(msca)
        elif psca is not None:
            colorbar(psca)


    def plot_ped_mod_on_pro(self, model, profiles, select=None, **kwargs):
        '''Plot potential energy deficit of model data correspponding to profiles position

        :Params: See: :func:`coloc_mod_on_pro`

        '''
        self.verbose('Plotting PED of colocalized %s on %s\nselect: %s\nkwargs: %s',
            model.__class__.__name__, profiles.__class__.__name__, select, kwargs)

        # =====
        # Lecture des donnees de stratification
        lons_mod, lats_mod, deps_mod, temp_mod, sal_mod, pres_mod, dens_mod = \
            self.coloc_strat_mod_on_pro(model, profiles, select, **kwargs)

        # =====
        # Calcul ped modele
        ddeps = meshweights(deps_mod, axis=0)
        # Densité moyenne
        dmean = MV2.average(dens_mod, axis=0, weights=ddeps)
        # Anomalie de densité
        danom = dens_mod-dmean
        # Énergie potentielle disponible
        ape = danom * g
        ape *= ddeps
        # Deficit
        ped_mod = MV2.average(ape, axis=0, weights=ddeps)
        ped_mod.units = 'J.m^{-2}'
        ped_mod.long_name = u"Definit d'energie potentielle"
        ped_mod.setAxisList(lons_mod.getAxisList())
        model.verbose('PED:  %s', ped_mod)

        # =====
        # Calcul ped profiles
        ped_pro, lats_pro, lons_pro = profiles.get_ped(select)
        profiles.verbose('PED:  %s', ped_pro)

        # =====
        # Tracés
        vmin = ped_pro.min()
        vmax = ped_pro.max()
        if ped_mod.count():
            vmin = min(ped_mod.min(), vmin)
            vmax = max(ped_mod.max(), vmax)
        else:
            self.warning('No model data')
        if 'latitude' in select:
            lat_min, lat_max = select['latitude'][:2]
        else:
            la = model.get_latitude()
            lat_min, lat_max = min(la), max(la)
        if 'longitude' in select:
            lon_min, lon_max = select['longitude'][:2]
        else:
            lo = model.get_longitude()
            lon_min, lon_max = min(lo), max(lo)

        levels = auto_scale((vmin, vmax))
        vmin = levels[0]
        vmax = levels[-1]
        # Création de la palette
        cmap = cmap_magic(levels)
        # On trace de champ du modèle moyené sur la période choisie
        m = map2(lon=(lon_min, lon_max), lat=(lat_min, lat_max), show=False)
        # On trace le modèle en fond
        msca = None
        try: msca = m.map.scatter(lons_mod, lats_mod, s=100, c=ped_mod, vmin=vmin, vmax=vmax, cmap=cmap, label='model')
        except: self.exception('Failed plotting model data')
        # On trace les observations avec des points plus petits
        psca = None
        try: psca = m.map.scatter(lons_pro, lats_pro, s=25, c=ped_pro, vmin=m.vmin, vmax=m.vmax, cmap=m.cmap, label='profiles')
        except: self.exception('Failed plotting profiles data')

        # Colorbar
        if msca is not None:
            colorbar(msca)
        elif psca is not None:
            colorbar(psca)




