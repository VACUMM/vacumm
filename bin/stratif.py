#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar (2010)
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
__doc__ = ''


import optparse

import pylab

from vacumm.data.misc.coloc import Colocator
from vacumm.data.misc.profile import ProfilesDataset
from vacumm.data.model.mars3d import MARS3D
from vacumm.misc.bases import Object
from vacumm.misc.atime import Intervals, strftime
from vacumm.misc.log import default as log


#
# ./stratif.py --cfgfile myconf.cfg -t 2004-01-01,2004-01-15,ccb,7,days
#


if __name__ == '__main__':
    
    parser = optparse.OptionParser(
        description='',
        usage='%prog [options]',
        version=__date__)
    parser.add_option('--cfgfile', action='store', dest='cfgfile', default=None, help='Configuration file')
    parser.add_option('-t', '--time', action='store', dest='time', default=None, metavar='min,max,[bb],step,unit', help='Time selection:'
        '\n- min,max: specify the time range to operate.'
        '\n- bb: optionnal, time bounds open/closed selection.'
        '\n- step,unit: period covered by each plot.'
        '\nEx:  "2001-01,2001-01-15T00,7,days"'
        '\n     "2001-06,2001-09,co,1,month"')
    parser.add_option('-b', '--bbox', action='store', dest='bbox', default=None, metavar='lonmin,latmin,lonmax,latmax',  help='Restrict processed zone to the specified bounding box')
    colocs = ('nearest','interp')#colocs = (None,'nearest','interp')
    parser.add_option('--coloc', action='store', dest='coloc', default=None, choices=colocs, help='Use colocation method, one of %s, default is %%default)'%(colocs,))
    parser.add_option('--deep', action='store_true', dest='deep', default=None, help='Deep water mode')
    plots = ('mld','ped')
    parser.add_option('--plot', action='append', choices=plots, dest='plots', default=None, help='Specify one or more figure in %s, default is %%default'%(plots,))
    parser.add_option('--show', action='store_true', dest='show', default=None, help='Show figures')
    parser.add_option('-o', '--output', action='store', dest='output', default='%(plot)s-%(tmin)s-%(tmax)s.png', metavar='pattern',  help='Output files pattern (default: %default)')
    options, args = parser.parse_args()
    
    if options.time is None: parser.error('Missing time parameter')
    if options.bbox:
        bbox = options.bbox.split(',')
        if len(bbox) != 4: parser.error('Invalid bbox parameter: %s', options.bbox)
        lonmin, latmin, lonmax, latmax = map(float, bbox)
    if options.plots is None:
        #parser.error('You must specify at least one plot command')
        options.plots = plots
    if options.cfgfile:
        Object.load_default_config(options.cfgfile, nested=True)
        MARS3D.load_default_config(Object.get_default_config())
        ProfilesDataset.load_default_config(Object.get_default_config())
    
    time = options.time.split(',')
    if len(time) == 4: tmin, tmax, tstep, ustep = time; tbb = 'co'
    elif len(time) == 5: tmin, tmax, tbb, tstep, ustep = time
    else: parser.error('Invalid time parameter - 4 (tmin, tmax, tstep, ustep) or 5 (tmin, tmax, tbb, tstep, ustep) parameters are needed.')
    tstep = int(tstep)
    dtfmt = '%Y%m%dT%H%M%S'
    select=dict()
    
    mars = MARS3D()
    profiles = ProfilesDataset()
    if options.cfgfile: 
        mars.load_default_config(options.cfgfile, nested=True)
        profiles.load_default_config(options.cfgfile, nested=True)
        mars.load_config(options.cfgfile, nested=True)
        profiles.load_config(options.cfgfile, nested=True)
    if not options.bbox:
        #mars_grid = mars.get_grid()
        lat, lon = mars.get_latitude(), mars.get_longitude()
        lonmin, latmin, lonmax, latmax = min(lon), min(lat), max(lon), max(lat)
    
    for itv in Intervals((tmin, tmax, tbb),(tstep, ustep)):
        log.notice('Interval %s', itv)
        
        select = dict(time=itv, latitude=(latmin,latmax,'ccb'), longitude=(lonmin,lonmax,'ccb'))
        plotkw = dict(plot_show=False)
        
        # ==========
        # Methode 1. Tracé des observations (points) sur fond 2D grillés de champs de modèles moyennés sur la période choisie.
        if not options.coloc:
            
            if 'mld' in options.plots:
                pylab.figure()
                kw = plotkw.copy()
                m = c = None
                try: m,c = mars.plot_mld(select, **kw)
                except: log.exception('Cannot plot mars mld')
                kw.update(map=m, cmap=c)
                try: profiles.plot_mld(select, **kw)
                except: log.exception('Cannot plot profiles mld')
                pylab.legend(loc='best')
                pylab.title('Mixed Layer Depth\ntime: %s'%(itv,))
                output = options.output%dict(plot='mld', tmin=strftime(dtfmt, itv[0]), tmax=strftime(dtfmt, itv[1]))
                log.notice('Saving %s', output)
                pylab.savefig(output)
            
            if 'ped' in options.plots:
                pylab.figure()
                kw = plotkw.copy()
                m = c = None
                try: m,c = mars.plot_ped(select, **kw)
                except: log.exception('Cannot plot mars ped')
                kw.update(map=m, cmap=c)
                try: profiles.plot_ped(select, **kw)
                except: log.exception('Cannot plot profiles ped')
                pylab.legend(loc='best')
                pylab.title('Potential Energy Deficit\ntime: %s'%(itv,))
                output = options.output%dict(plot='ped', tmin=strftime(dtfmt, itv[0]), tmax=strftime(dtfmt, itv[1]))
                log.notice('Saving %s', output)
                pylab.savefig(output)
    
        # ==========
        # Methode 2. Tracé des observations et des champs du modèles colocalisés avec les observations à l’aide d’un gros point
        #    (ou carré) pour le modèle et d’un point normal pour les observations.
        else:
            
            if 'mld' in options.plots:
                pylab.figure()
                try: Colocator().plot_mld_mod_on_pro(mars, profiles, select, method=options.coloc, deep=options.deep)
                except: log.exception('Cannot plot mld')
                output = options.output%dict(plot='mld', tmin=strftime(dtfmt, itv[0]), tmax=strftime(dtfmt, itv[1]))
                pylab.legend(loc='best')
                pylab.title('Mixed Layer Depth\ntime: %s'%(itv,))
                log.notice('Saving %s', output)
                pylab.savefig(output)
            
            if 'ped' in options.plots:
                pylab.figure()
                try: Colocator().plot_ped_mod_on_pro(mars, profiles, select, method=options.coloc)
                except: log.exception('Cannot plot ped')
                pylab.legend(loc='best')
                pylab.title('Potential Energy Deficit\ntime: %s'%(itv,))
                output = options.output%dict(plot='ped', tmin=strftime(dtfmt, itv[0]), tmax=strftime(dtfmt, itv[1]))
                log.notice('Saving %s', output)
                pylab.savefig(output)
                
    if options.show: pylab.show()
    


