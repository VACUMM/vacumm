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

from vacumm.misc.bases import Object
from vacumm.data.model.mars3d import MARS3D
from vacumm.data.misc.satellite import Satellite
from vacumm.misc.atime import Intervals, strftime
from vacumm.misc.log import default as log
from vacumm.misc.plot import xrotate

#
# ./localize.py --cfgfile myconf.cfg -v TEMP -t 2004-01-01,2004-01-15,ccb,7,days -c -7,-6,47
#


if __name__ == '__main__':
    
    parser = optparse.OptionParser(
        description='',
        usage='%prog [options]',
        version=__date__)
    parser.add_option('--cfgfile', action='store', dest='cfgfile', default=None, help='Configuration file')
    parser.add_option('-v', '--variables', action='append', dest='variables', default=None, metavar='varname[,alias]*', help='Variables to be processed. Use aliases when varname differ between datasets. This option may be repeated to produce figures for each variable definition.')
    parser.add_option('-t', '--time', action='store', dest='time', default=None, metavar='min,max,[bb],step,unit', help='Time selection:'
        '\n- min,max: specify the time range to operate.'
        '\n- bb: optionnal, time bounds open/closed selection.'
        '\n- step,unit: period covered by each plot.'
        '\nEx:  "2001-01,2001-01-15T00,7,days"'
        '\n     "2001-06,2001-09,co,1,month"')
    parser.add_option('-l', '--level', action='store', dest='level', default=-1, metavar='level', help='Level indice for model data (0 based, default: %default)')
    parser.add_option('-m', '--meridional', action='store_true', dest='meridional', default=False,  help='Specify the section is meridional, at a longitude and along given latgitude range (default is zonal)')
    parser.add_option('-c', '--coords', action='store', dest='coords', default=None, metavar='xy1,xy2,xy',  help='Hovmoller\'s section coordinates xy1 and xy2 are longitudes and xy a latitude if the section is zonal, latitudes and a longitude if the section is meridonal')
    parser.add_option('-b', '--bbox', action='store', dest='bbox', default=None, metavar='lonmin,latmin,lonmax,latmax',  help='Restrict processed zone to the specified bounding box')
    ops = ('min','max')
    parser.add_option('--op', action='append', choices=ops, dest='operations', default=None, help='Specify one or more operation in %s, default is %%default'%(ops,))
    plots = ('loc','val')
    parser.add_option('--plot', action='append', choices=plots, dest='plots', default=None, help='Specify one or more figure in %s, default is %%default'%(plots,))
    parser.add_option('-o', '--output', action='store', dest='output', default='%(plot)s-%(op)s-%(var)s-%(tmin)s-%(tmax)s.png', metavar='pattern',  help='Output files pattern  (default: %default)')
    parser.add_option('--show', action='store_true', dest='show', default=None, help='Show figures')
    options, args = parser.parse_args()
    
    if options.variables is None: parser.error('Missing variables parameter')
    if options.time is None: parser.error('Missing time parameter')
    if options.coords is None: parser.error('Missing coords parameter')
    if options.bbox:
        bbox = options.bbox.split(',')
        if len(bbox) != 4: parser.error('Invalid bbox parameter: %s', options.bbox)
        lonmin, latmin, lonmax, latmax = map(float, bbox)
    if options.operations is None:
        options.operations = ops
        #parser.error('You must specify at least one operation command')
    if options.plots is None:
        #parser.error('You must specify at least one plot command')
        options.plots = plots
    if options.cfgfile:
        Object.load_default_config(options.cfgfile, nested=True)
        MARS3D.load_default_config(Object.get_default_config())
        Satellite.load_default_config(Object.get_default_config())
    
    variables = map(lambda v:v.split(','), options.variables)
    time = options.time.split(',')
    if len(time) == 4: tmin, tmax, tstep, ustep = time; tbb = 'co'
    elif len(time) == 5: tmin, tmax, tbb, tstep, ustep = time
    else: parser.error('Invalid time parameter')
    tstep = int(tstep)
    coords = options.coords.split(',')
    if len(coords) != 3: parser.error('Invalid coords parameter: %s', options.coords)
    xymin, xymax, xy = map(float, coords)
    dtfmt = '%Y%m%dT%H%M%S'
    
    mars = MARS3D()
    sat = Satellite()
    if options.cfgfile:
        mars.load_default_config(options.cfgfile, nested=True)
        sat.load_default_config(options.cfgfile, nested=True)
        mars.load_config(options.cfgfile, nested=True)
        sat.load_config(options.cfgfile, nested=True)
    
    for varname in variables:
        log.notice('Variable %s', varname)
        
        for itv in Intervals((tmin, tmax, tbb),(tstep, ustep)):
            log.notice('Interval %s', itv)
            
            select = dict(time=itv[:2])
            if options.bbox: select.update(dict(latitude=(latmin,latmax), longitude=(lonmin,lonmax)))
            datakw = dict(
                varname=('temperature', 'TEMP'),
                select=select,
                xorymin=xymin, xorymax=xymax, xory=xy, meridonal=options.meridional,
                cur_axes_rect=[0.1, 0.1, 0.6, 0.7], map_axes_rect=[0.75, 0.1, 0.2, 0.3],
                plot_show=False)
            
            if options.meridional:
                lat = ([xymin], [xymax])
                lon = (xy, xy)
            else:
                lon = ([xymin], [xymax])
                lat = (xy, xy)
            
            if 'loc' in options.plots:
                for op in options.operations:
                    log.notice('loc %s', op)
                    pylab.figure()
                    xlkw = datakw.copy()
                    xlkw.update(extrema=op)
                    
                    xlkw.update(pmap=False, cur_color='b', cur_linestyle='', cur_marker='o')
                    try: cur, var, lat, lon = sat.plot_extrema_location(**xlkw)
                    except: log.exception('Cannot plot satellite')
                    
                    xlkw.update(pmap=False, cur_color='g', cur_linestyle='-', cur_marker='+')
                    xlkw['select'].update(level=slice(int(options.level),None))
                    try: cur, var, lat, lon = mars.plot_extrema_location(**xlkw)
                    except: log.exception('Cannot plot mars')
                    
                    xrotate(45)
                    pylab.legend(loc='best')
                    
                    try: mars.plot_trajectory_map(lon, lat, **xlkw)
                    except: log.exception('Cannot plot trajectory')
                    
                    output = options.output%dict(plot='loc', op=op, var=varname[0], tmin=strftime(dtfmt, itv[0]), tmax=strftime(dtfmt, itv[1]))
                    log.info('Saving %s', output)
                    pylab.savefig(output)
            
            if 'val' in options.plots:
                for op in options.operations:
                    log.notice('val %s', op)
                    pylab.figure()
                    lckw = datakw.copy()
                    lckw.update(operation=op, cur_vmin=0, cur_vmax=20)
                    
                    lckw.update(pmap=False, cur_color='b', cur_linestyle='', cur_marker='o')
                    try: cur, var, lat, lon = sat.plot_localized_computed_values(**lckw)
                    except: log.exception('Cannot plot satellite')
                    
                    lckw.update(pmap=False, cur_color='g', cur_linestyle='-', cur_marker='+')
                    lckw['select'].update(level=slice(int(options.level),None))
                    try: cur, var, lat, lon = mars.plot_localized_computed_values(**lckw)
                    except: log.exception('Cannot plot mars')
                    
                    xrotate(45)
                    pylab.legend(loc='best')
                    
                    try: mars.plot_trajectory_map(lon, lat, **lckw)
                    except: log.exception('Cannot plot trajectory')
                    
                    output = options.output%dict(plot='val', op=op, var=varname[0], tmin=strftime(dtfmt, itv[0]), tmax=strftime(dtfmt, itv[1]))
                    log.info('Saving %s', output)
                    pylab.savefig(output)
    
    if options.show: pylab.show()
    


