#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Plot a netcdf variable"""

# Arguments
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys, os, shutil
parser = ArgumentParser(description="Plot a netcdf variable.",
    formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('ncfile', help='netcdf file name or opendap url')
parser.add_argument('var', help="name of variable, up to three : var OR (u,v) OR (mod,u,v)", nargs='*')
parser.add_argument('-v', '--verbose', help='verbose mode', action='store_true')
gread = parser.add_argument_group('reading parameters')
gread.add_argument('-g', '--generic', help='name of a generic dataset class to use for reading')
gread.add_argument('-x', '--lon', help='longitude selection')
gread.add_argument('-y', '--lat', help='latitude selection')
gread.add_argument('-z', '--dep', help='depth selection')
gread.add_argument('-t', '--time', help='time selection')
gread.add_argument('-s', '--select', help='extra selection specifications', default='')
gread.add_argument('--dvar', help='depth variable name if not an axis')
goper = parser.add_argument_group('operators')
goper.add_argument('--add', help='add a constant', type=float)
goper.add_argument('--mul', help='multiply by a factor', type=float)
gplot = parser.add_argument_group('plotting parameters')
gplot.add_argument('-c', '--color', help='main color')
gplot.add_argument('-l', '--linewidth', help='main linewidth', type=float)
gplot.add_argument('-f', '--fill', help='filling mode for 2D', default='contourf',
    choices=['contourf', 'pcolormesh', 'pcolor', 'imshow', 'False'])
gplot.add_argument('--nocontour', help='do not plot contours for 2D', action='store_true')
gplot.add_argument('--nocb', help='do not add a colorbar', action='store_true')
gplot.add_argument('--title', help='change the title')
gplot.add_argument('--res', help='change the shoreline resolution',
    choices=['auto', 'c', 'l', 'i', 'h', 'f', 's', 'off'], default='auto')
gplot.add_argument('--cmap', help='colormap name')
gplot.add_argument('--vmin', help='min value')
gplot.add_argument('--vmax', help='max value')
gplot.add_argument('--qnorm', help='quiver norm', choices=[0, 1, 2, 3], default=0, type=int)
gplot.add_argument('--qonly', help='do not plot a scalar field with quiver', action='store_true')
gplot.add_argument('--clabel', help='add contour labels', action='store_true')
gplot.add_argument('--figsize', help='figure size')
gplot.add_argument('-p', '--plot', help='extra plot parameters', default='')
args = parser.parse_args()

# Imports
from vcmq import cdms2, curve2, map2, hov2, stick2, plot2d, DS, isdepthup, N, re
from vacumm.misc.atime import RE_MATCH_TIME

# Functions
def get_sel(arg):
    if arg is None: return
    arg = arg.strip()
#    if re.match(r'[\(\[].*', arg): # classic selector
#        return tuple(eval(arg))
    if RE_MATCH_TIME(arg) and not arg.isdigit(): # date
        return arg
    if re.match('[\w\s]*:[\w\s]*(:[\w\s]*)?$', arg): # args to slice
        return slice(*[(None if a=='' else eval(a)) for a in arg.split(':')])
    arg = eval(arg)
    if isinstance(arg, int): # slice at a single index
        return slice(arg, (arg+1) if arg!=-1 else None)
    if isinstance(arg, float): # selection of a single position
        return (arg, arg, 'cob')
    return arg

# Open reader or file
if args.generic:
    f = DS(args.ncfile, args.generic, logger_level='info' if args.verbose else 'error')
else:
    f = cdms2.open(args.ncfile)

# Find a variable name if not provided
if not args.var:
    if args.generic:
        allvars = f.dataset[0].listvariables()
    else:
        allvars = f.listvariables()
    for var in allvars:
        if 'bounds' not in var:
            if args.generic:
                var = '+'+var
            args.var = [var]
            break
else:
    args.var = args.var[:3]

# Reading parameters
kwread = eval('dict(%s)'%args.select)
kwread['squeeze'] = 1
for ss in 'lon', 'lat', 'dep', 'time':
    if getattr(args, ss) is not None:
        kwread[ss] = get_sel(getattr(args, ss))

# Read variable
vv = [f(vname, **kwread) for vname in args.var]
if args.dvar:
    depth = f(args.dvar, **kwread)
f.close()
if args.verbose:
    for v in vv:
        v.info()
order = vv[0].getOrder()
ndim = vv[0].ndim

# Reduce rank when needed
if ndim>2:
    for i in xrange(ndim-2):
        if order[i]=='z':
            islice = -1 if isdepthup(var[0].getAxis(0)) else 0
        else:
            islice = 0
        vv = [v[islice] for v in vv]
order = vv[0].getOrder()
ndim = vv[0].ndim
if args.dvar and depth.ndim>ndim:
    for i in xrange(depth.ndim-ndim):
        depth = depth[0]

# Operators
if args.add:
    for v in vv:
        v[:] += args.add
if args.mul:
    for v in vv:
        v[:] *= args.mul

# Plotting parameters
kwplot = eval('dict(%s)'%args.plot)
if args.color is not None and args.color.startswith('('):
    args.color = eval(args.color)
if args.fill=='False':
    args.fill = False
kwplot['contour'] = not args.nocontour
kwplot['colorbar'] = not args.nocb
if args.dvar:
    kwplot['yaxis'] = depth
for att in ('title', 'cmap', 'fill', 'clabel', 'vmin', 'vmax'):
    kwplot[att] = getattr(args, att)
for att in ('color', 'linewidth'):
    val = getattr(args, att)
    if val is not None:
        kwplot[att] = val
kwplot['quiver_norm'] = args.qnorm
if len(vv)==2 and args.qonly:
    kwplot['contour'] = kwplot['fill'] = False
if args.figsize:
    kwplot['figsize'] = tuple([float(v) for v in args.figsize.split(',')])
if args.res=='off':
    args.res = None
kwplot['res'] = args.res

#del kwplot['color']
#kwplot = {'color':None}

# Adaptative plot
if ndim==2:
    if order in ['yx', 'xy']:
        map2(vv, **kwplot)
    elif 't' in order:
        hov2(vv, **kwplot)
    else:
        plot2d(vv, **kwplot)
else:
    if len(vv)>1:
        stick2(*vv, **kwplot)
    else:
        curve2(vv, **kwplot)
