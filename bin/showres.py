#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Show the grid resolution of netcdf file"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] ncfile",
    description="Show the grid (rectangular or curvilinear) resolution of netcdf file.")
parser.add_option('-v', '--var', action='store', dest='vname',
    help='name of the variable from which to get grid')
parser.add_option('-m', '--meters', action='store_true', dest='meters',
    default=False,
    help='resolution in meters? [default: False -> degrees]')

# Parse
(options, args) = parser.parse_args()
if len(args)==0:
    parser.error('You must provide a valid netcdf file as first argument')

# Netcdf file
ncfile = args[0]
if not ncfile.startswith('http://') and not os.path.exists(ncfile):
    sys.exit('File not found: '+ncfile)
try:
    import cdms2
    f = cdms2.open(ncfile)
except:
    sys.exit("Can't open "+ncfile)

# Guess grid
grid = None
if options.vname is not None: # from specified variable
    if options.vname not in f.listvariables():
        sys.exit('Variable "%s" not found'%options.vname)
    grid = f[options.vname].getGrid()
    if grid is None:
        sys.exit('Variable "%s" has no grid'%options.vname)
else: # loop on all variables
    for vname in f.listvariables():
        grid = f[vname].getGrid()
        if grid is not None: break
    else:
        sys.exit('No variable with a grid found')
grid = (grid.getLongitude().clone(), grid.getLatitude().clone())

# VACUMM path
here = os.path.dirname(__file__)
for rpath in [os.path.join(here, '../lib/python')]:
    path = os.path.abspath(rpath)
    sys.path.append(path)

# Resolution
try:
    from vacumm.misc.grid import resol
except:
    sys.exit('Error when loading vacumm')
xres, yres = resol(grid, proj=options.meters)
if options.meters:
    units = 'km'
    xres /= 1000.
    yres /= 1000.
else:
    units = '°'
print 'dx=%(xres)g%(units)s, dy=%(yres)g%(units)s'%locals()
f.close()
