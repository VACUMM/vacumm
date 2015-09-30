#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Show grid info of netcdf file"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] ncfile",
    description="Show grid (rectangular or curvilinear) info of a netcdf file.")
parser.add_option('-v', '--var', action='store', dest='vname',
    help='Name of the variable from which to get grid')
parser.add_option('-p', '--plot', action='store_true', dest='plot',
    default=False,
    help='Plot the grid on a map [default: %default]')
parser.add_option('-o', '--out', action='store', dest='out',
    help='Name of an output file where to store the plot')
parser.add_option('-z', '--zoom', action='store', dest='zoom', type='float',
    help='Plot zoom [default: %default]', default=.8)
parser.add_option('-f', '--figsize', action='store', dest='figsize',
    help='Size of figure in inches, like "5,6" or simply "5"')
parser.add_option('-s', '--shoreline', action='store', dest='shoreline',  default='auto',
    help="GSHHS shoreline as one of 'c', 'l', 'i', 'h', 'f' [default: %default]")
parser.add_option('--xmin', action='store', dest='xmin', type='float',
    help="Min longitude of the map")
parser.add_option('--xmax', action='store', dest='xmax', type='float',
    help="Max longitude of the map")
parser.add_option('--ymin', action='store', dest='ymin', type='float',
    help="Min latitude of the map")
parser.add_option('--ymax', action='store', dest='ymax', type='float',
    help="Max latitude of the map")
parser.add_option('--proj', action='store', dest='proj', default='merc',
    help="Map projection [default: %default]")

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
except:
    sys.exit("Can't import python module: cdms2")
try:
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
lon = grid.getLongitude().clone()
lat = grid.getLatitude().clone()

# VACUMM path
here = os.path.dirname(__file__)
for rpath in [os.path.join(here, '../lib/python')]:
    path = os.path.abspath(rpath)
    sys.path.append(path)
try:
    from vacumm.misc.grid import resol, get_xy, create_grid
    from vacumm.misc import lonlab, latlab, scalebox
    from vacumm.misc.plot import map2, add_grid
    import cdms2
except:
    sys.exit('Error when loading vacumm')

# Dimension
gg = (lon, lat)
lonn, latn = get_xy(gg, num=True)
grid = create_grid(*gg)
print 'Dimensions           : nx=%i ny=%i'%grid.shape[::-1]
print 'Axes                 : lon="%s" lat="%s" (%iD)'%(lon.id, lat.id, lonn.ndim)


# Extension
print 'Zonal extent         : %s -> %s  [%s -> %s]'%(
    lonn.min(), lonn.max(), lonlab(lonn.min(), decimal=False), lonlab(lonn.max(), decimal=False))
print 'Meridional extent    : %s -> %s  [%s -> %s]'%(
    latn.min(), latn.max(), latlab(latn.min(), decimal=False), latlab(latn.max(), decimal=False))

# Resolution
lonres, latres = resol(gg, proj=False)
xres, yres = resol(gg, proj=True)
xres /= 1000.
yres /= 1000.
print 'Zonal resolution     : %g° / %gkm'%(lonres, xres)
print 'Meridional resolution: %g° / %gkm'%(latres, yres)

# Plot
if options.plot or options.out:

    xmin, ymin, xmax, ymax = scalebox(grid, 1/options.zoom)
    for att in 'xmin', 'xmax', 'ymin', 'ymax':
        if getattr(options, att) is not None: exec att+" = %s"%getattr(options, att)
    if options.figsize is not None:
        try:
            options.figsize = eval(options.figsize)
            if not isinstance(options.figsize, tuple):
                options.figsize = options.figsize, options.figsize
        except:
            options.figsize = None
    if options.proj.lower() == 'none': options.proj = None
    m = map2(lon=(xmin, xmax), lat=(ymin, ymax), show=False, res=options.shoreline,
        title='Grid from '+os.path.basename(ncfile), figsize=options.figsize,
        proj = options.proj)
    add_grid(grid, m=m, centers=True, borders=True, color='b', markersize=5)
    m.legend(zorder=200, alpha=.7, title='Grid: %ix%i'%grid.shape[::-1])
    if options.out:
        m.savefig(options.out)
        print 'Plot saved to '+ options.out
    if options.plot: m.show()
f.close()
