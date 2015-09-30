#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Show grid info of netcdf file"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] [ncfile]",
    description="Plot the bathymetry in an area. "
    "Value are considered negative under the sea and positive on land.")
parser.add_option('--xmin', '--x0', '--lonmin', action='store', dest='xmin', type='float',
    help="Min longitude")
parser.add_option('--xmax', '--x1', '--lonmax', action='store', dest='xmax', type='float',
    help="Max longitude")
parser.add_option('--ymin', '--y0', '--latmin', action='store', dest='ymin', type='float',
    help="Min latitude")
parser.add_option('--ymax', '--y1', '--latmax', action='store', dest='ymax', type='float',
    help="Max latitude")
parser.add_option('--ncfile', action='store', dest='ncfile',
    help='Name of netcdf file for a direct access [obsolete, use first argument instead]')
parser.add_option('-n', '--name', action='store', dest='name',
    help='Name of bathymetry as stored in configuration file')
parser.add_option('-c', '--cfgfile', action='store', dest='cfgfile',
    help='Name of bathymetry configuration file where to get info (file, names)')
parser.add_option('--varname', action='store', dest='varname',
    help='Name of the bathymetric netcdf variable')
parser.add_option('--lonname', action='store', dest='lonname',
    help='Name of longitude axis')
parser.add_option('--latname', action='store', dest='latname',
    help='Name of latitude axis')
parser.add_option('-t', '--title', action='store', dest='title',
    help='Title of the plot')
parser.add_option('--land', action='store_true', dest='land',
    help='Do not hide land topography [default: %default]', default=False)
parser.add_option('-r', '--reverse', action='store_true', dest='reverse',
    help='Reverse bathymetry [default: %default]', default=False)
parser.add_option('--zmin', '--z0', action='store', dest='zmin',
    type='float', help='Minimal altitude (<0 under the sea, >0 on land)')
parser.add_option('--zmax', '--z1', action='store', dest='zmax',
    type='float', help='Maximal altitude (<0 under the sea, >0 on land)')
parser.add_option('--maxalt', action='store', dest='maxalt',
    type='float', help='Max altitude to plot (positive float)')
parser.add_option('--noshadow', action='store_true', dest='noshadow',
    help="Do not plot shadow for the relief", default=False)
parser.add_option('--clsize', action='store', dest='clsize',
    type='float', help='Size of contour labels [default: %default, 0 to turn off]', default=8)
parser.add_option('--clglow', action='store_true', dest='clglow',
    help='Add glow to contour labels [default: %default]', default=False)
parser.add_option('-o', '--out', action='store', dest='figname',
    help='Name of an output file where to store the plot')
parser.add_option('-s', '--figsize', action='store', dest='figsize',
    help='Size of figure in inches, like "5,6" or simply "5"')
parser.add_option('--proj', action='store', dest='proj', default='merc',
    help="Map projection [default: %default]")
parser.add_option('--res', action='store', dest='res', default='auto',
    help="Shoreline resolution [default: %default]")
parser.add_option('--force', '-f', action='store_true', dest='force',
    help="Don't ask question", default=False)

# Parse
(options, args) = parser.parse_args()

# Target config or netcdf
if len(args)>0:
    ncfile = args[0]
else:
    ncfile = options.ncfile
if ncfile is not None:
    options.force = True

# Bounds
defaults = [('xmin', -720), ('xmax', 720), ('ymin', -90), ('ymax', 90)]
for bname, bdef in defaults:
    if getattr(options, bname) is None:
        desc = parser.get_option('--'+bname).help
        if options.force:
            value = bdef
        else:
            while True:
                try:
                    value = float(raw_input('Please set %s (or use -f or --%s options): '%(desc.lower(), bname)).strip())
                    break
                except ValueError:
                    print 'Please enter a float'
        if getattr(options, bname) is None:
            setattr(options, bname, value)
lon = (options.xmin, options.xmax)
lat = (options.ymin, options.ymax)

# Netcdf file direct access
if ncfile is not None:
    if not ncfile.startswith('http://') and not os.path.exists(ncfile):
        parser.error('Netcdf file not found: '+ncfile)
    options.cfgfile = None
    options.name = ncfile

# Import
try:
    from vacumm.bathy.bathy import NcGriddedBathy
except:
    sys.exit("Can't import vacumm bathymetry module")

# Get and plot bathy
kw = {'maxvalue':None} if options.land else {}
b = NcGriddedBathy(lon=lon, lat=lat, name=options.name,
    cfgfile=options.cfgfile, varname=options.varname, lonname=options.lonname,
    latname=options.latname, reverse=options.reverse, **kw)
kw = {}
if options.land:
    kw['fillcontinents'] = False
    kw['drawcoastlines_color'] = 'r'
    kw['drawcoastlines_linewidth'] = 1
if options.figsize is not None:
    options.figsize = eval(options.figsize)
if options.clsize is not None:
    if options.clsize:
        kw['clabel_size'] = options.clsize
    else:
        kw['clabel'] = options.clsize
b.plot(proj=options.proj, title=options.title, figsize=options.figsize,
    shadow=not options.noshadow,
    zmin=options.zmin, zmax=options.zmax, clabel_glow=options.clglow,
    savefig=options.figname, show=options.figname is None, res=options.res, **kw)
