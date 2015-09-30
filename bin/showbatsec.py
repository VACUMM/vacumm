#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Plot a section of bathymetry"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] [ncfile]",
    description="Plot a section of bathymetry.")
parser.add_option('--x0', '--lon0', action='store', dest='x0', type='float',
    help="Longitude of first point")
parser.add_option('--y0', '--lat0', action='store', dest='y0', type='float',
    help="Latitude of first point")
parser.add_option('--x1', '--lon1', action='store', dest='x1', type='float',
    help="Longitude of second point")
parser.add_option('--y1', '--lat1', action='store', dest='y1', type='float',
    help="Latitude of second point")
parser.add_option('--ncfile', action='store', dest='ncfile',
    help='Name of netcdf file for a direct access [obsolete, use first argument instead]')
parser.add_option('-n', '--name', action='store', dest='name',
    help='Name of bathymetry as stored in configuration file)')
parser.add_option('-c', '--cfgfile', action='store', dest='cfgfile',
    help='Name of bathymetry configuration file where to get info (file, names)')
parser.add_option('--varname', action='store', dest='varname',
    help='Name of the bathymetric netcdf variable')
parser.add_option('--lonname', action='store', dest='lonname',
    help='Name of longitude axis')
parser.add_option('--latname', action='store', dest='latname',
    help='Name of latitude axis')
parser.add_option('--along', action='store', dest='along', default='auto',
    help='Along wich coordinates to plot: "lon", "lat", "m" ("m" or "km" or "dist"), "auto" [default: %default]')
parser.add_option('-t', '--title', action='store', dest='title',
    help='Title of the plot')
parser.add_option('-r', '--reverse', action='store_true', dest='reverse',
    help='Reverse bathymetry [default: %default]', default=False)
parser.add_option('--zmin', action='store', dest='zmin',
    type='float', help='Min altitude to plot (negative under water)')
parser.add_option('--zmax', action='store', dest='zmax',
    type='float', help='Max altitude to plot (negative under water)')
parser.add_option('--nomap', action='store_false', dest='map',
    help='So add a small map of the transect situation', default=True)
parser.add_option('--mapscale', action='store', dest='mapscale', default=2,
    type='float', help='Rescale the map area [default: %default')
parser.add_option('-o', '--out', action='store', dest='figname',
    help='Name of an output file where to store the plot')
parser.add_option('-s', '--figsize', action='store', dest='figsize',  default="7,4",
    help='Size of figure in inches, like "5,6" or simply "5"')

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
for bname in 'x0', 'y0', 'x1', 'y1':
    if getattr(options, bname) is None:
        desc = parser.get_option('--'+bname).help
        while True:
            try:
                value = float(raw_input('Please set %s: '%desc.lower().strip()))
                break
            except ValueError:
                print 'Please enter a float'
        if getattr(options, bname) is None:
            setattr(options, bname, value)
xmin = min(options.x0, options.x1)
xmax = max(options.x0, options.x1)
ymin = min(options.y0, options.y1)
ymax = max(options.y0, options.y1)
lon = (xmin, xmax)
lat = (ymin, ymax)

# Netcdf file direct access
if options.ncfile is not None:
    if not options.ncfile.startswith('http://') and not os.path.exists(options.ncfile):
        parser.error('Netcdf file not found: '+options.ncfile)
    options.cfgfile = None
    options.name = options.ncfile

# Import
try:
    from vacumm.bathy.bathy import NcGriddedBathy
    from vacumm.misc.grid import transect_specs
    from vacumm.misc.grid.regridding import grid2xy
    from vacumm.misc.axes import create_lon, create_lat
    import cdms2, numpy as N
    from vacumm.misc.plot import curve2, minimap
    from vacumm.misc.color import cmap_custom
except:
    sys.exit("Can't import needed vacumm modules")

# Load bathy
b = NcGriddedBathy(lon=lon, lat=lat, name=options.name,
    cfgfile=options.cfgfile, varname=options.varname, lonname=options.lonname,
    latname=options.latname, imargin=2, jmargin=2, maxvalue=None,
    reverse=options.reverse)

# Interpolate on section
bathy = b.bathy()
lons, lats, xx, yy = transect_specs(bathy.getGrid(),
    options.x0, options.y0, options.x1, options.y1, getxy=True)
tbathy = grid2xy(bathy, lons, lats)
if not isinstance(options.along, basestring) or options.along not in ['lon', 'lat', 'm', 'km', 'dist', 'auto']:
    options.along = 'auto'
if options.along=='auto':
    options.along = 'lon' if xx.max()-xx.min() > yy.max()-yy.min() else 'lat'
if options.along=='lon':
    taxis = create_lon(lons)
elif options.along=='lat':
    taxis = create_lat(lats)
else:
    xy = N.sqrt((xx-xx[0])**2+(yy-yy[0])**2)
    if xy.max()-xy.min()>2500 or options.along=='km':
        units = 'km'
        xy *= 1e-3
    else:
        units = 'm'
    taxis = cdms2.createAxis(xy)
    taxis.units = units
    taxis.long_name = 'Distance'
tbathy.setAxis(0, taxis)

# Plot
c = curve2(tbathy, 'k', bgcolor=(.95, .95, 1), order='d-',
    show=False, yfmt='%gm', yunits=False, bottom=0.12, zorder=10,
    figsize=eval(options.figsize), title=options.title, ymin=options.zmin, ymax=options.zmax)
xylims = c.axes.axis()
xx = c.get_xdata(scalar=0)
yy = c.get_ydata(scalar=0)
if xylims[2]<0.: # ocean
#    c.axes.axhspan(xylims[2], 0, color=(.8, .8, 1), zorder=2, linewidth=0)
    c.axes.imshow(N.resize(N.arange(256), (2,256)).T,
        cmap=cmap_custom([((.3,0,.5),0),  ('#000055', .2), ('#007da9', .7), ('#7cbed3',1)]),#'ocean',
        extent=[xylims[0], xylims[1], 0., tbathy.min()], aspect='auto')
c.axes.fill_between(xx, yy, xylims[2], color='.5', zorder=9, linewidth=0) # earth
c.axes.axhline(zorder=8, color='b', linewidth=.8) # sea surface
if options.map:
    m = minimap((lons, lats), zoom=1./options.mapscale)
    xm, ym = m(lons, lats)
    m.axes.plot(xm, ym, color='r')
if options.figname:
    c.savefig(options.figname)
else:
    c.show()

