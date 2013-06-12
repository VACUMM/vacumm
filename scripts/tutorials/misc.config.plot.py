#!/bin/env python
"""Sample executable script using config"""

# Init commandline parser
from optparse import OptionParser, sys, os
parser = OptionParser(usage="Usage: %prog [options] ncfile",
    description="Script to plot a 2D netcdf variable", add_help_option=False )
parser.add_option('-h','--help', action='store_true', 
    dest="help", help='show a reduced help')
parser.add_option('--long-help', action='store_true', 
    dest="long_help", help='show an extended help')

# Define config manager 
from vacumm.misc.config import ConfigManager, print_short_help
cfgm = ConfigManager('config.ini')

# Get config from commandline options
cfgpatch = cfgm.opt_parse(parser)

# Help
opts = parser.values
args = parser.largs
if opts.long_help:
    parser.print_help()
    sys.exit()
elif opts.help:
    print_short_help(parser)
    sys.exit()

# Check mandatory arguments
import os
if len(args) < 1 or not os.path.exists(args[0]):
    parser.error('You must provide a valid netcdf file name as first argument')
ncfile = args[0]

# Complete config
# - load personal file and default values
cfg = cfgm.load(opts.cfgfile)
# - patch with commandeline options
cfgm.patch(cfg, cfgpatch)

# Check config
print 'Current configuration:', cfg

# Now use it
# - load zoom
lon = (cfg['zoom']['lon']['min'], cfg['zoom']['lon']['max'])
lat = (cfg['zoom']['lat']['min'], cfg['zoom']['lat']['max'])
# - read var
import cdms2
f = cdms2.open(ncfile)
var = f(cfg['var'], lon=lon, lat=lat, time=slice(0, 1))
f.close()
# - plot
from vacumm.misc.plot import map2
long_name = var.long_name
map2(var, title=cfg['title']%locals(), savefigs=__file__, show=False)

