#!/bin/env python
"""Sample executable script using config"""

# Init commandline parser
import sys, os
from argparse import ArgumentParser
parser = ArgumentParser(
    description="Script to plot a 2D netcdf variable", add_help=False) # 
parser.add_argument('ncfile', help='input netcdf file')
parser.add_argument('-h','--help', action='store_true', help='show a reduced help')
parser.add_argument('--long-help', action='store_true', help='show an extended help')

# Configuration
from vacumm.misc.config import cfgargparse
cfg, args = cfgargparse('misc.config.argparse.ini', parser)
print 'Current configuration:', cfg

# Now use it
# - load zoom
lon = (cfg['zoom']['lon']['min'], cfg['zoom']['lon']['max'])
lat = (cfg['zoom']['lat']['min'], cfg['zoom']['lat']['max'])
# - read var
import cdms2
f = cdms2.open(args.ncfile)
var = f(cfg['var'], lon=lon, lat=lat, time=slice(0, 1))
f.close()
# - plot
from vacumm.misc.plot import map2
long_name = var.long_name
map2(var, title=cfg['title']%locals(), savefigs=__file__, show=False)

