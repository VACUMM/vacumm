#!/bin/env python
"""Sample executable script using config"""
from __future__ import print_function

from vcmq import data_sample, code_file_name, cfgargparse, map2
from argparse import ArgumentParser
import cdms2

# FORCE COMMANDLINE ARGUMENTS FOR THE EXAMPLE
ARGUMENTS = "--var=temp --zoom-lat-min=48 "+data_sample('mars3d.xy.nc')
import sys
sys.argv = sys.argv[:1] + ARGUMENTS.split()

############################################################################
############################################################################
############################################################################

# Init commandline parser
parser = ArgumentParser(description="Script to plot a 2D netcdf variable")
parser.add_argument('ncfile', help='input netcdf file')

# Configuration
cfg, args = cfgargparse('misc.config.argparse.ini', parser)
print('Current configuration:', cfg)

# Now use it
# - load zoom
lon = (cfg['zoom']['lon']['min'], cfg['zoom']['lon']['max'])
lat = (cfg['zoom']['lat']['min'], cfg['zoom']['lat']['max'])
# - read var
f = cdms2.open(args.ncfile)
var = f(cfg['var'], lon=lon, lat=lat, time=slice(0, 1))
f.close()
# - plot
long_name = var.long_name
map2(var, title=cfg['title']%locals(), show=False, close=False)

