#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os, sys

parser = ArgumentParser(description='Mask variables of a NetCDF file')
parser.add_argument('inpfile', help='Netcdf input file')
parser.add_argument('outfile', help='Netcdf output file')
parser.add_argument('-v', '--var', dest='variables', help='Variable to process (by default, all gridded variables)')
args = parser.parse_args()

print parser.inpfile
print parser.outfile
print parser.variables