#!/usr/bin/env python
"""Show attributes of generic variable or axis"""

# Arguments
from argparse import ArgumentParser
import sys, os, shutil
parser = ArgumentParser(description=__doc__)
parser.add_argument('name', help='Name of axis or variable')
parser.add_argument('-s', '--separator', default=',', 
    help='separator between attribute key:value pairs [default: %(default)s]')
parser.add_argument('-v', '--varname', help='variable name')
parser.add_argument('-f', '--filename', help='file name')
args = parser.parse_args()


# Imports
from vacumm.data.cf import *

# Get attributes
atts = cf2atts(args.name)

# Output
# - ncatted output
varname = args.name if args.varname is None else args.varname
cmd = []
for attname, val in atts.items():
    atttype = 'c'
    mode = 'o'
    if atttype=='c':
        val = repr(val)
    cmd.append('-a %(attname)s,%(varname)s,%(mode)s,%(atttype)s,%(val)s'%vars())
cmd = ' '.join(cmd)
if args.filename is not None:
    cmd = 'ncatted -O %s %s'%(cmd, args.filename)
    
print cmd






