#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Print your configuration of VACUMM"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(#usage="Usage: %prog [options]", 
    description="Print your configuration of VACUMM")
parser.add_option('--system', '-s', action='store_true', dest='s',
    help="Print only system information", default=False)
parser.add_option('--direc', '-d', action='store_true', dest='d',
    help="Print only important VACUMM directories", default=False)
parser.add_option('--config', '-c', action='store_true', dest='c',
    help="Print only VACUMM config itself", default=False)
parser.add_option('--packages', '-p', action='store_true', dest='p',
    help="Print only version of important packages", default=False)
parser.add_option('--no-header', action='store_true', dest='noheader',
    help="Do not add headers", default=False)
    
# Parse
(options, args) = parser.parse_args()

# Import
try:
    from vacumm.config import print_config
except ImportError, e:
    parser.error("Error importing vacumm:\n"+e.message)
    
# Options
if not (options.s or options.d or options.c or options.p):
    options.s = options.d = options.c = options.p = True

# Print
print_config(system=options.s, direc=options.d, config=options.c, packages=options.p, 
    headers=not options.noheader)
