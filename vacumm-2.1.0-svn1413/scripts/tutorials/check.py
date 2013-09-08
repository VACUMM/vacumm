#!/usr/bin/env python
"""This script run all tutorial script to check if they work

Help: 

    check --help
"""

# Parse arguments
from optparse import OptionParser, sys, os
parser = OptionParser(usage="Usage: %prog [options] [pattern]",
    description="Check the execution of tutorials.\n\n[pattern] defaults to *.py.",)
parser.add_option('-e','--exclude', action='append', dest="exclude", 
    help='add a glob pattern to exclude')
parser.add_option('-l','--loglevel', action='store', dest="loglevel", 
    choices=['debug', 'info', 'error'], default='info', 
    help='logging level to console [%default]')
options, args = parser.parse_args()
if len(args)==0:
    args = ['*.py']

# Retreive the list of scripts
include = []
from glob import glob
# - include
for pat in args:
    include.extend(glob(pat))
include = [os.path.abspath(fn) for fn in include]
# - exclude
if options.exclude is not None:
    exclude = []
    for pat in options.exclude:
        exclude.extend(glob(pat))
    exclude = [os.path.abspath(fn) for fn in exclude] + [os.path.abspath(__file__)]
    files = []
    for fn in include:
        if fn not in exclude:
            files.append(fn)
else:
    files = include
files.sort()

# Setup the logger
from vacumm.misc.io import Logger
logger = Logger('CHECK', 
    logfile='.'.join([os.path.splitext(__file__)[0], 'log']), 
    cfmt='[%(name)s] %(message)s', 
    ffmt='%(asctime)s: [%(name)s] %(message)s', 
    full_line=True,
)
logger.set_loglevel(console=options.loglevel, file='debug')
        
# Run them all
import subprocess
igood = 0
ibad = 0
from vacumm.config import data_sample
for script in files:
    if os.path.realpath(script)==os.path.realpath(__file__): 
        logger.debug('%s: %s', os.path.basename(script), 'SKIPPED')
        continue
    out, err = subprocess.Popen(["python", script,  '--var=temp', '--zoom-lat-min=48',  data_sample('mars3d.xy.nc')], 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    script = os.path.basename(script)
    if len(err)==0 or ('warning' in err.lower() and 'raise' not in err):
        logger.info('%s: %s', script, 'OK')
        igood +=1
    else:
        logger.error('%s: %s', script, 'FAILED')
        ibad +=1
    if len(out): logger.debug(' OUT: '+out)
    if len(err): logger.debug(' ERR: '+err)
logger.set_loglevel(console='debug')
msg = 'Results (%i checks): '%(igood+ibad)
gb = []
if igood:
    gb.append('%i OK'%igood)
if ibad:
    gb.append('%i FAILED'%ibad)
msg = 'Results (%i checks): %s'%(igood+ibad, ', '.join(gb))
logger.debug(msg)
