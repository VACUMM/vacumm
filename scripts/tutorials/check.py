#!/usr/bin/env python
"""This script run all tutorial script to check if they work

Help:

    check --help
"""

# Parse arguments
from argparse import ArgumentParser
import sys, os
parser = ArgumentParser(
    description="Check the execution of tutorials.\n\n[pattern] defaults to *.py.",)
parser.add_argument('pattern', default='*.py', help='file pattern  [%(default)s]',
    nargs='*')
parser.add_argument('-e','--exclude', action='append',
    help='add a glob pattern to exclude')
parser.add_argument('-l','--loglevel',
    choices=['debug', 'info', 'error'], default='info',
    help='logging level to console [%(default)s]')
parser.add_argument('-c', '--cfgfile', default='check.cfg', help='configuration file')
args = parser.parse_args()

# Retreive the list of scripts
include = []
from glob import glob
# - include
if not isinstance(args.pattern, list):
    args.pattern = [args.pattern]
for pat in args.pattern:
    include.extend(glob(pat))
include = [os.path.abspath(fn) for fn in include]
# - exclude
if args.exclude is not None:
    exclude = []
    for pat in args.exclude:
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
import logging
logfile = '.'.join([os.path.splitext(__file__)[0], 'log'])
if os.path.exists(logfile):
    os.remove(logfile)
logger = logging.getLogger('CHECK')
file = logging.FileHandler(logfile)
file.setFormatter(
    logging.Formatter('%(asctime)s: %(name)s [%(levelname)-8s] %(message)s',
        '%Y-%m-%d %H:%M'))
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(name)s [%(levelname)-8s] %(message)s'))
file.setLevel('DEBUG')
console.setLevel(args.loglevel.upper())
logger.setLevel('DEBUG')
logger.addHandler(file)
logger.addHandler(console)

# Configuration: arguments to scripts
from ConfigParser import SafeConfigParser
from vacumm.config import data_sample
import re
cfg = SafeConfigParser()
cfg.read(args.cfgfile)
re_data_sample = re.compile(r'(data_sample\([^\)]+\))', re.I)
def convert_data_samples(argline, split=True):
    aargs = re_data_sample.split(argline)
    for i, arg in enumerate(aargs):
        m = re_data_sample.match(arg)
        if m:
            arg = eval(arg)
        aargs[i] = arg
    if split: return aargs
    return ''.join(aargs)
def get_args(scriptname, split=True):
    if scriptname.endswith('.py'):
        scriptname = scriptname[:-3]
    scriptname = os.path.basename(scriptname)
    if cfg.has_option('args', scriptname):
        argline = cfg.get('args', scriptname)
        return convert_data_samples(argline, split=split)


# Run them all
import subprocess
igood = 0
ibad = 0
for script in files:
    if os.path.realpath(script)==os.path.realpath(__file__):
        logger.debug('%s: %s', os.path.basename(script), 'SKIPPED')
        continue
    out, err = subprocess.Popen(["python", script] + (get_args(script) or []),
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    script = os.path.basename(script)
    if len(err)==0 or 'Traceback'  not in err:
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
