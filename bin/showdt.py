#!/usr/bin/env python
"""Show time axis dates of netcdf file"""

# Arguments
from argparse import ArgumentParser
import sys, os, shutil
parser = ArgumentParser(description="Show time step of a netcdf file.")
parser.add_argument('ncfile', help='Netcdf file name')
parser.add_argument('-t', '--time', help='Name of the time axis variable')
parser.add_argument('-v', '--var', dest='vname',
    help='Name of the variable from which to guess time')
parser.add_argument('-u', '--units', default="a", 
    help='Units as one of "s" for seconds, "m" for minutes, "h" for hours, "d" for days,'
        ' "a" for auto (the defaults), and time axis units if empty ""')
parser.add_argument('-s', '--stat', default='median', 
    help='Statistics to guess time step, as one of "min", "max" ,"mean" or "median" (default: %(default)s)')
parser.add_argument('-f', '--format', default='%(value)g %(units)s', 
    help='Output format (default: %(default)s)')
   
    
# Parse
args = parser.parse_args()

# Netcdf file
if '://' not in args.ncfile and not os.path.exists(args.ncfile):
    sys.exit('File not found: '+args.ncfile)
try:
    import cdms2, cdtime
    f = cdms2.open(args.ncfile)
except:
    sys.exit("Can't open "+args.ncfile)
    
# Guess time
taxis = None
# - from specified axis
if args.time is not None:
    if args.time in f.listdimension():
        taxis = f.getAxis(args.time)
    else:
        sys.exit('Time axis "%s" does not exist'%args.time)
# -from specified variable
if taxis is None and args.vname is not None:
    if args.vname in f.listvariables():
        var = f[args.vname]
        taxis = var.getTime()
        if taxis is None:
            sys.exit('Variable "%s" found but has no time axis'%args.vname)
    else:
        sys.exit('Variable "%s" not found'%args.vname)
# - guess through variables
if taxis is None:
    for vname in f.listvariables():
        taxis = f[vname].getTime()
        if taxis is not None: break
    else:
        sys.exit('No valid time axis found')
taxis = taxis.clone()
f.close()
taxis._data_ = taxis._data_.astype('d')

# Check length
if len(taxis)<2:
    sys.exit('Length of your time axis must be at least 2')
    
# Check units
try:
    ctimes = taxis.asComponentTime()
except:
    sys.exit('Time axis seems to have bad units: %s'%getattr(taxis, 'units', ''))
units = taxis.units.split()[0].lower()

# Numeric values
for uu in ['seconds', 'minutes', 'hours', 'days']:
    if args.units and args.units[0] != uu[0]: continue
    if uu != units:
        taxis.toRelativeTime(uu+' since 2000')
        units = uu
    break
tt = taxis.getValue()

# Time step
if not args.stat in ['min', 'max', 'mean', 'median']:
    args.stat = "median"
import numpy 
value = getattr(numpy, args.stat)(numpy.abs(numpy.diff(tt)))

# Auto units
if args.units[0].startswith('a') and units not in ['months', 'years'] and \
    (value<1 or value>100):
    testlist =  ['seconds', 'minutes', 'hours', 'days'][::1-2*int(value>100)]
    rt0 = taxis.asRelativeTime()[0]
    cdmsunits = getattr(cdtime, units.title())
    if testlist[0]!=units:
        for aunits in testlist:
            faunits = aunits+' since 2000'
            valuenew = rt0.add(value, cdmsunits).torel(faunits).value-rt0.torel(faunits).value
            if valuenew>1 and valuenew<=100:
                value = valuenew
                units = aunits
                break

# Display
if value<=1 and units.endswith('s'): units = units[:-1]
try:
    print args.format%locals()
except:
    parser.error("Can't format result with this pattern: %s"%args.format)
