#!/usr/bin/env python
"""Show time axis dates of netcdf file"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] ncfile", 
    description="Show time step of a netcdf file.")
parser.add_option('-t', '--time', action='store', dest='tname',
    help='Name of the time axis variable')
parser.add_option('-v', '--var', action='store', dest='vname',
    help='Name of the variable from which to guess time')
parser.add_option('-u', '--units', action='store', dest='units', default="a", 
    help='Units as one of "s" for seconds, "m" for minutes, "h" for hours, "d" for days,'
        ' "a" for auto (the defaults), and time axis units if empty ""')
parser.add_option('-s', '--stat', action='store', default='median', dest='stat',
    help='Statistics to guess time step, as one of "min", "max" ,"mean" or "median" (default: %default)')
parser.add_option('-f', '--format', action='store', default='%(value)g %(units)s', dest='format',
    help='Output format (default: %default)')
   
    
# Parse
(options, args) = parser.parse_args()
if len(args)==0:
    parser.error('you must provide a valid netcdf file as first argument')

# Netcdf file
ncfile = args[0]
if not os.path.exists(ncfile):
    sys.exit('File not found: '+ncfile)
try:
    import cdms2, cdtime
    f = cdms2.open(ncfile)
except:
    sys.exit("Can't open "+ncfile)
    
# Guess time
taxis = None
# - from specified axis
if options.tname is not None:
    if options.tname in f.listdimension():
        taxis = f.getAxis(options.tname)
    else:
        sys.exit('Time axis "%s" does not exist'%options.tname)
# -from specified variable
if taxis is None and options.vname is not None:
    if options.vname in f.listvariables():
        var = f[options.vname]
        taxis = var.getTime()
        if taxis is None:
            sys.exit('Variable "%s" found but has no time axis'%options.vname)
    else:
        sys.exit('Variable "%s" not found'%options.vname)
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
    if options.units and options.units[0] != uu[0]: continue
    if uu != units:
        taxis.toRelativeTime(uu+' since 2000')
        units = uu
    break
tt = taxis.getValue()

# Time step
if not options.stat in ['min', 'max', 'mean', 'median']:
    options.stat = "median"
import numpy 
value = getattr(numpy, options.stat)(numpy.abs(numpy.diff(tt)))

# Auto units
if options.units[0].startswith('a') and units not in ['months', 'years'] and \
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
    print options.format%locals()
except:
    parser.error("Can't format result with this pattern: %s"%options.format)
