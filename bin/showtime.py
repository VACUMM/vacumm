#!/usr/bin/env python
"""Show time axis dates of netcdf file"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] ncfile",
    description="Show time axis dates of a netcdf file.")
parser.add_option('-t', '--time', action='store', dest='tname',
    help='name of the time axis or variable')
parser.add_option('-v', '--var', action='store', dest='vname',
    help='name of the variable from which to guess time')
parser.add_option('-s', '--slice', action='store', dest='tslice',
    help='time slice (examples: "2", ":-2:4")')
parser.add_option('-m', '--minmax', action='store_true', default=False, dest='minmax',
    help='show min/max only (default: %default)')
parser.add_option('-n', '--ncol', action='store', dest='ncol', default=5, type='int',
    help='number of comlumns for output (default: %default)')
parser.add_option('-f', '--format', action='store', default='%Y-%m-%d %H:%M:%S', dest='format',
    help='date format (default: %default)')


# Parse
(options, args) = parser.parse_args()
if len(args)==0:
    parser.error('you must provide a valid netcdf file as first argument')

# Netcdf file
ncfile = args[0]
if not ncfile.startswith('http://') and not os.path.exists(ncfile):
    sys.exit('File not found: '+ncfile)
try:
    import cdms2
    f = cdms2.open(ncfile)
except:
    sys.exit("Can't open "+ncfile)

# Guess time
taxis = None
# - from specified axis
if options.tname is not None:
    if options.tname in f.listdimension():
        taxis = f.getAxis(options.tname)
    elif options.tname in f.listvariables():
        v = f.getVariable(options.tname)
        taxis = cdms2.createAxis(v[:], id='time')
        taxis.units = v.units
        print 'Created time axis from variable "%s" with units "%s"'%(options.tname, taxis.units)
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
ctimes = taxis.asComponentTime()

# Slice
if options.tslice is not None:
    try:
        ss = []
        for s in options.tslice.split(':')[:3]:
            ss.append(None if s=='' else int(s))
        if len(ss)>1:
            tslice = slice(*ss)
        else:
            tslice = ss[0]
    except:
        parser.error('Incorrect slice option: '+options.tslice)
    try:
        ctimes = ctimes[tslice]
        if not isinstance(ctimes, list): ctimes = [ctimes]
    except:
        sys.exit('Error when slicing dates (slice arg="" and nt=%i)'%(ss, len(ctimes)))
if len(ctimes)==0:
    sys.exit('No valid time to show')
# Display
from datetime import datetime
def ct2str(ct):
    ssec = int(round(ct.second))%60
    msec = int(round(ct.second-ssec))*1000
    try:
        return datetime(ct.year, ct.month, ct.day, ct.hour, ct.minute, ssec, msec).strftime(options.format)
    except Exception, e:
        sys.exit('Error while formatting time "%s" with format "%s": %s'%(str(ct), options.format, e))
if options.minmax: # min/max only
    if len(ctimes)==1:
        print 'min/max:', ct2str(ctimes[0])
    else:
        print 'min:', ct2str(min(ctimes))
        print 'max:', ct2str(max(ctimes))

else: # all
    s = 's' if len(ctimes) else ''
    print '%i value%s:'%(len(ctimes), s)
    nt = len(ctimes)
    for l in xrange(nt/options.ncol+1):
        i0 = l*options.ncol
        i1 = i0+options.ncol
        print ', '.join([ct2str(ct) for ct in ctimes[i0:i1]])
