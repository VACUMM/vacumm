#!/usr/bin/env python
"""Plot mixed layer depth from netcdf file"""

# Arguments
import argparse, os
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("ncfile", help="model file", nargs="?",
    default="menor.nc")
parser.add_argument("-m", "--model", help="model type [default: %(default)s]",
    choices = ['ocean', 'mars', 'nemo'], default='mars')
parser.add_argument('-t', '--time', help='time selection, like \'("2000","2001","co")\'')
parser.add_argument('-x', '--lon', help='longitude selection, like \'(40,42)\'')
parser.add_argument('-y', '--lat', help='latitude selection, like \'(40,42)\'')
parser.add_argument('--title', help='alternative title')
parser.add_argument('--fillmode', help='fill mode [default: %(default)s]',
    choices=['nofill', 'pcolormesh', 'contourf'], default='contourf')
parser.add_argument('--disp', help='display the figure', action='store_true')
parser.add_argument('-o', '--outputfig', help='file name for output figure')
options = parser.parse_args()

# Check arguments
select = {}
for selname in 'lon', 'lat', 'time':
    val = getattr(options, selname)
    if val is not None:
        try:
            lon = eval(options.lon)
        except:
            parser.error('--%s: wrong argument'%selname)
    select[selname] = val

# Imports
from vcmq import DS, map2, MV2, data_sample

# The file
ncfile = options.ncfile
if not os.path.exists(ncfile):
    ncfile = data_sample(ncfile)
if not os.path.exists(ncfile):
    parser.error('File not found: '+ncfile)

# Read temperature and depth
ds = DS(ncfile, options.model, **select)
mld = ds.get_mld(mode='deltatemp', squeeze=1)

# Time average (if needed)
if mld.ndim==3:
    mld = MV2.average(mld, axis=0)

# Plot
map2(mld, proj='merc', figsize=(6, 6), autoresize=0,
    title=options.title, colorbar_shrink=0.7, right=1,
    show=options.disp, savefig=options.outputfig, savefig_verbose=True,
    close=True)


