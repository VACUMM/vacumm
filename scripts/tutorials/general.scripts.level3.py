#!/usr/bin/env python
"""Plot mixed layer depth from netcdf file"""
import argparse, os
from vcmq import DS, map2, MV2, data_sample

# Arguments
parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("ncfile", help="model file", nargs="?",
                    default="menor.nc")
parser.add_argument(
        "-m", "--model",
        help="model type [default: %(default)s]",
        choices=['ocean', 'mars', 'nemo'], default='mars')
parser.add_argument('-t', '--time',
                    help='time selection, like \'("2000","2001","co")\'')
parser.add_argument('-x', '--lon', default="(3.,4.5)",
                    help='longitude selection, like \'(40,42)\'')
parser.add_argument('-y', '--lat', default="(42,42.8)",
                    help='latitude selection, like \'(40,42)\'')
parser.add_argument('--title', help='alternative title')
parser.add_argument(
        '--fillmode', help='fill mode [default: %(default)s]',
        choices=['nofill', 'pcolormesh', 'contourf'], default='contourf')
parser.add_argument('--disp', help='display the figure', action='store_true')
parser.add_argument('-o', '--outputfig', help='file name for output figure')
args = parser.parse_args()

# Check arguments
select = {}
for selname in 'lon', 'lat', 'time':
    val = getattr(args, selname)
    if val is not None:
        try:
            val = eval(val)
        except Exception:
            parser.error('--{}: wrong argument'.format(selname))
    select[selname] = val

# The file
ncfile = args.ncfile
if not os.path.exists(ncfile):
    ncfile = data_sample(ncfile)
if not os.path.exists(ncfile):
    parser.error('File not found: '+ncfile)

# Read temperature and depth
ds = DS(ncfile, args.model, **select)
mld = ds.get_mld(mode='deltatemp', deltatemp=0.1, squeeze=1)

# Time average (if needed)
if mld.ndim == 3:
    mld = MV2.average(mld, axis=0)

# Plot
map2(mld, proj='merc', figsize=(6, 4), autoresize=0, fill='contourf',
     title=args.title, colorbar_shrink=0.7, right=1,
     show=args.disp, savefig=args.outputfig, savefig_verbose=True,
     close=False, cmap='cmocean_dense')
