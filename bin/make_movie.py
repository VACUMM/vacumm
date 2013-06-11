#!/usr/bin/env python
"""Helper to make a movie using :func:`vacumm.misc.plot.make_movie`"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(usage="Usage: %prog [options] <pattern or files>", 
    description="Make a .gif or .mpg animation from image files.\n"
    "You need programs 'convert' to create a .gif and 'ffmpeg' to create a .mpg.")
parser.add_option('-o', '--output', action='store', dest='output',
    help='output file name [default: anim.gif or anim.mpg]')
parser.add_option('-d', '--delay', action='store', dest='delay', default=1., 
    help='delay in seconds between two images [default: %default]', type='float')
parser.add_option('-w', '--windows', action='store_true', dest='windows',
    help='use the windows codec [default: %default]',  default=True)
parser.add_option('-c', '--clean', action='store_true', dest='clean',
    help='remove image files once the animation is created [default: %default]',  
    default=False)
parser.add_option('-v', '--verbose', action='store_true', dest='verbose',
    help='be verbose with system calls [default: %default]',  
    default=False)
   
# Parse
(options, args) = parser.parse_args()
if len(args)==0:
    figpats = ['*.png']
#    parser.error('you must provide a valid file pattern as first argument')
else:
    figpats = args

# Check files
import glob
files = [glob.glob(figpat) for figpat in figpats]
if not files:
    sys.exit('No file found with this pattern: %s'%figpat)
elif options.verbose:
    print 'Found %i files'%len(files)
if options.output is None:
    if options.windows:
        options.output = 'anim.mpg'
    else:
        options.output = 'anim.gif'
    
# Call to make_movie
from vacumm.misc.plot import make_movie
make_movie(figpat, options.output, delay=options.delay, windows=options.windows, 
    clean=options.clean, verbose=options.verbose)
if options.verbose: print 'Created animation file: '+options.output

    
