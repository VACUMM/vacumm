#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Get or update scripts used for the VACUMMM courses"""


# Arguments
import argparse, os, sys
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("workdir", help="working directory where to put files [default: %(default)s]", nargs="?", 
    default=os.getcwd())
parser.add_argument("-r", "--replace", action='store_true', 
    help="do not backup (.bak) a file if it already exists")
args = parser.parse_args()

# Get the list of scripts
from vacumm.config import get_scripts_dir
try:
    cdir = get_scripts_dir('courses')
except:
    sys.exit("Can't find the directory hosting the courses script files")
from glob import glob
cfiles = glob(os.path.join(cdir, "courses_*.py"))
if not cfiles:
    print "No script files found in directory: "+cdir
    sys.exit()
    
# Working directory
if not os.path.exists(args.workdir):
    os.mkdir(args.workdir)
    print 'Created working directory: '+os.path.realpath(args.workdir)
else:
    print "Working directory: "+os.path.realpath(args.workdir)

# Copy files
import shutil
for cfile in cfiles:
    
    # Local file
    basefile = os.path.basename(cfile)
    cwfile = os.path.join(args.workdir, basefile)
    
    # Backup ?
    msg = 'Installed '+basefile
    if os.path.exists(cwfile):
        if args.replace:
            os.remove(cwfile)
            msg = "Replaced "+basefile
        else:
            bfile = cwfile+'.bak'
            bmsg = "Backuped %s to %s"%(basefile, basefile+'.bak')
            if os.path.exists(bfile):
                os.remove(bfile)
                bmsg += ' (old backup file removed)'
            os.rename(cwfile, bfile)
            print bmsg
           
    # Get it
    shutil.copy(cfile, cwfile)
    print msg
