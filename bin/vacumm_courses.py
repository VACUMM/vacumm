#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Get/update scripts used for the VACUMMM courses or show the web page"""

targets = ["courses/courses_*.py", "courses/courses_*.txt", "courses/myfile.f90", 
    "tutorials/general.scripts.*.py", "tutorials/misc.config.argparse.*", 
    "courses/courses_advanced_cfgm.ini", "courses/courses_advanced_cfgm.cfg"]
url = "/courses/index.html"

# Init function
def init(args):
    
    # Get the list of scripts
    from vacumm.config import get_scripts_dir
    try:
        sdir = get_scripts_dir()
    except:
        sys.exit("Can't find the directory hosting the script files")
    print 'Source directory for scripts: '+sdir
    from glob import glob
    cfiles = []
    for pattern in targets:
        cfiles.extend(glob(os.path.join(sdir, pattern)))
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
    import shutil, stat
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
        os.chmod(cwfile, os.stat(cwfile).st_mode|stat.S_IXUSR)
    

# Show function
def show(args):
    
    from vacumm import help
    help(url, recent=not args.ifremer)

# Arguments
import argparse, os, sys
parser = argparse.ArgumentParser(description=__doc__)
subparsers = parser.add_subparsers()
# - init
parser_init = subparsers.add_parser('init', help='initialize the courses by copying all needed files')
parser_init.add_argument("workdir", help="working directory where to put files [default: %(default)s]", nargs="?", 
    default=os.getcwd())
parser_init.add_argument("-r", "--replace", action='store_true', 
    help="do not backup (.bak) a file if it already exists")
parser_init.set_defaults(func=init)
# - show 
parser_show = subparsers.add_parser('show', help='open the VACUMM website to courses page')
parser_show.add_argument("-i", "--ifremer", action='store_true', 
    help="use the IFREMER VACUMM website")
parser_show.set_defaults(func=show)
# - parse
args = parser.parse_args()
args.func(args)

