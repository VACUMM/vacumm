#!/usr/bin/env python

# Imports
import argparse
from glob import glob
import os

# Commandline arguments
parser = argparse.ArgumentParser(description="Generate the rst files for the test scripts")
parser.add_argument("pattern", help="file patern", nargs="*", default="test_*.py")
parser.add_argument("-d", "--testdir", 
    help="Directory of test scripts [default: %(default)s]",
    default="../../../../scripts/test")
#parser.add_argument("-c", "--noclean", action='store_true', 
#    help='Do not remove unneeded rst files')
parser.add_argument("-n", "--new", action='store_true',
    help='Only add new files (do not modify existing files)')
parser.add_argument("-v", "--verbose", action='store_true',
    help='Be more verbose')
parser.add_argument("-q", "--quiet", action='store_true',
    help='Be quiet!')
options = parser.parse_args()

if options.quiet: options.verbose = False

def check_and_write(fname, content):
    """Check that content has changed and write it if needed"""
    # Read
    if os.path.exists(fname):
        status = 'Updated'
        f = open(fname)
        oldcontent = f.read()
        f.close()
        if oldcontent==content:
            if options.verbose: print 'No update needed for '+fname
            return 0
    else:
         status = 'Created'
            
    # Write
    f = open(fname, 'w')
    f.write(content)
    f.close()
    if not options.quiet: 
        print status+' '+fname
    return 1 if status == 'Updated' else 2

# Loop on test files
script_basenames = []
scfiles = glob(os.path.join(options.testdir, options.pattern))
scfiles.sort()
rstfiles = []
skipped = []
svnadd = []
svnrm = []
svnci = []
for i, script_path in enumerate(scfiles):
    
    # Basenames
    script_name = os.path.basename(script_path)
    script_basename, ext = os.path.splitext(script_name)
    script_basenames.append(script_basename)
    
    # Output file
    rstfile = script_basename+'.rst'
    if options.new and os.path.exists(rstfile):
        skipped.append(rstfile)
        continue
    
    # Title
    fi = open(script_path)
    title = fi.readline()
    if title.startswith('#'):
        title = fi.readline()
    fi.close()
    if len(title)>7:
        title = title[3:-4]
    title = ':program:`%(script_name)s` -- %(title)s'%locals()
    content = title+'\n'
    content += '='*len(title)+'\n\n\n'
    
    # Figure?
    figpaths = []
    for figpat in '', '_[0-9]', '_[0-9][0-9]':
        figbases = [script_basename+figpat]
        figbases.append(figbases[0].replace('.', '-'))
        for figbase in figbases:
            figpath = figbase+'.png'
            figpaths.extend(glob(os.path.join(options.testdir, figpath)))
    figpaths = list(set(figpaths))
    figpaths.sort()
    for figpath in figpaths:
        content += '.. figure:: %(figpath)s\n\n'%locals()
        
    # Code content
    content += ".. literalinclude:: %(script_path)s\n\n"%locals()
    
    
    # Check and write
    status = check_and_write(rstfile, content)
    rstfiles.append(rstfile)
    if status:
        svnci.append(rstfile)
        if status==2:
            svnadd.append(rstfile)

# Check unused rst files
for rstfile in glob('test_*.rst'):
    if rstfile not in rstfiles+skipped:
        if options.verbose: print 'Unused: '+rstfile
        svnrm.append(rstfile)
svnci.extend(svnrm)

# Subversion commands
if not options.quiet and svnci:
    print 'Subversion commands:\n'
    if svnadd: print "svn add "+' '.join(svnadd)
    if svnrm: print "svn rm "+' '.join(svnrm)
    print "svn ci "+' '.join(svnci)
    

