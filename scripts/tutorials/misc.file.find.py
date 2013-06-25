#!/usr/bin/env python

from os.path import abspath, dirname, join, realpath

import vacumm.misc.file as F

# This file
# realpath() is not really required in this example, this is to advice that in case
# you were executing a link to this file, the code would still work (it can be a good
# practice to do so in your scripts if you want that kind of support)
thisfile = realpath(__file__)

# The directory containing this file
thisdir = dirname(__file__)

# The vacumm's doc/sphinx/source directory
# abspath() is also not really required in this example, but this would make paths
# look nicer (without '..' in them)
# docsrcdir = abspath(join(dirname(__file__), '..', '..'))
docsrcdir = join(dirname(__file__), '..', '..')

fkw = dict(abspath=False)

# Using wildcard expressions
# ==========================

# A simple example
print 'Searching any file starting with "misc.file." in this directory:\n'
print ' ','\n  '.join(F.find('misc.file.*', thisdir, **fkw))
print


# Controlling the search depth
print 'Searching any file starting with "misc.file." in the doc source directory:\n'
print ' ','\n  '.join(F.find('misc.file.*', docsrcdir, depth=None, **fkw))
print
# we've set the depth to None to supress the limit in recursion (which defaults
# to 0 meaning 0 levels from the search directory)


# Multiple patterns
print 'Searching any file starting with "misc.file" and ending with ".py" or ".rst" in the doc source directory:\n'
print ' ','\n  '.join(F.find(('misc.file*.py', 'misc.file*.rst',), docsrcdir, depth=None, **fkw))
print


# How to exclude
print 'Searching any file starting with "misc.file." in the doc source directory, excluding "library" and ".svn" directories and "*.rst" files:\n'
print ' ','\n  '.join(F.find('*/misc.file.*', docsrcdir, depth=None, exclude=('*/library/*', '*/.svn/*', '*.rst'), matchall=True, **fkw))
print
# by default only file names are evaluated
# we've set matchall to True because .svn is part of the whole paths to be evaluated
# also note that we had to change the pattern no match whole paths


# Search directories only
# By default only files are looked up
print 'Searching any directory named python in the doc source directory:\n'
print ' ','\n  '.join(F.find('python', docsrcdir, depth=None, files=False, dirs=True, **fkw))
print


# Using callbacks
# ===============

# Controlling the search depth
print 'Showing search results with callbacks:\n'
def ondir(e):
    # we want to limit the outputs of this tutorial...
    if '/.svn' not in e:
        print '  + evaluating directory:', e
def onfile(e):
    # we want to limit the outputs of this tutorial...
    if '/.svn' not in e and 'misc.file.' in e:
        print '  - evaluating file:', e
def onmatch(e):
    print '  * matching entry:', e
F.find('*/misc.file.*', docsrcdir, depth=None, exclude=('*/library/*', '*/.svn/*', '*.rst'), matchall=True, ondir=ondir, onfile=onfile, onmatch=onmatch, **fkw)
print
# note that ondir and onfile callbacks are called even if the path does not matches


# Using regular expressions
# =========================

# A simple example
print 'Searching any file starting with "misc.file." in this directory:\n'
print ' ','\n  '.join(F.efind('misc\.file\..*', thisdir, **fkw))
print

# An advanced example which uses regex match results
print 'Searching any file starting with "misc.file." in this directory and show matching groups:\n'
for filepath, matchobj in F.xefind('misc\.file\.(.*)\.(?:py|rst)', thisdir, getmatch=True, **fkw):
    print '  file path: %s, match groups: %s'%(filepath, matchobj.groups())
print
# grouping is done with parenthesis (), we could also have used named group (?P<groupname1>groupexp) and matchobj.groupdict()
# note that the second group is a non capturing one, we needed this to include both '.py' and '.rst' files
# also note that we use xefind instead of efind, this way we could avoid buildnig a huge list of (filepath,matchobj) couples
# this is not the case here, just an example


