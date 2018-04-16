#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import vacumm.misc.file as F


print 'Changing the working directory to vacumm\'s root directory'
print
PWD = os.getcwd()
os.chdir(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
searchdir = '.'
fkw = dict(abspath=False)

# Using wildcard expressions
# ==========================

# A simple example
print 'Searching any file starting with "setup." (depth=0):\n'
print ' ','\n  '.join(F.find('setup*', searchdir, depth=0, **fkw))
print


# Controlling the search depth
print 'Searching any file starting with "misc.file." (depth=None):\n'
print ' ','\n  '.join(F.find('misc.file.*', searchdir, depth=None,
    exclude=('*/build/*', ), **fkw))
print
# we've set the depth to None to supress the limit in recursion (which defaults
# to 0 meaning 0 levels from the search directory)


# Multiple patterns
print 'Searching any file starting with "misc.config*" and ending with ".py", ".ini" or ".cfg":\n'
print ' ','\n  '.join(F.find(('misc.config*.py', 'misc.config*.ini', 'misc.config*.cfg',), searchdir, depth=None, **fkw))
print


# How to exclude
print 'Searching any file starting with "misc.file.", excluding "library" and ".svn" directories and "*.rst" files:\n'
print ' ','\n  '.join(F.find('*/misc.file.*', searchdir, depth=None,
    exclude=('*/library/*', '*/.svn/*', '*.rst'), matchall=True, **fkw))
print
# by default only file names are evaluated
# we've set matchall to True because .svn is part of the whole paths to be evaluated
# also note that we had to change the pattern no match whole paths


# Search directories only
# By default only files are looked up
print 'Searching any directory named misc:\n'
print ' ','\n  '.join(F.find('misc', searchdir, depth=None, files=False, dirs=True, **fkw))
print


# Using callbacks
# ===============

# Controlling the search depth
print 'Showing search results with callbacks:\n'
def ondir(e):
    # we want to limit the outputs of this tutorial...
    if '/.svn' not in e and '/.git' not in e and '/build' not in e:
        print '  + evaluating directory:', e
def onfile(e):
    # we want to limit the outputs of this tutorial...
    if '/.svn' not in e and 'misc.file.' in e and '/.git' not in e and '/build' not in e:
        print '  - evaluating file:', e
def onmatch(e):
    print '  * matching entry:', e
F.find('*/misc.file.*', searchdir, depth=None,
    exclude=('*/library/*', '*/.svn/*', '*.rst', '*/.??*', '*/build/*'),
    matchall=True, ondir=ondir, onfile=onfile, onmatch=onmatch, **fkw)
print
# note that ondir and onfile callbacks are called even if the path does not matches


# Using regular expressions
# =========================

# A simple example
print 'Searching any file starting with "misc.file.":\n'
print ' ','\n  '.join(F.efind('misc\.file\..*', searchdir, **fkw))
print


# An advanced example which uses regex match results
print 'Searching any file starting with "misc.file." and show matching groups:\n'
for filepath, matchobj in F.xefind('misc\.file\.(.*)\.(?:py|rst)', searchdir, getmatch=True, **fkw):
    print '  file path: %s, match groups: %s'%(filepath, matchobj.groups())
print
# grouping is done with parenthesis (), we could also have used named group (?P<groupname1>groupexp) and matchobj.groupdict()
# note that the second group is a non capturing one, we needed this to include both '.py' and '.rst' files
# also note that we use xefind instead of efind, this way we could avoid buildnig a huge list of (filepath,matchobj) couples
# this is not the case here, just an example


os.chdir(PWD)
