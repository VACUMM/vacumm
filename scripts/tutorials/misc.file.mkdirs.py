#!/usr/bin/env python

from os.path import abspath, dirname, join, realpath
from shutil import rmtree
from tempfile import mkdtemp

import vacumm.misc.file as F

# This file
# realpath() is not really required in this example, this is to advice that in case
# you were executing a link to this file, the code would still work (it can be a good
# practice to do so in your scripts if you want that kind of support)
thisfile = realpath(__file__)

# The directory containing this file
thisdir = dirname(__file__)

tmpdir = mkdtemp(dir=thisdir)

adir = join(tmpdir, 'path/to/adir')
afile = join(adir, 'afile')

anotherdir = join(tmpdir, 'path/to/anotherdir')
anotherfile = join(anotherdir, 'anotherfile')

try:
    
    print 'Created directory:', F.mkdirs(adir)
    F.mkfdirs(afile) # no effect, already created
    print 'Created directories:',  F.mkfdirs((afile, anotherfile))
    F.mkdirs((adir, anotherdir)) # no effect, already created
    
finally:
    rmtree(tmpdir)

