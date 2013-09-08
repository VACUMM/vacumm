#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os.path import abspath, dirname, join, realpath
from shutil import rmtree
from tempfile import mkdtemp

import vacumm.misc.file as F

tmpdir = mkdtemp(dir='.')

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

