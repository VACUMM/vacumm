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

afile = join(tmpdir, 'afile')
count = 2
suffix = '.backup%d'

try:
    
    F.rollover(afile, count=count, suffix=suffix)
    # nothing done, afile does not exists
    with file(afile, 'w') as f: f.write('0')
    # we created afile
    # afile contains 0
    F.rollover(afile, count=count, suffix=suffix)
    # did a copy of afile to afile.backup1
    with file(afile, 'w') as f: f.write('1')
    # afile contains 1
    # afile.backup1 contains 0
    F.rollover(afile, count=count, suffix=suffix)
    # did a copy of afile.backup1 to afile.backup2, and a copy of afile to afile.backup2
    with file(afile, 'w') as f: f.write('2')
    # afile contains 2
    # afile.backup1 contains 1
    # afile.backup2 contains 0
    F.rollover(afile, count=count, suffix=suffix)
    # afile.backup2 removed, copy of afile.backup1 to afile.backup2, copy of afile to afile.backup2
    # afile contains 2
    # afile.backup1 contains 2
    # afile.backup2 contains 1
    
finally:
    # cleaning
    rmtree(tmpdir)