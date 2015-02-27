#!/usr/bin/env python
# -*- coding: utf-8 -*-

from shutil import rmtree
from tempfile import mkdtemp
from os.path import join

import vacumm.misc.file as F

tmpdir = mkdtemp(dir='.')

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
print 'Empty output'
