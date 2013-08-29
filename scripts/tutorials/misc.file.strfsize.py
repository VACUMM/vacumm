#!/usr/bin/env python
# -*- coding: utf-8 -*-

import vacumm.misc.file as F

sizes = (
    1 * 10**3,
    1 * 2**10,
    1 * 10**6,
    1 * 2**20,
    1 * 10**9,
    1 * 2**30,
    1 * 10**12,
    1 * 2**40,
    1 * 10**15,
    1 * 2**50,
)

for size in sizes:
    fsize = F.strfsize(size, si=False)
    fsisize = F.strfsize(size, si=True)
    print 'size: %(size)14d, formatted: %(fsize)8s (CEI, SI: %(fsisize)8s)'%vars()

