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

