#!/usr/bin/env python
# encoding: utf-8
""" UV-CDAT - Time (time, datetime, cdtime, ...) """

# Imports
import matplotlib
matplotlib.use('qt4agg')

import time, datetime, cdtime

# -----------------------------------------------------------------------------------------------------------
# ---- Basic class: time
print 10*'-'+' time '+10*'-'
mytime = time.localtime()
print 'Heure locale: ', time.asctime()

# ---- Class as 
print 10*'-'+' datetime '+10*'-'


print 10*'-'+' cdtime '+10*'-'

# sys.exit() # End of the execution
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

