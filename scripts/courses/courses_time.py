#!/usr/bin/env python
# -*- coding: utf8 -*- 
""" UV-CDAT - Time (time, datetime, cdtime, ...) """

# Imports
import matplotlib
matplotlib.use('qt4agg')

import time, datetime, cdtime
from vacumm.misc.atime import lindates, IterDates, Intervals,strftime,strptime, create_time, ch_units
import numpy as N

# -----------------------------------------------------------------------------------------------------------
# ---- Basic class: time (from Python) ----
print 10*'-'+' time '+10*'-'
mytime = time.localtime()
print 'Heure locale: ', time.asctime()

# ---- Datetime (from Matplotlib) ---- 
print 10*'-'+' datetime '+10*'-'
# !! do not work for dates before 1900 !!
mytime = datetime.datetime(2000,10,1,2)
print mytime.day,mytime.second
# - increment
mytime2 = mytime + datetime.timedelta(1,1) # (day,second)
print mytime2.day,mytime.second

# ---- Cdtime (from UV-CDAT) ----
print 10*'-'+' cdtime '+10*'-'
# ---- Create an object 'comptime' : absolute time
# - everything specified (year, month, day, hour, minute, second)
ctime = cdtime.comptime(2000,1,1,0,0,0)
# => Practice: try to create time only with a reduced number of arguments
# - From a string
ctime = cdtime.s2c('2000-1-1 0')
# - we check
print ctime
# - we check the numerical values
print ctime.year,ctime.day

# ---- Create an object 'reltime' : relative time
# - arg given explicitly (value, CF units)
rtime = cdtime.reltime(50, 'years since 1950')
print '%s | %s | %s' %(rtime,rtime.value,rtime.units)
# - from sting and units
rtime = cdtime.s2r('2000-1-1 0','years since 1950')
print rtime.value

# ---- Operations
# - add/subtract
print ctime.add(1,cdtime.Year),'|',rtime.add(-1,cdtime.Year)
# - convert
rtime2 = ctime.torel('days since 2000')
ctime2 = rtime.tocomp().add(1,cdtime.Year)
# - compare
print rtime2 == rtime
print ctime2 <= ctime


# => Practice: Create a cdtime array and print the most recent date.

# ---- Check types
from vacumm.misc.atime import is_comptime,is_reltime,is_cdtime
print is_comptime(ctime),is_reltime(rtime),is_cdtime(ctime)


# ---- VACUMM Bonus ----
# Read from a string and a format
# => Practice: check strptime in google
mytime = strptime('1950-01-01 07:00:00','%Y-%m-%d %H:%M:%S') # => Practice: Try different formats
# - Check
print mytime.year,mytime.minute

# We choose the french language 
import locale
locale.setlocale(locale.LC_ALL,'fr_FR')

# Write in a different format
print strftime('%e %B %Y a %Hh%M',mytime) # => Practice: Try different formats

# ---- Time axes ---- 
print 10*'-'+' Time axes '+10*'-'
units = 'hours since 2000-01-15 06:00'
taxis = create_time(N.arange(6.)*48, units)
print taxis.getValue()
print 'Units before change: ',taxis.units
# - change units
ch_units(taxis, 'days since 2000-1-15 06:00', copy=0)
print taxis.getValue()
print 'Units after change: ',taxis.units

# ------------------------------------------------------------------------------------------------------------
# ---- lindates ----
print 10*'-'+' lindates '+10*'-'
dates = lindates('2000', '2002-05', 3, 'months')
print dates
# => Practice: Try different arguments.

# ---- Intervals ----
print 10*'-'+' Intervals '+10*'-'
for itv in Intervals(('2000','2001','co'),(2,'month')): print itv # => Practice: try different Intervals (2 days, 2 hours)
print '-- List ...'
print Intervals(('2000','2001','co'),12).tolist()

print '-- Reverse ...'
for itv in Intervals((cdtime.comptime(2000), '2001'), 
    'month',reverse=True): print itv

# ---- IterDates ----
print 10*'-'+' IterDates '+10*'-'
print '-- Example 1 --'
for date in IterDates(('2000','2001'),(1,'month')): print date
# => Practice: Iteration on every days in January 2012.




# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

