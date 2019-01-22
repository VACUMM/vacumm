"""Manipulate time objects"""

from __future__ import print_function
from datetime import datetime as Datetime
import cdtime
from vcmqm import (create_time, date2num, strtime, comptime, numtime, add,
    adatetime, utc_to_paris, paris_to_utc)

# Comversions
ct = comptime('2000-01')
print(ct.day)
# -> 1
print(comptime(['2000', '2001']))
# -> [2000-1-1 0:0:0.0, 2001-1-1 0:0:0.0]
print(numtime(ct))
# -> 730120.0
print(date2num(Datetime(2000, 1, 1)))
# -> 730120.0
print(strtime([730120.0, cdtime.reltime(1,  'years since 1999')]))
['2000-1-1 0:0:0.0', '2000-1-1 0:0:0.0']

# Time axis
taxis = create_time((0, 3.), 'days since 2000-01-01')
print(adatetime(taxis))
# -> [datetime.datetime(2000, 1, 1, 0, 0),
#     datetime.datetime(2000, 1, 2, 0, 0),
#     datetime.datetime(2000, 1, 3, 0, 0)]

# Additions and substractions
print(add(ct, 2, 'days'))
# -> 2000-1-3 0:0:0.0
print(add(taxis, 1, 'month')[:])
# -> [ 31.  32.  33.]

# Time zones
print(utc_to_paris('2000-01-01 12:00'))
# -> 2000-1-1 13:0:0.0
print(utc_to_paris('2000-06-01 12:00'))
# -> 2000-6-1 14:0:0.0
taxis = create_time((0, 2.), 'hours since 2000-01-01')
print(paris_to_utc(taxis)[:]-taxis[:])
# -> [-1. -1.]
