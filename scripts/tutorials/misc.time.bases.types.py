"""Discover basic time types"""
from __future__ import print_function
import time, datetime, cdtime
from vcmqm import is_datetime, is_cdtime

# Basic time type: time
# - local time
mytime = time.localtime()
print(mytime)
#  -> (2007, 11, 12, 17, 5, 20, 0, 316, 0)
year = mytime[0]
# - string
print(time.asctime())
#  -> Mon Nov 12 17:08:13 2007

# Another time time type: datetime.time
mytime = datetime.datetime(2000, 10, 1, 2)
print(mytime.day, mytime.second)
#  -> 1 0
# - increment
mytime2 = mytime + datetime.timedelta(1, 1)  # (day,second)
print(mytime2.day, mytime.second)
#  -> 2 1

# CDAT time: cdtime
ctime = cdtime.comptime(2000, 10)
print(mytime.year, mytime.month)
#  -> 2000 10

# Check types
print(is_datetime(mytime), is_cdtime(ctime))
# -> True True
print(is_datetime(ctime), is_cdtime(mytime))
# -> False False
