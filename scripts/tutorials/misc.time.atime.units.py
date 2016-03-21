import numpy as N
from vacumm.misc.atime import (are_valid_units,  ch_units,  are_same_units, mpl,
    strftime, tz_to_tz, now, to_utc, utc_to_paris)
from vacumm.misc.axes import create_time


# We define units
units = 'hours since 2000-01-15 06:00'

# Is is well formatted?
print are_valid_units(units)
#  -> True

# Same units?
print are_same_units('hours since 2000-1-15 06', units)
#  -> True

# Change axis time units
taxis = create_time(N.arange(6.)*48, units)
# - before
print taxis.units, taxis[0:2]
print taxis.asComponentTime()[0:2]
#  -> hours since 2000-01-15 06:00 [  0.  48.]
#  -> [2000-1-15 6:0:0.0, 2000-1-17 6:0:0.0]
# - change
ch_units(taxis, 'days since 2000-1-15 06', copy=0)
# - after
print taxis.units, taxis[0:2]
print taxis.asComponentTime()[0:2]
#  -> days since 2000-1-15 06 [ 0.  2.]
#  -> [2000-1-15 6:0:0.0, 2000-1-17 6:0:0.0]

# Matplotlib times
taxis_mpl = mpl(taxis)
print taxis_mpl[0], taxis_mpl.units
#  -> 730134.25 days since 0001

# Change the time units of a variable
import MV2
var = MV2.array(MV2.arange(len(taxis)), dtype='f', axes=[taxis])
ch_units(var, 'hours since 2000-01-15 06')
print var.getTime()[0:2]
#  -> [  0.  48.]


# Change time zone
# - UTC time now
t_utc = now(True)
print strftime('%H:%M', t_utc), t_utc.hour
#  -> 15:54 15
# - Paris time
t_paris = utc_to_paris(t_utc)
print strftime('%H:%M', t_paris), t_paris.hour
#  -> 17:54 17
# - back to UTC
print tz_to_tz(t_paris,'Europe/Paris','UTC').hour
#  -> 15
# - we can work on a time string!
print to_utc('2000-10-01 10:20', 'Europe/Paris')
#  -> '2000-10-1 8:20:0.0'
