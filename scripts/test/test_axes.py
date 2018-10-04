"""Test the :class:`~vacumm.misc.axes` module"""
from vacumm.misc.axes import (islon, islat, isdep, isforecast, istime,
                              create_lon, create_lat, create_dep,
                              create_time, create_forecast)


#%% Create
lon = create_lon(range(5), long_name='My longitudes')
assert lon.id == 'lon'
assert lon.long_name == 'My longitudes'
assert islon(lon)

