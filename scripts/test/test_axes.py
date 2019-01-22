"""Test the :class:`~vacumm.misc.axes` module"""
from vacumm.misc.axes import (islon, islat, isdep, isforecast, istime,
                              create_lon, create_lat, create_dep,
                              create_alt, isalt,
                              create_time, create_forecast, create_axis)


# %% Longitude
lon = create_lon(range(5), long_name='My longitudes')
assert lon.id == 'lon'
assert lon.long_name == 'My longitudes'
assert islon(lon)
lon = create_axis(range(2))
lon.units = 'degrees_east'
assert islon(lon)

# %% Latitude
lat = create_lat(range(5))
assert lat.axis == 'Y'
assert islat(lat)

# %% Depth
dep = create_dep(range(5))
assert dep.units == 'm'
assert isdep(dep)

# %% Altitude
alt = create_alt(range(5))
assert alt.units == 'm'
assert isalt(alt)

# %% Time
time = create_time(['2000', '2001'], units='years since 1990')
assert istime(time)
assert time[0] == 10

# %% Forecast
forecast = create_forecast(range(5), units='days')
assert isforecast(forecast)

