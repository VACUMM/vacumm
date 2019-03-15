"""Test the :class:`~vacumm.misc.cf.SGLocator` class"""
from vacumm.misc.cf import SGLocator
import xarray as xr
import numpy as np

# %% Minimal
sgl = SGLocator()
assert sgl.loc == ''
assert sgl == ''

# %% Explicit
sgl = SGLocator('u')
assert sgl == sgl.loc == 'u'

# %% From attributes
sgln = SGLocator.from_name('sst_u')
assert sgln == 'u'
assert sgln.attrs['name'] == 'sst'
sgls = SGLocator.from_standard_name('sst_at_v_location')
assert sgls == 'v'
assert sgls.attrs['standard_name'] == 'sst'
sgll = SGLocator.from_long_name('SST at W location')
assert sgll == 'w'
assert sgll.attrs['long_name'] == 'SST'

# %% From data array
da = xr.DataArray(np.ones((2, 3, 4))*10, dims=('time', 'y', 'x'), name='sst_u')
assert SGLocator.from_dataarray(da) == 'u'
da = xr.DataArray(np.ones((2, 3, 4))*10, dims=('time', 'y', 'x'),
                  attrs={'standard_name': 'sst_at_v_location'})
sgl = SGLocator.from_dataarray(da)
assert sgl == 'v'
assert sgl.attrs['standard_name'] == 'sst'
da = xr.DataArray(np.ones((2, 3, 4))*10, dims=('time', 'y', 'x'),
                  attrs={'long_name': 'SST at W location'})
assert SGLocator.from_dataarray(da) == 'w'

# %% To
sglu = SGLocator('u', attrs={'name': 'temp',
                             'standard_name': 'sea_water_temperature',
                             'long_name': 'Temperature',
                             'foo': 'bar'})
sglv = sglu.to('v')
assert sglv == 'v'
assert sglv.attrs['name'] == 'temp'
assert sglv.attrs['foo'] == 'bar'

# %% Format attributes
attrs = sglv.format_attrs()
assert attrs['name'] == 'temp_v'
assert attrs['standard_name'] == 'sea_water_temperature_at_v_location'
assert attrs['long_name'] == 'Temperature at V location'

# %% Format data array
lon = xr.DataArray(np.zeros((3, 4)), dims=('y', 'x'),
                   attrs={'long_name': 'Longitude'})
lat = xr.DataArray(np.zeros((3, 4)), dims=('y', 'x'))
da = xr.DataArray(np.ones((2, 3, 4))*10,
                  dims=('time', 'y', 'x'),
                  coords={'lon': lon, 'lat': lat},
                  attrs={'name': 'temp',
                         'standard_name': 'sea_water_temperature',
                         'long_name': 'Temperature',
                         'foo': 'bar'})
da = SGLocator('v').format_dataarray(da, coords=('lat', 'lon'),
                                     dims=('x', 'y'))
assert da.name == 'temp_v'
assert da.standard_name == 'sea_water_temperature_at_v_location'
assert 'lon_v' in da.coords
assert da.lon_v.long_name == 'Longitude at V location'
assert 'time' in da.dims
assert 'y_v' in da.dims

# %% Format cf specs
specs = {'name': ['temp', 'temperature', 'TEMP'],
         'standard_name': ['sea_water_temperature',
                           'sea_water_potential_temperature'],
         'long_name': ['Temperature'],
         'units': ['degrees_celsius'],
         'cmap': "thermal",
         'coords': {
                 'x': 'lon',
                 'y': 'lat',
                 'z': {'ocean': 'depth'},
                 },
         }
new_specs = SGLocator('v').format_cf_specs(specs)
assert new_specs == {
        'name': ['temp_v', 'temperature_v', 'TEMP_v'],
        'standard_name': ['sea_water_temperature_at_v_location',
                          'sea_water_potential_temperature_at_v_location'],
        'long_name': ['Temperature at V location'],
        'units': ['degrees_celsius'],
        'cmap': 'thermal',
        'coords': {'x': 'lon_v', 'y': 'lat_v', 'z': {'ocean': 'depth'}}}
