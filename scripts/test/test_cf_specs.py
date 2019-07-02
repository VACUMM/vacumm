"""Test the :class:`~vacumm.misc.cf.CFSpecs` class"""
from vacumm.misc.bases import code_file_name
from vacumm.misc.cf import get_cf_specs, get_cf_cmap, CFSpecs, CFContext
import xarray as xr

# %% Get current specs
cfs = get_cf_specs(cache=False)

# %% Content
assert 'temp' in cfs.variables
ts = cfs.variables['temp']
assert isinstance(ts['name'], list)
assert ts['name'][0] == 'temp'
assert ts['long_name'][0] == 'Temperature'


# %% Register new variable from direct specs
cfs.variables.register('banana', long_name='Good bananas')
assert 'banana' in cfs.variables
assert 'Good bananas' in cfs.variables['banana']['long_name']

# %% Change existing variable
cfs.variables.register('temp', long_name='Another temperature')
assert 'TEMP' in cfs.variables['temp']['name']  # keep default
assert 'Another temperature' in cfs.variables['temp']['long_name']
assert 'banana' in cfs.variables  # don't erase

# %% Register from a config dict
cfgdict = {'variables': {'ssb': {'standard_name': 'sea_surface_banana'}}}
cfs.variables.register_from_cfg(cfgdict)
assert "ssb" in cfs.variables
assert cfs.variables['ssb']['standard_name'][0] == 'sea_surface_banana'
assert cfs.variables['ssb']['long_name'][0] == 'Sea surface banana'

# %% Register from a config file
cfgfile = code_file_name(ext='cfg')
cfs.variables.register_from_cfg(cfgfile)
assert "funnyvar" in cfs.variables
assert "nicetemp" in cfs.variables
assert cfs.variables["nicetemp"]['standard_name'][0] == "nice_temperature"
assert "temp" in cfs.variables["nicetemp"]['name']

# %% Register from a config string
cfgstring = """[variables]
[[nicetemp]]
name=verynicetemp

[[temp]]
name=mytemp
"""
cfs.variables.register_from_cfg(cfgstring)
assert "verynicetemp" in cfs.variables['nicetemp']['name']
assert "mytemp" in cfs.variables['temp']['name']

# %% Calling get_cf_specs
assert get_cf_specs("nicetemp") is cfs.variables["nicetemp"]
assert get_cf_specs("nicetemp", category='coords') is None
assert get_cf_specs("lon", category='coords') is cfs.coords["lon"]

# %% Search
lon = xr.DataArray([1, 2], name='X',
                   attrs={'standard_name': 'longitude_at_u_location'},
                   dims='X')
da = xr.DataArray([10, 20], coords=[lon], dims='X', name="mytemp")
da2 = da.copy()
da2.attrs['long_name'] = 'My salinity'
ds = xr.Dataset({'mytemp': da, 'psal': da2})
assert cfs.variables.is_matching(da, 'temp')
assert not cfs.variables.is_matching(da, 'banana')
assert cfs.variables.is_matching_any(da) == 'temp'
assert cfs.variables.is_matching(ds['psal'], 'sal')
assert cfs.variables.search(ds, 'sal').name == 'psal'
assert cfs.coords.search(ds, 'lon', staggering='uv').name == 'X'

# %% Misc functions
assert get_cf_cmap('temp').name == 'thermal'
assert get_cf_cmap(da).name == 'thermal'

# %% Contexts
cfs_a = CFSpecs({'variables': {'alpha': {}}})
cfs_b = CFSpecs({'variables': {'beta': {}}})
old_cfs = get_cf_specs()
with CFContext(cfs_a) as cfs_c:
    cfs = get_cf_specs()
    assert cfs is not old_cfs
    assert cfs is cfs_c
    assert cfs is cfs_a
    assert 'temp' not in cfs.variables.names
    assert 'alpha' in cfs.variables.names

    with CFContext(cfs_b):
        cfs = get_cf_specs()
        assert cfs is cfs_b
        assert 'beta' in cfs.variables.names
        assert 'alpha' not in cfs.variables.names

    assert get_cf_specs() is cfs_a
assert get_cf_specs() is old_cfs
