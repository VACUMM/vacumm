"""Test the :class:`~vacumm.misc.cf.VarSpecs` and :class:`~vacumm.misc.cf.AxisSpecs` classes"""
from vacumm.misc.bases import code_file_name
from vacumm.misc.cf import get_cf_specs

# %% Get current specs
cf_specs = get_cf_specs()

# %% Content
cf_var_specs = cf_specs.variables
assert 'temp' in cf_var_specs
ts = cf_var_specs['temp']
assert isinstance(ts['name'], list)
assert ts['name'][0] == 'temp'
assert ts['long_name'][0] == 'Temperature'

# %% Register new variable from direct specs
cf_var_specs.register('banana', long_name='Good bananas')
assert 'banana' in cf_var_specs
assert 'Good bananas' in cf_var_specs['banana']['long_name']

# %% Change existing variable
cf_var_specs.register('temp', long_name='Another temperature')
assert 'TEMP' in cf_var_specs['temp']['name']  # keep default
assert 'Another temperature' in cf_var_specs['temp']['long_name']
assert 'banana' in cf_var_specs  # don't erase

# %% Register from a config dict
cfgdict = {'variables': {'ssb': {'standard_name': 'sea_surface_banana'}}}
cf_var_specs.register_from_cfg(cfgdict)
assert "ssb" in cf_var_specs
assert cf_var_specs['ssb']['standard_name'][0] == 'sea_surface_banana'
assert cf_var_specs['ssb']['long_name'][0] == 'Sea surface banana'

# %% Register from a config file
cfgfile = code_file_name(ext='cfg')
cf_var_specs.register_from_cfg(cfgfile)
assert "funnyvar" in cf_var_specs
assert "nicetemp" in cf_var_specs
assert cf_var_specs["nicetemp"]['standard_name'][0] == "nice_temperature"
assert "temp" in cf_var_specs["nicetemp"]['name']

# %% Register from a config string
cfgstring = """[variables]
[[nicetemp]]
name=verynicetemp
"""
cf_var_specs.register_from_cfg(cfgstring)
assert "verynicetemp" in cf_var_specs['nicetemp']['name']

# %% Other checks
assert get_cf_specs("nicetemp") is cf_var_specs["nicetemp"]
assert get_cf_specs("nicetemp", category='coords') is None

##%% Copies
#new_var_specs = cf_var_specs.copy()
#assert "nicetemp" in new_var_specs
#cfg_patch = """[variables]
#[[nicetemp]]
#name=verynicetemp
#"""
#CF_SPECS_PATCHED = {}
#very_new_var_specs = new_var_specs.copy_and_update(cfg=cfg_patch,
#                                                   parent=CF_SPECS_PATCHED)
#assert "nicetemp" in very_new_var_specs
#assert very_new_var_specs["nicetemp"] is not new_var_specs["nicetemp"]
#assert "verynicetemp" in very_new_var_specs["nicetemp"]["name"]
#
