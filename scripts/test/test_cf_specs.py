"""Test the :class:`~vacumm.misc.cf.VarSpecs` and :class:`~vacumm.misc.cf.AxisSpecs` classes"""
from vcmq import code_file_name, MV2
from vacumm.misc.cf import CF_VAR_SPECS, get_cf_specs

# %% Content
assert 'temp' in CF_VAR_SPECS
ts = CF_VAR_SPECS['temp']
assert isinstance(ts['name'], list)
assert ts['name'][0] == 'temp'
assert ts['long_name'][0] == 'Temperature'

# %% Register new variable from direct specs
CF_VAR_SPECS.register('banana', long_name='Good bananas')
assert 'banana' in CF_VAR_SPECS
assert 'Good bananas' in CF_VAR_SPECS['banana']['long_name']

# %% Change existing variable
CF_VAR_SPECS.register('temp', long_name='Another temperature')
assert 'TEMP' in CF_VAR_SPECS['temp']['name']  # keep default
assert 'Another temperature' in CF_VAR_SPECS['temp']['long_name']
assert 'banana' in CF_VAR_SPECS  # don't erase

# %% Register from a config file
cfgfile = code_file_name(ext='cfg')
CF_VAR_SPECS.register_from_cfg(cfgfile)
assert "funnyvar" in CF_VAR_SPECS
assert "nicetemp" in CF_VAR_SPECS
assert CF_VAR_SPECS["nicetemp"]['standard_name'][0] == "nice_temperature"
assert "temp" in CF_VAR_SPECS["nicetemp"]['name']

# %% Other checks
assert get_cf_specs("nicetemp") is CF_VAR_SPECS["nicetemp"]
assert get_cf_specs("nicetemp", category='axes') is None

#%% Copies
new_var_specs = CF_VAR_SPECS.copy()
assert "nicetemp" in new_var_specs
cfg_patch = """[variables]
[[nicetemp]]
name=verynicetemp
"""
CF_SPECS_PATCHED = {}
very_new_var_specs = new_var_specs.copy_and_update(cfg=cfg_patch,
                                                   parent=CF_SPECS_PATCHED)
assert "nicetemp" in very_new_var_specs
assert very_new_var_specs["nicetemp"] is not new_var_specs["nicetemp"]
assert "verynicetemp" in very_new_var_specs["nicetemp"]["name"]

