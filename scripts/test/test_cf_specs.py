"""Test the :class:`~vacumm.data.cf.VarSpecs` and :class:`~vacumm.data.cf.AxisSpecs` classes"""
from vcmq import code_file_name
from vacumm.data.cf import CF_VAR_SPECS, CF_AXIS_SPECS, register_cf_variable


# Content
assert 'temp' in CF_VAR_SPECS
ts = CF_VAR_SPECS['temp']
assert isinstance(ts ['id'], list)
assert ts['id'][0] == 'temp'
assert ts['long_name'][0] == 'Temperature'


# Register new variable from direct specs
CF_VAR_SPECS.register('banana', long_name='Good bananas')
assert 'banana' in CF_VAR_SPECS
assert 'Good bananas' in CF_VAR_SPECS['banana']['long_name']

# Change existing variable
CF_VAR_SPECS.register('temp', long_name='Another temperature')
assert 'TEMP' in CF_VAR_SPECS['temp']['id'] # keep default
assert 'Another temperature' in CF_VAR_SPECS['temp']['long_name']
assert 'banana' in CF_VAR_SPECS # don't erase

# Register from a config file
cfgfile = code_file_name(ext='cfg')
CF_VAR_SPECS.register_from_cfg(cfgfile)
assert "funnyvar" in CF_VAR_SPECS
assert "nicetemp" in CF_VAR_SPECS
assert "temp" in CF_VAR_SPECS["nicetemp"]['id']

