"""Test :func:`vacumm.misc.config.get_spec`"""
from configobj import ConfigObj
from vcmq import get_validator, cfg2rst


# Basic
cfg = """
[secA] # Section A

    [[secAA]] # Section AA

    ages = 2, 5
""".splitlines()
cfg = ConfigObj(cfg)
basic = cfg2rst(cfg, mode='basic')

# Values
values = cfg2rst(cfg, mode='values')

# Specs
cfgspecs = """
[secA] # Section A

    [[secAA]] # Section AA

    ages = integers(default=0, min=0, max=12) # age of children
""".splitlines()
cfgspecs = ConfigObj(cfgspecs, list_values=False,
                 interpolation=False, raise_errors=True)
specs =  cfg2rst(cfgspecs, mode='specs')

