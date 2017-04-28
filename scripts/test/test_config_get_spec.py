"""Test :func:`vacumm.misc.config.get_spec`"""
from vcmq import get_spec


sspec = 'integers(default=0, min=0, max=10)'

spec = get_spec(sspec)

assert spec['default'] == 0
assert spec['args'] == []
assert spec['funcname'] == 'integers'
assert spec['opttype'] == 'integers'
assert spec['iterable'] == True
assert spec['kwargs'] == {'max': '10', 'min': '0'}
