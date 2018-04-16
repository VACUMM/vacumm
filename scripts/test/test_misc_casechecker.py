"""Test :class:`~vacumm.misc.misc.CaseChecker`"""
from vcmq import CaseChecker, VACUMMError, assert_raises

c1 = CaseChecker(['mode1', 'mode2'], casename='mode')
assert c1.isvalid('mode1') is True
assert c1.isvalid('mode3') is False
assert c1.isvalid(None) is True


c2 = CaseChecker(['-mode1'])
assert c2.isvalid('mode1') is False
assert c2.isvalid('mode3') is True
assert c2.isvalid(None) is True

#assert_raises(VACUMMError, c1.check, 'mode5')
try:
    c1.check('mode5')
except Exception, e: #VACUMMError, e:
    assert e.message == "Invalid mode 'mode5': it must be one of ['mode1', 'mode2']"

c3 = CaseChecker(['mode?', '-mode[56]'])
assert c3.isvalid('mode3') is True
assert c3.isvalid('mode55') is False
assert c3.isvalid('mode5') is False

