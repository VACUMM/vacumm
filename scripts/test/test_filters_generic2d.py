"""Test function :func:`~vacumm.misc.filters.generic2d`"""
from vcmq import generic2d, N

# Constant no mask
var0 = N.ones((10, 10))
wei0 = N.ones((3, 3))
N.testing.assert_allclose(generic2d(var0, wei0), 1.)

# Constant with mask
var1 = N.ma.array(var0)
var1[0:3, 0:3] = N.ma.masked
wei1 = wei0
#N.testing.assert_allclose(generic2d(var1, wei1).compressed(), 1.)

# Constant with variable weight
var2 = var1
wei2 = wei1.copy()
wei2[1, 1] = 2.
wei2[[0, 2, 2, 0], [0, 0, 2, 2]] = 0.
N.testing.assert_allclose(generic2d(var2, wei2, mask='minimal').compressed(), 1.)


