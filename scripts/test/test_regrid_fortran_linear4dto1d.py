"""Test the fortran function :f:func:`linear4dto1d`"""
from utils import np
from vacumm.fortran.interp import linear4dto1d

nxi = 7
nyi = 6
nzi = 5
nti = 4
no = 100

vfunc = lambda t, z, y, x: (1*x + 2.35*y + 3.24*z -0.65*t)

tti, zzi, yyi, xxi = np.mgrid[0:nti-1:nti*1j, 0:nzi-1:nzi*1j,
                              0:nyi-1:nyi*1j, 0:nxi-1:nxi*1j]
xi = xxi[0, 0, 0, :]
yi = yyi[0, 0, :, 0]
zi = zzi[0, :, 0, 0]
ti = tti[:, 0, 0, 0]
vi = vfunc(tti, zzi, yyi, xxi)

np.random.seed(0)
xyztomin = -0.5
xo = np.random.uniform(xyztomin, nxi-1.5, no)
yo = np.random.uniform(xyztomin, nyi-1.5, no)
zo = np.random.uniform(xyztomin, nzi-1.5, no)
to = np.random.uniform(xyztomin, nti-1.5, no)
vo_truth = vfunc(to, zo, yo, xo)
out = (xo < 0) | (yo < 0) | (zo < 0) | (to < 0)

vo_interp = linear4dto1d(xi, yi, zi, ti, vi, xo, yo, zo, to)

vo_truth[out] = np.nan
np.testing.assert_almost_equal(vo_interp.ravel(), vo_truth.ravel())
