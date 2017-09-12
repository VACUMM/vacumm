"""Test the fortran function :f:func:`linear4dto1d`"""
from vcmq import N
from vacumm.misc.grid._interp_ import linear4dto1d

nxi = 7
nyi = 6
nzi = 5
nti = 4
mv = 1.e20
no = 100

vfunc = lambda t, z, y, x: (1*x + 2.35*y + 3.24*z -0.65*t)

tti, zzi, yyi, xxi = N.mgrid[0:nti-1:nti*1j, 0:nzi-1:nzi*1j,
    0:nyi-1:nyi*1j, 0:nxi-1:nxi*1j]
xi = xxi[0, 0, 0, :]
yi = yyi[0, 0, :, 0]
zi = zzi[0, :, 0, 0]
ti = tti[:, 0, 0, 0]
vi = vfunc(tti, zzi, yyi, xxi)

N.random.seed(0)
xyztomin = -0.5
xo = N.random.uniform(xyztomin, nxi-1.5, no)
yo = N.random.uniform(xyztomin, nyi-1.5, no)
zo = N.random.uniform(xyztomin, nzi-1.5, no)
to = N.random.uniform(xyztomin, nti-1.5, no)
vo_truth = N.ma.array(vfunc(to, zo, yo, xo))
out = (xo<0)|(yo<0)|(zo<0)|(to<0)

mv = 1e20
vo_interp = linear4dto1d(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp = N.ma.masked_values(vo_interp, mv)

vo_truth = N.ma.masked_where(out, vo_truth)
N.testing.assert_almost_equal(vo_interp.compressed(), vo_truth.compressed())
