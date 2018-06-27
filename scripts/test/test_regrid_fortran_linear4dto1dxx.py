"""Test the fortran function :f:func:`linear4dto1dxx`"""
from vcmq import N, P, meshcells, minmax, code_file_name, os
from vacumm.misc.grid._interp_ import linear4dto1dxx
# vo = linear4dto1dxx(xxi,yyi,zzi,ti,vi,xo,yo,zo,to,mv,[nxi,nyi,nyix,nxiy,nyiz,nxiz,nzi,nti,ntiz,no,nex])

# %% Bases

nex = 4
nexz = 2
nxi = 7
nyi = 6
nzi = 5
nti = 4
mv = 1.e20
no = 1
#no = 1

vfunc = lambda t, z, y, x: (1*x + 2.35*y + 3.24*z -0.65*t)

tti, zzi, yyi, xxi = N.mgrid[0:nti-1:nti*1j, 0:nzi-1:nzi*1j,
    0:nyi-1:nyi*1j, 0:nxi-1:nxi*1j]

zzi = zzi[None]
zzi = N.repeat(zzi, nexz, axis=0)


# %% Pure 1D axes

xi = xxi[0, 0, 0:1, :] # (nyix=1,nxi)
yi = yyi[0, 0, :, 0:1] # (nyi,nxiy=1)
zi = zzi[0:1, 0:1, :, 0:1, 0:1] # (nexz=1,ntiz=1,nzi,nyiz=1,nxiz=1)
ti = tti[:, 0, 0, 0] # (nti)
vi = vfunc(tti, zzi, yyi, xxi)
vi = N.resize(vi, (nex, )+vi.shape[1:]) # (nex,nti,nzi,nyi,nxi)

N.random.seed(0)
xyztomin = -0.5
xo = N.random.uniform(xyztomin, nxi-1.5, no)
yo = N.random.uniform(xyztomin, nyi-1.5, no)
zo = N.random.uniform(xyztomin, nzi-1.5, no)
to = N.random.uniform(xyztomin, nti-1.5, no)
vo_truth = N.ma.array(vfunc(to, zo, yo, xo))

mv = 1e20
vo_interp = linear4dto1dxx(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp = N.ma.masked_values(vo_interp, mv)

N.testing.assert_almost_equal(vo_interp[0], vo_truth)


# %% Single point in space

xi = xxi[0, 0, 0:1, :1] # (nyix=1,nxi)
yi = yyi[0, 0, :1, 0:1] # (nyi,nxiy=1)
zi = zzi[0:1, 0:1, :, 0:1, 0:1] # (nexz=1,ntiz=1,nzi,nyiz=1,nxiz=1)
ti = tti[:, 0, 0, 0] # (nti)
vi = vfunc(tti, zzi, yyi, xxi)[:, :, :, :1, :1]
vi = N.resize(vi, (nex, )+vi.shape[1:]) # (nex,nti,nzi,1,1)

N.random.seed(0)
xyztomin = -0.5
xo = N.random.uniform(xyztomin, nxi-1.5, no)
yo = N.random.uniform(xyztomin, nyi-1.5, no)
zo = N.random.uniform(xyztomin, nzi-1.5, no)
to = N.random.uniform(xyztomin, nti-1.5, no)
vo_truth = N.ma.array(vfunc(to, zo, yi[0], xi[0]))

mv = 1e20
vo_interp = linear4dto1dxx(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp = N.ma.masked_values(vo_interp, mv)

N.testing.assert_almost_equal(vo_interp[0], vo_truth)


# %% Constant time

xi = xxi[0, 0, 0:1, :] # (nyix=1,nxi)
yi = yyi[0, 0, :, 0:1] # (nyi,nxiy=1)
zi = zzi[0:1, 0:1, :, 0:1, 0:1] # (ntiz=1,nzi,nyiz=1,nxiz=1)
ti = tti[:1, 0, 0, 0] # (1)
vi = vfunc(tti, zzi, yyi, xxi)[:, :1, :, :, :]
vi = N.resize(vi, (nex, )+vi.shape[1:]) # (nex,1,nzi,nyi,nxi)

N.random.seed(0)
xyztomin = -0.5
xo = N.random.uniform(xyztomin, nxi-1.5, no)
yo = N.random.uniform(xyztomin, nyi-1.5, no)
zo = N.random.uniform(xyztomin, nzi-1.5, no)
to = N.random.uniform(xyztomin, nti-1.5, no)
vo_truth = N.ma.array(vfunc(ti, zo, yo, xo))

mv = 1e20
vo_interp = linear4dto1dxx(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp = N.ma.masked_values(vo_interp, mv)

N.testing.assert_almost_equal(vo_interp[0], vo_truth)


# %% Variable depth with 1D X/Y + T

xi = xxi[0, 0, 0:1, :] # (nyix=1,nxi)
yi = yyi[0, 0, :, 0:1] # (nyi,nxiy=1)
zi = zzi[:, :, :, :, :] # (nexz,ntiz=nti,nzi,nyiz=nyi,nxiz=nxi)
ti = tti[:, 0, 0, 0] # (nti)
vi = vfunc(tti, zzi, yyi, xxi)
vi = N.resize(vi, (nex, )+vi.shape[1:]) # (nex,nti,nzi,nyi,nxi)

N.random.seed(0)
xyztomin = -0.5
xo = N.random.uniform(xyztomin, nxi-1.5, no)
yo = N.random.uniform(xyztomin, nyi-1.5, no)
zo = N.random.uniform(xyztomin, nzi-1.5, no)
to = N.random.uniform(xyztomin, nti-1.5, no)
vo_truth = N.ma.array(vfunc(to, zo, yo, xo))

mv = 1e20
vo_interp = linear4dto1dxx(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp = N.ma.masked_values(vo_interp, mv)

N.testing.assert_almost_equal(vo_interp[0], vo_truth)


# %% 2D X/Y with no other axes (pure curvilinear)

xi = xxi[0, 0] # (nyix=nyi,nxi)
yi = yyi[0, 0]# (nyi,nxiy=nxi)
zi = zzi[0:1, 0:1, 0:1, 0:1, 0:1] # (nexz=1,ntiz=1,1,nyiz=1,nxiz=1)
ti = tti[:1, 0, 0, 0] # (1)

vi = vfunc(tti, zzi, yyi, xxi)[:, :1, :1, :, :]
vi = N.resize(vi, (nex, )+vi.shape[1:]) # (nex,1,1,nyi,nxi)

N.random.seed(0)
xyztomin = -0.5
xo = N.random.uniform(xyztomin, nxi-1.5, no)
yo = N.random.uniform(xyztomin, nyi-1.5, no)
zo = N.random.uniform(xyztomin, nzi-1.5, no)
to = N.random.uniform(xyztomin, nti-1.5, no)
vo_truth = N.ma.array(vfunc(ti, zi.ravel()[0], yo, xo))

mv = 1e20
vo_interp = linear4dto1dxx(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp_rect = linear4dto1dxx(xi[:1],yi[:, :1],zi,ti,vi,xo,yo,zo,to,mv=mv)

N.testing.assert_almost_equal(vo_interp[0], vo_truth)
N.testing.assert_almost_equal(vo_interp_rect[0], vo_truth)

# %% Same coordinates

xi = xxi[0, 0, 0:1, :] # (nyix=1,nxi)
yi = yyi[0, 0, :, 0:1] # (nyi,nxiy=1)
zi = zzi[0:1, 0:1, :, 0:1, 0:1] # (nexz=1,ntiz=1,nzi,nyiz=1,nxiz=1)
ti = tti[:, 0, 0, 0] # (nti)
vi = vfunc(tti, zzi, yyi, xxi)
vi = N.resize(vi, (nex, )+vi.shape[1:]) # (nex,nti,nzi,nyi,nxi)

xyzo = N.meshgrid(xi, yi, zi, ti, indexing='ij')
xo = xyzo[0].ravel()
yo = xyzo[1].ravel()
zo = xyzo[2].ravel()
to = xyzo[3].ravel()
vo_truth = N.ma.array(vfunc(to, zo, yo, xo))

mv = 1e20
vo_interp = linear4dto1dxx(xi,yi,zi,ti,vi,xo,yo,zo,to,mv=mv)
vo_interp = N.ma.masked_values(vo_interp, mv)

#N.testing.assert_almost_equal(vo_interp[0], vo_truth)




