"""Test :func:`~vacumm.misc.plot.hov2` with a TZ variable"""

# Imports
from vcmq import N, MV2, create_dep, create_time, hov2, rc

# Init data with z 1D
nt = 10
nz = 8
var = N.dot(N.hanning(nt).reshape(nt, 1), N.hanning(nz).reshape(1, nz))
var = MV2.array(var)
time = create_time((0., nt), units="days since 2000")
z1d = create_dep((-nz+1, 1.))
var.setAxis(0, time)
var.setAxis(1, z1d)
z2d = N.resize(z1d, var.shape)
z2d *= N.resize((N.arange(1., nt+1)/nt).reshape(1, nt), (nz, nt)).T

# Plot with z 1D
rc('font', size=8)
kw = dict(bgcolor='0.5', date_fmt="%a", show=False)
hov2(var, subplot=211, **kw)

# Plot with z 2D
hov2(var, xaxis=z2d, subplot=212, twin='x', **kw)





