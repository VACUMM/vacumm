"""Test :func:`~vacumm.misc.plot.section2` with a Z- variable"""

# Imports
from vcmq import N, MV2, cdms2, create_dep, rc, section2, code_file_name, os

# Init data with z 1D
nz = 8
nd = 10
var = N.dot(N.hanning(nz).reshape(nz, 1), N.hanning(nd).reshape(1, nd))
var = MV2.array(var)
d = cdms2.createAxis(N.arange(nd))
d.units='km'
d.long_name='Distance'
z1d = create_dep((-nz+1, 1.))
var.setAxis(0, z1d)
var.setAxis(1, d)
z2d = N.resize(z1d[:].reshape(1, nz), (nd, nz)).T
z2d *= N.arange(1., nd+1)/nd

# Plot with z 1D
rc('font', size=8)
kw = dict(show=False, bgcolor='0.5')
section2(var, subplot=211, **kw)

# Plot with z 2D
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
section2(var, yaxis=z2d, subplot=212, savefig=figfile, close=True, **kw)

# Result
result = dict(files=figfile)


