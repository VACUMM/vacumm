"""Test :func:`~vacumm.misc.plot.section2` in quiver mode"""

# Imports
from vcmq import create_lon, N, MV2, create_dep, os, code_file_name, section2
from vacumm.misc.phys.units import deg2m

# Init data with z 1D
nz = 11
nx = 11
x = create_lon(N.arange(nx))
xm = deg2m(x[:],lat=45.) # meters
dx = xm[:].ptp()
z = create_dep((-nz+1, 1.), units='m', long_name='Depth')
dz = z[:].ptp()
scale = dz/dx
u = MV2.ones((nz,nx)) # 1 m/s
w = u*scale           # 1 m/s * scale
for var in u,w:
    var.setAxis(0, z)
    var.setAxis(1, x)
    var.units = 'm/s'

# Plot
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
s = section2((u,w), quiver_norm=1, fill=False, show=False,
    savefig=figfile, close=True, axes_aspect=1)



