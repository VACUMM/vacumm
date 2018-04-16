"""Basic stick plot"""
from vacumm.config import data_sample
from vacumm.misc.plot import stick2

# Read speed
import cdms2
f = cdms2.open(data_sample('mars3d.t.nc'))
u = f('u')
v = f('v')
f.close()

# Plot
stick2(u, v,  title='Current in the Iroise Sea', units='m/s',
    bottom=.2, top=.85, quiver_headwidth=2, close=True,
    quiverkey_value=.5, quiver_width=0.002, quiver_scale=5.,
    savefigs=__file__, figsize=(5.5,  3),  show=False)

