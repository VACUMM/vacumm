"""Basic stick plot"""
from vcmq import data_sample, stick2, cdms2

# Read speed
f = cdms2.open(data_sample('mars3d.t.nc'))
u = f('u')
v = f('v')
f.close()

# Plot
stick2(u, v,  title='Current in the Iroise Sea', units='m/s',
    bottom=.2, top=.85, quiver_headwidth=2, close=False,
    quiverkey_value=.5, quiver_width=0.002, quiver_scale=5.,
    figsize=(5.5,  3), show=False)

