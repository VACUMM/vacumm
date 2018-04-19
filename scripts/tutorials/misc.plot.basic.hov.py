from __future__ import print_function
import cdms2
from vcmq import hov2

# Lecture de l'elevation de la surface
from vacumm.config import data_sample
select = dict(time=slice(0, 13), squeeze=1)
f = cdms2.open(data_sample('mars3d.xt.xe.nc'))
xe = f('xe', **select)
f.close()


# Plot
h = hov2(xe, fmt='%.1f', left=.15, right=.98, linewidths=1.5,
    clabel_glow=2, nmax=20, figsize=(6, 6), close=False,
    fill='contour', show=False)

