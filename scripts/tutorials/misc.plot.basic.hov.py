# Lecture de l'elevation de la surface
import cdms2, numpy as N, MV2
from vacumm.config import data_sample
select = dict(time=slice(0, 13), squeeze=1)
f = cdms2.open(data_sample('mars3d.xt.xe.nc'))
xe = f('xe', **select)
f.close()


# Plot
from vacumm.misc.plot import hov2
import vacumm.misc.plot
print vacumm.misc.plot
h = hov2(xe, fmt='%.1f', left=.15, right=.98, linewidths=1.5,
    clabel_glow=2, nmax=20, figsize=(6, 6), close=True,
    fill='contour', show=False, savefigs=__file__, savefigs_pdf=True)

