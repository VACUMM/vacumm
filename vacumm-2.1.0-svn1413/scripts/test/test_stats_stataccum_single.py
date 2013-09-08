"""Test :class:`~vacumm.misc.stats.StatAccum` for a single variable"""

# Imports
from vcmq import MV2, N
from vacumm.misc.stats import StatAccum

# Setup masked data
nt, ny, nx = 20, 15, 10
var = MV2.arange(1.*nt*ny*nx).reshape(nt, ny, nx)
var[3:13, 5:9, 3:7] = MV2.masked
var.getAxis(0).designateTime()
var.getAxis(0).units = 'months since 2000'

# Direct
var2d = var.asma().reshape(nt, -1)
dstats = dict(
    tavail = 1.*(~var.mask).sum(axis=0)/var.shape[0], 
    tmean = var.mean(axis=0), 
    tstd = var.std(axis=0), 
    savail = 1.*(~var2d.mask).sum(axis=1)/var2d.shape[1], 
    smean = var2d.mean(axis=1), 
    sstd = var2d.std(axis=1), 
)

# Indirect using StatAccum
sa = StatAccum(
#    tmean=True, tstd=True, tavail=True, 
    smean=True)
sa += var[:10]
sa += var[10:]
istats = sa.get_stats()

# Unitest
result = []
for key in sorted(dstats.keys()):
    value = N.ma.allclose(dstats[key], istats[key])
    print key, value
    result.append(('assertTrue', value))
