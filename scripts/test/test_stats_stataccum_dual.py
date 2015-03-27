"""Test :class:`~vacumm.misc.stats.StatAccum` for a pair of variables"""

# Imports
from vcmq import MV2, N
from vacumm.misc.stats import StatAccum

# Setup masked data
nt, ny, nx = 20, 15, 10
var1 = MV2.arange(1.*nt*ny*nx).reshape(nt, ny, nx)
var1.getAxis(0).designateTime()
var1.getAxis(0).units = 'months since 2000'
var2 = var1.clone()
var2[:] += N.ma.sin(var1)*100
var1[3:13, 5:9, 3:7] = MV2.masked
var2[5:15, 7:12, 1:5] = MV2.masked
mask = var1.mask|var2.mask # common mask

# Direct
varm1 = var1.asma().copy()
varm2 = var2.asma().copy()
varm1[mask] = N.ma.masked
varm2[mask] = N.ma.masked
var2dm1 = varm1.reshape(nt, -1)
var2dm2 = varm2.reshape(nt, -1)
dstats = dict(
    tavail = 1.*(~varm1.mask).sum(axis=0)/nt,
    tmean = (varm1.mean(axis=0), varm2.mean(axis=0)),
    tbias = varm2.mean(axis=0)-varm1.mean(axis=0),
    tstd = (varm1.std(axis=0), varm2.std(axis=0)),
    trms = N.ma.sqrt(((varm1-varm2)**2).mean(axis=0)),
    tcrms = (varm1-varm2).std(axis=0),
    tcorr = ((varm1-varm1.mean(axis=0))*(varm2-varm2.mean(axis=0))).mean(axis=0)/ \
        (varm1.std(axis=0)*varm2.std(axis=0)),
    savail = 1.*(~var2dm1.mask).sum(axis=1)/var2dm1.shape[1],
    smean = (var2dm1.mean(axis=1), var2dm2.mean(axis=1)),
    sbias = var2dm2.mean(axis=1)-var2dm1.mean(axis=1),
    sstd = (var2dm1.std(axis=1), var2dm2.std(axis=1)),
    srms = N.ma.sqrt(((var2dm1-var2dm2)**2).mean(axis=1)),
    scrms = (var2dm1-var2dm2).std(axis=1),
    scorr = ((var2dm1.T-var2dm1.T.mean(axis=0))*(var2dm2.T-var2dm2.T.mean(axis=0))).mean(axis=0)/ \
        (var2dm1.std(axis=1)*var2dm2.std(axis=1)),
)

# Indirect using StatAccum
sa = StatAccum(tall=True, sall=True)
sa += var1[:10], var2[:10]
sa += var1[10:], var2[10:]
istats = sa.get_stats()

# Unitest
result = []
for key in sorted(dstats.keys()):
    if isinstance(dstats[key], tuple):
        value = ()
        for i in 0, 1:
            value += N.ma.allclose(dstats[key][i], istats[key][i]),
    else:
        value = N.ma.allclose(dstats[key], istats[key])
    result.append(('assertTrue', value))
