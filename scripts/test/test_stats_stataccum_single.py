"""Test :class:`~vacumm.misc.stats.StatAccum` for a single variable"""

# Imports
from vcmq import cdms2, MV2, N, StatAccum

# Setup masked data
nt, ny, nx = 20, 15, 10
var = MV2.arange(1.*nt*ny*nx).reshape(nt, ny, nx)
var[3:13, 5:9, 3:7] = MV2.masked
var.getAxis(0).designateTime()
var.getAxis(0).units = 'months since 2000'
vmax = var.max()
bins = N.linspace(-0.1*vmax, 0.9*vmax, 14)
nbins = len(bins)
maskyx = var.mask.all(axis=0)
maskt = var.mask.reshape((nt, -1)).all(axis=1)

# Direct
varm = var.asma()
var2d = var.asma().reshape(nt, -1)
dstats = dict(
    tavail = 1.*(~var.mask).sum(axis=0)/var.shape[0],
    tmean = var.mean(axis=0),
    tstd = var.std(axis=0),
    tmin = var.min(axis=0),
    tmax = var.max(axis=0),
    savail = 1.*(~var2d.mask).sum(axis=1)/var2d.shape[1],
    smean = var2d.mean(axis=1),
    sstd = var2d.std(axis=1),
    smin = var2d.min(axis=1),
    smax = var2d.max(axis=1),
)
dstats['thist'] = N.ma.zeros((nbins-1, ny, nx), 'l')
dstats['shist'] = N.ma.zeros((nt, nbins-1), 'l')
for ibin in range(nbins-1):
    valid = (varm>=bins[ibin])&(varm<bins[ibin+1])
    dstats['thist'][ibin] = valid.filled(False).sum(axis=0)
    dstats['shist'][:, ibin] = valid.filled(False).reshape((nt, -1)).sum(axis=1)
dstats['thist'][:, maskyx] = N.ma.masked
dstats['shist'][maskt] = N.ma.masked


# Indirect using StatAccum
sa = StatAccum(tall=True, sall=True, bins=bins)
sa += var[:7]
sa += var[7:]
istats = sa.get_stats()

# Unitest
result = []
for key in sorted(dstats.keys()):
    if cdms2.isVariable(dstats[key]):
        dstats[key] = dstats[key].asma()
    N.testing.assert_allclose(dstats[key], istats[key].asma())
