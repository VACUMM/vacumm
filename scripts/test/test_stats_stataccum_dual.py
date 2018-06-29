"""Test :class:`~vacumm.misc.stats.StatAccum` for a pair of variables"""

# Imports
from vcmq import MV2, N
from vacumm.misc.stats import StatAccum
import warnings
warnings.filterwarnings('error', "boolean index did not match", N.VisibleDeprecationWarning)

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
vmax = var2.max()
bins = N.linspace(-0.1*vmax, 0.9*vmax, 14)
nbins = len(bins)
maskyx = mask.all(axis=0)
maskt = mask.reshape((nt, -1)).all(axis=1)

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
    tmin = (varm1.min(axis=0), varm2.min(axis=0)),
    tmax = (varm1.max(axis=0), varm2.max(axis=0)),
    savail = 1.*(~var2dm1.mask).sum(axis=1)/var2dm1.shape[1],
    smean = (var2dm1.mean(axis=1), var2dm2.mean(axis=1)),
    sbias = var2dm2.mean(axis=1)-var2dm1.mean(axis=1),
    sstd = (var2dm1.std(axis=1), var2dm2.std(axis=1)),
    srms = N.ma.sqrt(((var2dm1-var2dm2)**2).mean(axis=1)),
    scrms = (var2dm1-var2dm2).std(axis=1),
    scorr = ((var2dm1.T-var2dm1.T.mean(axis=0))*(var2dm2.T-var2dm2.T.mean(axis=0))).mean(axis=0)/ \
        (var2dm1.std(axis=1)*var2dm2.std(axis=1)),
    smin = (var2dm1.min(axis=1), var2dm2.min(axis=1)),
    smax = (var2dm1.max(axis=1), var2dm2.max(axis=1)),
)
dstats['thist'] = N.ma.zeros((nbins-1, ny, nx), 'l'), N.ma.zeros((nbins-1, ny, nx), 'l')
dstats['shist'] = N.ma.zeros((nt, nbins-1), 'l'), N.ma.zeros((nt, nbins-1), 'l')
for ivar, varm in enumerate([varm1, varm2]):
    for ibin in xrange(nbins-1):
        valid = (varm>=bins[ibin])&(varm<bins[ibin+1])
        dstats['thist'][ivar][ibin] = valid.filled(False).sum(axis=0)
        dstats['shist'][ivar][:, ibin] = valid.filled(False).reshape((nt, -1)).sum(axis=1)
    dstats['thist'][ivar][:, maskyx] = N.ma.masked
    dstats['shist'][ivar][maskt, :] = N.ma.masked

# Indirect using StatAccum
sa = StatAccum(tall=True, sall=True, bins=bins)
sa += var1[:7], var2[:7]
sa += var1[7:], var2[7:]
istats = sa.get_stats()

# Unitest
result = []
for key in sorted(dstats.keys()):
    if isinstance(dstats[key], tuple):
        for i in 0, 1:
            value = N.ma.allclose(dstats[key][i], istats[key][i])
            result.append(('assertTrue', value))
    else:
        value = N.ma.allclose(dstats[key], istats[key])
        result.append(('assertTrue', value))
