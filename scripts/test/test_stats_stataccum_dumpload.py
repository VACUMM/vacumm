"""Test :class:`~vacumm.misc.stats.StatAccum` with dumping and loading """

# %% Imports
from vcmq import MV2, N, code_file_name, StatAccum, cdms2
from numpy.testing import assert_array_almost_equal
from numpy.random import seed

cdms2.setAutoBounds('off')

# %% Setup masked data
seed(0)
nt, ny, nx = 20, 3, 3
var1 = MV2.array(N.random.random((nt, ny, nx)))
var1.getAxis(0).designateTime()
var1.getAxis(0).units = 'months since 2000'
var1.getAxis(0).id = 'time'
var1.getAxis(1).id = 'lat'
var1.getAxis(1).designateLatitude()
var1.getAxis(2).id = 'lon'
var1.getAxis(2).designateLongitude()
var1.units = 'm'
var1.id = 'ssh'
var2 = var1.clone()
var2[:] = N.random.random((nt, ny, nx))
var1[3:13, :1, :1] = MV2.masked
var2[5:15, -1:, -1:] = MV2.masked
var2.long_name = 'Sea level'
var2.id = 'sla'
mask = var1.mask | var2.mask  # common mask
vmax = var2.max()
bins = N.linspace(-0.1*vmax, 0.9*vmax, 14)
nbins = len(bins)
restart_file5 = code_file_name(ext='5.nc')
restart_file7 = code_file_name(ext='7.nc')
print restart_file5

# %% Normal
sab = StatAccum(tall=True, sall=True, bins=bins)
sab += var1[:5], var2[:5]

# %% Dump
sab.dump(restart_file5)
sab += var1[5:7], var2[5:7]
sab.dump(restart_file7)

# %% Load from scratch
sa5 = StatAccum(restart=True, restart_file=restart_file5)
sa7 = StatAccum(restart=True, restart_file=restart_file7)

# %% Intermediate loads
sai = StatAccum(tall=True, sall=True, bins=bins)
sai += var1[:5], var2[:5]
print sai.load(restart_file7)  # erase
saw = StatAccum(tall=True, sall=True, bins=bins)
saw += var1[:7], var2[:7]
print saw.load(restart_file5)  # not loaded

# %% Finish stats
sab += var1[7:], var2[7:]
sa5 += var1[5:], var2[5:]
sa7 += var1[7:], var2[7:]
sai += var1[7:], var2[7:]

# %% Get stats
sab_stats = sab.get_stats()
sa5_stats = sa5.get_stats()
sa7_stats = sa7.get_stats()
sai_stats = sai.get_stats()

# %% Result
result = []
for sname in sab_stats.keys():
    for sat_stats in sa5_stats, sa7_stats, sai_stats:
        vrefs = sab_stats[sname]
        vtests = sat_stats[sname]
        if not isinstance(sab_stats[sname], tuple):
            vrefs = vrefs,
            vtests = vtests,
    for vref, vtest in zip(vrefs, vtests):
        assert_array_almost_equal(vref.asma(), vtest.asma())
