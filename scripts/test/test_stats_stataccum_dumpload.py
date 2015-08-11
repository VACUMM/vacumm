"""Test :class:`~vacumm.misc.stats.StatAccum` with dumping and loading """

# Imports
from vcmq import MV2, N, code_file_name
from vacumm.misc.stats import StatAccum

# Setup masked data
nt, ny, nx = 20, 15, 10
var1 = MV2.arange(1.*nt*ny*nx).reshape(nt, ny, nx)
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
var2[:] += N.ma.sin(var1)*100
var1[3:13, 5:9, 3:7] = MV2.masked
var2[5:15, 7:12, 1:5] = MV2.masked
var2.long_name = 'Sea level'
var2.id = 'sla'
mask = var1.mask|var2.mask # common mask
vmax = var2.max()
bins = N.linspace(-0.1*vmax, 0.9*vmax, 14)
nbins = len(bins)
restart_file = code_file_name(ext='nc')
print restart_file

# Normal
sa0 = StatAccum(tall=True, sall=True, bins=bins)
sa0 += var1[:7], var2[:7]

# Dump
sa0.dump(restart_file)

# Load from scratach
sa1 = StatAccum(restart=True, restart_file=restart_file)

## Intermediate load
#sa1 = StatAccum(tall=True, sall=True, bins=bins)
#sa1 += var1[:5], var2[:5]
#sa1.load(restart_file)
#
# Finish stats
sa0 += var1[7:], var2[7:]
sa1 += var1[7:], var2[7:]
#sa2 += var1[7:], var2[7:]

# Check reloading
iterindex = sa1.load(restart_file)
print iterindex

# Get stats
sa0_stats = sa0.get_stats()
sa1_stats = sa1.get_stats()

# Result
result = []
for sname in sa0_stats.keys():
    result.append(('assertTrue', N.ma.allclose(sa0_stats[sname], sa1_stats[sname])))
    print sname, N.ma.allclose(sa0_stats[sname], sa1_stats[sname])

print 'Done', restart_file
