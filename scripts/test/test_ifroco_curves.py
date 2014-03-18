# -*- coding: utf8 -*-
"""Plot a PREVIMER curve"""

import matplotlib
matplotlib.use('agg')

from matplotlib.dates import HourLocator

# Some params
field = 'sst'
ncpat = "/home/coriolis_exp/spool/co01/co0123/co012302/co01230207_v2/d1/PREVIMER_D1-MARS3D-R2/best_estimate/2013/PREVIMER_D1-MARS3D-R2_%Y%m%dT%H%MZ.nc"
ndayfore=3
logofile = "logo2_previmer_cetmef.png"
copyright = u"© Previmer"

# Imports
from vacumm.report.ifroco.curves import round_date, P, DS, paris_to_utc, now, \
    cfgget, N, grid2xy, plot_curves, cdtime, utc_to_paris, Day12hFormatter

# Time
t0 = round_date(now(utc=False), 'day', 'floor')
t1 = t0.add(ndayfore, cdtime.Day)
t0 = paris_to_utc(t0)
t1 = paris_to_utc(t1)

# Read
ds = DS(ncpat, 'mars', time=(t0, t1), logger_level='critical')
data = ds.get(field)
utc_to_paris(data)

# Read some config values
long_name = cfgget('long_name', field)
long_long_name = cfgget('long_name_model', field)
units = cfgget('units', field)
vmin = cfgget('vmin', field)
vmax = cfgget('vmax', field)
points = cfgget('points', field)

# Interpolate to points
lons = N.array([d[0] for d in points])
lats = N.array([d[1] for d in points])
ptdata = grid2xy(data, lons, lats)
datas = []
for i in xrange(ptdata.shape[1]):
    datas.append(ptdata[:, i])
    datas[i].long_name = points[i][2].decode('utf8')

# Plot
x = plot_curves(datas, vmin=vmin, vmax=vmax, units=units, long_name=long_name, 
    title=long_long_name+u' aux bouées', logos=logofile, copyright=copyright, 
    date_locator=HourLocator(byhour=[0, 12]), date_formatter=Day12hFormatter('%a %d/%m'), 
    hlitvs_units='day', savefig=__file__,
#    ylocator=MaxNlocator(integer=True), yminor_locator=MaxNlocator(integer=True, steps=2), 
    close=False
    )
#P.show()
