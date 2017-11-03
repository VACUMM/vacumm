# -*- coding: utf8 -*-
import cdms2
from vacumm.misc.plot import curve2
from vacumm.config import data_sample

# Trace de la moyenne spatiale
# - lecture
f = cdms2.open(data_sample('mars3d.t.nc'))
tsst = f('temp')
f.close()
# - plot
c = curve2(tsst, title=u'Série temporelle de SST',units = r'$^{\circ}C$',
    subplot=211, show=False, )

# Trace de la moyenne meridienne
# - lecture
f = cdms2.open(data_sample('mars3d.xy.nc'))
zsst = cdms2.MV2.average(f('temp'), axis=0)
f.close()
# - plot
curve2(zsst, title=u'SST zonale', color='r', show=False, close=True,
    subplot=212,top=.9,hspace=.4,left=.15,bottom=.07,
    savefigs=__file__)




