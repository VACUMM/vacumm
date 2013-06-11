# -*- coding: utf8 -*-
# Read sea level at Brest
from vacumm.tide.sonel_mareg import get_slv_shom
import MV2
sea_level = get_slv_shom('BREST', time_range=('2006-08', '2006-09', 'co'))
sea_level[:] -= MV2.average(sea_level)

# Filtering
from vacumm.tide.filters import demerliac, godin
cotes, tide = demerliac(sea_level, get_tide=True)

# Plots
from vacumm.misc.plot import curve, savefigs
import pylab as P
kwplot = dict(date_fmt='%d/%m', show=False, date_rotation=0)
# - tidal signal
curve(sea_level, 'k', subplot=211, label='Original', title='Sea level at Brest',**kwplot)
curve(tide, 'b', label='Tidal signal', **kwplot)
P.legend().legendPatch.set_alpha(.7)
# - surcotes/decotes
curve(cotes, 'r', subplot=212, hspace=.3, label='Demerliac', **kwplot)
P.legend().legendPatch.set_alpha(.7)
savefigs(__file__, savefigs_pdf=True)
