# -*- coding: utf8 -*-
# Read sea level at Brest
from vcmq import cdms2, P, curve2, savefigs, data_sample
f = cdms2.open(data_sample("tide.sealevel.BREST.mars.nc"))
sea_level = f('sea_level')
f.close()

# Filtering
from vacumm.tide.filters import demerliac
cotes, tide = demerliac(sea_level, get_tide=True)

# Plots
kwplot = dict(date_fmt='%d/%m', show=False, date_rotation=0)
# - tidal signal
curve2(sea_level, 'k', subplot=211, label='Original', title='Sea level at Brest',**kwplot)
curve2(tide, 'b', label='Tidal signal', **kwplot)
P.legend().legendPatch.set_alpha(.7)
# - surcotes/decotes
curve2(cotes, 'r', subplot=212, hspace=.3, label='Demerliac', **kwplot)
P.legend().legendPatch.set_alpha(.7)
savefigs(__file__, savefigs_pdf=True)
P.close()
