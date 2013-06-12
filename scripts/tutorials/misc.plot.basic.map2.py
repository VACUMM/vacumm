# Lecture et masquage de la SST
import cdms2
from vacumm.config import data_sample
f = cdms2.open(data_sample('mars3d.xy.nc'))
sst = f('temp', time=slice(0, 1),  lat=(47.8, 48.6), lon=(-5.2, -4.25),  squeeze=1) 
f.close()
import pylab as P
# Trace de la carte
from vacumm.misc.plot import map2 as map
m = map(sst, title='SST en iroise', 
    vmin=9, top=.9, figsize=(6, 5), clabel_glow=True, contour_alpha=1, 
    show=False, savefigs=__file__)
