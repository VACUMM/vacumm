# Lecture et masquage de la vitesse
import cdms2
from vacumm.config import data_sample
zoom = dict(lat=(48.3, 48.55), lon=(-5.2, -4.8))
f = cdms2.open(data_sample('mars3d.xy.nc'))
u = f('u', **zoom) 
v = f('v', **zoom) 
f.close()
for var in u, v:
    var[:]*100.
    var[:] = cdms2.MV.masked_object(var, 0., copy=0)

# Trace des vecteurs
from vacumm.misc.plot import map2
import pylab as P ; P.rc('font', size=9)
map2((u, v),  title='Vitesses en surface', show=False, 
    units='$cm s^{-1}$', figsize=(6, 5.), 
    quiver_alpha=.5, quiver_samp=2, 
    quiver_scale=15,  quiver_width=0.006, #quiver_headaxislength=2, 
    quiverkey_value = 2., nofill=True, 
    left=.08, right=1., top=.9, proj='merc', quiver_norm=1, 
    savefigs=__file__, savefigs_pdf=True,)

