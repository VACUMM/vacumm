# Lecture et masquage de la vitesse
import cdms2
from vacumm.config import data_sample
f = cdms2.open(data_sample('mars3d.t.nc'))
u = f('u') 
v = f('v') 
f.close()


# Trace de la carte
from vacumm.misc.plot import stick2
stick2(u, v,  title='Courants en Iroise', units='m/s', 
    bottom=.2, top=.85, quiver_headwidth=2, 
    quiverkey_value=.5, quiver_width=0.002, quiver_scale=5., 
    savefigs=__file__, figsize=(5.5,  3),  show=False)

