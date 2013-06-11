# Lecture de la meteo
import cdms2, numpy as N, MV2
from vacumm.config import data_sample
select = dict(time=('2008-10-12', '2008-10-15', 'cc'), squeeze=1)
f = cdms2.open(data_sample('wrf_1d.nc'))
u = f('u10m', **select)
v = f('v10m', **select)
p = f('psfc', **select)
r = f('rain', **select)
f.close()
nt = len(u)


# Plots
from vacumm.misc.plot import curve2, stick2, bar2
import matplotlib.pyplot as P
P.rc('font', size=9)
P.figure(figsize=(5.5, 6))
    
# - vecteurs
colors = N.array(["#000000"]*nt) # noir a la base
colors[MV2.sqrt(u**2+v**2).filled(0.)>5.] = '#ff0000'         # rouge si module > 4.
qv = stick2(u, v, color=colors.tolist(), xticklabels=False, subplot=211, 
    shadow=True, right=.85, hspace=.2, left=.14, hldays=True, 
    key=1, key_size=11, key_color='.3', ylabel_color='r', units='m.s^{-1}', 
    mod=True, mod_color='.4', mod_linewidth=2, title='Vent', quiver_scale = 50., 
    quiver_headwidth=3,quiver_headlength=3,quiver_headaxislength=2,
    quiverkey_pos=(.05, 1.05), quiverkey_color='k', 
    quiver_width=.006, alpha=.8, show=False)
fbox = qv.axes.get_position(True) # positions pour ajouter Nord et Sud
qv.figtext(fbox.xmax,fbox.ymin,' Sud',fontweight='heavy', ha='left', size=11)
qv.figtext(fbox.xmax,fbox.ymax,' Nord',va='top', ha='left', fontweight='heavy', size=11)


# - precipitations
print 'precip'
dr = r.clone() # variation != accumulation
dr[:-1] = N.diff(r) ; dr[-1] = r[-1]-r[-2] 
lr = bar2(dr, width=.8, subplot=212, hldays=True, 
    long_name='Precipitations', shadow=True, ylabel = '%(long_name)s [$%(units)s$]', 
    title=False, dayhl=True, ylabel_color='#008888',zorder=100, 
    key=2, key_size=11, key_color='.3', log=True, ymin=0.01,
    color='#00ffff', linewidth=.2, edgecolor='#888888', show=False)
lr.axes.yaxis.grid(False)

print 'pres'
# - pression
p[:] /= 100.
lp = curve2(p, color='g', vminmax=1015, ymaxmin=1025, twin='x', 
    zorder=150, title=False, ylabel = '%(long_name)s [$%(units)s$]', 
    ylabel_color='g', shadow=True, show=False)
    
lp.savefigs(__file__, pdf=True)

