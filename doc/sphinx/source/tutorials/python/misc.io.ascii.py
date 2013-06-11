# On recupere une variable cdms
import cdms2 as cdms
from vacumm.config import data_sample
f = cdms.open(data_sample('mars3d.t.nc'))
temp = f('temp')
u =  f('u')
v =  f('v')
f.close()

# On recupre le temps
time =  temp.getTime() # axe
ctime = time.asComponentTime() # temps cdtime.comptime()

# Creation a la main au format : YYYY/MM/DDZHH:MM TEMP U V
from vacumm.misc.atime import strftime, ch_units, strptime
f = open('misc.io.ascii.1.dat', 'w')
f.write('# Ligne de commentaire\n')
for it in xrange(len(temp)):
    t = strftime('%Y/%m/%dZ%H:%M',  ctime[it])
    f.write('%s %.4f %f %f\n' % (t, temp[it], u[it], v[it]))
f.close()

# Verification rapide (deux premieres lignes)
f = open('misc.io.ascii.1.dat')
print ''.join(f.readlines()[:3])
f.close()
#  -> # Ligne de commentaire
# -> 2008/01/07Z00:00 12.0672 0.359631 0.156422
# -> 2008/01/07Z01:00 11.9421 0.493174 0.158244

# Ecriture rapide via numpy
# - creation
import numpy as N
time_units = 'hours since %s'%ctime[0]
newtime = ch_units(time, time_units)[:]
data = N.array([newtime,temp.filled(999.), 
    u.filled(999.), v.filled(999.)],copy=0)
f = open('misc.io.ascii.2.dat', 'w')
f.write('# Ligne de commentaire\n')
N.savetxt(f, data.transpose(), fmt='%.3f', delimiter='\t')
#   note : on peut donner f ou le nom du fichier a N.savetxt()
# - verification
f = open('misc.io.ascii.2.dat')
print ''.join(f.readlines()[:2])
f.close()
#  -> # Ligne de commentaire
#  -> 0.000 12.067  0.360   0.156

# Lecture avancee
# - convertisseur pour le temps
def convtime(s):
    return strptime(s, '%Y/%m/%dZ%H:%M').torel(time_units).value
# - chargement partiel
tt, uu, vv = N.loadtxt('misc.io.ascii.1.dat',  comments='#', 
    usecols=[0, 2, 3], converters={0:convtime}, unpack=True)
# - verification
print tt[0], uu[0], vv[0]
#  -> 0.0 0.359631 0.156422
