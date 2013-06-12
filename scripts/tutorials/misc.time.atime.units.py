# On definit des unites
units = 'hours since 2000-01-15 06:00'

# Le format est-il bon ?
from vacumm.misc.atime import *
print are_good_units(units)
#  -> True

# Memes unites ?
print are_same_units('hours since 2000-1-15 06', units)
#  -> True

# Changer les unites d'un axe de temps
from vacumm.misc.axes import create_time
import numpy as N
taxis = create_time(N.arange(6.)*48, units)
# - avant changement
print taxis.units, taxis[0:2]
print taxis.asComponentTime()[0:2]
#  -> hours since 2000-01-15 06:00 [  0.  48.]
#  -> [2000-1-15 6:0:0.0, 2000-1-17 6:0:0.0]
# - changement
ch_units(taxis, 'days since 2000-1-15 06', copy=0)
# - apres changement
print taxis.units, taxis[0:2]
print taxis.asComponentTime()[0:2]
#  -> days since 2000-1-15 06 [ 0.  2.]
#  -> [2000-1-15 6:0:0.0, 2000-1-17 6:0:0.0]

# Le temps matplotlib
taxis_mpl = mpl(taxis)
print taxis_mpl[0], taxis_mpl.units
#  -> 730134.25 days since 0001

# Changer les unites de temps d'une variable
import MV2
var = MV2.array(MV2.arange(len(taxis)), dtype='f', axes=[taxis])
ch_units(var, 'hours since 2000-01-15 06')
print var.getTime()[0:2]
#  -> [  0.  48.]


# Changements de fuseau horaire
# - maintenant a l'heure UTC
t_utc = now(True)
print strftime('%H:%M', t_utc), t_utc.hour
#  -> 15:54 15
# - heure de paris
t_paris = utc_to_paris(t_utc)
print strftime('%H:%M', t_paris), t_paris.hour
#  -> 17:54 17
# - retour en UTC
print tz_to_tz(t_paris,'Europe/Paris','UTC').hour
#  -> 15
# - travail sur une chaine de caracteres !
print to_utc('2000-10-01 10:20', 'Europe/Paris')
#  -> '2000-10-1 8:20:0.0'
